#' @title Get gene expression by featureCounts
#' @description save featureCounts output to work_path/expr
#'
#' @param paras list fromJSON(paras_file)
#' Default bam_path, gtf, ref, script, featureCounts_path,
#' dataset, and core from paras
#' @param bam_path directory to single cell bam file
#' @param gtf the gene annotation in gtf format
#' @param ref fasta file
#' @param dataset the name of dataset
#' @param pair type of read: either "paired" for paired-end or "single" for single-end
#' @param core the number of threads

#' @return featureCounts output path
#' @export

getGeneExpression <- function(
    paras, bam_path=paras$Basic$bam_path,
    gtf=paras$Basic$refgenome$gtf_path,
    ref=paras$Basic$refgenome$ref_path,
    dataset=paras$DataSet, pair=paras$Basic$paired,
    core=paras$Basic$core) {
    #script
    dir_shell <- system.file("shell", package = "SCSES")
    filename <- "run_featurecounts.sh"
    script <- file.path(dir_shell, filename)
    featureCounts_path=paras$Basic$featureCounts_path
    # output
    work_path = paras$Basic$work_path
    featurecounts.work.path <- paste0(work_path, "/expr/")
    log_file <- paste0(work_path, "/runfeatureCounts.log")
    cmd <- paste(
        "bash", script,
        featurecounts.work.path,
        ref,
        gtf,
        bam_path,
        core,
        dataset,
        pair,
        featureCounts_path, ">>", log_file, "2>&1"
    )
    msg <- paste0("[", Sys.time(), "] ", "Detect gene expression: ", cmd)
    print(msg)
    res <- system(command = cmd, intern = T, wait = T)
    if (!is.null(attributes(res)) && attributes(res)$status == 1) {
        msg <- paste0("[", Sys.time(), "] ", "Run featureCounts Error.")
        stop(msg)
    } else {
        msg <- paste0("[", Sys.time(), "] ", "Detect gene expression Finish.")
        print(msg)
        return(featurecounts.work.path)
    }
}

#' @title Get gene expression matrix
#' @description save gene expression count and TPM matrix to work_path/rds/
#'
#'
#' @param paras list fromJSON(paras_file)
#' Default dataset, filter.mt, filter.rp from paras
#' @param expr_path directory to the featureCounts output
#' @param dataset the name of dataset
#' @param filter.mt filter out mitochondrial genes
#' @param filter.rp filter out ribosomal genes

#' @return gene expression matrix path
#' @export
#'
#' @importFrom stats median

getEXPmatrix <- function(
    paras, expr_path = paste0(paras$Basic$work_path, "/expr/"),
    dataset = paras$DataSet,
    filter.mt = paras$Basic$filter_sc$filter.mt,
    filter.rp = paras$Basic$filter_sc$filter.rp) {
    fea_file <- paste0(expr_path, "/", dataset, "_count.txt")
    output_path <- paste0(dirname(expr_path), "/rds/")
    if (!dir.exists(output_path)) {
        dir.create(output_path)
    }
    fea <- read.table(fea_file, comment.char = "#", sep = "\t", header = T, check.names = F)
    row.names(fea) <- fea$Geneid
    # filter Mitochondrial genes
    if (filter.mt) {
        fea <- fea[!grepl("^MT-|^mt-", row.names(fea)), ]
    }
    # filter ribosomal genes
    if (filter.rp) {
        fea <- fea[!grepl("^RP[SL]|^Rp[sl]", row.names(fea)), ]
    }
    fea_kb <- fea$Length / 1000
    fea <- fea[, -c(1:6), drop = F]
    fea_rpk <- fea / fea_kb
    fea_tpm <- t(t(fea_rpk) / colSums(fea_rpk) * 10^6)
    fea_tpm <- as.data.frame(fea_tpm)
    colnames(fea) <- sapply(colnames(fea), function(x) {
        unlist(strsplit(x, "/"))[length(unlist(strsplit(x, "/")))]
    })
    colnames(fea_tpm) <- sapply(colnames(fea_tpm), function(x) {
        unlist(strsplit(x, "/"))[length(unlist(strsplit(x, "/")))]
    })
    fea_tpm <- log2(fea_tpm + 1)
    saveRDS(fea, paste0(output_path, "/count.rds"))
    saveRDS(fea_tpm, paste0(output_path, "/TPM.rds"))
    return(output_path)
}


#' @title Get gene expression matrix from Cell Ranger
#' @description Read count matrix from 10X CellRanger hdf5 file
#' save raw and normalized UMI counts sparse matrix to work_path/rds/
#'
#'
#' @param paras list fromJSON(paras_file)
#' Default dataset, filter.mt, filter.rp from paras
#' @param expr_path directory to the cellranger output
#' @param dataset the name of dataset
#' @param filter.mt filter out mitochondrial genes
#' @param filter.rp filter out ribosomal genes
#' @param sample_name the name of samples in the dataset
#' If the dataset contains multiple samples, the Cell Ranger outputs for
#' different samples need to be placed in different folders under the
#' expr_path directory, named after the sample_name.

#' @return gene expression matrix path
#' @export
#'
#' @importFrom stats median

get10XEXPmatrix <- function(
    paras, expr_path,
    dataset = paras$DataSet,
    filter.mt = paras$Basic$filter_sc$filter.mt,
    filter.rp = paras$Basic$filter_sc$filter.rp,sample_name) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Please install Seurat to read files")
    }
    # input
    scelist = lapply(sample_name, function(x) {
        print(x)
        file_name = list.files(paste0(expr_path,'/',x,'/'),
            pattern = "filtered_feature_bc_matrix.h5",
            recursive = TRUE, full.names = TRUE
        )
        if (length(file_name)==0) {
            stop("File not found")
        }
        sce = Seurat::Read10X_h5(file_name)
        colnames(sce) = paste0(x, "_", colnames(sce))
        return(sce)
    })
    expr = do.call(what = cbind, args = scelist)
    # output
    output_path <- paste0(dirname(expr_path), "/rds/")
    if (!dir.exists(output_path)) {
        dir.create(output_path)
    }
    # filter Mitochondrial genes
    if (filter.mt) {
        expr <- expr[!grepl("^MT-|^mt-", row.names(expr)), ]
    }
    # filter ribosomal genes
    if (filter.rp) {
        expr <- expr[!grepl("^RP[SL]|^Rp[sl]", row.names(expr)), ]
    }
    library_size <- colSums(expr, na.rm = T)
    median_transcript_count <- median(library_size)
    expr_norm <- t(t(expr) * (median_transcript_count / library_size))
    expr_norm <- log2(expr_norm + 1)
    saveRDS(expr, paste0(output_path, "/count.rds"))
    saveRDS(expr_norm, paste0(output_path, "/count_norm.rds"))
    return(output_path)
}