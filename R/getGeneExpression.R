#' @title Get gene expression by featureCounts
#'
#' @description  This function runs featureCounts to quantify gene
#' expression from BAM files using a shell script. It processes
#' single-cell RNA-seq BAM files against a reference genome annotation
#' to generate count matrices, and saves the featureCounts output
#' to the specified working directory.
#'
#' @param paras paras A list object parsed from SCSES JSON parameter file using
#'   \code{readSCSESconfig()}.
#' @param bam_path Character string specifying the directory path containing
#'   single-cell BAM files. Default is taken from \code{paras$Basic$bam_path}.
#' @param gtf Character string specifying the path to the gene annotation
#'   file in GTF format. Default is taken from
#'   \code{paras$Basic$refgenome$gtf_path}.
#' @param ref Character string specifying the path to the reference genome
#'   FASTA file. Default is taken from \code{paras$Basic$refgenome$ref_path}.
#' @param dataset Character string specifying the dataset name, used in
#'   output file naming. Default is taken from \code{paras$DataSet}.
#' @param pair Character string indicating the sequencing type. Should be
#'   either "paired" for paired-end reads or "single" for single-end reads.
#'   Default is taken from \code{paras$Basic$paired}.
#' @param core Integer specifying the number of CPU threads to use for
#'   parallel processing. Default is taken from \code{paras$Basic$core}.
#'
#' @return Character string of the featureCounts output directory path
#'   (\code{work_path/expr/}) where the count files are saved.
#'
#' The output from this function can be used as input for \code{\link{getEXPmatrix}}
#' to generate processed expression matrices.
#'
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

#' @title Get gene expression matrix from featureCounts output
#' @description This function reads gene expression count data
#' from featureCounts output file, performs optional gene filtering
#' (mitochondrial and ribosomal genes), calculates TPM (Transcripts
#' Per Million) normalization, applies log2 transformation, and saves
#' both raw count and normalized TPM matrices to the specified output directory.
#'
#' @param paras paras A list object parsed from SCSES JSON parameter file using
#'   \code{readSCSESconfig()}.
#' @param expr_path Character string specifying the directory path to the
#'   featureCounts output files. Default is constructed from
#'   \code{paras$Basic$work_path/expr/}.
#' @param dataset Character string specifying the dataset name, used to
#'   construct the input filename as \code{dataset_count.txt}.
#'   Default is taken from \code{paras$DataSet}.
#' @param filter.mt Logical value indicating whether to filter out
#'   mitochondrial genes (genes with names starting with "MT-" or "mt-").
#'   Default: TRUE (from \code{paras$Basic$filter_sc$filter.mt}).
#' @param filter.rp Logical value indicating whether to filter out
#'   ribosomal genes (genes with names starting with "RPS", "RPL", "Rps",
#'   or "Rpl"). Default: TRUE (from \code{paras$Basic$filter_sc$filter.rp}).
#'
#' @return Character string of the output directory path where the RDS files
#'   are saved. The directory contains two files: \code{count.rds} (raw counts)
#'   and \code{count_norm.rds} (log2-transformed TPM values).
#' @note
#' The input file should be a standard featureCounts output with tab-separated
#' values containing columns: Geneid, Chr, Start, End, Strand, Length, followed
#' by sample count columns. The required input files can be generated using
#' \code{\link{getGeneExpression}} function. The output directory \code{rds/}
#' is created in \code{paras$Basic$work_path} if it doesn't exist.
#'
#' @export
#'
#' @importFrom stats median

getEXPmatrix <- function(
    paras, expr_path = paste0(paras$Basic$work_path, "/expr/"),
    dataset = paras$DataSet,
    filter.mt = paras$Basic$filter_sc$filter.mt,
    filter.rp = paras$Basic$filter_sc$filter.rp) {
    fea_file <- paste0(expr_path, "/", dataset, "_count.txt")
    output_path <- paste0(paras$Basic$work_path, "/rds/")
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
    saveRDS(fea_tpm, paste0(output_path, "/count_norm.rds"))
    return(output_path)
}


#' @title Get gene expression matrix from Cell Ranger
#' @description This function reads count matrix from 10X CellRanger HDF5 file, performs
#' optional gene filtering (mitochondrial and ribosomal genes), applies median
#' normalization and log2 transformation, then saves both raw and normalized
#' UMI count sparse matrices to the specified output directory.
#'
#' @param paras A list object parsed from SCSES JSON parameter file using
#'   \code{readSCSESconfig()}.
#' @param expr_path Character string specifying the directory path to the
#'   Cell Ranger output folder containing the HDF5 file.
#' @param filter.mt Logical value indicating whether to filter out
#'   mitochondrial genes (genes with names starting with "MT-" or "mt-").
#'   Default: TRUE (from \code{paras$Basic$filter_sc$filter.mt}).
#' @param filter.rp Logical value indicating whether to filter out
#'   ribosomal genes (genes with names starting with "RPS", "RPL", "Rps",
#'   or "Rpl"). Default: TRUE (from \code{paras$Basic$filter_sc$filter.rp}).
#' @param cells Character vector of cell barcodes to subset. If provided,
#'   only cells with barcodes present in this vector will be retained.
#'   If \code{NULL} (default), all cells in the matrix are kept.
#'
#' @return Character string of the output directory path where the RDS files
#'   are saved. The directory contains two files: \code{count.rds} (raw counts)
#'   and \code{count_norm.rds} (normalized and log2-transformed counts).
#'
#' @note
#' This function requires the Seurat package to be installed. The output
#' directory \code{rds/} is created in the \code{paras$Basic$work_path}
#' @export
#'
#' @importFrom stats median

get10XEXPmatrix <- function(
    paras, expr_path,
    filter.mt = paras$Basic$filter_sc$filter.mt,
    filter.rp = paras$Basic$filter_sc$filter.rp,
    cells = NULL) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Please install Seurat to read files")
  }
  # input
  file_name = list.files(expr_path,
                         pattern = "filtered_feature_bc_matrix.h5",
                         recursive = TRUE, full.names = TRUE
  )
  if (length(file_name)==0) {
    stop("filtered_feature_bc_matrix.h5 File not found")
  }
  expr = Seurat::Read10X_h5(file_name)
  if(!is.null(cells)){
    expr=expr[,which(colnames(expr)%in%cells)]
  }
  # output
  output_path <- paste0(paras$Basic$work_path, "/rds/")
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



#' @title Get gene expression matrix from Cell Ranger of multiple samples
#' @description Read count matrix from 10X CellRanger hdf5 file
#' save raw and normalized UMI counts sparse matrix to work_path/rds/
#'
#'
#' @param paras A list object parsed from SCSES JSON parameter file using
#'   \code{readSCSESconfig()}.
#' @param expr_path Character string specifying the directory path to the
#'   Cell Ranger output folder containing the HDF5 file.
#' @param filter.mt Logical value indicating whether to filter out
#'   mitochondrial genes (genes with names starting with "MT-" or "mt-").
#'   Default: TRUE (from \code{paras$Basic$filter_sc$filter.mt}).
#' @param filter.rp Logical value indicating whether to filter out
#'   ribosomal genes (genes with names starting with "RPS", "RPL", "Rps",
#'   or "Rpl"). Default: TRUE (from \code{paras$Basic$filter_sc$filter.rp}).
#' @param sample_name the name of directory for different samples in the expr_path,
#' different samples seperated by ";"
#' If the dataset contains multiple samples, the Cell Ranger outputs for
#' different samples need to be placed in different folders under the
#' expr_path directory, named after the sample_name.

#' @return gene expression matrix path
#' @export
#'
#' @importFrom stats median

getmulti10XEXPmatrix <- function(
    paras, expr_path,
    filter.mt = paras$Basic$filter_sc$filter.mt,
    filter.rp = paras$Basic$filter_sc$filter.rp,sample_name) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Please install Seurat to read files")
    }
    # input
    sample_name=unlist(strsplit(sample_name,";"))
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
