#' @include utilities.R
#' @include eventFeatureFunction.R
#'

#' @title PCA
#' @description Perform a principal components analysis on the given data
#' and return the results
#' @param data a matrix or data frame
#' @param id_select used features
#'
#' @return principal components retain ~90% of the variation in the data
#'
#' @export
#' @importFrom stats prcomp var

PCA_D_reduct <- function(data, id_select) {
    data <- data[id_select, ]
    D_Reduct_res <- prcomp(t(data[which(apply(data, 1, var) != 0), ]), center = TRUE, scale. = TRUE)
    sdev <- D_Reduct_res$sdev
    var_prop <- sdev^2 / sum(sdev^2)
    cumulative_variance <- cumsum(var_prop)
    n_components <- which(cumulative_variance >= 0.9)[1]
    D_Reduct <- D_Reduct_res$x[, 1:n_components]
    return(D_Reduct)
}

#' @title PCA for sparse matrix
#' @description Perform a principal components analysis on the given sparse
#' matrix and return the results
#' @param data a sparse matrix
#' @param id_select used features
#'
#' @return principal components retain ~90% of the variation in the data
#'
#' @export
#' @importFrom irlba irlba
#' @importFrom stats var

PCA_D_reduct_sparse <- function(data, id_select) {
    data <- data[id_select, ]
    data <- data[which(apply(data, 1, var) != 0), ]
    data <- t(scale(t(data)))
    pca.results <- irlba(A = t(data), nv = min(100, nrow(data) - 1, ncol(data) - 1))
    u <- pca.results$u
    s <- pca.results$d
    v <- pca.results$v

    feature.loadings <- pca.results$v # 右矩阵v是 gene 的loading
    sdev <- s / sqrt(max(1, ncol(data) - 1)) # 特征根 / 根号(列数 - 1)
    # 左奇异矩阵 u 是 cell embeddings，乘以 var 权重
    cell.embeddings <- u %*% diag(s)
    explained_variance <- (sdev^2) / sum(sdev^2)
    cumulative_variance <- cumsum(explained_variance)
    n_components <- which(cumulative_variance >= 0.9)[1]
    D_reduct_res <- cell.embeddings[, 1:n_components]
    return(D_reduct_res)
}

#' @title Create R package for a reference sequence file
#' @description BSgenome used for sequence extraction
#' @param ref_path reference fasta file
#' @param out_path path to save the package
#' @param pkg the name of pkg
#' @param install_lib library path that save the package
#'
#' @return package path
#'
#' @export
#' @importFrom Biostrings readDNAStringSet writeXStringSet
#' @importFrom BSgenome forgeBSgenomeDataPkg
#'
createBSgenome <- function(ref_path,out_path,pkg, install_lib = NULL) {
    out_path = paste0(out_path, "/", pkg, "/")
    dir.create(out_path, recursive = T)
    # 0. Create seperated reference sequence file----
    fa <- readDNAStringSet(ref_path, format = "fasta", nrec = -1)
    chrs <- names(fa)
    chrs <- unlist(lapply(X = as.list(chrs), FUN = function(x) {
        return(unlist(strsplit(x = x, split = " "))[1])
    }))
    names(fa) <- chrs
    for (chr in chrs)
    {
        writeXStringSet(x = fa[chr], filepath = paste0(out_path, "/", chr, ".fa"), format = "fasta")
    }
    # 1. Prepare Seed Files----
    content <- ""
    content <- paste0(content, "Package: ", pkg, "\n")
    content <- paste0(content, "Title: ", "Created reference BSgenome object", "\n")
    content <- paste0(content, "Description: ", "This is a self created file", "\n")
    content <- paste0(content, "Version: ", "1.0.0", "\n")
    content <- paste0(content, "organism: ", "self defination", "\n")
    content <- paste0(content, "common_name: ", "NA", "\n")
    content <- paste0(content, "genome: ", "NA", "\n")
    content <- paste0(content, "provider: ", "self", "\n")
    content <- paste0(content, "release_date: ", Sys.Date(), "\n")
    content <- paste0(content, "source_url: ", normalizePath(out_path, "/"), "\n")
    content <- paste0(content, "organism_biocview: ", "NA", "\n")
    content <- paste0(content, "BSgenomeObjname: ", pkg, "\n")
    content <- paste0(content, "seqnames: ", paste0("c(", paste(paste0("\"", chrs, "\""), collapse = ","), ")"), "\n")
    content <- paste0(content, "circ_seqs: ", "character(0)", "\n")
    content <- paste0(content, "seqs_srcdir:", out_path, "\n")
    content <- paste0(content, "seqfiles_suffix: ", ".fa", "\n")
    write(x = content, file = paste0(out_path, "/seed.txt"))

    # 2. Build and Install package----
    msg = paste0("[", Sys.time(), "] ", "Build ", pkg, " package...", "")
    print(msg)

    forgeBSgenomeDataPkg(paste0(out_path, "/seed.txt"), destdir = out_path)

    log_file <- paste0(out_path, "/buildpkg.log")
    cmd <- paste("R CMD build", paste0(out_path, "/", pkg),
                 "--no-build-vignettes --no-manual", ">>", log_file, "2>&1")
    system(command = cmd, wait = T)

    cmd <- paste("R CMD check", paste0(out_path, "/", pkg),
                 "--no-vignettes --no-manual", ">>", log_file, "2>&1")
    system(command = cmd, wait = T)

    cmd <- paste("R CMD REMOVE", pkg, ">>", log_file, "2>&1")
    system(command = cmd, wait = T)

    msg = paste0("[", Sys.time(), "] ", "Install ", pkg, " package...", "")
    print(msg)

    log_file <- paste0(out_path, "/installpkg.log")

    if (!is.null(install_lib)) {
      cmd <- paste("R CMD INSTALL", paste0(pkg, "_1.0.0.tar.gz"),
                   "-l", install_lib, ">>", log_file, "2>&1")
    } else {
      cmd <- paste("R CMD INSTALL", paste0(pkg, "_1.0.0.tar.gz"),
                   ">>", log_file, "2>&1")
    }

    # msg = paste0("[", Sys.time(), "] ", "Install ", pkg, " package...", "")
    # print(msg)
    #
    # log_file <- paste0(out_path, "/installpkg.log")
    # cmd <- paste("R CMD INSTALL", paste0(pkg, "_1.0.0.tar.gz"),
    #              ">>", log_file, "2>&1")
    system(command = cmd, wait = T)
    file.remove(paste0(pkg, "_1.0.0.tar.gz"))
    unlink(paste0(pkg, ".Rcheck"),recursive=TRUE)
    return(out_path)
}

#' @title Construct KNN Graph for Single Cells
#'
#' @description This function constructs k-nearest neighbor (KNN) graphs for single cells
#'   based on multiple data modalities including RNA-binding protein (RBP) expression,
#'   raw junction read counts, and PSI values. It computes cell-to-cell distances
#'   using specified distance metrics and applies random walk-based similarity refinement
#'   to capture higher-order relationships in the cellular neighborhood structure.
#'   The resulting similarity matrices are used for downstream imputation.
#'
#' @param paras A list object containing SCSES configuration parameters, typically
#'   loaded using \code{readSCSESconfig(paras_file)}.
#'
#' @param rds_path Character string specifying the path to processed RDS data files.
#'   This directory should contain the preprocessed single-cell data matrices
#'   (expr.rds, rc.rds, psi.rds) generated by upstream processing functions.
#'   Default: extracted from \code{paras$Basic$work_path/rds_processed/}.
#'
#' @param output_path Character string specifying the output directory for cell
#'   similarity results. Default: \code{work_path/imputation/cell_similarity/}.
#'
#' @param feature_num Integer specifying the number of highly variable features
#'   to select. Default: 1000 (extracted from \code{paras$Task$impute$feature_num}).
#'
#' @param rbp Character string or file path specifying RNA-binding protein (RBP)
#'   information. Can be either:
#'   \itemize{
#'     \item File path to RBP annotation file
#'     \item Character string containing RBP gene symbols
#'     \item \code{NULL} to use default RBP set from configuration
#'   }
#'   Default: \code{paras$Task$impute$rbp}.
#'
#' @param cell_similarity_data Character string specifying features for
#'   cell similarity calculation, separated by semicolons.Valid options include:
#'   \itemize{
#'     \item \code{"EXP_RBP"}: TPM/normalized UMI counts for RNA-binding proteins
#'       (source: rds_processed/expr.rds)
#'     \item \code{"RC"}: Raw junction read counts for splicing events
#'       (source: rds_processed/rc.rds)
#'     \item \code{"PSI"}: PSI values for splicing events
#'       (source: rds_processed/psi.rds)
#'   }
#'   Example: \code{"EXP_RBP;RC;PSI"}
#'
#' @param distance_method Character string specifying the distance metric for
#'   cell-to-cell similarity calculation. Supported methods:
#'   \itemize{
#'     \item \code{"euclidean"}
#'     \item \code{"cosine"}
#'   }
#'   Default: \code{"euclidean"}.
#'
#' @param alpha_cell Numeric value between 0 and 1 specifying the random walk
#' probability (1-restart probability).
#' Default: 0.8 (from \code{paras$Task$impute$KNN$cell$alpha}).
#'
#' @param decay_cell Numeric threshold of change in the similarity matrix.
#' The algorithm stops when the change falls below this threshold between
#' consecutive iterations. Default: 0.05 (from \code{paras$Task$impute$KNN$cell$decay}).
#'
#' @param kcell_max Integer specifying the maximum number of nearest neighbors
#'   to consider for each cell.
#'   Default: 50 (from \code{paras$Task$impute$KNN$cell$kmax}).
#'
#' @param kcell_min Integer specifying the minimum number of nearest neighbors
#'   to maintain for each cell.
#'   Default: 5 (from \code{paras$Task$impute$KNN$cell$kmin}).
#'
#' @param cell.select Character vector specifying cell identifiers to include
#'   in the analysis. Cell names should match the single-cell BAM file names
#'   excluding the .bam suffix. If \code{NULL}, all available cells are processed.
#'   Default: \code{NULL} (process all cells).
#'
#' @return Character string specifying the path to the cell similarity output
#'   directory (\code{work_path/imputation/cell_similarity/}). The function saves
#'   two main output files:
#'   \itemize{
#'     \item \code{cell.similars.rds}: List containing similarity matrices for
#'       different features
#'     \item \code{dyk.cell.rds}: List containing the number of neighbors
#'       determined for each cell
#'   }
#'
#' @export
#' @import parallel
#' @import Matrix
#' @importFrom reticulate source_python py_module_available use_python
#' @importFrom stats var
#'
#'
getCellSimilarity <- function(
    paras, rds_path = NULL,output_path=NULL,feature_num = paras$Task$impute$feature_num,
    rbp = paras$Task$impute$rbp,cell_similarity_data = paras$Task$impute$cell_similarity_data,
    distance_method = paras$Task$impute$KNN$cell$distance_method,
    alpha_cell = paras$Task$impute$KNN$cell$alpha,
    decay_cell = paras$Task$impute$KNN$cell$decay,
    kcell_max = paras$Task$impute$KNN$cell$kmax,
    kcell_min = paras$Task$impute$KNN$cell$kmin,cell.select = NULL) {
    msg <- paste0("[", Sys.time(), "] ", "Calculate cell similarity...")
    print(msg)
    #script
    srcpath = system.file("python", package = "SCSES")
    py_path <- paras$Basic$python_path
    Sys.setenv(RETICULATE_PYTHON=py_path)
    # Input
    if (is.null(rds_path)) {
        rds_path <- paste0(paras$Basic$work_path, "/rds_processed/")
    }
    print(paste0("Input: ", rds_path))
    # Output
    if (is.null(output_path)) {
        output_path <- paste0(paras$Basic$work_path, "/imputation/cell_similarity/")
    }
    if(dir.exists(output_path)){
        unlink(output_path, recursive = TRUE)
    }
    dir.create(output_path, recursive = T)
    print(paste0("Output: ", output_path))

    # validate parameters
    feature_num <- check.int.or.null(x = feature_num, default = 1000)
    print(paste0("feature_num=", feature_num, "  checked"))
    cell_similarity_data <- unlist(strsplit(cell_similarity_data, ";"))
    cell_similarity_data <- check.valid(x = cell_similarity_data, select = c("EXP_RBP", "RC", "PSI"))
    print(paste0("cell_similarity_data=", paste(cell_similarity_data, collapse = ";"), "  checked"))
    distance_method <- check.valid(x = distance_method, select = c("euclidean", "cosine"))
    print(paste0("distance_method=", paste(distance_method, collapse = ";"), "  checked"))
    alpha_cell <- check.double.or.null(x = alpha_cell, default = 0.8)
    print(paste0("alpha_cell=", paste(alpha_cell, collapse = ";"), "  checked"))
    kcell_max <- check.int.or.null(x = kcell_max, default = 50)
    print(paste0("kcell_max=", paste(kcell_max, collapse = ";"), "  checked"))
    kcell_min <- check.int.or.null(x = kcell_min, default = 5)
    print(paste0("kcell_min=", paste(kcell_min, collapse = ";"), "  checked"))
    decay_cell <- check.double.or.null(x = decay_cell, default = 0.05)
    print(paste0("decay_cell=", paste(decay_cell, collapse = ";"), "  checked"))
    # validate input
    rds_files = list.files(path = rds_path, pattern = "*rds")
    for(type in cell_similarity_data){
        msg = paste("Checking data:", type)
        print(msg)
        if (type == "PSI") {
            if ("psi.rds" %in% rds_files) {
                psi <- readRDS(file = paste0(rds_path, "/psi.rds"))
            } else {
                stop(paste("Preprocessed psi must be saved in", rds_path))
            }
            if (nrow(psi) < 3) {
                stop("It is unreliable to calculate cell similarity based on less than 3 events.")
            }else{
                print(paste(
                    nrow(psi), "events", "are used to calculate cell similarity"
                ))
            }
        } else if (type == "RC") {
            if ("rc.rds" %in% rds_files) {
                rc <- readRDS(file = paste0(rds_path, "/rc.rds"))
            } else {
                stop(paste("Preprocessed read count must be saved in", rds_path))
            }
            if (nrow(rc) < 6) {
                stop("It is unreliable to calculate cell similarity based on less than 6 events associated reads.")
            } else {
                print(paste(
                    nrow(rc), "reads", "are used to calculate cell similarity"
                ))
            }

        } else if (type == "EXP_RBP") {
            if ("expr.rds" %in% rds_files) {
                expr <- readRDS(file = paste0(rds_path, "/expr.rds"))
            } else {
                stop(paste("Preprocessed gene expression must be saved in", rds_path))
            }
            # validate rbp
            if (check.path.or.character(rbp) == "path") {
                rbp <- readLines(rbp)
            } else if (check.path.or.character(rbp) == "character") {
                rbp <- rbp
            } else {
                stop("rbp must be a path or character string")
            }
            rbp <- intersect(rbp, rownames(expr))
            if (length(rbp) < 3) {
                stop("It is unreliable to calculate cell similarity based on less than 3 rbps.")
            }else{
                print(paste(
                    length(rbp), "rbps", "are used to calculate cell similarity"
                ))
            }
        }
    }
    # calculate cell similarity
    if (is(expr,"Matrix")){
        fun_pca = get("PCA_D_reduct_sparse")
    }else {
       fun_pca = get("PCA_D_reduct")
    }
    cell.similars.res <- list()
    cell.similars <- list()
    dyk_cell <- list()
    param_cell_list <- list(similar_method = distance_method, kmax = kcell_max, kmin = kcell_min, alpha = alpha_cell, decay = decay_cell)
    for (type in cell_similarity_data) {
        msg = paste(paste0("[", Sys.time(), "]"), "Computing cell similarity based on", type)
        print(msg)
        if (type == "PSI") {
            data <- psi
        } else if (type == "RC") {
            data <- rc
        } else if (type == "EXP_RBP") {
            data <- expr[rbp, ]
        }
        if (!is.null(cell.select)) {
            cell.select <- intersect(cell.select, colnames(data))
        } else {
            cell.select <- colnames(data)
        }
        if (length(cell.select) == 0) {
            stop("The number of cells is 0.")
        }else {
           print(paste(
                "Calculate similarity among",length(cell.select), "cells."
           ))
        }
        if (feature_num > nrow(data)) {
            print("The number of features is greater than the number of rows in the input data.")
            print(paste(
              "Total",nrow(data), "features will be used"
            ))
            feature_num = nrow(data)
        }
        v <- apply(data, 1, var)
        id_select <- row.names(data)[order(v, decreasing = T)[1:feature_num]]
        D_reduct_res <- fun_pca(data = data, id_select = id_select)
        if (!py_module_available(module = "sklearn")) {
            stop(
                "Cannot find sklearn, please install through pip (e.g. pip install scikit-learn)."
            )
        }
        if (!py_module_available(module = "numpy")) {
            stop(
                "Cannot find numpy, please install through pip (e.g. pip install numpy)."
            )
        }
        if (!py_module_available(module = "pandas")) {
            stop(
                "Cannot find pandas, please install through pip (e.g. pip install pandas)."
            )
        }
        if(nrow(D_reduct_res)>5000){
            if (!py_module_available(module = "scipy")) {
                stop(
                    "Cannot find scipy, please install through pip (e.g. pip install scipy)."
                )
            }
            script_cell_similarity = paste0(srcpath, "/Dynamic_Kcell2.py")
            source_python(script_cell_similarity)
            cell.similars.res[[type]] <- cell_similarity(D_reduct_res, param_cell_list)
        }else if(nrow(D_reduct_res)>3){
            script_cell_similarity = paste0(srcpath, "/Dynamic_Kcell.py")
            source_python(script_cell_similarity)
            cell.similars.res[[type]] <- cell_similarity(D_reduct_res, param_cell_list)
        }else{
            script_cell_similarity = paste0(srcpath, "/Cell_similarity.py")
            source_python(script_cell_similarity)
            cell.similars.res[[type]] <- cell_similarity(D_reduct_res, list(similar_method = distance_method))
        }
        cell.similars[[type]] <- cell.similars.res[[type]][[1]]
        colnames(cell.similars[[type]]) = colnames(data)
        rownames(cell.similars[[type]]) = colnames(data)
        dyk_cell[[type]] <- cell.similars.res[[type]][[2]]
        names(cell.similars.res[[type]][[2]]) = colnames(data)
    }
    saveRDS(cell.similars, paste0(output_path, "/cell.similars.rds"))
    saveRDS(dyk_cell, paste0(output_path, "/dyk.cell.rds"))
    msg <- paste0("[", Sys.time(), "] ", "Calculate cell similarity Finish.")
    print(msg)
    return(output_path)
}


#' @title Construct KNN Graph for Splicing Events
#' @description
#' Computes splicing event features and constructs k-nearest neighbors graph
#' for given dataset by event types. This function integrates sequence features,
#' and RNA-binding protein information to build similarity networks for
#' splicing event imputation.
#'
#' @param paras A list containing configuration parameters, typically loaded via
#'   \code{readSCSESconfig(paras_file)}.
#' @param rds_path Character string specifying the path to processed RDS data files.
#'   This directory should contain the preprocessed single-cell data matrices
#'   (expr.rds, rc.rds, psi.rds) generated by upstream processing functions.
#'   Default: extracted from \code{paras$Basic$work_path/rds_processed/}.
#' @param output_path Character string specifying the output directory for event
#'   similarity results. Default: \code{work_path/imputation/event_similarity/}.
#' @param core Integer specifying the number of CPU threads to use for parallel
#'   processing. Default: \code{paras$Basic$core}.
#' @param pkg Character string specifying the genome package name.
#'   Default: \code{paras$Basic$refgenome$genome_name}.
#' @param ref_path Character string specifying the path to reference FASTA file.
#'   Default: \code{paras$Basic$refgenome$ref_path}.
#' @param phast.path Character string specifying the path to phastCons.bw file
#'   for conservation scores. Default: \code{paras$Task$impute$event_features$phast_path}.
#' @param chr.prefix Character string specifying the prefix added to chromosome IDs.
#'   Should be either "chr" or "". Default: "" (from \code{paras$Task$impute$event_features$chr_prefix}).
#' @param ae.para List containing parameters for encoding sequence features using
#'   autoencoder methods. Default: \code{paras$Task$impute$event_features$AE}.
#' @param rbp Character string or file path specifying RNA-binding protein (RBP)
#'   information. Can be either:
#'   \itemize{
#'     \item File path to RBP annotation file
#'     \item Character string containing RBP gene symbols
#'     \item \code{NULL} to use default RBP set from configuration
#'   }
#'   Default: \code{paras$Task$impute$rbp}.
#' @param kevent Integer specifying the number of nearest neighbors to consider
#'   in KNN graph construction. Default: 10 (from \code{paras$Task$impute$KNN$event$k}).
#' @param alpha_event Numeric value between 0 and 1 specifying the random walk
#' probability (1-restart probability).
#' Default: 0.8 (from \code{paras$Task$impute$KNN$event$alpha}).
#' @param decay_event Numeric threshold of change in the similarity matrix.
#' The algorithm stops when the change falls below this threshold between
#' consecutive iterations. Default: 0.05 (from \code{paras$Task$impute$KNN$event$decay}).
#'
#' @return Character string representing the path to saved event similarity results.
#'   The function saves a list of event similarity matrices for different types
#'   of splicing events to: \code{work_path/imputation/event_similarity/event.similars.rds}
#'
#' @export
#' @importFrom reshape2 melt
#' @import parallel
#' @importFrom reticulate source_python py_module_available use_python
#' @import rtracklayer
#' @importFrom stats rbinom
#' @import BSgenome
#' @import Biostrings
#' @import R.matlab
#' @import raveio
#' @import rhdf5
#' @import hdf5r
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom Biostrings readDNAStringSet writeXStringSet
#' @importFrom BSgenome forgeBSgenomeDataPkg
#'
#'
getEventSimilarity <- function(
    paras,rds_path = NULL,output_path=NULL,core = paras$Basic$core,
    pkg = paras$Basic$refgenome$genome_name,ref_path = paras$Basic$refgenome$ref_path,
    phast.path = paras$Task$impute$event_features$phast_path,
    chr.prefix = paras$Task$impute$event_features$chr_prefix,
    ae.para = paras$Task$impute$event_features$AE,
    rbp = paras$Task$impute$rbp,
    kevent = paras$Task$impute$KNN$event$k,
    alpha_event = paras$Task$impute$KNN$event$alpha,
    decay_event = paras$Task$impute$KNN$event$decay) {
    msg <- paste0("[", Sys.time(), "] ", "Calculate events KNN...")
    print(msg)
    # Input
    if (is.null(rds_path)) {
        rds_path <- paste0(paras$Basic$work_path, "/rds_processed/")
    }
    print(paste0("Input: ", rds_path))
    # Output
    if (is.null(output_path)) {
        output_path <- paste0(paras$Basic$work_path, "/imputation/event_similarity/")
    }
    if (dir.exists(output_path)) {
        unlink(output_path, recursive = TRUE)
    }
    dir.create(output_path, recursive = T)
    print(paste0("Output: ", output_path))
    # script
    matlab_path <- system.file("matlab", package = "SCSES")
    mcr_path = paras$Basic$mcr_path
    mat_combineDistance = paste0(matlab_path, "/combineDistance/run_combineDistance.sh")
    mat_knn_similarity = paste0(matlab_path, "/knn_similarity_from_path/run_knn_similarity_from_path.sh")
    py_srcpath <- system.file("python", package = "SCSES")
    py_path <- paras$Basic$python_path
    Sys.setenv(RETICULATE_PYTHON=py_path)
    source_python(paste0(py_srcpath, "/AE.py"))
    # validate parameters----
    alpha_event <- check.double.or.null(x = alpha_event, default = 0.8)
    print(paste0("alpha_event=", paste(alpha_event, collapse = ";"), "  checked"))
    kevent <- check.int.or.null(x = kevent, default = 10)
    print(paste0("kevent=", paste(kevent, collapse = ";"), "  checked"))
    decay_event <- check.double.or.null(x = decay_event, default = 0.05)
    print(paste0("decay_event=", paste(decay_event, collapse = ";"), "  checked"))

    # validate input----
    rds_files = list.files(path = rds_path, pattern = "*rds")
    msg <- paste("Checking data...")
    print(msg)
    if (!all(c("psi.rds", "expr.rds", "event.rds") %in% rds_files)) {
        stop(paste("Preprocessed psi, gene expression, and events annotation must be saved in", rds_path))
    }
    psi <- readRDS(file = paste0(rds_path, "/psi.rds"))
    expr <- readRDS(file = paste0(rds_path, "/expr.rds"))
    event <- readRDS(file = paste0(rds_path, "/event.rds"))
    event.info <- readRDS(file = paste0(rds_path, "/event.info.list.rds"))
    event_ids <- row.names(psi)
    if (!all(event_ids %in% event$event)) {
        stop("The event information in event.rds dose not match event id in psi matrix")
    }
    if (nrow(psi) < 3) {
        stop("The number of splicing events is less than 3.")
    }
    print("Checking events...")
    event_types = names(event.info)
    event_types <- check.valid(x = event_types, select = c("A3SS", "A5SS", "AL", "SE", "MXE", "RI"))
    print(paste0("event_type=", paste(event_types, collapse = ";"), "  checked"))

    # get events feature----
    msg <- paste0("[", Sys.time(), "] ", "Calculate events feature...")
    print(msg)
    msg = paste0("[", Sys.time(), "] ", "step1 Creating BSgenome for ", pkg, "=======", "")
    print(msg)
    if (!requireNamespace(pkg, quietly = TRUE)) {
      libpath=.libPaths()
      libpath=unlist(strsplit(libpath," "))[1]
      createBSgenome(ref_path = ref_path,
                     out_path = output_path, pkg = pkg,install_lib = libpath
      )
    }
    library(pkg, character.only = T, quietly = T)
    bs.genome = get(pkg)
    msg = paste0("[", Sys.time(), "] ", "step2 Extracting features =======")
    print(msg)
    for (type in event_types)
    {
        events = event.info[[type]]$event
        if(length(events)>1){
            msg = paste(paste0("[", Sys.time(), "]"), "Extracting", type, "features...")
            print(msg)
            # 0. region extraction----
            msg = paste0("[", Sys.time(), "] ", "Loading events...", "")
            print(msg)
            cluster = makeCluster(core)
            fun = get(paste0(type, ".info.seperation"))
            rr = clusterEvalQ(cl = cluster, {
              library(rtracklayer)
              NULL
            })
            clusterExport(cl = cluster, varlist = c("GRanges"), envir = environment())
            events.info = parLapply(cluster, X = as.list(events), fun = fun)
            events.info = do.call(what = rbind, args = events.info)
            stopCluster(cluster)

            msg = paste0("[", Sys.time(), "] ", "Parsing events region...", "")
            print(msg)
            splice.region.fun = get(paste0(type, ".splice.seq.extraction"))
            splice.region = splice.region.fun(events.info, bs.genome = bs.genome,
                                              core = core)
            exon.region.fun = get(paste0(type, ".exon.seq.extraction"))
            exon.region = exon.region.fun(events.info, bs.genome = bs.genome, core = core)
            # 1. length feature----
            msg = paste0("[", Sys.time(), "] ", "Extracting length features.", "")
            print(msg)
            length.feature = lapply(X = as.list(seq(1, nrow(events.info))), FUN = function(x) {
                lf.fun = get(paste0(type, ".length.feature"))
                result = lf.fun(events.info[x, ])
                return(result)
            })
            length.feature = do.call(what = rbind, args = length.feature)
            for (i in seq(1, ncol(length.feature) - 1))
            {
                len = length.feature[, i]
                len.quanl = quantile(x = len, probs = seq(0, 1, 0.05))
                len[len < len.quanl["5%"]] = len.quanl["5%"]
                len[len > len.quanl["95%"]] = len.quanl["95%"]
                length.feature[, i] = len
            }
            # 2. motif feature----
            msg = paste0("[", Sys.time(), "] ", "Extracting motif features.", "")
            print(msg)
            splice.region = lapply(X = splice.region, FUN = function(region) {
                seqs = strsplit(x = region$seq, split = "")
                seqs = do.call(what = rbind, args = seqs)
                colnames(seqs) = paste0("p", seq(1, ncol(seqs)))

                seq2 = melt(seqs)
                count = as.data.frame.matrix(x = table(seq2[, c("value", "Var2")]))

                PWM = t(t(count) / colSums(count))
                colnames(PWM) = paste0("p", seq(1, ncol(PWM)))

                score = lapply(X = as.list(x = seq(1, ncol(PWM))), FUN = function(index) {
                    return(PWM[seqs[, index], index])
                })
                score = do.call(what = cbind, args = score)
                region$motif.score = rowSums(score)
                region$motif.score = scale(region$motif.score)[, 1]
                return(region)
            })
            # 3. conservation feature----
            msg = paste0("[", Sys.time(), "] ", "Extracting conservation features.")
            print(msg)
            bwf = BigWigFile(phast.path)
            print("Checking chromosome prefix...")
            chr.prefix <- check.valid(x = chr.prefix, select = c("", "chr"))
            splice.region = phastScore.extraction(splice.region, bwf = bwf, chr.prefix = chr.prefix, core = core)
            exon.region = phastScore.extraction(exon.region, bwf = bwf, chr.prefix = chr.prefix, core = core)

            # 4. Kmer feature----
            msg = paste0("[", Sys.time(), "] ", "Extracting kmer features.", "")
            print(msg)
            kmer.fun = get(paste0(type, ".kmer.extraction"))
            kmer = kmer.fun(splice.region, exon.region)
            # 5. A ratio feature----
            msg = paste0("[", Sys.time(), "] ", "Extracting A Ratio features.", "")
            print(msg)
            save.obj = c("length.feature", "splice.region", "exon.region", "kmer")
            if (type == "RI") {
                A.aRatio = aPercentage(region = exon.region$junction, core = core)
                save.obj = c(save.obj, "A.aRatio")
            }
            if (type == "A3SS" | type == "A5SS" | type == "AL") {
                I.region = GRanges(sub(pattern = "junction:", replacement = "", x = events.info$junction1))
                I.region = getSeq(x = bs.genome, I.region)
                I.region = data.frame(event = events.info$event, seq = as.character(I.region))
                I.aRatio = aPercentage(region = I.region, seperate = 100, core = core)

                A.aRatio = aPercentage(region = exon.region$A, seperate = 10, core = core)

                save.obj = c(save.obj, "A.aRatio", "I.aRatio")
            }

            # 6. Saving Result----
            msg = paste0("[", Sys.time(), "] ", "Saving Result", "")
            print(msg)
            save(list = save.obj, file = paste0(output_path, "/", type, "_feature.RData"))

            msg = paste(paste0("[", Sys.time(), "]"), "Extracting", type, "features Finished")
            print(msg)
        }
    }
    # combine events feature----
    msg = paste0("[", Sys.time(), "] ", "step3 Combining events feature =======")
    print(msg)
    for (type in event_types)
    {
        inpath_type = paste0(output_path, type, "_feature.RData")
        if(file.exists(inpath_type)){
            outpath_type = paste0(output_path, type, "_feature.txt")
            msg = paste(paste0("[", Sys.time(), "]"), "Parsing", type, "features...")
            print(msg)
            load(file = inpath_type)
            kmer.f = do.call(what = cbind, args = kmer)
            exon.f = lapply(names(exon.region), function(x) {
                out = data.frame(value = exon.region[[x]]$phastscore, row.names = exon.region[[x]]$event)
                colnames(out) = x
                return(out)
            })
            exon.f = do.call(what = cbind, args = exon.f)
            splice.f = lapply(names(splice.region), function(x) {
                out = data.frame(
                    value1 = splice.region[[x]]$phastscore,
                    value2 = splice.region[[x]]$motif.score,
                    row.names = splice.region[[x]]$event
                )
                colnames(out) = paste0(x, c("phastscore", "motif.score"))
                return(out)
            })
            splice.f = do.call(what = cbind, args = splice.f)
            length.f = length.feature[, 1:(ncol(length.feature) - 1)]
            row.names(length.f) = length.feature$event
            length.f = length.f[row.names(kmer.f), ]
            feature_m = cbind(length.f, kmer.f, exon.f, splice.f)
            feature_m[is.na(feature_m)] <- 0

            if (type %in% c("A3SS", "A5SS", "AL")) {
                feature_m = cbind(feature_m, A.aRatio[rownames(feature_m), ], I.aRatio[rownames(feature_m), ])
            }
            if (type == "RI") {
                feature_m = cbind(feature_m, A.aRatio[rownames(feature_m), ])
            }
            write.table(feature_m, outpath_type, quote = F, sep = "\t")
            msg = paste(paste0("[", Sys.time(), "]"), "Parsing", type, "features Finished")
            print(msg)
        }
    }
    # AE embedding----
    if (!py_module_available(module = "keras")) {
        stop(
            "Cannot find keras, please install through conda (e.g. conda install tensorflow keras)."
        )
    }
    if (!py_module_available(module = "pandas")) {
        stop(
            "Cannot find pandas, please install through pip (e.g. pip install pandas)."
        )
    }
    if (!py_module_available(module = "numpy")) {
        stop(
            "Cannot find numpy, please install through pip (e.g. pip install numpy)."
        )
    }
    if (!py_module_available(module = "scipy")) {
        stop(
            "Cannot find numpy, please install through pip (e.g. pip install scipy)."
        )
    }
    msg = paste0("[", Sys.time(), "] ", "step4 Encoding events feature =======")
    print(msg)
    event.features = list()
    for (type in event_types)
    {
        inpath = paste0(output_path, "/", type, "_feature.txt")
        log_file <- paste0(output_path, "/py_AE_", type, ".log")
        if(file.exists(inpath)){
            print(paste(paste0("[", Sys.time(), "]"), type, "event encoding..."))
            feature = read.table(inpath)
            feature = as.matrix(feature)
            layers = gsub(pattern = "\\[|\\]", replacement = "", ae.para[[type]]$layer)
            layers = as.numeric(unlist(strsplit(x = layers, split = ",")))
            flag = T
            while (flag) {
                auto.feature = model_training_parameter(
                    feature,
                    ae.para[[type]]$embedding,
                    layers,
                    ae.para[[type]]$epoch,
                    log_file
                )
                vars <- apply(auto.feature, 2, var)
                valid <- which(vars != 0)
                auto.feature <- auto.feature[, valid]
                flag <- length(valid) <= 1
            }
            rownames(auto.feature) = rownames(feature)
            event.features[[type]] = as.matrix(auto.feature)
            print(paste(paste0("[", Sys.time(), "]"), type, "event encoding Finish."))
        }
    }
    msg = paste0("[", Sys.time(), "] ", "step5 Calculate splicing regulation distance and Combine distance =======")
    print(msg)
    # Calculate splicing regulation distance and Combine distance----
    # validate rbp
    if (check.path.or.character(rbp) == "path") {
        rbp <- readLines(rbp)
    } else if (check.path.or.character(rbp) == "character") {
        rbp <- rbp
    } else{
        stop("rbp must be a path or character string")
    }
    rbp <- intersect(rbp, rownames(expr))
    rbp.var = apply(X = expr[rbp, ], MARGIN = 1, FUN = var)
    rbp = rbp[which(rbp.var != 0)]
    if (length(rbp) < 3) {
        print("It is unreliable to calculate event similarity based on less than 3 rbps.")
    } else {
        print(paste(
            length(rbp), "rbps are used to calculate splicing regulation information"
        ))
    }
    event.psi = list()
    for (type in names(event.features))
    {
        event.psi[[type]] = as.matrix(psi[rownames(event.features[[type]]), colnames(expr)])
    }
    token = paste(Sys.getpid(), rbinom(1, size = 1000000000, prob = 0.5), sep = "-")
    datapath = paste0(output_path, "/combine_feature_data", token, ".h5")
    distance_path = paste0(output_path, "/feature_combine.h5")
    print("Save data")
    saveHdf5File(datapath, list(
        psi = event.psi, rbp = as.matrix(expr[rbp, ]),
        feature = event.features, method = "seuclidean"
    ))
    print("Save data Finished")
    log_file <- paste0(output_path, "/mat_combineEventFeatureDistance.log")
    if (file.exists(log_file)) {
      file.remove(log_file)
    }
    cmd = paste("bash", mat_combineDistance, mcr_path, datapath,
                distance_path, ">>", log_file, "2>&1")
    system(cmd, wait = T)
    file.remove(datapath)

    # Get event similarity----
    msg = paste0("[", Sys.time(), "] ", "step6 Calculate combined event similarity =======")
    print(msg)
    h5createGroup(file = distance_path, group = "event_names")
    for (type in names(event.features))
    {
        h5write(obj = event.info[[type]]$event, file = distance_path, name = paste0("/event_names/", type), )
    }

    log_file <- paste0(output_path, "/mat_calculateEventSimilarity.log")
    if (file.exists(log_file)) {
      file.remove(log_file)
    }
    event.similars = list()
    for (type in event_types)
    {
        gc()
        msg <- paste(paste0("[", Sys.time(), "]"), "Calculate ", type, "event Similarity")
        print(msg)
        events = event.info[[type]]
        if(nrow(events)==1){
            event_similar = matrix(1)
            dimnames(event_similar) = list(events$event, events$event)
        }else{
          if (nrow(events) < kevent) {
            stop(sprintf(
              "The number of events is %d, which is smaller than the specified K-event (%d).
              Please set 'kevent' to a value less than %d and greater than 1.",
              nrow(events), kevent, nrow(events)
            ))
          }
          token <- paste(Sys.getpid(), rbinom(1, size = 1000000000, prob = 0.5), sep = "-")
          inpath <- paste0(dirname(output_path), "/knn_data_", token, ".h5")
          outpath <- paste0(dirname(output_path), "/knn_result_", token, ".mat")
          paras <- list()
          paras[["k"]] <- ifelse(kevent>=nrow(events),nrow(events)-1,kevent)
          paras[["alpha"]] <- alpha_event
          paras[["decay"]] <- decay_event

          h5createFile(inpath)
          h5write(obj = distance_path, inpath, "dist_path")
          h5write(obj = type, inpath, "type")
          h5write(paras, inpath, "paras")

          cmd <- paste(
            "bash", mat_knn_similarity,
            mcr_path, inpath, outpath, ">>", log_file, "2>&1"
          )
          system(command = cmd, wait = T)

          similar <- read_mat(outpath)$similar
          event_ids <- h5read(distance_path, paste0("/event_names/", type))
          rownames(similar) <- event_ids
          colnames(similar) <- event_ids
          file.remove(inpath, outpath)
        }
        msg <- paste(paste0("[", Sys.time(), "]"), "Calculate", type, "event Similarity Finished")
        print(msg)
        event.similars[[type]] = similar
    }
    msg = paste0("[", Sys.time(), "] ", "Calculate event similarity Finished.")
    print(msg)
    saveRDS(event.similars, paste0(output_path, "/event.similars.rds"))
    gc()
    return(output_path)
}



