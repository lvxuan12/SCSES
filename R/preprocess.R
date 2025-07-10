#' @include utilities.R

#' @title MXE Filter
#' @description MXE need strict filter of upstream and downstream junction reads
#' @param event data frame of event information
#' @param rc data frame of reads associated with splicing events
#' @param min.Cell the minimum number of cells that are required to have at least
#'     one reads of both upstream and downstream junction supporting exon

#' @return data frame of filtered event information
#'
#' @export

filterMXE <- function(event, rc, min.Cell) {
    event_mxe <- event[which(event$type == "MXE"), ]
    tmp.array.ex <- array(0, dim = c(nrow(event_mxe), ncol(rc), 2))
    for (i in c(2, 3))
    {
        tmp <- rc[event_mxe[, i], ]
        tmp <- as.matrix(tmp)
        tmp.array.ex[, , i - 1] <- tmp > 1
    }
    tmp.mean.ex <- rowSums(tmp.array.ex, dims = 2, na.rm = T)
    row.names(tmp.mean.ex) <- event_mxe$event
    colnames(tmp.mean.ex) <- colnames(rc)
    tmp.array.re <- array(0, dim = c(nrow(event_mxe), ncol(rc), 2))
    for (i in c(4, 5))
    {
        tmp <- rc[event_mxe[, i], ]
        tmp <- as.matrix(tmp)
        tmp.array.re[, , i - 3] <- tmp > 1
    }
    tmp.mean.re <- rowSums(tmp.array.re, dims = 2, na.rm = T)
    row.names(tmp.mean.re) <- event_mxe$event
    colnames(tmp.mean.re) <- colnames(rc)
    idx <- which(rowSums(tmp.mean.ex == 2) >= min.Cell | rowSums(tmp.mean.re == 2) >= min.Cell)
    event_new <- rbind(event[which(event$type != "MXE"), ], event_mxe[idx, ])
    return(event_new)
}

#' @title Filter and Process Gene Expression and Splicing Data
#' @description
#' Performs quality control and filtering of raw single-cell RNA-seq
#' and splicing events. This function processes raw data from
#' \code{work_path/rds/} and saves filtered, normalized data to
#' \code{work_path/rds_processed/} for downstream analysis.
#'
#' The preprocessing pipeline includes:
#' \itemize{
#'   \item Cell quality control based on library size, gene detection, and mitochondrial content
#'   \item Gene filtering based on expression across cells
#'   \item Splicing event filtering based on PSI value and read count support
#'   \item Read count normalization and log transformation
#' }
#'
#' @param paras A list of parameters obtained from \code{readSCSESconfig(paras_file)}.
#' @param min.percentCells.gene Numeric value (0-1) specifying the minimum
#'   percentage of cells required to have at least 1 read for a gene to be retained.
#'   Default: 0.9 (\code{paras$Basic$filter_sc$min.percentCells.gene}).
#' @param min.percentCells.event Numeric value (0-1) specifying the maximum
#'   percentage of cells allowed to have extreme PSI values (0, 1, or NA)
#'   for an event to be retained.
#'   Default: 0.9 (\code{paras$Basic$filter_sc$min.percentCells.event}).
#' @param min.Cell Integer specifying the minimum number of cells required to have
#'   at least \code{min.RC} inclusion/exclusion reads for an event to be retained.
#'   Default: 10 (\code{paras$Basic$filter_sc$minCell}).
#' @param min.RC Integer specifying the minimum number of inclusion/exclusion reads
#'   required in at least \code{min.Cell} cells for an event to be retained.
#'   Default: 5 (\code{paras$Basic$filter_sc$minRC}).
#' @param min.nCount Integer specifying the minimum total read count required per cell.
#'   Cells below this threshold are filtered out.
#'   Default: 1 (\code{paras$Basic$filter_sc$min.nCount}).
#' @param min.nFeatures Integer specifying the minimum number of detected genes
#'   required per cell. Cells below this threshold are filtered out.
#'   Default: 1 (\code{paras$Basic$filter_sc$min.nFeatures}).
#' @param max.percentMT Numeric value (0-1) specifying the maximum percentage of
#'   mitochondrial gene expression allowed per cell. Cells above this threshold
#'   are filtered out. Default: 1 (\code{paras$Basic$filter_sc$max.percentMT}).
#' @param cell.select Character vector of cell IDs to include in analysis.
#'   Should match single-cell BAM file names excluding the \code{.bam} suffix.
#'   If \code{NULL}, all common cells across datasets are used. Default: \code{NULL}.
#'
#' @return
#' Character string containing the path to the processed data directory
#' (\code{work_path/rds_processed/}). The following files are saved:
#' \describe{
#'   \item{\code{expr.rds}}{Filtered and quality-controlled gene expression matrix}
#'   \item{\code{psi.rds}}{Filtered PSI values matrix}
#'   \item{\code{rc.rds}}{Normalized and log-transformed read count matrix for junctions}
#'   \item{\code{event.rds}}{Filtered splicing event annotation data frame}
#'   \item{\code{event.info.list.rds}}{Structured event information organized by
#'         splicing type (A3SS, A5SS, SE, MXE, RI)}
#' }
#'
#' @export

preprocessEvent <- function(
    paras,
    min.percentCells.gene=paras$Basic$filter_sc$min.percentCells.gene,
    min.percentCells.event=paras$Basic$filter_sc$min.percentCells.event,
    min.Cell=paras$Basic$filter_sc$minCell,
    min.RC=paras$Basic$filter_sc$minRC,
    min.nCount=paras$Basic$filter_sc$min.nCount,
    min.nFeatures=paras$Basic$filter_sc$min.nFeatures,
    max.percentMT=paras$Basic$filter_sc$max.percentMT, cell.select = NULL) {
    msg <- paste0("[", Sys.time(), "] ", "Processing raw data...")
    print(msg)
    rds_path <- paste0(paras$Basic$work_path, "/rds/")
    print(paste0("Input: ", rds_path))
    output_path <- paste0(paras$Basic$work_path, "/rds_processed/")
    print(paste0("Output: ", output_path))
    if (!dir.exists(output_path)) {
        dir.create(output_path)
    }
    # psi input
    psi <- readRDS(file = paste0(rds_path, "/psi.rds"))
    # reads count input
    rc <- readRDS(file = paste0(rds_path, "/rc.rds"))
    event <- readRDS(file = paste0(rds_path, "/event.rds"))
    # expr input
    expr <- readRDS(file = paste0(rds_path, "/count_norm.rds"))
    print("Before filtering")
    print(paste0("expr: ", nrow(expr), "*", ncol(expr)))
    print(paste0("psi: ", nrow(psi), "*", ncol(psi)))
    print(paste0("rc: ", nrow(rc), "*", ncol(rc)))
    if (!((all(colnames(psi) %in% colnames(expr)) & ncol(psi) == ncol(expr)) | (all(colnames(rc) %in% colnames(expr)) & ncol(rc) == ncol(expr)))) {
        print("Cell ids of raw data is mismatch! The insection of cells will be used to further analysis.")
    }
    cell_id <- Reduce(intersect, list(colnames(expr), colnames(psi), colnames(rc)))
    if (!is.null(cell.select)) {
        cell.select <- intersect(cell.select, cell_id)
    } else {
        cell.select <- cell_id
    }
    psi <- psi[, cell.select, drop = F]
    rc <- rc[, cell.select, drop = F]
    expr <- expr[, cell.select, drop = F]
    # filter invalid cells and genes
    anno_qc <- data.frame(
        Library_size = colSums(expr),
        Detected_gene_num = colSums(expr != 0),
        MT_gene_Fraction = colSums(expr[grepl("^MT-|^mt-", row.names(expr)), ]) / colSums(expr),
        RP_gene_Fraction = colSums(expr[grepl("^RP[SL]|^Rp[sl]", row.names(expr)), ]) / colSums(expr),
        row.names = colnames(expr)
    )
    keep_cols <- row.names(anno_qc)[which(anno_qc$Library_size >= min.nCount &
        anno_qc$Detected_gene_num >= min.nFeatures &
        anno_qc$MT_gene_Fraction <= max.percentMT)]
    expr <- expr[, keep_cols, drop = F]
    expr <- expr[which(apply(expr, 1, function(x) {
        return((length(which(is.na(x))) + length(which(x == 0))) <= ncol(expr) * min.percentCells.gene)
    })), , drop = F]
    saveRDS(as.matrix(expr), paste0(output_path, "/expr.rds"))
    psi <- psi[, keep_cols, drop = F]
    rc <- rc[, keep_cols, drop = F]
    # filter invalid chr
    psi <- psi[which(sapply(row.names(psi), function(x) {
        chr <- strsplit(strsplit(x, "@")[[1]][[1]], ":")[[1]][[2]]
        grepl("^[1-9,X,Y]|chr[1-9,X,Y]", chr)
    })), , drop = F]
    events_geneid <- unlist(lapply(as.list(row.names(psi)), function(x) {
        unlist(strsplit(x, "[|]"))[length(unlist(strsplit(x, "[|]"))) - 1]
    }))
    psi_id <- match(events_geneid, row.names(expr))
    psi <- psi[which(!is.na(psi_id)), , drop = F]
    # filter invalid events
    psi <- psi[which(rowSums(psi == 0) <= ncol(psi) * min.percentCells.event), , drop = F]
    psi <- psi[which(rowSums(psi == 1) <= ncol(psi) * min.percentCells.event), , drop = F]
    psi <- psi[which(rowSums(is.na(psi)) <= ncol(psi) * min.percentCells.event), , drop = F]
    event <- event[match(row.names(psi), event$event), ]
    rc_ex <- calcu_ex_rc(rc, event)
    rc_in <- calcu_in_rc(rc, event)
    psi <- psi[which(rowSums(rc_ex > min.RC) >= min.Cell | rowSums(rc_in > min.RC) >= min.Cell), ]
    # filter MXE events
    event_new <- filterMXE(event = event, rc = rc, min.Cell = min.Cell)
    psi <- psi[intersect(row.names(psi), event_new$event), ]
    event <- event[match(row.names(psi), event$event), ]
    # match rc with psi
    rc_name <- unique(c(event[, 2], event[, 3], event[, 4], event[, 5]))
    rc_name <- rc_name[which(!is.na(rc_name))]
    rc <- rc[match(rc_name, row.names(rc)), , drop = F]
    library_size <- colSums(rc, na.rm = T)
    median_transcript_count <- median(library_size)
    rc_norm <- sweep(rc, MARGIN = 2, median_transcript_count / library_size, FUN = "*")
    rc_norm <- log2(rc_norm + 1)
    saveRDS(event, paste0(output_path, "/event.rds"))
    saveRDS(as.matrix(psi), paste0(output_path, "/psi.rds"))
    saveRDS(as.matrix(rc_norm), paste0(output_path, "/rc.rds"))
    print("After filtering")
    print(paste0("expr: ", nrow(expr), "*", ncol(expr)))
    print(paste0("psi: ", nrow(psi), "*", ncol(psi)))
    print(paste0("rc: ", nrow(rc), "*", ncol(rc)))
    event_types = unique(event$type)
    event.info = lapply(X = as.list(event_types), FUN = function(type) {
        cur.event = event[event$type == type, ]
        if (type == "RI" | type == "SE") {
            info = cur.event[, c("event", "exclusion1", "retention1", "retention2")]
        } else if (type == "MXE") {
            info = cur.event[, c("event", "exclusion1", "exclusion2", "retention1", "retention2")]
        } else if (type == "A3SS" | type == "A5SS" | type == "AL") {
            info = cur.event[, c("event", "exclusion1", "retention1")]
        }
        return(info)
    })
    names(event.info) = event_types
    saveRDS(event.info, paste0(output_path, "/event.info.list.rds"))
    msg <- paste0("[", Sys.time(), "] ", "Successfully processed data.")
    print(msg)
    return(output_path)
}
