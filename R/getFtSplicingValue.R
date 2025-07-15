#' @title Quantify Raw Reads of Splicing Events for Classifier Fine-Tuning
#' @description This function quantifies raw read counts for predefined splicing
#'   events specifically selected for classifier fine-tuning.
#'
#' @param paras A list object containing SCSES configuration parameters, typically
#'   loaded using \code{readSCSESconfig(paras_file)}.
#' @param bam_path Character string specifying the directory containing single-cell
#'   BAM files. Each BAM file should represent one cell and be properly indexed.
#'   Default: \code{paras$Basic$bam_path}.
#' @param core Integer specifying the number of CPU cores to use for parallel
#'   processing. \code{paras$Basic$core}.
#' @param sequence Character string specifying the sequencing protocol type.
#'   Must be either "full-length" (e.g., Smart-seq2) or "UMI" (e.g., 10x Genomics).
#'   This determines which Java quantification tool to use.
#'   Default: \code{paras$Basic$sequence}.
#' @param genome_name Character string specifying the reference genome version.
#'   Must be either "hg19" or "hg38". This determines which predefined event
#'   set to use for quantification. Default: \code{paras$Basic$refgenome$genome_name}.
#'
#' @return Character string specifying the path to the output directory
#'   (\code{work_path/splicing_value_ft/}) containing read count
#'   RDS files. This path can be used as input for \code{\link{FtClassifier}}
#'
#' @export
#' @import parallel

getFtRawRC <- function(
    paras, bam_path = paras$Basic$bam_path,
    core = paras$Basic$core, sequence = paras$Basic$sequence,
    genome_name = paras$Basic$refgenome$genome_name) {
    options("scipen" = 100)
    if (!genome_name %in% c("hg19", "hg38")) {
      stop("Error: genome_name must be either 'hg19' or 'hg38'")
    }
    # script
    jar_path <- system.file("java", package = "SCSES")
    java_path <- paras$Basic$java_path
    # input
    # Built-in data
    print("Loading splicing events for classifer fine tune...")
    extdata_path <- system.file("extdata", package = "SCSES")
    event_path <- paste0(extdata_path, "/ftevents/", genome_name,"/")
    event_types <- gsub(".txt", "", list.files(event_path, "*.txt"))
    print("Checking cells...")
    cells <- list.files(paras$Basic$bam_path, "*bam$", recursive = F)
    if (length(cells) == 0) {
        msg <- paste("There is no bam file in", bam_path, "!")
        stop(msg)
    }
    print(paste(length(cells), "cells are considered."))
    # output
    rc_path <- paste0(paras$Basic$work_path, "/splicing_value_ft/")
    print(paste0("Output: ", rc_path))
    if (!dir.exists(rc_path)) {
        dir.create(rc_path)
    }
    for (type in event_types) {
        event_file <- paste0(event_path, "/", type, ".txt")
        if (!file.exists(event_file)) {
            msg <- paste("There is no", type, "events in", event_path, "!")
            stop(msg)
        }
        msg <- paste0("[", Sys.time(), "] ", "Counting reads of ", type, " events...")
        print(msg)
        outpath_per_cell <- paste0(rc_path, "/", type, "_rjm")
        dir.create(outpath_per_cell)
        log_file <- paste0(event_path, "/java_getFtRC_", type, ".log")
        if (file.exists(log_file)) {
          file.remove(log_file)
        }
        if (type == "RI") {
            cluster <- makeCluster(spec = core)
            clusterExport(cl = cluster, varlist = c("jar_path", "event_file", "outpath_per_cell", "bam_path", "sequence"), envir = environment())
            l <- parLapply(cl = cluster, X = as.list(cells), fun = function(cell) {
                if (sequence=="UMI") {
                    jar.file <- paste0(jar_path, "/ujm3.RI.jar")
                } else {
                    jar.file <- paste0(jar_path, "/rjm3.RI.jar")
                }
                cmd <- paste(
                    java_path, "-Xmx5120m -cp", paste0(jar_path, "/lib"), "-jar", jar.file,
                    bam_path,
                    cell,
                    event_file,
                    outpath_per_cell,
                    "0,16,99,147,83,163", ">>", log_file, "2>&1"
                )
                print(cmd)
                system(cmd, wait = T, show.output.on.console = T)
            })
            stopCluster(cl = cluster)
            files <- list.files(outpath_per_cell, pattern = "bam")
            tmp <- read.table(file = paste0(outpath_per_cell, "/", files[1]), header = F, sep = "\t")
            events <- tmp$V1
            data <- matrix(data = NA, nrow = length(events), ncol = length(files))
            rownames(data) <- events
            colnames(data) <- files
            progress_rate <- 0
            for (i in seq(1, length(files)))
            {
                f <- files[i]
                if (i > (progress_rate * length(files))) {
                    print(paste("Reading RJC File Progress:", paste0(progress_rate * 100, "%")))
                    progress_rate <- round(i / length(files), digits = 1)
                }
                tmp <- read.table(file = paste0(outpath_per_cell, "/", f), header = F, sep = "\t")
                rownames(tmp) <- tmp$V1
                data[events, f] <- tmp[events, 2]
            }
            print(paste("Reading RJC File Progress:100%"))
            colnames(data) <- sub("\\.rjm\\.txt$", "", colnames(data))
            events <- readLines(event_file)
            events <- lapply(as.list(events), function(x) {
                info <- unlist(strsplit(x, split = "isoform[0-9]=|@|\\|"))
                return(data.frame(
                    event = x, exclusion1 = info[grepl(pattern = "junction", x = info)],
                    exclusion2 = NA,
                    retention1 = info[grepl(pattern = "retention", x = info)][1],
                    retention2 = info[grepl(pattern = "retention", x = info)][2]
                ))
            })
            events <- do.call(rbind, events)
            rc_exclusion <- data[match(events$exclusion1, row.names(data)), , drop = F]
            rc_retention1 <- data[match(events$retention1, row.names(data)), , drop = F]
            rc_retention2 <- data[match(events$retention2, row.names(data)), , drop = F]
            saveRDS(
                list(
                    rc_exclusion = rc_exclusion, rc_retention1 = rc_retention1,
                    rc_retention2 = rc_retention2, events = events
                ),
                paste0(rc_path, "/", type, "_rc.rds")
            )
        } else {
            cluster <- makeCluster(spec = core)
            clusterExport(cl = cluster, varlist = c("jar_path", "event_file", "outpath_per_cell", "bam_path", "sequence"), envir = environment())
            l <- parLapply(cl = cluster, X = as.list(cells), fun = function(cell) {
                if (sequence=="UMI") {
                    jar.file <- paste0(jar_path, "/ujm3.jar")
                } else {
                    jar.file <- paste0(jar_path, "/rjm3.jar")
                }
                cmd <- paste(
                    java_path,"-Xmx5120m -cp", paste0(jar_path, "/lib"), "-jar", jar.file,
                    bam_path,
                    cell,
                    event_file,
                    outpath_per_cell,
                    "0,16,99,147,83,163", ">>", log_file, "2>&1"
                )
                print(cmd)
                system(cmd, wait = T)
            })
            stopCluster(cl = cluster)
            files <- list.files(outpath_per_cell, pattern = "bam")
            tmp <- read.table(file = paste0(outpath_per_cell, "/", files[1]), header = F, sep = "\t")
            events <- tmp$V1
            data <- matrix(data = NA, nrow = length(events), ncol = length(files))
            rownames(data) <- events
            colnames(data) <- files
            progress_rate <- 0
            for (i in seq(1, length(files)))
            {
                f <- files[i]
                if (i > (progress_rate * length(files))) {
                    print(paste("Reading RJC File Progress:", paste0(progress_rate * 100, "%")))
                    progress_rate <- progress_rate + 0.1
                }
                tmp <- read.table(file = paste0(outpath_per_cell, "/", f), header = F, sep = "\t")
                rownames(tmp) <- tmp$V1
                data[events, f] <- tmp[events, 2]
            }
            colnames(data) <- sub("\\.rjm\\.txt$", "", colnames(data))
            events <- readLines(event_file)
            events <- lapply(as.list(events), function(x) {
                info <- unlist(strsplit(x, split = "\\|"))
                info <- gsub(pattern = "isoform[0-9]=", replacement = "", x = info)
                iso1 <- unlist(strsplit(x = info[1], split = "@"))
                iso2 <- unlist(strsplit(x = info[2], split = "@"))
                iso1 <- iso1[grepl(pattern = "junction", x = iso1)]
                iso2 <- iso2[grepl(pattern = "junction", x = iso2)]
                if (length(iso1) > 1) {
                    iso1_1 <- iso1[1]
                    iso1_2 <- iso1[2]
                } else {
                    iso1_1 <- iso1[1]
                    iso1_2 <- NA
                }
                if (length(iso2) > 1) {
                    iso2_1 <- iso2[1]
                    iso2_2 <- iso2[2]
                } else {
                    iso2_1 <- iso2[1]
                    iso2_2 <- NA
                }
                return(data.frame(
                    event = x, exclusion1 = iso1_1, exclusion2 = iso1_2,
                    retention1 = iso2_1, retention2 = iso2_2
                ))
            })
            events <- do.call(rbind, events)
            if (type %in% c("A3SS", "A5SS", "AL")) {
                rc_exclusion <- data[match(events$exclusion1, row.names(data)), , drop = F]
                rc_retention <- data[match(events$retention1, row.names(data)), , drop = F]
                saveRDS(
                    list(rc_exclusion = rc_exclusion, rc_retention = rc_retention, events = events),
                    paste0(rc_path, "/", type, "_rc.rds")
                )
            } else if (type == "SE") {
                rc_exclusion <- data[match(events$exclusion1, row.names(data)), , drop = F]
                rc_retention1 <- data[match(events$retention1, row.names(data)), , drop = F]
                rc_retention2 <- data[match(events$retention2, row.names(data)), , drop = F]
                saveRDS(
                    list(
                        rc_exclusion = rc_exclusion, rc_retention1 = rc_retention1,
                        rc_retention2 = rc_retention2, events = events
                    ),
                    paste0(rc_path, "/", type, "_rc.rds")
                )
            } else if (type == "MXE") {
                rc_exclusion1 <- data[match(events$exclusion1, row.names(data)), , drop = F]
                rc_exclusion2 <- data[match(events$exclusion2, row.names(data)), , drop = F]
                rc_retention1 <- data[match(events$retention1, row.names(data)), , drop = F]
                rc_retention2 <- data[match(events$retention2, row.names(data)), , drop = F]
                saveRDS(
                    list(
                        rc_exclusion1 = rc_exclusion1, rc_exclusion2 = rc_exclusion2,
                        rc_retention1 = rc_retention1, rc_retention2 = rc_retention2,
                        events = events
                    ),
                    paste0(rc_path, "/", type, "_rc.rds")
                )
            }
        }
        msg <- paste0("[", Sys.time(), "] ", "Counting reads of ", type, " events Finish.")
        print(msg)
    }
    return(rc_path)
}

#' @title Calculate PSI Values of Splicing Events for Classifier Fine-Tuning
#' @description This function computes PSI values from raw
#'   read counts for predefined splicing events used in classifier fine-tuning.
#'   The function processes read count data generated
#'   by \code{\link{getFtRawRC}} and applies event-specific PSI calculation formulas
#'   to produce training-ready splicing quantification data.
#'
#' @param paras A list object containing SCSES configuration parameters, typically
#'   loaded using \code{readSCSESconfig(paras_file)}.
#' @return Character string specifying the path to the output directory
#'   (\code{work_path/splicing_value_ft/}) containing PSI RDS files.
#'   This path can be used as input for \code{\link{FtClassifier}}
#'
#' @export
#'
getFtRawPSI <- function(paras) {
    # input
    print("Checking raw reads...")
    rc_path <- paste0(paras$Basic$work_path, "/splicing_value_ft/")
    if (length(list.files(rc_path, "*_rc.rds")) == 0) {
        stop(paste("No raw reads data in", rc_path))
    }
    print("Loading splicing events for classifer fine tune...")
    # Built-in data
    extdata_path <- system.file("extdata", package = "SCSES")
    genome_name = paras$Basic$refgenome$genome_name
    event_path <- paste0(extdata_path, "/ftevents/", genome_name, "/")
    event_types <- gsub(".txt", "", list.files(event_path, "*.txt"))
    for (type in event_types) {
        msg <- paste0("[", Sys.time(), "] ", "Calculating PSI value of ", type, " events...")
        print(msg)
        data <- readRDS(file = paste0(rc_path, "/", type, "_rc.rds"))
        events <- data$events
        if (type == "RI" | type == "SE") {
            rc_exclusion <- data$rc_exclusion
            rc_retention1 <- data$rc_retention1
            rc_retention2 <- data$rc_retention2
            iso1 <- rc_exclusion
            rc_retention <- array(0, dim = c(nrow(rc_retention1), ncol(rc_retention1), 2))
            rc_retention[, , 1] <- as.matrix(rc_retention1)
            rc_retention[, , 2] <- as.matrix(rc_retention2)
            iso2 <- rowMeans(rc_retention, dims = 2, na.rm = T)
            psi <- iso2 / (iso1 + iso2)
            psi <- as.data.frame(psi)
            row.names(psi) <- events$event
            colnames(psi) <- colnames(rc_exclusion)
        } else if (type == "MXE") {
            rc_exclusion1 <- data$rc_exclusion1
            rc_exclusion2 <- data$rc_exclusion2
            rc_retention1 <- data$rc_retention1
            rc_retention2 <- data$rc_retention2
            rc_exclusion <- array(0, dim = c(nrow(rc_exclusion1), ncol(rc_exclusion2), 2))
            rc_exclusion[, , 1] <- as.matrix(rc_exclusion1)
            rc_exclusion[, , 2] <- as.matrix(rc_exclusion2)
            iso1 <- rowMeans(rc_exclusion, dims = 2, na.rm = T)
            rc_retention <- array(0, dim = c(nrow(rc_retention1), ncol(rc_retention1), 2))
            rc_retention[, , 1] <- as.matrix(rc_retention1)
            rc_retention[, , 2] <- as.matrix(rc_retention2)
            iso2 <- rowMeans(rc_retention, dims = 2, na.rm = T)
            psi <- iso2 / (iso1 + iso2)
            psi <- as.data.frame(psi)
            row.names(psi) <- events$event
            colnames(psi) <- colnames(rc_exclusion1)
        } else if (type == "A3SS" | type == "A5SS" | type == "AL") {
            rc_exclusion <- data$rc_exclusion
            rc_retention <- data$rc_retention
            iso1 <- rc_exclusion
            iso2 <- rc_retention
            psi <- iso2 / (iso1 + iso2)
            psi <- as.data.frame(psi)
            row.names(psi) <- events$event
            colnames(psi) <- colnames(rc_exclusion)
        }
        saveRDS(psi, paste0(rc_path, "/", type, "_psi.rds"))
        msg <- paste0("[", Sys.time(), "] ", "Calculating PSI value of ", type, " events Finish.")
        print(msg)
    }
    return(rc_path)
}
