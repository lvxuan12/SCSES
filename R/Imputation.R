#' @include utilities.R

#' @title Perform SCSES Imputation
#' @description
#' Performs PSI imputation using three distinct strategies
#' that leverage cell similarity and splicing event similarity:
#' \itemize{
#'   \item \strong{Strategy 1}: Direct imputation of raw PSI values using cell similarities
#'   \item \strong{Strategy 2}: Imputation of raw read counts using cell similarities,
#'         followed by PSI calculation from imputed counts
#'   \item \strong{Strategy 3}: Sequential imputation combining Strategy 2 results with
#'         event similarity-based imputation
#' }
#'
#' @param paras A list containing configuration parameters, typically loaded via
#'   \code{readSCSESconfig(paras_file)}.
#' @param output_path Character string specifying the output directory path.
#'   If \code{NULL}, defaults to \code{work_path/imputation/}.
#' @param decay_impute Numeric value (0-1) representing the convergence threshold for
#'   imputation iterations. The algorithm stops when changes between consecutive
#'   iterations fall below this threshold. Default: 0.05 (from \code{paras$Task$impute$decay_impute}).
#' @param rc Matrix of normalized splicing event-associated read counts, or path to
#'   an RDS file. Default: \code{rds_processed/rc.rds}.
#' @param psi Matrix of PSI values, or path to an RDS file.
#'   Default: \code{rds_processed/psi.rds}.
#' @param event.info List containing splicing event information, or path to an RDS file.
#'   Default: \code{rds_processed/event.info.list.rds}.
#' @param cell_similarity Named list of cell similarity matrices for different features
#'   (EXP_RBP, RC, PSI), or path to an RDS file.
#'   Default: \code{imputation/cell_similarity/cell.similars.rds}.
#' @param event_similarity Named list of event similarity matrices for different splicing
#'   event types, or path to an RDS file.
#'   Default: \code{imputation/event_similarity/event.similars.rds}.
#'
#' @return
#' Character string containing the file path to the saved imputation results
#' (\code{work_path/imputation/Imputed_seperated_*.rds}). The saved object contains
#' a list with two main components:
#' \describe{
#'   \item{\code{cell}}{Results from Strategy 1 and Strategy 2. Matrix names follow
#'     the pattern \code{[CellSimilarityType]_[Strategy]}, where:
#'     \itemize{
#'       \item Cell similarity types: EXP_RBP, RC, PSI
#'       \item Strategy indicators: PSI (Strategy 1), RC (Strategy 2)
#'     }
#'   }
#'   \item{\code{cell_event}}{Results from Strategy 3. Matrix names follow the pattern
#'     \code{[CellSimilarityType]_PSI}, representing the results after event similarity
#'     imputation
#'   }
#' }
#'
#' @export
#' @import R.matlab
#' @import raveio
#' @import rhdf5
#' @import hdf5r
#' @importFrom stats rbinom
#'
ImputationAll <- function(
    paras, output_path = NULL,
    decay_impute = paras$Task$impute$decay_impute,
    rc = NULL, psi = NULL, event.info = NULL,
    cell_similarity = NULL, event_similarity = NULL) {
    # script----
    matlab_path <- system.file("matlab", package = "SCSES")
    mat_scses <- paste0(matlab_path, "/scses/run_scses.sh")
    mcr_path <- paras$Basic$mcr_path
    msg <- paste0("[", Sys.time(), "] ", "Get imputed result using cell similarity and event similarity.")
    print(msg)
    # validate input----
    print(paste("Checking data..."))
    rc <- check.readRDS(rc, default = paste0(paras$Basic$work_path, "/rds_processed/rc.rds"))
    print("rc checked")
    psi <- check.readRDS(psi, default = paste0(paras$Basic$work_path, "/rds_processed/psi.rds"))
    print("psi checked")
    event.info <- check.readRDS2(event.info, default = paste0(paras$Basic$work_path, "/rds_processed/event.info.list.rds"))
    print("event.info checked")
    cell_similarity <- check.readRDS2(cell_similarity, default = paste0(paras$Basic$work_path, "/imputation/cell_similarity/cell.similars.rds"))
    print("cell similarity checked")
    event_similarity <- check.readRDS2(event_similarity, default = paste0(paras$Basic$work_path, "/imputation/event_similarity/event.similars.rds"))
    print("event similarity checked")
    decay_impute <- check.double.or.null(x = decay_impute, default = 0.05)
    print(paste0("decay=", paste(decay_impute, collapse = ";"), "  checked"))
    # validate event type----
    print(paste("Checking event type"))
    event_type1 <- names(event.info)
    event_type2 <- names(event_similarity)
    event_type1 <- check.valid(x = event_type1, select = c("A3SS", "A5SS", "AL", "SE", "MXE", "RI"))
    event_type2 <- check.valid(x = event_type2, select = c("A3SS", "A5SS", "AL", "SE", "MXE", "RI"))
    event_type <- intersect(event_type1, event_type2)
    print(paste0("event_type=", paste(event_type, collapse = ";"), "  checked"))
    # validate cell similarity type----
    print(paste("Checking cell similarity type"))
    cell_similarity_data <- names(cell_similarity)
    cell_similarity_data <- check.valid(x = cell_similarity_data, select = c("EXP_RBP", "RC", "PSI"))
    print(paste0("cell_similarity_data=", paste(cell_similarity_data, collapse = ";"), "  checked"))

    # Output
    if (is.null(output_path)) {
        output_path <- paste0(paras$Basic$work_path, "/imputation/")
    }
    dir.create(output_path)
    print(paste0("Output: ", output_path))

    # log
    log_file <- paste0(output_path, "/mat_Imputation.log")
    if (file.exists(log_file)) {
      file.remove(log_file)
    }

    psi.psi.imputation.cell <- list()
    rc.psi.imputation.cell <- list()
    psi.imputation.event <- list()
    for (type in event_type)
    {
        gc()
        events <- event.info[[type]]
        psi.psi.results.cell <- list()
        rc.psi.results.cell <- list()
        psi.results.event <- list()
        for (data_type in cell_similarity_data)
        {
            msg <- paste0("[", Sys.time(), "] ", "Running ", paste0("Event_type=", type, ";cell_similarity_feature=", data_type))
            print(msg)
            # similarity data
            cell_similarity_type <- cell_similarity[[data_type]]
            event_similarity_type <- event_similarity[[type]]
            similarity_type <- list(cell = as.matrix(cell_similarity_type),
                                    event = as.matrix(event_similarity_type))

            token <- paste(Sys.getpid(), rbinom(1, size = 1000000000, prob = 0.5), sep = "-")
            datapath <- paste0(output_path, "/imputation_data_", token, ".h5")
            resultpath <- paste0(output_path, "/imputation_result_", token, ".mat")
            # rc data
            rc_type <- rc[do.call(what = c, args = events[, -1]), , drop = F]
            rc_data <- list()
            for (position in colnames(events[, -1]))
            {
                tmp <- list(as.matrix(rc_type[events[[position]], , drop = F]))
                rc_data <- c(rc_data, tmp)
            }
            names(rc_data) <- colnames(events[, -1])
            # psi data
            psi_type <- psi[events$event, , drop = F]
            psi_data <- list(as.matrix(psi_type))
            names(psi_data) <- "all"
            all_data <- list(PSI = psi_data, RC = rc_data)
            msg <- paste0("[", Sys.time(), "] ", "Save data")
            print(msg)
            saveHdf5File(datapath, list(
                similar = similarity_type, data = all_data,
                parameter = list(decay = decay_impute)
            ))
            msg <- paste0("[", Sys.time(), "] ", "Save data Finished")
            print(msg)
            # run scses
            cmd <- paste("bash", mat_scses, mcr_path, datapath,
                         resultpath, ">>", log_file, "2>&1")
            print(cmd)
            system(cmd, wait = T)
            # Output formatting
            result <- read_mat(resultpath)
            result[[1]] <- NULL
            names(result) <- gsub("data/", "", names(result))
            result <- lapply(result, function(x) {
                if (is.vector(x)) {
                    x <- t(as.matrix(x))
                }
                row.names(x) <- rownames(psi_type)
                colnames(x) <- colnames(psi_type)
                return(x)
            })
            file.remove(datapath, resultpath)
            psi.psi.results.cell[[data_type]] <- result[[1]]
            rc.psi.results.cell[[data_type]] <- result[[2]]
            psi.results.event[[data_type]] <- result[[3]]
        }
        psi.psi.imputation.cell[[type]] <- psi.psi.results.cell
        rc.psi.imputation.cell[[type]] <- rc.psi.results.cell
        psi.imputation.event[[type]] <- psi.results.event
    }
    res_list1 <- list()
    for (data_type in cell_similarity_data) {
        PSI_imputed_res <- do.call(what = rbind, args = lapply(psi.psi.imputation.cell, FUN = function(x) {
            return(x[[data_type]])
        }))
        RC_PSI_imputed_res <- do.call(what = rbind, args = lapply(rc.psi.imputation.cell, FUN = function(x) {
            return(x[[data_type]])
        }))
        res_list1[[paste0(data_type, "_PSI")]] <- PSI_imputed_res
        res_list1[[paste0(data_type, "_RC")]] <- RC_PSI_imputed_res
    }
    res_list2 <- list()
    for (data_type in cell_similarity_data) {
        PSI_imputed_res <- do.call(what = rbind, args = lapply(psi.imputation.event, FUN = function(x) {
            return(x[[data_type]])
        }))
        res_list2[[paste0(data_type, "_PSI")]] <- PSI_imputed_res
    }
    psi_imputed_seperated <- list(cell = res_list1, cell_event = res_list2)
    msg <- paste0("[", Sys.time(), "] ", "Get imputed result using cell similarity and event similarity Finish.")
    print(msg)
    output_name <- paste0(
        output_path, "/Imputed_seperated_",
        rbinom(1, size = 1000000000, prob = 0.5),
        ".rds"
    )
    saveRDS(psi_imputed_seperated, output_name)
    return(output_name)
}

#' @title Perform SCSES Imputation for Specific Cell Similarity Type
#' @description
#' Executes three PSI imputation strategies using a single
#' specified type of cell similarity feature. This function provides
#' imputation analysis for one cell similarity type at a time.
#'
#' The three imputation strategies are:
#' \itemize{
#'   \item \strong{Strategy 1}: Direct imputation of raw PSI values using cell similarities
#'   \item \strong{Strategy 2}: Imputation of raw read counts using cell similarities,
#'         followed by PSI calculation from imputed counts
#'   \item \strong{Strategy 3}: Sequential imputation combining Strategy 2 results with
#'         event similarity-based imputation
#' }
#'
#' @param paras A list containing configuration parameters, typically loaded via
#'   \code{readSCSESconfig(paras_file)}.
#' @param output_path Character string specifying the output directory path.
#'   If \code{NULL}, defaults to \code{work_path/imputation/}.
#' @param decay_impute Numeric value (0-1) representing the convergence threshold for
#'   imputation iterations. The algorithm stops when changes between consecutive
#'   iterations fall below this threshold. Default: 0.05 (from \code{paras$Task$impute$decay_impute}).
#' @param rc Matrix of normalized splicing event-associated read counts, or path to
#'   an RDS file. Default: \code{rds_processed/rc.rds}.
#' @param psi Matrix of PSI values, or path to an RDS file.
#'   Default: \code{rds_processed/psi.rds}.
#' @param event.info List containing splicing event information, or path to an RDS file.
#'   Default: \code{rds_processed/event.info.list.rds}.
#' @param cell_similarity Named list of cell similarity matrices for different features
#'   (EXP_RBP, RC, PSI), or path to an RDS file.
#'   Default: \code{imputation/cell_similarity/cell.similars.rds}.
#' @param event_similarity Named list of event similarity matrices for different splicing
#'   event types, or path to an RDS file.
#'   Default: \code{imputation/event_similarity/event.similars.rds}.
#' @param cell_similarity_type Character string specifying the single cell similarity
#'   type to use for imputation. Must be one of: \code{"EXP_RBP"}, \code{"RC"}, or \code{"PSI"}.
#'
#' @return
#' Character string containing the file path to the saved imputation results
#' (\code{work_path/imputation/Imputed_seperated_*.rds}). The saved object contains
#' a list with two main components:
#' \describe{
#'   \item{\code{cell}}{Results from Strategy 1 and Strategy 2. Matrix names follow
#'     the pattern \code{[cell_similarity_type]_[Strategy]}, where:
#'     \itemize{
#'       \item Cell similarity types: EXP_RBP, RC, PSI
#'       \item Strategy indicators: PSI (Strategy 1), RC (Strategy 2)
#'     }
#'   }
#'   \item{\code{cell_event}}{Results from Strategy 3. Matrix names follow the pattern
#'     \code{[cell_similarity_type]_PSI}, representing the results after event similarity
#'     imputation
#'   }
#' }
#'
#'
#' @export
#' @import R.matlab
#' @import raveio
#' @import rhdf5
#' @import hdf5r
#' @importFrom stats rbinom
#'
Imputation <- function(
    paras, output_path = NULL,
    decay_impute = paras$Task$impute$decay_impute,
    rc = NULL, psi = NULL, event.info = NULL,
    cell_similarity = NULL, event_similarity = NULL,
    cell_similarity_type) {
    # script----
    matlab_path <- system.file("matlab", package = "SCSES")
    mat_scses <- paste0(matlab_path, "/scses/run_scses.sh")
    mcr_path <- paras$Basic$mcr_path
    msg <- paste0("[", Sys.time(), "] ", "Get imputed result using cell similarity and event similarity.")
    print(msg)
    # validate input----
    print(paste("Checking data..."))
    rc <- check.readRDS(rc, default = paste0(paras$Basic$work_path, "/rds_processed/rc.rds"))
    print("rc checked")
    psi <- check.readRDS(psi, default = paste0(paras$Basic$work_path, "/rds_processed/psi.rds"))
    print("psi checked")
    event.info <- check.readRDS2(event.info, default = paste0(paras$Basic$work_path, "/rds_processed/event.info.list.rds"))
    print("event.info checked")
    cell_similarity <- check.readRDS2(cell_similarity, default = paste0(paras$Basic$work_path, "/imputation/cell_similarity/cell.similars.rds"))
    print("cell similarity checked")
    event_similarity <- check.readRDS2(event_similarity, default = paste0(paras$Basic$work_path, "/imputation/event_similarity/event.similars.rds"))
    print("event similarity checked")
    decay_impute <- check.double.or.null(x = decay_impute, default = 0.05)
    print(paste0("decay=", paste(decay_impute, collapse = ";"), "  checked"))
    # validate event type----
    print(paste("Checking event type"))
    event_type1 <- names(event.info)
    event_type2 <- names(event_similarity)
    event_type1 <- check.valid(x = event_type1, select = c("A3SS", "A5SS", "AL", "SE", "MXE", "RI"))
    event_type2 <- check.valid(x = event_type2, select = c("A3SS", "A5SS", "AL", "SE", "MXE", "RI"))
    event_type <- intersect(event_type1, event_type2)
    print(paste0("event_type=", paste(event_type, collapse = ";"), "  checked"))
    # validate cell similarity type----
    print(paste("Checking cell similarity type"))
    cell_similarity_data <- names(cell_similarity)
    cell_similarity_data <- intersect(cell_similarity_data, cell_similarity_type)
    cell_similarity_data <- check.valid(x = cell_similarity_data, select = c("EXP_RBP", "RC", "PSI"))
    print(paste0("cell_similarity_data=", paste(cell_similarity_data, collapse = ";"), "  checked"))

    # Output
    if (is.null(output_path)) {
        output_path <- paste0(paras$Basic$work_path, "/imputation/")
    }
    dir.create(output_path)
    print(paste0("Output: ", output_path))

    # log
    log_file <- paste0(output_path, "/mat_Imputation_",cell_similarity_type,".log")
    if (file.exists(log_file)) {
      file.remove(log_file)
    }

    psi.psi.imputation.cell <- list()
    rc.psi.imputation.cell <- list()
    psi.imputation.event <- list()
    for (type in event_type)
    {
        gc()
        events <- event.info[[type]]
        msg <- paste0("[", Sys.time(), "] ", "Running ", paste0("Event_type=", type, ";cell_similarity_feature=", cell_similarity_type))
        print(msg)
        similarity_type <- list(
            cell = cell_similarity[[cell_similarity_type]],
            event = event_similarity[[type]]
        )

        token <- paste(Sys.getpid(), rbinom(1, size = 1000000000, prob = 0.5), sep = "-")
        datapath <- paste0(output_path, "/imputation_data_", token, ".h5")
        resultpath <- paste0(output_path, "/imputation_result_", token, ".mat")
        # rc data
        rc_type <- rc[do.call(what = c, args = events[, -1]), , drop = F]
        rc_data <- list()
        for (position in colnames(events[, -1]))
        {
            tmp <- list(as.matrix(rc_type[events[[position]], , drop = F]))
            rc_data <- c(rc_data, tmp)
        }
        names(rc_data) <- colnames(events[, -1])
        # psi data
        psi_type <- psi[events$event, , drop = F]
        psi_data <- list(as.matrix(psi_type))
        names(psi_data) <- "all"
        all_data <- list(PSI = psi_data, RC = rc_data)
        print("Save data")
        saveHdf5File(datapath, list(
            similar = similarity_type, data = all_data,
            parameter = list(decay = decay_impute)
        ))
        print("Save data Finished")
        # run scses
        cmd <- paste("bash", mat_scses, mcr_path, datapath,
                     resultpath, ">>", log_file, "2>&1")
        print(cmd)
        system(cmd, wait = T)
        # Output formatting
        result <- read_mat(resultpath)
        result[[1]] <- NULL
        names(result) <- gsub("data/", "", names(result))
        result <- lapply(result, function(x) {
            if (is.vector(x)) {
                x <- t(as.matrix(x))
            }
            row.names(x) <- rownames(psi_type)
            colnames(x) <- colnames(psi_type)
            return(x)
        })
        file.remove(datapath, resultpath)
        psi.psi.imputation.cell[[type]] <- result[[1]]
        rc.psi.imputation.cell[[type]] <- result[[2]]
        psi.imputation.event[[type]] <- result[[3]]
    }
    res_list1 <- list()
    PSI_imputed_res <- do.call(what = rbind, args = psi.psi.imputation.cell)
    RC_PSI_imputed_res <- do.call(what = rbind, args = rc.psi.imputation.cell)
    res_list1[[paste0(cell_similarity_type, "_PSI")]] <- PSI_imputed_res
    res_list1[[paste0(cell_similarity_type, "_RC")]] <- RC_PSI_imputed_res
    res_list2 <- list()
    PSI_imputed_res <- do.call(what = rbind, args = psi.imputation.event)
    res_list2[[paste0(cell_similarity_type, "_PSI")]] <- PSI_imputed_res
    psi_imputed_seperated <- list(cell = res_list1, cell_event = res_list2)
    msg <- paste0("[", Sys.time(), "] ", "Get imputed result using cell similarity and event similarity Finish.")
    print(msg)
    output_name <- paste0(
        output_path, "/Imputed_seperated_",cell_similarity_type,"_",
        rbinom(1, size = 1000000000, prob = 0.5),
        ".rds"
    )
    saveRDS(psi_imputed_seperated, output_name)
    return(output_name)
}

