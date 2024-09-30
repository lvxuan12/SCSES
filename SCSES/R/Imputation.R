#' @include utilities.R

#' @title Perform SCSES imputation
#' @description run three PSI imputation strategies used
#' All types of cell similarity and splicing events are used
#' first: strategy1 impute the raw PSI with cell similarities
#' second: strategy2 first impute raw read counts with cell similarities
#' and calculate the imputed PSI
#' third: strategy3 first impute the PSI using strategy 2, then further impute
#'  the results using event similarities

#'
#' @param paras list fromJSON(paras_file)
#' Default decay_impute from paras
#' @param output_path path to event similarity
#' @param decay_impute threshold of change in the similarity matrix
#' @param rc matrix of normalized splicing events associated read counts or path
#' to a rds file, default: rds_processed/rc.rds
#' @param psi matrix of psi value or path to a rds file,
#' default: rds_processed/psi.rds
#' @param event.info list of events information or path to a rds file,
#' default: rds_processed/event.info.list.rds
#' @param cell_similarity a list of different types of cell similarity named by the data
#' type (EXP_RBP, RC, and PSI) or path to a rds file, 
#' default: imputation/cell_similarity/cell.similars.rds
#' @param event_similarity a list of event similarity named by different types of
#' splicing events or path to a rds file, 
#' default: imputation/event_similarity/event.similars.rds
#' 
#' @return path to the result of three imputation strategies
#' save a list of three imputation strategies result to
#' work_path/imputation/Imputed_seperated_*.rds.
#' The result contains two lists, each of which contains multiple imputation results.
#' Different event types are combined.
#' The first list is the result of strategy1 and strategy2, where the name of each matrix 
#' is divided by _, the first half represents different cell similarities, and the second 
#' half represents PSI for strategy1 and RC for strategy2; 
#' The second list is the result of strategy3, again where the names of each matrix are 
#' divided by _, with the first half representing different cell similarities
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
    decay_impute <- check.int.or.null(x = decay_impute, default = 0.05)
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
            msg <- paste0("[", Sys.time(), "] ", "Running ", paste0("Event_type=", type, ";similarity_type=", data_type))
            print(msg)
            # similarity data
            cell_similarity_type <- cell_similarity[[data_type]]
            event_similarity_type <- event_similarity[[type]]
            similarity_type <- list(cell = cell_similarity_type, event = event_similarity_type)

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
                similarity = similarity_type, data = all_data,
                parameter = list(decay = decay_impute)
            ))
            msg <- paste0("[", Sys.time(), "] ", "Save data Finished")
            print(msg)
            # run scses
            cmd <- paste("bash", mat_scses, mcr_path, datapath, resultpath)
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
                row.names(x) <- rownames(psi)
                colnames(x) <- colnames(psi)
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
    return(output_path)
}

#' @title Perform SCSES imputation
#' @description run three PSI imputation strategies
#' for a sepcific type of cell similarity
#'
#' @param paras list fromJSON(paras_file)
#' Default decay_impute from paras
#' @param output_path path to event similarity
#' @param decay_impute threshold of change in the similarity matrix
#' @param rc matrix of normalized splicing events associated read counts or path
#' to a rds file, default: rds_processed/rc.rds
#' @param psi matrix of psi value or path to a rds file,
#' default: rds_processed/psi.rds
#' @param event.info list of events information or path to a rds file,
#' default: rds_processed/event.info.list.rds
#' @param cell_similarity a list of different types of cell similarity named by the data
#' type (EXP_RBP, RC, and PSI) or path to a rds file,
#' default: imputation/cell_similarity/cell.similars.rds
#' @param event_similarity a list of event similarity named by different types of
#' splicing events or path to a rds file,
#' default: imputation/event_similarity/event.similars.rds
#' @param cell_similarity_type one data type used to calculate cell similarity (EXP_RBP, RC, PSI)
#' 
#' 
#' @return path to the result of three imputation strategies
#' save a list of result of three imputation strategies to
#' work_path/imputation/Imputed_seperated_event_${cell_similarity_type}_*.rds
#' The result contains two lists, each of which contains multiple imputation results.
#' The first list is the result of strategy1 and strategy2, 
#' The second list is the result of strategy3
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
    decay_impute <- check.int.or.null(x = decay_impute, default = 0.05)
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

    psi.psi.imputation.cell <- list()
    rc.psi.imputation.cell <- list()
    psi.imputation.event <- list()
    for (type in event_type)
    {
        gc()
        events <- event.info[[type]]
        msg <- paste0("[", Sys.time(), "] ", "Running ", paste0("Event_type=", type, ";similarity_type=", cell_similarity_type))
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
        msg <- paste0("[", Sys.time(), "] ", "Save data")
        print(msg)
        saveHdf5File(datapath, list(
            similarity = similarity_type, data = all_data,
            parameter = list(decay = decay_impute)
        ))
        msg <- paste0("[", Sys.time(), "] ", "Save data Finished")
        print(msg)
        # run scses
        cmd <- paste("bash", mat_scses, mcr_path, datapath, resultpath)
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
            row.names(x) <- rownames(psi)
            colnames(x) <- colnames(psi)
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
    return(output_path)
}

