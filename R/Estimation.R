#' @include utilities.R

#' @title Calculate features for classifer
#'
#' @param cell_similarity matrix of cell similarity
#' @param psi matrix of psi value
#' @param rc matrix of normalized splicing events associated read counts
#' @param expr matrix of gene expression (TPM for full-length protocols
#' and normalized UMI counts for droplet-based protocols)
#' @param event data.frame of events information
#' @param k the number of neighbor cells
#' @param rc_imputed matrix of psi value calculated using reads after
#' imputation with information from similar cells  (Strategy 2)
#' @param psi_imputed matrix of imputed psi with information from similar cells (Strategy 1)
#'
#'
#' @return a data.frame of features
#'
#' @export

getClassifierFeature <- function(
    cell_similarity, psi, rc, expr, event,
    k, rc_imputed, psi_imputed) {
    k <- as.numeric(k)
    cell_similarity <- as.matrix(cell_similarity)
    rc_ex <- calcu_ex_rc(rc, event)
    rc_in <- calcu_in_rc(rc, event)
    idx_na <- which(rc_in == 0 & rc_ex == 0)
    idx_0 <- which(rc_in == 0 & rc_ex != 0)
    idx_1 <- which(rc_in != 0 & rc_ex == 0)
    idx_other <- which(rc_in != 0 & rc_ex != 0)
    m <- matrix(data = NA, nrow = nrow(psi), ncol = ncol(psi))
    m[idx_na] <- -1
    m[idx_0] <- 0
    m[idx_1] <- 1
    m[idx_other] <- 2
    cs <- cell_similarity
    diag(cs) <- 0
    cs[which(cs != 0)] <- 1
    data_knn <- t(t(as.matrix(psi) %*% t(cs)) / k)
    cs <- cell_similarity
    cs[which(cs != 0)] <- 1
    deltapsi_self_knn <- abs(psi - data_knn)
    tmp.array <- array(0, dim = c(nrow(event), ncol(rc), 4))
    for (i in seq(2, 5))
    {
        tmp <- rc[match(event[, i], row.names(rc)), ]
        tmp.array[, , i - 1] <- as.matrix(tmp)
    }
    tmp.mean <- rowMeans(tmp.array, dims = 2, na.rm = T)
    rc_self <- tmp.mean
    cs <- cell_similarity
    diag(cs) <- 0
    cs[which(cs != 0)] <- 1
    rc_knn <- t(t(tmp.mean %*% t(cs)) / k)
    knn_vc <- knn_cv(psi, cs)
    cs <- cell_similarity
    diag(cs) <- 0
    cs[which(cs != 0)] <- 1
    rc_in_knn <- t(t(rc_in %*% t(cs)) / k)
    rc_ex_knn <- t(t(rc_ex %*% t(cs)) / k)
    rc_in_0 <- rc_in == 0
    rc_ex_0 <- rc_ex == 0
    rc_in_non_0 <- rc_in != 0
    rc_ex_non_0 <- rc_ex != 0
    # na num
    rc_total_na <- rc_in_0 + rc_ex_0
    cs <- cell_similarity
    diag(cs) <- 0
    cs[which(cs != 0)] <- 1
    knn_rc_na <- t(t((rc_total_na == 2) %*% t(cs)) / k)
    # 0 num
    rc_total_0 <- rc_in_0 + rc_ex_non_0
    cs <- cell_similarity
    diag(cs) <- 0
    cs[which(cs != 0)] <- 1
    knn_rc_0 <- t(t((rc_total_0 == 2) %*% t(cs)) / k)
    # 1 num
    rc_total_1 <- rc_in_non_0 + rc_ex_0
    cs <- cell_similarity
    diag(cs) <- 0
    cs[which(cs != 0)] <- 1
    knn_rc_1 <- t(t((rc_total_1 == 2) %*% t(cs)) / k)
    deltapsi_impute_train_rc <- abs(rc_imputed - psi)
    deltapsi_impute_train_psi <- abs(psi_imputed - psi)
    expr <- expr[, colnames(psi)]
    gene <- unlist(lapply(row.names(psi), function(x) {
        unlist(strsplit(x, "\\|"))[3]
    }))

    expr2 <- rbind(expr[intersect(gene, rownames(expr)), , drop = FALSE],
                   matrix(0,
                          nrow = length(setdiff(gene, rownames(expr))),
                          ncol = ncol(expr),
                          dimnames = list(setdiff(gene, rownames(expr)),
                                          colnames(expr))))[gene, ]
    feature_df <- data.frame(
        knn_vc = unlist(as.vector(knn_vc)),
        rc_knn = unlist(as.vector(rc_knn)),
        rc_self = unlist(as.vector(rc_self)),
        rc_in_self = unlist(as.vector(rc_in)),
        rc_ex_self = unlist(as.vector(rc_ex)),
        rc_in_knn = unlist(as.vector(rc_in_knn)),
        rc_ex_knn = unlist(as.vector(rc_ex_knn)),
        deltapsi_self_knn = unlist(as.vector(deltapsi_self_knn)),
        psi_self = unlist(as.vector(psi)),
        psi_knn = unlist(as.vector(data_knn)),
        knn_rc_na = unlist(as.vector(knn_rc_na)),
        knn_rc_0 = unlist(as.vector(knn_rc_0)),
        knn_rc_1 = unlist(as.vector(knn_rc_1)),
        deltapsi_impute_train_rc = unlist(as.vector(deltapsi_impute_train_rc)),
        deltapsi_impute_train_psi = unlist(as.vector(deltapsi_impute_train_psi)),
        expr = unlist(as.vector(expr2)),
        psi_self_type = unlist(as.vector(m))
    )
    # feature preprocess
    feature_df[is.na(feature_df)] <- 0
    # normalize input features
    normalize <- function(x) {
        return((x - min(x)) / (max(x) - min(x)))
    }
    feature_df$rc_knn <- normalize(feature_df$rc_knn)
    feature_df$rc_self <- normalize(feature_df$rc_self)
    feature_df$rc_in_self <- normalize(feature_df$rc_in_self)
    feature_df$rc_ex_self <- normalize(feature_df$rc_ex_self)
    feature_df$rc_in_knn <- normalize(feature_df$rc_in_knn)
    feature_df$rc_ex_knn <- normalize(feature_df$rc_ex_knn)
    feature_df$expr <- normalize(feature_df$expr)
    feature_df$psi_self_type.0 <- ifelse(feature_df$psi_self_type == 0, 1, 0)
    feature_df$psi_self_type.1 <- ifelse(feature_df$psi_self_type == 1, 1, 0)
    feature_df$psi_self_type.2 <- ifelse(feature_df$psi_self_type == 2, 1, 0)
    feature_df$`psi_self_type.-1` <- ifelse(feature_df$psi_self_type == -1, 1, 0)
    feature_df <- feature_df[, -which(colnames(feature_df) == "psi_self_type")]
    return(feature_df)
}

#' @title Fine-Tune the Classifier for Event-Cell Group Assignment
#' @description This function fine-tunes a two-stage classifier to determine
#'   the groups that event-cell pairs belong to. Model1 predicts probabilities
#'   for BD vs. TD, while Model2 distinguishes between TD+Info and TD-Info.
#'
#' @param paras A list object containing SCSES configuration parameters, typically
#'   loaded using \code{readSCSESconfig(paras_file)}.
#' @param rds_path Character string specifying the path to processed RDS data
#'   containing gene expression matrices. Default: \code{work_path/rds/}.
#' @param rds_ft_path Character string specifying the path to RDS data of
#'   fine-tuning splicing events (PSI, RC, event annotations).
#'   Default: \code{work_path/rds_ft/}.
#' @param rds_cell_similarity_path Character string specifying the path to cell
#'   similarity data.
#'   Default: \code{work_path/imputation/cell_similarity/}.
#' @param output_path Character string specifying the output directory for
#'   fine-tuned classifiers. Default: \code{work_path/classifer/}.
#' @param decay_impute Numeric value specifying the convergence threshold
#' for imputation. Default: 0.05 (from \code{paras$Task$impute$decay_impute}).
#' @param genome_name Character string specifying the genome version ("hg19" or "hg38").
#'   Default: from \code{paras$Basic$refgenome$genome_name}.
#'
#' @return Character string specifying the path to the directory containing
#'   the fine-tuned classifier models.
#'
#' @export
#' @importFrom stats rbinom coef predict quantile
#' @importFrom caret downSample
#' @importFrom glmnet cv.glmnet
#' @import R.matlab
#' @import raveio
#' @import rhdf5
#' @import hdf5r
#'
#'
FtClassifier <- function(
    paras, rds_path = NULL, rds_ft_path = NULL,
    rds_cell_similarity_path = NULL,
    output_path = NULL,
    decay_impute = paras$Task$impute$decay_impute,
    genome_name = paras$Basic$refgenome$genome_name) {
    # script----
    matlab_path <- system.file("matlab", package = "SCSES")
    mat_impute_v1 <- paste0(matlab_path, "/imputation1/run_imputation1.sh")
    mcr_path <- paras$Basic$mcr_path
    # build in data----
    extdata_path <- system.file("extdata", package = "SCSES")
    psi_gtex_path <- paste0(extdata_path, "/ftevents/", genome_name, "/psi/psi_gtex_select.txt")
    print("Reading true Ft PSI...")
    psi_gtex <- read.table(psi_gtex_path, header = T)
    psi_gtex <- rowMeans(psi_gtex)
    print("Loading Pre-training classifer...")
    model1_path = paste0(extdata_path, "/model/model_change_nonchange.rdata")
    model2_path = paste0(extdata_path, "/model/model_change_01_change_other.rdata")
    model1 <- get(load(model1_path))
    model2 <- get(load(model2_path))
    msg <- paste0("[", Sys.time(), "] ", "Classifer fine tune")
    print(msg)
    msg <- paste0("[", Sys.time(), "] ", "Processing raw Ft data...")
    print(msg)
    # input rds----
    if (is.null(rds_ft_path)) {
        rds_ft_path <- paste0(paras$Basic$work_path, "/rds_ft/")
    }
    if (is.null(rds_path)) {
        rds_path <- paste0(paras$Basic$work_path, "/rds/")
    }
    if (is.null(rds_cell_similarity_path)) {
        rds_cell_similarity_path <- paste0(paras$Basic$work_path, "/imputation/cell_similarity/")
    }
    print("Checking data...")
    rds_ft_files = list.files(path = rds_ft_path, pattern = "*rds")
    if (!all(c("psi.rds", "rc.rds", "event.rds") %in% rds_ft_files)) {
        stop(paste("Preprocessed psi, reads count, and events annotation must be saved in", rds_ft_path))
    }
    rds_files = list.files(path = rds_path, pattern = "*rds")
    if (!all(c("count_norm.rds") %in% rds_files)) {
        stop(paste("Gene expression must be saved in", rds_path))
    }
    # output----
    if (is.null(output_path)) {
        output_path <- paste0(paras$Basic$work_path, "/classifer/")
    }
    if (!dir.exists(output_path)) {
        dir.create(output_path)
    }
    # psi input
    psi <- readRDS(file = paste0(rds_ft_path, "/psi.rds"))
    # reads count input
    rc <- readRDS(file = paste0(rds_ft_path, "/rc.rds"))
    event <- readRDS(file = paste0(rds_ft_path, "/event.rds"))
    # expr input
    expr <- readRDS(file = paste0(rds_path, "/count_norm.rds"))
    cell_id <- Reduce(intersect, list(colnames(expr), colnames(psi), colnames(rc)))
    psi <- psi[, cell_id, drop = F]
    rc <- rc[, cell_id, drop = F]
    expr <- expr[, cell_id, drop = F]
    # cell similarity input
    cell_similarity <- readRDS(file = paste0(rds_cell_similarity_path, "/cell.similars.rds"))
    dyk_cell <- readRDS(file = paste0(rds_cell_similarity_path, "/dyk.cell.rds"))
    # validate cell similarity type----
    print(paste("Checking cell similarity type"))
    cell_similarity_data1 = names(cell_similarity)
    cell_similarity_data1 <- check.valid(x = cell_similarity_data1, select = c("EXP_RBP", "RC", "PSI"))
    cell_similarity_data2 = names(dyk_cell)
    cell_similarity_data2 <- check.valid(x = cell_similarity_data2, select = c("EXP_RBP", "RC", "PSI"))
    cell_similarity_data <- intersect(cell_similarity_data1, cell_similarity_data2)
    print(paste0("cell_similarity_data=", paste(cell_similarity_data, collapse = ";"), "  checked"))
    # validate parameters----
    decay_impute <- check.double.or.null(x = decay_impute, default = 0.05)

    psi <- psi[intersect(names(psi_gtex), row.names(psi)), ]
    event <- event[match(row.names(psi), event$event), ]
    event_types <- unique(event$type)
    event.info <- lapply(X = as.list(event_types), FUN = function(type) {
        cur.event <- event[event$type == type, ]
        if (type == "RI" | type == "SE") {
            info <- cur.event[, c("event", "exclusion1", "retention1", "retention2")]
        } else if (type == "MXE") {
            info <- cur.event[, c("event", "exclusion1", "exclusion2", "retention1", "retention2")]
        } else if (type == "A3SS" | type == "A5SS" | type == "AL") {
            info <- cur.event[, c("event", "exclusion1", "retention1")]
        }
        return(info)
    })
    names(event.info) <- event_types

    print("Calculating Classifier Features...")
    log_file <- paste0(output_path, "/mat_calculateFtFeature.log")
    if (file.exists(log_file)) {
      file.remove(log_file)
    }
    rc.imputation.cell <- list()
    rc.psi.imputation.cell <- list()
    psi.imputation.cell <- list()
    for (type in event_types)
    {
        gc()
        events <- event.info[[type]]
        rc.results.cell <- list()
        rc.psi.results.cell <- list()
        psi.results.cell <- list()
        for (data_type in cell_similarity_data)
        {
            cell_similarity_type <- cell_similarity[[data_type]]
            dyk_cell_type <- dyk_cell[[data_type]]
            # rc impute
            token <- paste(Sys.getpid(), rbinom(1, size = 1000000000, prob = 0.5), sep = "-")
            datapath <- paste0(output_path, "/imputation1_data_", token, ".h5")
            resultpath <- paste0(output_path, "/imputation1_result_", token, ".mat")
            cell_similarity_type <- as.matrix(cell_similarity_type)
            data = rc[do.call(what = c, args = events[, -1]), ]
            data <- as.matrix(data)
            all_data <- list()
            for (position in colnames(events[, -1]))
            {
                tmp <- list(data[events[[position]], , drop = F])
                all_data <- c(all_data, tmp)
            }
            names(all_data) <- colnames(events[, -1])
            msg <- paste0("[", Sys.time(), "] ", "Save data")
            print(msg)
            saveHdf5File(datapath, list(
                similar = cell_similarity_type, data = all_data,
                parameter = list(decay = decay_impute), data_type = "RC",
                similar_type = "cell"
            ))
            msg <- paste0("[", Sys.time(), "] ", "Save data Finished")
            print(msg)
            cmd <- paste("bash", mat_impute_v1, mcr_path, datapath,
                         resultpath, ">>", log_file, "2>&1")
            print(cmd)
            system(cmd, wait = T)
            rc_imputed_res_cell <- read_mat(resultpath)
            rc_imputed_res_cell[[1]] <- NULL
            names(rc_imputed_res_cell) <- names(all_data)
            for (position in names(rc_imputed_res_cell))
            {
                tmp <- data[events[[position]], , drop = F]
                if (!is.matrix(rc_imputed_res_cell[[position]])) {
                    rc_imputed_res_cell[[position]] <- matrix(rc_imputed_res_cell[[position]],
                        ncol = length(rc_imputed_res_cell[[position]]), nrow = 1
                    )
                }
                rownames(rc_imputed_res_cell[[position]]) <- rownames(tmp)
                colnames(rc_imputed_res_cell[[position]]) <- colnames(tmp)
            }
            file.remove(datapath, resultpath)
            rc_psi_imputed_res_cell <- rc_to_psi(
                data.frame(events, type = type),
                rc_imputed_res_cell
            )
            # psi impute
            data = psi[events$event, , drop = F]
            data <- as.matrix(data)
            token <- paste(Sys.getpid(), rbinom(1, size = 1000000000, prob = 0.5), sep = "-")
            datapath <- paste0(output_path, "/imputation1_data_", token, ".h5")
            resultpath <- paste0(output_path, "/imputation1_result_", token, ".mat")
            all_data <- list(data)
            names(all_data) <- "all"
            saveHdf5File(datapath, list(
                similar = cell_similarity_type,
                data = all_data, parameter = list(decay = decay_impute),
                data_type = "PSI", similar_type = "cell"
            ))

            cmd <- paste("bash", mat_impute_v1, mcr_path, datapath,
                         resultpath, ">>", log_file, "2>&1")
            system(cmd, wait = T)
            psi_imputed_res_cell <- read_mat(resultpath)
            psi_imputed_res_cell[[1]] <- NULL
            names(psi_imputed_res_cell) <- names(all_data)
            if (!is.matrix(psi_imputed_res_cell[[1]])) {
                psi_imputed_res_cell[[1]] <- matrix(psi_imputed_res_cell[[1]],
                    ncol = length(psi_imputed_res_cell[[1]]), nrow = 1
                )
            }
            rownames(psi_imputed_res_cell[[1]]) <- rownames(data)
            colnames(psi_imputed_res_cell[[1]]) <- colnames(data)
            psi_imputed_res_cell <- do.call(what = rbind, args = psi_imputed_res_cell)
            file.remove(datapath, resultpath)
            psi.results.cell[[data_type]] <- psi_imputed_res_cell
            rc.psi.results.cell[[data_type]] <- rc_psi_imputed_res_cell
            rc.results.cell[[data_type]] <- rc_imputed_res_cell
        }
        psi.imputation.cell[[type]] <- psi.results.cell
        rc.imputation.cell[[type]] <- rc.results.cell
        rc.psi.imputation.cell[[type]] <- rc.psi.results.cell
    }
    res_list <- list()
    for (data_type in cell_similarity_data) {
        PSI_imputed_res <- do.call(what = rbind, args = lapply(psi.imputation.cell, FUN = function(x) {
            return(x[[data_type]])
        }))
        RC_imputed_res <- lapply(rc.imputation.cell, FUN = function(x) {
            return(x[[data_type]])
        })
        RC_PSI_imputed_res <- do.call(what = rbind, args = lapply(rc.psi.imputation.cell, FUN = function(x) {
            return(x[[data_type]])
        }))
        res_list[[paste0(data_type, "_PSI")]] <- PSI_imputed_res
        res_list[[paste0(data_type, "_RC")]] <- list(rc = RC_imputed_res, psi = RC_PSI_imputed_res)
    }
    for (name in names(res_list)) {
        if (unlist(strsplit(name, "_"))[length(unlist(strsplit(name, "_")))] == "RC") {
            data_imputated <- res_list[[name]]$psi
            rc_list <- res_list[[name]]$rc
            event_types <- unique(events$type)
            rc.mean <- lapply(X = as.list(event_types), FUN = function(type) {
                cur.event <- events[events$type == type, ]
                cur_rc_list <- rc_list[[type]]
                tmp.array <- array(0, dim = c(nrow(cur_rc_list[[1]]), ncol(cur_rc_list[[1]]), length(cur_rc_list)))
                for (i in 1:length(cur_rc_list)) {
                    tmp.array[, , i] <- as.matrix(cur_rc_list[[i]])
                }
                tmp.array.mean <- rowMeans(tmp.array, dims = 2, na.rm = T)
                rownames(tmp.array.mean) <- cur.event[, 1]
                return(tmp.array.mean)
            })
            rc.mean <- do.call(what = rbind, args = rc.mean)
            rc.mean <- rc.mean[events[, 1], ]
            rc_v <- unlist(as.vector(rc.mean))
            min_reads <- quantile(rc_v[which(rc_v != 0)], probs = seq(0, 1, 0.1))[3]
            rm(rc_v)
            data_imputated[which(rc.mean < min_reads)] <- 0
            res_list[[name]] <- data_imputated
        }
    }
    output_name = paste0(
        output_path, "/Imputed_fine_tune_",
        rbinom(1, size = 1000000000, prob = 0.5),
        ".rds"
    )
    saveRDS(res_list, output_name)

    ## feature
    model_list <- list()
    for (type in cell_similarity_data) {
        msg <- paste0("[", Sys.time(), "] ", " Model training", paste0(";similarity_type=", type))
        print(msg)
        cell_similarity_type <- cell_similarity[[type]]
        dyk_cell_type <- dyk_cell[[type]]
        rc_imputed <- res_list[[paste0(type, "_RC")]]
        psi_imputed <- res_list[[paste0(type, "_PSI")]]
        feature_df <- getClassifierFeature(
            cell_similarity_type, psi, rc, expr,
            event, dyk_cell_type, rc_imputed, psi_imputed
        )
        feature_df <- as.data.frame(feature_df)
        psi_gtex_tmp <- psi_gtex[row.names(psi)]

        idx_change_other <- which(abs(as.matrix(psi) - psi_gtex_tmp) >= 0.05 & abs(rc_imputed - psi_gtex_tmp) < 0.9)
        idx_change_big <- which(abs(as.matrix(psi) - psi_gtex_tmp) >= 0.05 & abs(rc_imputed - psi_gtex_tmp) >= 0.9)
        idx_nonchange <- which(abs(as.matrix(psi) - psi_gtex_tmp) < 0.05)
        df <- rbind(
            cbind(feature_df[idx_change_other, ], group = "change_other"),
            cbind(feature_df[idx_change_big, ], group = "change_big"),
            cbind(feature_df[idx_nonchange, ], group = "nonchange")
        )
        saveRDS(df, paste0(output_path, "/", type, "_fine_tune_feature_df", ".rds"))
        X <- df[which(df$psi_self_type.2 == 0), ]
        X <- X[, -which(colnames(X) %in% c("deltapsi_impute_train_rc", "deltapsi_impute_train_psi"))]
        Y <- X[, "group"]
        Y <- ifelse(grepl("nonchange", Y), 1, 0)
        Y <- factor(Y, levels = c(0, 1), labels = c("idx_change", "idx_nonchange"))
        X <- X[, -which(colnames(X) %in% c("group"))]
        X <- as.matrix(X)
        train_ds <- downSample(x = X, y = Y)
        start_coef1 <- as.matrix(coef(model1[[type]], s = "lambda.min"))
        start_x_name <- rownames(start_coef1)[which(start_coef1[, 1] != 0)]
        start_x_name <- start_x_name[-1]
        train_ds <- train_ds[, c(start_x_name, "Class")]
        train_ds_X <- train_ds[, -ncol(train_ds)]
        train_ds_Y <- train_ds[, ncol(train_ds)]
        cv_fit1 <- cv.glmnet(as.matrix(train_ds_X), as.matrix(train_ds_Y),
            family = "binomial",
            type.measure = "class", standardize = T, alpha = 0,
            start_coef = coef(model1[[type]], s = "lambda.min"), intercept = T
        )
        X <- df[which(df$psi_self_type.2 == 0), ]
        X <- X[which(grepl("change_", X$group)), ]
        Y <- X[, "group"]
        Y <- ifelse(grepl("change_big", Y), 0, 1)
        Y <- factor(Y, levels = c(0, 1), labels = c("change_big", "change_other"))
        X <- X[, -which(colnames(X) %in% c("group", "dataset"))]
        X <- as.matrix(X)
        train_ds <- downSample(x = X, y = Y)
        start_coef2 <- as.matrix(coef(model2[[type]], s = "lambda.min"))
        start_x_name <- rownames(start_coef2)[which(start_coef2[, 1] != 0)]
        start_x_name <- start_x_name[-1]
        train_ds <- train_ds[, c(start_x_name, "Class")]
        train_ds_X <- train_ds[, -ncol(train_ds)]
        train_ds_Y <- train_ds[, ncol(train_ds)]
        cv_fit2 <- cv.glmnet(as.matrix(train_ds_X), as.matrix(train_ds_Y),
            family = "binomial",
            type.measure = "class", standardize = T, alpha = 0,
            start_coef = coef(model2[[type]], s = "lambda.min"), intercept = T
        )
        model_list[[type]] <- list(model1 = cv_fit1, model2 = cv_fit2)
    }
    saveRDS(model_list, paste0(output_path, "/model_ft.rds"))
    msg = paste0("[", Sys.time(), "] ", "Classifer fine tune Finish.")
    print(msg)
    return(output_path)
}


#' @title Predict the probabilities of WD, BD, and TD-Info for each event-cell pair.
#' first: probabilities of being WD but not ND for each event-cell pair;
#' second: probabilities of being BD but not TD for each event-cell pair;
#' third: probabilities of being TD+Info but not TD-Info for each event-cell pair.
#'
#' @param feature_df a data.frame of features for classifer
#' @param model1 logistical model to distinguish BD and TD
#' @param model2 logistical model to distinguish TD+Info and TD-Info
#' @param psi matrix of psi value before imputation
#'
#' @return a list of probability matrix
#'
#' @export
#' @importFrom stats predict
#'
getPredictProb <- function(feature_df, model1, model2, psi) {
    num_r = nrow(psi)
    num_c = ncol(psi)
    start_coef1 <- as.matrix(coef(model1, s = "lambda.min"))
    start_x_name1 <- rownames(start_coef1)[-1]
    start_coef2 <- as.matrix(coef(model2, s = "lambda.min"))
    start_x_name2 <- rownames(start_coef2)[-1]
    cutoff <- 0.5
    idx_01 <- which(psi == 0 | psi == 1)
    prob_m1 <- matrix(0, ncol = num_c, nrow = num_r)
    prob_m1[idx_01] <- 1

    prob <- c()
    block_size <- 1000000
    rows_feature <- nrow(feature_df)
    num_block <- ceiling(rows_feature / block_size)
    for (i in 1:num_block) {
        start_row <- (i - 1) * block_size + 1
        end_row <- min(i * block_size, rows_feature)
        current_block <- feature_df[start_row:end_row, ]
        current_prob <- predict(model1, as.matrix(current_block[, start_x_name1]), type = "response")
        prob <- c(prob, current_prob)
    }
    prob[which(prob < 0.01)] <- 0
    prob[which(prob > 0.99)] <- 1
    prob_m2 <- matrix(prob, byrow = F, ncol = num_c, nrow = num_r)
    feature_df_2 <- feature_df[which(prob_m2 < cutoff), ]
    prob2 <- c()
    block_size <- 1000000
    rows_feature <- nrow(feature_df_2)
    num_block <- ceiling(rows_feature / block_size)
    for (i in 1:num_block) {
        start_row <- (i - 1) * block_size + 1
        end_row <- min(i * block_size, rows_feature)
        current_block <- feature_df_2[start_row:end_row, ]
        current_prob <- predict(model2, as.matrix(current_block[, start_x_name2]), type = "response")
        prob2 <- c(prob2, current_prob)
    }

    prob2[which(prob2 < 0.01)] <- 0
    prob2[which(prob2 > 0.99)] <- 1
    prob_m3 <- matrix(1, byrow = F, ncol = num_c, nrow = num_r)
    prob_m3[which(prob < cutoff)] <- prob2
    return(list(prob_m1, prob_m2, prob_m3))
}

#' @title Combine Multiple Imputation Strategies Using Probabilistic Weighting
#'
#' @description This function combines results from different imputation strategies
#'   using a probabilistic weighting scheme based on event-cell pair classifications.
#'   Strategy 1 is applied to pairs in ND, Strategy 2 to pairs in BD and TD+Info,
#'   and Strategy 3 to pairs in TD-Info. The final imputation combines these strategies
#'   weighted by classification probabilities.
#'
#' @param paras A list object containing SCSES configuration parameters, typically
#'   loaded using \code{readSCSESconfig(paras_file)}.
#' @param rds_imputed_file Character string specifying the path to RDS file containing
#'   imputation results from multiple strategies. This should be the output file path
#'   returned by the \code{\link{ImputationAll}} function.
#' @param output_path Character string specifying the output directory for final
#'   combined imputation results. Default: \code{work_path/imputation/}.
#' @param cell_similarity Either a list of cell similarity matrices named by data
#'   types (EXP_RBP, RC, PSI) or a character string path to an RDS file.
#'   Default: \code{imputation/cell_similarity/cell.similars.rds}.
#' @param dyk_cell Either a list of neighbor numbers named by data types
#'   (EXP_RBP, RC, PSI) or a character string path to an RDS file.
#'   Default: \code{imputation/cell_similarity/dyk.cell.rds}.
#' @param rc Either a matrix of normalized splicing event read counts or a
#'   character string path to an RDS file. Default: \code{rds_processed/rc.rds}.
#' @param psi Either a matrix of PSI values or a character string path to an
#'   RDS file. Default: \code{rds_processed/psi.rds}.
#' @param expr Either a matrix of gene expression (TPM for full-length protocols,
#'   normalized UMI counts for droplet-based protocols) or a character string
#'   path to an RDS file. Default: \code{rds_processed/expr.rds}.
#' @param event Either a data.frame of event information or a character string
#'   path to an RDS file. Default: \code{rds_processed/event.rds}.
#'
#' @return Character string specifying the path to the RDS file containing
#'   the final combined imputation results. The returned list is named by
#'   similarity types, with each element containing the combined PSI matrix
#'   for that similarity approach.
#'
#' @export
#' @importFrom stats predict
#'
Estimation <- function(
    paras, rds_imputed_file, output_path = NULL, cell_similarity = NULL,
    dyk_cell = NULL, rc = NULL, psi = NULL, expr = NULL, event = NULL) {
    msg <- paste0("[", Sys.time(), "] ", "Combine imputed psi.")
    print(msg)
    # validate input----
    msg <- paste0("[", Sys.time(), "] ", "Loading data...")
    print(msg)
    print(paste0("Input: ", rds_imputed_file))
    psi_imputed_seperated <- readRDS(file = rds_imputed_file)
    print(paste("Checking data..."))
    rc <- check.readRDS(rc, default = paste0(paras$Basic$work_path, "/rds_processed/rc.rds"))
    print("rc checked")
    psi <- check.readRDS(psi, default = paste0(paras$Basic$work_path, "/rds_processed/psi.rds"))
    print("psi checked")
    expr <- check.readRDS(expr, default = paste0(paras$Basic$work_path, "/rds_processed/expr.rds"))
    print("expr checked")
    event <- check.readRDS2(event, default = paste0(paras$Basic$work_path, "/rds_processed/event.rds"))
    print("event checked")
    cell_similarity <- check.readRDS2(cell_similarity, default = paste0(paras$Basic$work_path, "/imputation/cell_similarity/cell.similars.rds"))
    print("cell similarity checked")
    dyk_cell <- check.readRDS2(dyk_cell, default = paste0(paras$Basic$work_path, "/imputation/cell_similarity/dyk.cell.rds"))
    print("dynamic cell knn checked")
    if(file.exists(paste0(paras$Basic$work_path, "/classifer/model_ft.rds"))){
      print("Fine tune model will be used.")
      model = readRDS(paste0(paras$Basic$work_path, "/classifer/model_ft.rds"))
    }else{
      print("Pre-trained model will be used.")
      extdata_path <- system.file("extdata", package = "SCSES")
      model1_path = paste0(extdata_path, "/model/model_change_nonchange.rdata")
      model2_path = paste0(extdata_path, "/model/model_change_01_change_other.rdata")
      model1 <- get(load(model1_path))
      model2 <- get(load(model2_path))
      model <- lapply(names(model1), function(x){
        return(list(model1=model1[[x]],
                    model2=model2[[x]]))
      })
      names(model) <- names(model1)
    }
    print("classifer checked")
    if (is.null(output_path)) {
        output_path <- paste0(paras$Basic$work_path, "/imputation/")
    }
    if (!dir.exists(output_path)) {
        dir.create(output_path)
    }
    print(paste0("Output: ", output_path))

    # validate cell similarity type----
    print(paste("Checking cell similarity type"))
    cell_similarity_data1 <- names(cell_similarity)
    cell_similarity_data1 <- check.valid(x = cell_similarity_data1, select = c("EXP_RBP", "RC", "PSI"))
    cell_similarity_data2 <- names(dyk_cell)
    cell_similarity_data2 <- check.valid(x = cell_similarity_data2, select = c("EXP_RBP", "RC", "PSI"))
    cell_similarity_data <- intersect(cell_similarity_data1, cell_similarity_data2)
    print(paste0("cell_similarity_data=", paste(cell_similarity_data, collapse = ";"), "  checked"))

    psi_imputed_final <- list()
    for (type in cell_similarity_data) {
        gc()
        cell_similar_type <- cell_similarity[[type]]
        dyk_cell_type <- dyk_cell[[type]]
        rc_imputed <- psi_imputed_seperated[["cell"]][[paste0(type, "_RC")]]
        psi_imputed <- psi_imputed_seperated[["cell"]][[paste0(type, "_PSI")]]
        feature_df <- getClassifierFeature(
            cell_similar_type, psi, rc, expr,
            event, dyk_cell_type, rc_imputed, psi_imputed
        )
        feature_df <- as.data.frame(feature_df)
        save(feature_df, file = paste0(output_path, "/", type, "_feature_df", ".rdata"))
        prob_m <- getPredictProb(feature_df,
            model1 = model[[type]][["model1"]],
            model2 = model[[type]][["model2"]],
            psi = psi
        )
        impute_rc <- psi_imputed_seperated[["cell"]][[paste0(type, "_RC")]]
        impute_psi <- psi_imputed_seperated[["cell"]][[paste0(type, "_PSI")]]
        impute_psi_v2 <- psi_imputed_seperated[["cell_event"]][[paste0(type, "_PSI")]]
        impute_res <- impute_psi * (1 - prob_m[[1]]) + prob_m[[1]] * prob_m[[2]] * impute_rc + prob_m[[1]] * (1 - prob_m[[2]]) * prob_m[[3]] * impute_rc + prob_m[[1]] * (1 - prob_m[[2]]) * (1 - prob_m[[3]]) * impute_psi_v2
        psi_imputed_final[[type]] <- impute_res
    }
    output_name <- paste0(
        output_path, "/Imputed_combined_",
        rbinom(1, size = 1000000000, prob = 0.5),
        ".rds"
    )
    saveRDS(psi_imputed_final, output_name)
    msg <- paste0("[", Sys.time(), "] ", "Combine imputed psi Finish.")
    print(msg)
    return(output_name)
}

