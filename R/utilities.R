######
# Parameter validation
######
check.int.or.null <- function(x, default) {
    if (is.numeric(x = x)) {
        x <- as.integer(x = x)
    } else if (is.integer(default)){
        if (is.null(x = x)) {
            x <- default
        } else if (x == "") {
            x <- default
        } else if (is.na(x = x)) {
            x <- default
        } else {
            stop(paste("Parameter must be a integer"))
        }
    } else {
        stop(paste("Parameter must be a integer"))
    }
    return(x)
}

check.double.or.null <- function(x, default) {
  if (is.double(x = x)) {
    x <- as.numeric(x = x)
  } else if (is.double(default)){
    if (is.null(x = x)) {
      x <- default
    } else if (x == "") {
      x <- default
    } else if (is.na(x = x)) {
      x <- default
    } else {
      stop(paste("Parameter must be a double"))
    }
  } else {
    stop(paste("Parameter must be a double"))
  }
  return(x)
}

check.path.or.matrix <- function(x) {
    if (is.matrix(x) | is.data.frame(x) | inherits(x, "sparseMatrix")) {
        return("variable")
    } else if (length(x) <= 1){
        if (is.null(x = x)) {
            return("")
        } else if ( x == "") {
            return("")
        }else if (is.na(x = x)) {
            return("")
        }else if (is.character(x)) {
            if(file.exists(x)){
                return("path")
            }else {
                stop(paste("Parameter read error"))
            }
        } else {
            stop(paste("Parameter read error"))
        }
    }else {
        stop(paste("Parameter must be a matrix or dataframe"))
    }
}
check.path.or.list <- function(x) {
    if (is.list(x)) {
        return("variable")
    } else if (length(x) <= 1) {
        if (is.null(x = x)) {
            return("")
        } else if (x == "") {
            return("")
        } else if (is.na(x = x)) {
            return("")
        } else if (is.character(x)) {
            if (file.exists(x)) {
                return("path")
            } else {
                stop(paste("Parameter read error"))
            }
        } else {
            stop(paste("Parameter read error"))
        }
    } else {
        stop(paste("Parameter must be a list"))
    }
}
check.readRDS2 <- function(x, default = NULL) {
    flag <- check.path.or.list(x)
    if (flag == "path") {
        x <- readRDS(file = x)
    } else if (flag == "") {
        if (file.exists(default)) {
            x <- readRDS(file = default)
        } else {
            stop(paste("Default path read error"))
        }
    } else if (flag == "variable") {
        x <- x
    } else {
        stop(paste("Parameter read error"))
    }
    if (is.list(x)) {
        return(x)
    } else {
        stop(paste("Parameter must be a list"))
    }
}
check.readRDS <- function(x,default=NULL){
    flag = check.path.or.matrix(x)
    if(flag == "path"){
        x <- readRDS(file = x)
    } else if(flag==""){
        if(file.exists(default)){
            x <- readRDS(file = default)
        }else{
            stop(paste("Default path read error"))
        }
    } else if (flag == "variable") {
        x <- x
    } else {
        stop(paste("Parameter read error"))
    }
    if (is.matrix(x) | is.data.frame(x) | inherits(x, "sparseMatrix")) {
        return(x)
    }else {
        stop(paste("Parameter must be a matrix or dataframe"))
    }
}

check.valid <- function(x, select) {
    invalid_x <- setdiff(x, select)
    valid_x <- intersect(x, select)
    if (length(valid_x) > 0 & length(invalid_x) > 0) {
        x <- valid_x
        print(paste(c("Invalid:", invalid_x), collapse = " "))
    } else if (length(valid_x) > 0 & length(invalid_x) == 0) {
        x <- valid_x
    } else {
        stop(paste(c("Invalid:", invalid_x), collapse = " "))
    }
    x
}

check.path.or.character <- function(x) {
    if (is.character(x)) {
        if (length(x) == 1) {
            if (file.exists(x)) {
                return("path")
            } else {
                return("character")
            }
        } else {
            return("character")
        }
    } else {
        stop("Parameter must be a character string")
    }
}

#' @title Save R variable to hdf5 file
#' @param con path to hdf5 file
#' @param dataset can be character, numeric, matrix, data.frame, list
#' @param con.type connection type
#' @return NULL
#'
#' @keywords internal
#' @export
#'
saveHdf5File <- function(con, dataset, con.type = "F") {
    if (class(con) == "character") {
        con <- H5Fcreate(con)
    } else if (class(con) != "H5IdComponent") {
        stop("Invalid conncetion variable!")
    }
    obj.names <- names(dataset)
    for (obj in obj.names)
    {
        if (class(dataset[[obj]])[1] == "character") {
            space <- H5Screate_simple(dims = 1)
            tid <- H5Tcopy(dtype_id = "H5T_C_S1")
            H5Tset_size(tid, size = nchar(dataset[[obj]]))
            # space=H5Screate("H5S_SCALAR")
            ds <- H5Dcreate(con, obj, dtype_id = tid, h5space = space)
            H5Dwrite(h5dataset = ds, buf = dataset[[obj]])
            H5Dclose(ds)
            H5Sclose(space)
        }
        if (class(dataset[[obj]])[1] == "numeric") {
            space <- H5Screate_simple(dims = c(length(dataset[[obj]]), 1))
            ds <- H5Dcreate(con, obj, dtype_id = "H5T_IEEE_F64BE", h5space = space)
            H5Dwrite(h5dataset = ds, buf = dataset[[obj]])
            H5Dclose(ds)
            H5Sclose(space)
        }
        if (class(dataset[[obj]])[1] %in% c("matrix", "dgCMatrix")) {
            space <- H5Screate_simple(dims = dim(dataset[[obj]]))
            ds <- H5Dcreate(con, obj, dtype_id = "H5T_IEEE_F64BE", h5space = space)
            H5Dwrite(h5dataset = ds, buf = dataset[[obj]])
            H5Dclose(ds)
            H5Sclose(space)
        }
        if (class(dataset[[obj]])[1] == "data.frame") {
            space <- H5Screate("H5S_SIMPLE")
            ds <- H5Dcreate(con, obj, dtype_id = "H5T_IEEE_F64BE", h5space = space)
            H5Dwrite(h5dataset = ds, buf = as.matrix(dataset[[obj]]))
            H5Dclose(ds)
            H5Sclose(space)
        }
        if (class(dataset[[obj]])[1] == "list") {
            group <- H5Gcreate(con, obj)
            saveHdf5File(group, dataset[[obj]], "G")
        }
    }
    get(paste0("H5", con.type, "close"))(con)
}

#' @title Calculate the mean of exclusion read counts
#' @description exclusion reads recorded in event
#' @param event data frame of event information
#' @param rc matrix, reads associated with splicing events

#' @return matrix of the mean of exclusion read counts
#'
#' @keywords internal
#' @export

calcu_ex_rc <- function(rc, event) {
    rc <- as.matrix(rc)
    tmp.array <- array(0, dim = c(nrow(event), ncol(rc), 2))
    for (i in seq(2, 3))
    {
        tmp <- rc[match(event[, i], row.names(rc)), , drop = F]
        tmp.array[, , i - 1] <- as.matrix(tmp)
    }
    tmp.mean <- rowMeans(tmp.array, dims = 2, na.rm = T)
    row.names(tmp.mean) <- event$event
    colnames(tmp.mean) <- colnames(rc)
    return(tmp.mean)
}

#' @title Calculate the mean of inclusion read counts
#' @description inclusion reads recorded in event
#' @param event data frame of event information
#' @param rc matrix, reads associated with splicing events

#' @return matrix of the mean of inclusion read counts
#'
#' @keywords internal
#' @export

calcu_in_rc <- function(rc, event) {
    rc <- as.matrix(rc)
    tmp.array <- array(0, dim = c(nrow(event), ncol(rc), 2))
    for (i in seq(4, 5))
    {
        tmp <- rc[match(event[, i], row.names(rc)), ]
        tmp.array[, , i - 3] <- as.matrix(tmp)
    }
    tmp.mean <- rowMeans(tmp.array, dims = 2, na.rm = T)
    row.names(tmp.mean) <- event$event
    colnames(tmp.mean) <- colnames(rc)
    return(tmp.mean)
}

#' @title Calculate the mean of read counts
#' @description associated reads name recorded in event
#' @param event data frame of event information
#' @param rc matrix, reads associated with splicing events

#' @return matrix of the mean of read counts
#'
#' @keywords internal
#' @export

calcu_rc <- function(rc, event) {
    tmp.array <- array(0, dim = c(nrow(event), ncol(rc), 4))
    for (i in seq(2, 5))
    {
        tmp <- rc[match(event[, i], row.names(rc)), ]
        tmp.array[, , i - 1] <- as.matrix(tmp)
    }
    tmp.mean <- rowMeans(tmp.array, dims = 2, na.rm = T)
    return(tmp.mean)
}

#' @title Calculate the coefficient of variation for psi of neighbor cells
#' @description associated reads name recorded in event
#' @param psi  matrix, psi associated with splicing events
#' @param cs matrix, cell similarity
#'
#' @return matrix
#'
#' @keywords internal
#' @export
knn_cv <- function(psi, cs) {
    psi <- as.matrix(psi)
    cs <- as.matrix(cs)
    cs[which(cs != 0)] <- 1
    sq_psi_e <- (psi**2) %*% t(cs)
    n_k <- rowSums(cs)
    sq_psi_e <- t(t(sq_psi_e) / n_k)
    e_psi_sq <- psi %*% t(cs)
    m <- t(t(e_psi_sq) / n_k)
    e_psi_sq <- (t(t(e_psi_sq) / n_k))^2
    v <- sq_psi_e - e_psi_sq
    return(v / m)
}

#' @title Calculate psi
#' @param rc matrix, reads associated with splicing events
#' @param event data frame of event information,
#' event types: A3SS, A5SS, AL, SE, MXE, RI
#'
#' @return data.frame of psi value
#'
#' @keywords internal
#' @export
#'
rc_to_psi <- function(event, rc) {
    event_type <- unique(event$type)
    psi_all <- data.frame()
    for (t in event_type) {
        df <- subset(event, type == t)
        if (t == "RI" | t == "SE") {
            iso1 <- rc[["exclusion1"]]
            rc_retention1 <- rc[["retention1"]]
            rc_retention2 <- rc[["retention2"]]
            rc_retention <- array(0, dim = c(nrow(rc_retention1), ncol(rc_retention1), 2))
            rc_retention[, , 1] <- as.matrix(rc_retention1)
            rc_retention[, , 2] <- as.matrix(rc_retention2)
            iso2 <- rowMeans(rc_retention, dims = 2)
            psi <- iso2 / (iso1 + iso2)
            psi <- as.data.frame(psi)
            row.names(psi) <- df$event
            psi[is.na(psi)] <- 0
        } else if (t == "MXE") {
            rc_exclusion1 <- rc[["exclusion1"]]
            rc_exclusion2 <- rc[["exclusion2"]]
            rc_retention1 <- rc[["retention1"]]
            rc_retention2 <- rc[["retention2"]]
            rc_exclusion <- array(0, dim = c(nrow(rc_exclusion1), ncol(rc_exclusion1), 2))
            rc_exclusion[, , 1] <- as.matrix(rc_exclusion1)
            rc_exclusion[, , 2] <- as.matrix(rc_exclusion2)
            iso1 <- rowMeans(rc_exclusion, dims = 2)
            rc_retention <- array(0, dim = c(nrow(rc_retention1), ncol(rc_retention1), 2))
            rc_retention[, , 1] <- as.matrix(rc_retention1)
            rc_retention[, , 2] <- as.matrix(rc_retention2)
            iso2 <- rowMeans(rc_retention, dims = 2)
            psi <- iso2 / (iso1 + iso2)
            psi <- as.data.frame(psi)
            row.names(psi) <- df$event
            psi[is.na(psi)] <- 0
        } else if (t == "A3SS" | t == "A5SS" | t == "AL") {
            iso1 <- rc[["exclusion1"]]
            iso2 <- rc[["retention1"]]
            psi <- iso2 / (iso1 + iso2)
            psi <- as.data.frame(psi)
            row.names(psi) <- df$event
            psi[is.na(psi)] <- 0
        }
        psi_all <- rbind(psi_all, psi)
    }
    colnames(psi_all) <- colnames(rc[[1]])
    return(as.matrix(psi_all))
}
