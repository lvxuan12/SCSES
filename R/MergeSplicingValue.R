#' @title Merge splicing events
#' @description Merge different types of splicing events
#' in work_path/splicing_value, and save merged splicing value to work_path/rds
#' @param paras A list object containing SCSES configuration parameters, typically
#'   loaded using \code{readSCSESconfig(paras_file)}.
#'
#' @return Merged splicing value path
#'
#' @export
#'
mergeSplicingValue <- function(paras){
  psi_path <- paste0(paras$Basic$work_path, "/splicing_value/")
  rc_path <- paste0(paras$Basic$work_path, "/splicing_value/")
  event_path <- paste0(paras$Basic$work_path, "/splicing_value/")
  output_path <- paste0(paras$Basic$work_path, "/rds/")
  if (!dir.exists(output_path)) {
    dir.create(output_path)
  }
  psi_rds <- list.files(psi_path, pattern = "_psi.rds", full.names = T)
  psi_all <- data.frame()
  for (f in psi_rds)
  {
    psi <- readRDS(file = f)
    psi_all <- rbind(psi_all, psi)
  }
  psi_all[is.na(psi_all)] <- 0
  saveRDS(psi_all, paste0(output_path, "/psi.rds"))
  rc_rds <- list.files(rc_path, pattern = "_rc.rds", full.names = T)
  rc_all <- data.frame()
  for (f in rc_rds)
  {
    data <- readRDS(file = f)
    name <- names(data)
    name <- name[grep("^rc", name)]
    for (x in name) {
      tmp <- data[[x]]
      tmp <- tmp[!duplicated(rownames(tmp)), , drop = F]
      tmp <- tmp[which(!row.names(tmp) %in% row.names(rc_all)), , drop = F]
      rc_all <- rbind(rc_all, tmp)
    }
  }
  rc_all[is.na(rc_all)] <- 0
  saveRDS(rc_all, paste0(output_path, "/rc.rds"))
  event_rds <- list.files(event_path, pattern = "_rc.rds", full.names = T)
  event_all <- data.frame()
  for (f in event_rds)
  {
    data <- readRDS(file = f)
    events <- data$events
    type <- gsub("_rc.rds", "", unlist(strsplit(f, "/"))[length(unlist(strsplit(f, "/")))])
    events$type <- type
    event_all <- rbind(event_all, events)
  }
  saveRDS(event_all, paste0(output_path, "/event.rds"))
  return(output_path)
}


#' @title Merge PSI
#' @description Merge PSI for different types of splicing events
#' in splicing_value_path, and save merged PSI to output_path
#' @param splicing_value_path The path to *_psi.rds
#' * represents different types of splicing events
#' @param output_path The path to save merged PSI

#' @return Merged PSI path
#'
#' @export
#'

mgergePSI <- function(splicing_value_path, output_path) {
  if (!dir.exists(output_path)) {
    dir.create(output_path)
  }
  psi_rds <- list.files(splicing_value_path, pattern = "_psi.rds", full.names = T)
  psi_all <- data.frame()
  for (f in psi_rds)
  {
    psi <- readRDS(file = f)
    psi_all <- rbind(psi_all, psi)
  }
  psi_all[is.na(psi_all)] <- 0
  saveRDS(psi_all, paste0(output_path, "/psi.rds"))
  return(output_path)
}

#' @title Merge Read counts
#' @description Merge Read counts for different types of splicing events
#' in splicing_value_path, and save merged Read counts to output_path
#' @param splicing_value_path The path to *_rc.rds
#' * represents different types of splicing events
#' @param output_path The path to save merged Read counts

#' @return Merged Read counts path
#'
#' @export
#'

mgergeRC <- function(splicing_value_path, output_path) {
  if (!dir.exists(output_path)) {
    dir.create(output_path)
  }
  rc_rds <- list.files(splicing_value_path, pattern = "_rc.rds", full.names = T)
  rc_all <- data.frame()
  for (f in rc_rds)
  {
    data <- readRDS(file = f)
    name <- names(data)
    name <- name[grep("^rc", name)]
    for (x in name) {
      tmp <- data[[x]]
      tmp <- tmp[!duplicated(rownames(tmp)), , drop = F]
      tmp <- tmp[which(!row.names(tmp) %in% row.names(rc_all)), , drop = F]
      rc_all <- rbind(rc_all, tmp)
    }
  }
  rc_all[is.na(rc_all)] <- 0
  saveRDS(rc_all, paste0(output_path, "/rc.rds"))
  return(output_path)
}

#' @title Merge event information
#' @description Merge event information for different types of splicing events
#' in splicing_value_path, and save merged event information to output_path
#' @param splicing_value_path The path to *_rc.rds
#' * represents different types of splicing events
#' @param output_path The path to save merged event information

#' @return Merged event information path
#'
#' @export
#'

mgergeEvent <- function(splicing_value_path, output_path) {
  if (!dir.exists(output_path)) {
    dir.create(output_path)
  }
  event_rds <- list.files(splicing_value_path, pattern = "_rc.rds", full.names = T)
  event_all <- data.frame()
  for (f in event_rds)
  {
    data <- readRDS(file = f)
    events <- data$events
    type <- gsub("_rc.rds", "", unlist(strsplit(f, "/"))[length(unlist(strsplit(f, "/")))])
    events$type <- type
    event_all <- rbind(event_all, events)
  }
  saveRDS(event_all, paste0(output_path, "/event.rds"))
  return(output_path)
}
