#' @title Merge splicing events for classifer fine tune
#' @description Merge different types of splicing events in
#' work_path/splicing_value_ft, and save to work_path/rds_ft
#' @param paras A list object containing SCSES configuration parameters, typically
#'   loaded using \code{readSCSESconfig(paras_file)}.
#' @return Merged splicing value path
#'
#' @export
#'
mergeFtSplicingValue <- function(paras) {
  psi_path <- paste0(paras$Basic$work_path, "/splicing_value_ft/")
  rc_path <- paste0(paras$Basic$work_path, "/splicing_value_ft/")
  event_path <- paste0(paras$Basic$work_path, "/splicing_value_ft/")
  output_path <- paste0(paras$Basic$work_path, "/rds_ft/")
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
