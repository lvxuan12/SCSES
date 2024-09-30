#' @title Create Pseudobulk bam by merging single cell bam files and index
#'
#' @param paras list fromJSON(paras_file)
#' Default bam_path, samtools_path, core from paras
#' @param bam_path directory to single cell bam file
#' @param core the number of threads
#'
#' @return Pseudobulk path
#' @export
#'
createPseudobulk <- function(
    paras, bam_path = paras$Basic$bam_path, core = paras$Basic$core) {
  #script
  samtools_path = paras$Basic$samtools_path
  # input
  print(paste0("Input: ", bam_path))
  # output
  pseudobulk.path <- paste0(paras$Basic$work_path, "/data/")
  print(paste0("Output: ", pseudobulk.path))
  msg <- paste0("[", Sys.time(), "] ", "Creating Pseudobulk directory...")
  print(msg)
  dir.create(path = pseudobulk.path, recursive = T)
  cmd <- paste(samtools_path, "merge", paste("-@", core), paste0(pseudobulk.path, "/all.bam"), paste0(bam_path, "/*.bam"), "--no-PG")
  msg <- paste0("[", Sys.time(), "] ", "Merge Bam Files: ", cmd)
  print(msg)
  system(command = cmd, wait = T)
  msg <- paste0("[", Sys.time(), "] ", "Merge Bam Files Finish.")
  print(msg)
  cmd <- paste(samtools_path, "index", paste0(pseudobulk.path, "/all.bam"))
  msg <- paste0("[", Sys.time(), "] ", "Bam File Index: ", cmd)
  print(msg)
  system(command = cmd, wait = T)
  msg <- paste0("[", Sys.time(), "] ", "Bam File Index Finish.")
  print(msg)
  return(pseudobulk.path)
}
