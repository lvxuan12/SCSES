#' @title Create Pseudobulk BAM by Merging Single Cell BAM Files
#'
#' @description This function merges multiple single-cell BAM files into a single
#' pseudobulk BAM file and creates an index for downstream analysis. It uses
#' samtools merge and index commands to perform the operations.
#'
#' @param paras A list object containing configuration parameters, typically
#'   loaded from \code{readSCSESconfig(paras_file)}.
#' @param bam_path Character string specifying the directory path containing
#'   single-cell BAM files. Default is extracted from \code{paras$Basic$bam_path}.
#' @param core Integer specifying the number of threads to use for parallel
#'   processing. Default is extracted from \code{paras$Basic$core}.
#' @param overwrite Logical value indicating whether to overwrite existing
#'   pseudobulk BAM file if it exists. Default is \code{TRUE}.
#'
#'
#' @return Character string containing the path to the directory where the
#'   pseudobulk BAM file and its index are created.
#' @export
#'
createPseudobulk <- function(
    paras, bam_path = paras$Basic$bam_path, core = paras$Basic$core,
    overwrite = TRUE) {
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
    if(file.exists(paste0(pseudobulk.path, "/all.bam"))){
        print("Pseudobulk bam file exists.")
        if (overwrite) {
            cmd <- paste(samtools_path, "merge", paste("-f -@", core), paste0(pseudobulk.path, "/all.bam"), paste0(bam_path, "/*.bam"), "--no-PG")
        }
    }else {
        cmd <- paste(samtools_path, "merge", paste("-@", core), paste0(pseudobulk.path, "/all.bam"), paste0(bam_path, "/*.bam"), "--no-PG")
    }
    msg <- paste0("[", Sys.time(), "] ", "Merge Bam Files: ", cmd)
    print(msg)
    system(command = cmd, wait = T)
    msg <- paste0("[", Sys.time(), "] ", "Merge Bam Files Finish.")
    print(msg)
    cmd <- paste(samtools_path, "index",paste("-@", core), paste0(pseudobulk.path, "/all.bam"))
    msg <- paste0("[", Sys.time(), "] ", "Bam File Index: ", cmd)
    print(msg)
    system(command = cmd, wait = T)
    msg <- paste0("[", Sys.time(), "] ", "Bam File Index Finish.")
    print(msg)
    return(pseudobulk.path)
}
