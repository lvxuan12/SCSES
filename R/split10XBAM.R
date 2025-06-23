#' @title get single cell bam file from 10X CellRanger output
#' @description this function will split possorted_genome_bam.bam
#' from 10X CellRanger output to single cell bam based bam tag "CB:"
#'
#' @param CellRanger_path directory to CellRanger output for one sample
#' @param out_path directory to save single cell bam
#' @param java_path directory to java
#' @param core the number of threads
#' @param times Process the cells in times iterations, default: 50

#' @return single cell bam directory
#'
#' @export

split10XBAM <- function(CellRanger_path,out_path,java_path,core,times=50) {
    options("scipen" = 100)
    # script
    jar_path <- system.file("java", package = "SCSES")
    dir_shell <- system.file("shell", package = "SCSES")
    script_split <- paste0(dir_shell, "/split_10XBAM.sh")
    # input
    bam_file = list.files(CellRanger_path,
        pattern = "possorted_genome_bam.bam$",
        recursive = TRUE, full.names = TRUE
    )
    if (length(bam_file) == 0) {
        stop("possorted_genome_bam.bam not found")
    }
    bc_file = list.files(CellRanger_path,
        pattern = "barcodes.tsv.gz",
        recursive = TRUE, full.names = TRUE
    )
    if (length(bc_file) == 0) {
        stop("barcodes.tsv.gz not found")
    }
    if (length(bc_file) > 1) {
        bc_file = bc_file[grep("filtered_feature_bc", bc_file)]
    }

    # output
    if (!dir.exists(out_path)) {
        dir.create(out_path, recursive = T)
    }

    msg <- paste0("[", Sys.time(), "] ", "Processing cell barcodes...")
    print(msg)
    cmd <- paste("cp", bc_file, paste0(out_path, "/barcodes.tsv.gz"))
    system(command = cmd, wait = T)
    cmd <- paste("gzip -d", paste0(out_path, "/barcodes.tsv.gz"))
    system(command = cmd, wait = T)

    # run split
    msg <- paste0("[", Sys.time(), "] ", "Running split...")
    print(msg)
    bc_file_new <- paste0(out_path, "/barcodes.tsv")
    jar_file <- paste0(jar_path, "/splitCellBam3.jar")
    lib_path <- paste0(jar_path, "/lib")
    log_file <- paste0(out_path, "/splitcell_", rbinom(1, size = 1000000000, prob = 0.5), ".log")
    if (file.exists(log_file)) {
        file.remove(log_file)
    }
    cmd <- paste(
        "bash", script_split,
        bam_file,
        out_path,
        jar_file,
        lib_path,
        java_path,
        core,
        bc_file_new,times, ">>", log_file, "2>&1"
    )
    res <- system(command = cmd, intern = T, wait = T)
    if (!is.null(attributes(res)) && attributes(res)$status == 1) {
        msg <- paste0("[", Sys.time(), "] ", "Run split Error.")
        stop(msg)
    } else {
        msg <- paste0("[", Sys.time(), "] ", "Run split Finish.")
        print(msg)
        return(paste0(out_path, "/bam/"))
    }
}
