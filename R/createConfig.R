#' @title Launch SCSES Configuration Shiny Application
#'
#' This function starts a Shiny web application that provides an interactive
#' interface for configuring SCSES parameters. The application runs on a 
#' local web server and can be accessed through a web browser.
#'
#' @param host A character string specifying the IPv4 address that the 
#'   application should listen on. Default is "127.0.0.1" (localhost only).
#' @param port An integer specifying the TCP port that the application should 
#'   listen on. Must be between 1 and 65535. Default is 9999.
#' @param launch.browser A logical value indicating whether to automatically 
#'   launch the system's default web browser to open the application. 
#'   Default is FALSE.
#'
#' @return Returns \code{NULL} invisibly. The function is called for its side 
#'   effect of launching the Shiny application.
#' @details
#' The function locates the Shiny application directory within the SCSES package
#' and launches it using \code{shiny::runApp()}. The application will continue 
#' running until manually stopped (Ctrl+C in R console).
#' 
#' \strong{Important Notes:}
#' \itemize{
#'   \item Setting \code{launch.browser = TRUE} may cause errors in headless 
#'         environments (servers without GUI) or when no default browser is configured
#'   \item Ensure the specified port is not already in use by another application
#' }
#' 
#' \strong{Server Environment Usage:}
#' 
#' For server environments, it is recommended to:
#' \enumerate{
#'   \item Set the host to the server's IP address
#'   \item Set \code{launch.browser = FALSE} to avoid browser launch errors
#'   \item Manually access the application URL shown in the console
#' }
#' 
#' After starting the application, you will see console output similar to:
#' \preformatted{Listening on http://123.678.112.78:9999}
#' 
#' Copy this URL and paste it into your web browser to access the application.
#'
#' @export
#' @importFrom shiny runApp
#' @import shinydashboard
#' @import shinyFiles
#' @import fs
#' @import shinyWidgets


createConfigshiny <-function(host="127.0.0.1",port = 9999,launch.browser=F,...){
    shiny_path <- system.file("shiny", package = "SCSES")
    shiny::runApp(shiny_path,
        host = host, port = port,launch.browser, ...
    )
    return(NULL)
}


#' @title read configure
#' @param paras_file path to configure file

#' @return list of parameters
#' @export
#' @importFrom jsonlite fromJSON toJSON

readSCSESconfig <- function(paras_file) {
    paras <- fromJSON(paras_file)
    bam_path <- paras$Basic$bam_path
    bam_file <- list.files(path = bam_path,pattern = "*bam",full.names = T)[1]
    #script
    samtools_path = paras$Basic$samtools_path
    cmd <- paste(samtools_path, "view",bam_file ,
                 "|sed -n '1p'| awk -F'\t' '{print $3}'")
    chr_bam=system(command = cmd, intern = TRUE)

    gtf <- paras$Basic$refgenome$gtf_path
    gtf <- read.table(gtf,header = F,sep="\t",nrows=1)
    gff <- paras$Basic$refgenome$gff_path
    gff <- read.table(gff,header = F,sep="\t",nrows=1)
    ref <- paras$Basic$refgenome$ref_path
    ref <- read.table(ref,header = F,sep="\t",nrows=1)

    if(grepl("^chr",gtf[1,1]) & grepl("^chr",gff[1,1]) &
       grepl("^>chr",ref[1,1]) & grepl("^chr",chr_bam)){
      phast.path <- paras$Task$impute$event_features$phast_path
      chromosomes <- seqnames(import(phast.path, which = GRanges("chr1:1-1")))
      chromosomes <- levels(chromosomes)
      paras$Task$event$remove_chr = "false"
      if(grepl("^chr",chromosomes[1])){
        paras$Task$impute$event_features$chr_prefix = ""
      }else{
        stop("The chromosome number prefixes in the bam file and the phast conservation file do not match.")
      }
    }else if(!grepl("^chr",gtf[1,1]) & !grepl("^chr",gff[1,1]) &
             !grepl("^>chr",ref[1,1]) & !grepl("^chr",chr_bam)){
      paras$Task$event$remove_chr = "true"
      if(grepl("^chr",chromosomes[1])){
        paras$Task$impute$event_features$chr_prefix = "chr"
        print("The chromosome number prefixes in the bam file and the phast conservation file do not match.")
      }else{
        paras$Task$impute$event_features$chr_prefix = ""
      }
    }else{
      stop("The chromosome prefixes in the bam files and reference/annotation file do not match.")
    }

    paras_new <- toJSON(paras, pretty = T, auto_unbox = T)
    paras_file_new <- paste0(paras$Basic$work_path,'/',paras$DataSet,'_new.json')
    cat(paras_new, file = paras_file_new, fill = FALSE, labels = NULL, append = FALSE)
    return(paras)
}
