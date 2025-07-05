#' @title Launch SCSES Configuration Shiny Application
#' @description This function starts a Shiny web application that provides an interactive
#' interface for configuring SCSES parameters. The application runs on a
#' local web server and can be accessed through a web browser to generate
#' configuration JSON files.
#'
#' @param host A character string specifying the IPv4 address that the
#'   application should listen on. Default is "127.0.0.1" (localhost only).
#' @param port An integer specifying the TCP port that the application should
#'   listen on. Default is 9999. Choose an available port if the default is occupied.
#' @param launch.browser A logical value indicating whether to automatically
#'   launch the system's default web browser to open the application.
#'   Default is \code{FALSE}. Setting \code{launch.browser = TRUE} may cause
#'   errors in headless environments (servers without GUI) or when no default
#'   browser is configured.
#'
#' @return Returns \code{NULL} invisibly. The function is called for its side
#'   effect of launching the Shiny application.
#'
#' @section Usage Scenarios:
#'
#' \strong{SCSES Docker Container:}
#'
#' If running within the SCSES docker container:
#' \preformatted{
#' createConfigshiny(host = "localhost", launch.browser = TRUE)
#' }
#'
#' \strong{Non-Docker Environment:}
#'
#' \preformatted{
#' createConfigshiny(host, port, launch.browser = FALSE)
#' }
#'
#' @section Workflow:
#' \enumerate{
#'   \item Run the function to start the Shiny server
#'   \item Copy the displayed URL and open it in your web browser (if \code{launch.browser = FALSE})
#'   \item Fill in the required parameters using the interactive interface
#'   \item Hover over widgets to see parameter descriptions
#'   \item Click "Create Config" to generate the JSON configuration file
#'   \item The configuration file will be saved to your specified work path
#' }
#'
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

#' @title Launch SCSES Configuration Shiny Application for test data
#' @description This function starts a Shiny web application that provides an interactive
#' interface for configuring SCSES parameters. The application runs on a
#' local web server and can be accessed through a web browser.
#'
#' @param host A character string specifying the IPv4 address that the
#'   application should listen on. Default is "127.0.0.1" (localhost only).
#' @param port An integer specifying the TCP port that the application should
#'   listen on. Default is 9999.
#' @param launch.browser A logical value indicating whether to automatically
#'   launch the system's default web browser to open the application.
#'   Default is FALSE. Setting launch.browser = TRUE may cause errors in headless
#'   environments (servers without GUI) or when no default browser is configured
#'
#' @return Returns \code{NULL} invisibly. The function is called for its side
#'   effect of launching the Shiny application.
#'
#' \strong{SCSES docker container Usage:}
#'
#' If you use it in the SCSES docker container, the command should be :
#'
#' \code{createConfigshiny(host = "localhost",launch.browser=TRUE)}
#'
#' \strong{non-docker Usage:}
#'
#' For non-docker users, the full command should be:
#'
#' \code{createConfigshiny(host, port, launch.browser=FALSE)}
#'
#' After running \code{createConfigshiny}, you will see a URL appear in the console.
#'
#' Copy this URL and paste it into your web browser to access the application.
#' An interactive window will popup, which allow you to fill some parameters,
#' such as Bam File Path, and Work Path.
#'
#' The meaning of each parameter can be found by hovering the mouse over the widget.
#'
#' Finally, you can click Create Config button and a json file will be generated
#' in the work_path you provided if successful.
#'
#' @export
#' @importFrom shiny runApp
#' @import shinydashboard
#' @import shinyFiles
#' @import fs
#' @import shinyWidgets


createDemoConfigshiny <-function(host="127.0.0.1",port = 9999,launch.browser=F,...){
  shiny_path <- system.file("shiny_demo", package = "SCSES")
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
