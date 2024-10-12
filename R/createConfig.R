#' @title Create configure file
#' @param host The IPv4 address that the application should listen on, default: 127.0.0.1
#' @param port The TCP port that the application should listen on, default: 9999
#' @param browser The HTML browser to be used by shiny
#' @return NULL
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
    gtf <- paras$Basic$refgenome$gtf_path
    gtf <- read.table(gtf,header = F,sep="\t",nrows=1)
    ref <- paras$Basic$refgenome$ref_path
    ref <- read.table(ref,header = F,sep="\t",nrows=1)
    phast.path <- paras$Task$impute$event_features$phast_path
    chromosomes <- seqnames(import(phast.path, which = GRanges("chr1:1-1")))
    chromosomes <- levels(chromosomes)
    if(grepl("^chr",gtf[1,1])){
      if(grepl("^>chr",ref[1,1])){
        paras$Task$event$remove_chr = "false"
        if(grepl("^chr",chromosomes[1])){
          paras$Task$impute$event_features$chr_prefix = ""
        }else{
          stop("The chromosome number prefixes in the reference genome file and the phast conservation file do not match.")
        }
      }else{
        stop("The chromosome number prefixes in the reference genome file and the annotation file do not match.")
      }
    }else{
      if(grepl("^>chr",ref[1,1])){
        stop("The chromosome number prefixes in the reference genome file and the annotation file do not match.")
      }else{
        paras$Task$event$remove_chr = "true"
        if(grepl("^chr",chromosomes[1])){
          paras$Task$impute$event_features$chr_prefix = "chr"
          print("The chromosome number prefixes in the reference genome file and the phast conservation file do not match.")
        }else{
          paras$Task$impute$event_features$chr_prefix = ""
        }
      }
    }
    paras_new <- toJSON(paras, pretty = T, auto_unbox = T)
    paras_file_new <- paste0(paras$Basic$work_path,'/',paras$DataSet,'_new.json')
    cat(paras_new, file = paras_file_new, fill = FALSE, labels = NULL, append = FALSE)
    return(paras)
}
