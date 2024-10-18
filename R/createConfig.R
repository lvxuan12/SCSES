#' @title Create configure file
#' @param host The IPv4 address that the application should listen on, default: 127.0.0.1
#' @param port The TCP port that the application should listen on, default: 9999
#' @param launch.browser If true, the system's default web browser will be launched automatically after the app is started.
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
