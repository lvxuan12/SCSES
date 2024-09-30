#' @title Create configure file
#' @param host The IPv4 address that the application should listen on, default: 127.0.0.1
#' @param port The TCP port that the application should listen on, default: 9999
 
#' @return NULL 
#' @export
#' @importFrom shiny runApp

createConfigshiny <-function(host="127.0.0.1",port = 9999,...){
    shiny::runApp("/disk/lvxuan/Single-Splicing/src/app/",
        host = host, port = port, ...
    )
    return(NULL)
}


#' @title read configure
#' @param paras_file path to configure file

#' @return list of parameters
#' @export
#' @importFrom jsonlite fromJSON

readSCSESconfig <- function(paras_file) {
    paras <- fromJSON(paras_file)
    return(paras)
}
