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

createConfigshiny <-function(host="127.0.0.1",port = 9999,browser = 'firefox',...){
    options(browser = browser)
    shiny_path <- system.file("shiny", package = "SCSES")
    shiny::runApp(shiny_path,
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
