% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createConfig.R
\name{createConfigshiny}
\alias{createConfigshiny}
\title{Launch SCSES Configuration Shiny Application}
\usage{
createConfigshiny(host = "127.0.0.1", port = 9999, launch.browser = F, ...)
}
\arguments{
\item{host}{A character string specifying the IPv4 address that the
application should listen on. Default is "127.0.0.1" (localhost only).}

\item{port}{An integer specifying the TCP port that the application should
listen on. Default is 9999. Choose an available port if the default is occupied.}

\item{launch.browser}{A logical value indicating whether to automatically
launch the system's default web browser to open the application.
Default is \code{FALSE}. Setting \code{launch.browser = TRUE} may cause
errors in headless environments (servers without GUI) or when no default
browser is configured.}
}
\value{
Returns \code{NULL} invisibly. The function is called for its side
  effect of launching the Shiny application.
}
\description{
This function starts a Shiny web application that provides an interactive
interface for configuring SCSES parameters. The application runs on a
local web server and can be accessed through a web browser to generate
configuration JSON files.
}
\section{Usage Scenarios}{


\strong{SCSES Docker Container:}

If running within the SCSES docker container:
\preformatted{
createConfigshiny(host = "localhost", launch.browser = TRUE)
}

\strong{Non-Docker Environment:}

\preformatted{
createConfigshiny(host, port, launch.browser = FALSE)
}
}

\section{Workflow}{

\enumerate{
  \item Run the function to start the Shiny server
  \item Copy the displayed URL and open it in your web browser (if \code{launch.browser = FALSE})
  \item Fill in the required parameters using the interactive interface
  \item Hover over widgets to see parameter descriptions
  \item Click "Create Config" to generate the JSON configuration file
  \item The configuration file will be saved to your specified work path
}
}

