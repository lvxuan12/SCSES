% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detectEvent.R
\name{getA5SSid}
\alias{getA5SSid}
\title{Generate A5SS event id}
\usage{
getA5SSid(file, outfile, remove.chr = F, junctionReads)
}
\arguments{
\item{file}{directory to MAJIQ output}

\item{outfile}{directory to output file A5SS.txt}

\item{remove.chr}{if TRUE,remove the "chr" character of event id}

\item{junctionReads}{the minimum total number of reads for any junction}
}
\value{
The number of A5SS events
}
\description{
Generate A5SS event id
}
