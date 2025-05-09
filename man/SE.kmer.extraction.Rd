% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eventFeatureFunction.R
\name{SE.kmer.extraction}
\alias{SE.kmer.extraction}
\title{get kmer of exons and splice region associated with SE event}
\usage{
SE.kmer.extraction(splice.region, exon.region)
}
\arguments{
\item{splice.region}{data.frame with sequence of splice region, \code{SE.splice.seq.extraction}}

\item{exon.region}{data.frame with sequence of exon region, \code{SE.exon.seq.extraction}}
}
\value{
vector of sequence kmer
}
\description{
get kmer of exons and splice region associated with SE event
}
\keyword{internal}
