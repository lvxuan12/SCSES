% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detectEvent.R
\name{getRIevent}
\alias{getRIevent}
\title{Detect RI events}
\usage{
getRIevent(
  work_path,
  bam_path,
  gtf,
  ref,
  core,
  readlength,
  script_irfinder,
  irfinder_path,
  samtools_path,
  star_path,
  log_file,
  star_ref_path = NULL
)
}
\arguments{
\item{work_path}{directory to deposite all SCSES results in}

\item{bam_path}{directory to Pseudobulk bam file}

\item{gtf}{the gene annotation in gtf format}

\item{ref}{fasta file}

\item{core}{the number of threads}

\item{readlength}{the length of each read}

\item{script_irfinder}{directory to script for running IRFinder}

\item{irfinder_path}{directory to executable file of IRFinder}

\item{samtools_path}{directory to samtools}

\item{star_path}{directory to STAR}

\item{log_file}{file saving stdout and stderr information}
}
\value{
Splicing event directory
}
\description{
Detect RI events
}
