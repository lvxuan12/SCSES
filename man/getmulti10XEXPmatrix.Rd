% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getGeneExpression.R
\name{getmulti10XEXPmatrix}
\alias{getmulti10XEXPmatrix}
\title{Get gene expression matrix from Cell Ranger of multiple samples}
\usage{
getmulti10XEXPmatrix(
  paras,
  expr_path,
  filter.mt = paras$Basic$filter_sc$filter.mt,
  filter.rp = paras$Basic$filter_sc$filter.rp,
  sample_name
)
}
\arguments{
\item{paras}{A list object parsed from SCSES JSON parameter file using
\code{readSCSESconfig()}.}

\item{expr_path}{Character string specifying the directory path to the
Cell Ranger output folder containing the HDF5 file.}

\item{filter.mt}{Logical value indicating whether to filter out
mitochondrial genes (genes with names starting with "MT-" or "mt-").
Default: TRUE (from \code{paras$Basic$filter_sc$filter.mt}).}

\item{filter.rp}{Logical value indicating whether to filter out
ribosomal genes (genes with names starting with "RPS", "RPL", "Rps",
or "Rpl"). Default: TRUE (from \code{paras$Basic$filter_sc$filter.rp}).}

\item{sample_name}{the name of directory for different samples in the expr_path,
different samples seperated by ";"
If the dataset contains multiple samples, the Cell Ranger outputs for
different samples need to be placed in different folders under the
expr_path directory, named after the sample_name.}
}
\value{
gene expression matrix path
}
\description{
Read count matrix from 10X CellRanger hdf5 file
save raw and normalized UMI counts sparse matrix to work_path/rds/
}
