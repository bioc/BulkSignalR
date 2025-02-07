% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataPrepare.R
\name{BSRDataModel}
\alias{BSRDataModel}
\title{Prepare a BSRDataModel object from expression data}
\usage{
BSRDataModel(
  counts,
  normalize = TRUE,
  symbol.col = NULL,
  min.count = 10,
  prop = 0.1,
  method = c("UQ", "TC"),
  log.transformed = FALSE,
  min.LR.found = 80,
  species = "hsapiens",
  conversion.dict = NULL,
  UQ.pc = 0.75,
  x.col = NULL,
  y.col = NULL,
  barcodeID.col = NULL
)
}
\arguments{
\item{counts}{A table or matrix of read counts.}

\item{normalize}{A logical indicating whether \code{counts} should be
normalized according to \code{method} or if it was normalized beforehand.}

\item{symbol.col}{The index of the column containing the gene symbols in case
those are not the row names of \code{counts} already.}

\item{min.count}{The minimum read count of a gene to be considered expressed
in a sample.}

\item{prop}{The minimum proportion of samples where a gene must be
expressed higher than \code{min.count} to keep that gene.}

\item{method}{The normalization method ('UQ' for upper quartile or 'TC'
for total count). If \code{normalize==FALSE}, then method must be
used to document the name of the normalization method applied by the user.}

\item{log.transformed}{A logical indicating whether expression data were
already log2-transformed, e.g., some microarray data.}

\item{min.LR.found}{The minimum number of ligands or receptors found in
\code{count} row names after eliminating the rows containing too many
zeros according to \code{min.count} and \code{prop}.}

\item{species}{Data were obtained for this organism.}

\item{conversion.dict}{Correspondence table of HUGO gene symbols
human/nonhuman. Not used unless the organism is different from human.}

\item{UQ.pc}{Percentile for upper-quartile normalization, number
between 0 and 1 (in case the default 0.75 - hence the name - is not
appropriate).}

\item{x.col}{In a SpatialExperiment object, the index of the column
containing x coordinates in the dafaframe returned by rowData(), usually 
named array_row}

\item{y.col}{In a SpatialExperiment object, the index of the column
containing y coordinates in the dafaframe returned by rowData(), usually 
named array_col}

\item{barcodeID.col}{In a SpatialExperiment object, the index of the column
containing barcodeID in the dafaframe returned by colData(), usually named
barcode_id}
}
\value{
A BSRModelData object with empty model parameters.
}
\description{
Take a matrix or data frame containing RNA sequencing,
microarray, or expression proteomics data and return a BSRDataModel
object ready for subsequent training. Normally, BSRDataModel objects
are not instantiated directly, but through this function.
}
\details{
The \code{counts} matrix or table should be provided with expression
  levels of protein coding genes in each samples (column) and
  \code{rownames(counts)} set to HUGO official gene symbols.
  For commodity, it is also possible 
  to provide \code{counts} with the
  gene symbols stored in one of its columns. This column must be specified
  with \code{symbol.col}. In such a case, \code{BSRDataModel} will extract
  this column and use it to set the row names. Because row names must be
  unique, \code{BSRDataModel} will eliminate rows with duplicated gene
  symbols by keeping the rows with maximum average expression. Gene symbol
  duplication may occur in protein coding genes after genome alignment
  due to errors in genome feature annotation files (GTF/GFF), where a handful
  of deprecated gene annotations might remain, or
  some genes are not given their fully specific symbols. If your read count
  extraction pipeline does not take care of this phenomenon, the maximum mean
  expression selection strategy implemented here should solve this difficulty
  for the sake of inferring ligand-receptor interactions.

  If \code{normalize} is \code{TRUE} then normalization is performed
  according to \code{method}. If those two simple methods are not satisfying,
  then it is possible to provide a pre-normalized matrix setting
  \code{normalize} to \code{FALSE}. In such a case, the parameter
  \code{method} must be used to document the name of the normalization
  algorithm used.

  In case proteomic or microarray data are provided, \code{min.count} must be
  understood as its equivalent with respect to those data.
}
\examples{
data(sdc, package = "BulkSignalR")
idx <- sample(nrow(sdc), 4000)
bsrdm <- BSRDataModel(sdc[idx, c("N22","SDC17")],
normalize = FALSE,method="UQ")
}
