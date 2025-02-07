% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BSRSignature-class.R
\name{BSRSignature}
\alias{BSRSignature}
\title{Extract gene signatures of LR pair activity}
\usage{
BSRSignature(obj, pval.thres = NULL, qval.thres = NULL, with.pw.id = FALSE)
}
\arguments{
\item{obj}{BSRinference object.}

\item{pval.thres}{P-value threshold.}

\item{qval.thres}{Q-value threshold.}

\item{with.pw.id}{A logical indicating whether the ID of a pathway
should be concatenated to its name.}
}
\value{
A BSRSignature object containing a gene signature for each triple
ligand-receptor pair. A reduction to the best pathway
for each pair is automatically performed and the gene signature is
comprised of the ligand, the receptor,
and all the target genes with rank equal or superior to \code{pairs$rank}.
}
\description{
Obtains gene signatures reflecting ligand-receptor as well as
receptor downstream activity to
score ligand-receptor pairs across samples subsequently with
\code{"\link[=BSRInference-class]{scoreLRGeneSignatures}"}
}
\examples{
data(bsrinf, package = "BulkSignalR")

bsrinf.redP <- reduceToPathway(bsrinf)
bsrsig.redP <- BSRSignature(bsrinf, qval.thres = 0.001)

}
