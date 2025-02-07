% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modelTraining.R
\name{.getEmpiricalNull}
\alias{.getEmpiricalNull}
\title{Sampling of correlations downstream the receptors null distribution}
\usage{
.getEmpiricalNull(obj)
}
\arguments{
\item{obj}{A BSRDatamodel without learned paramaters.}
}
\value{
A list of \code{n.rand} tables such as output by
  \code{.checkReceptorSignaling}. Each table is computed from a randomized
  expression matrix (randomized \code{ncounts}).
}
\description{
Perform receptor downstream analysis with
\code{.checkReceptorSignaling} based on randomized expression data and
ligand-receptor pairs selected from the same randomized data.
}
\details{
A large number of correlations 
  (ligand-receptor and receptor-downstream target genes) 
  is reported in each randomized matrix. Therefore,
  \code{n.rand} should be
  given a modest value to avoid unnecessarily long computations.

  See \code{\link{.checkReceptorSignaling}} for more details about the
  parameters.
}
\keyword{internal}
