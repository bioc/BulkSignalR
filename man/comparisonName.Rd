% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BSRInferenceComp-class.R,
%   R/BSRSignatureComp-class.R
\name{comparisonName}
\alias{comparisonName}
\alias{comparisonName,BSRInferenceComp-method}
\alias{comparisonName,BSRSignatureComp-method}
\title{Comparison name accessor}
\usage{
\S4method{comparisonName}{BSRInferenceComp}(x)

\S4method{comparisonName}{BSRSignatureComp}(x)
}
\arguments{
\item{x}{BSRSignatureComp object}
}
\value{
cmp.name

cmp.name
}
\description{
Comparison name accessor

Comparison name accessor
}
\examples{
bsrinf <- new("BSRInferenceComp")
comparisonName(bsrinf)
}
