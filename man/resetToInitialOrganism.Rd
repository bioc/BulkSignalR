% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BSRInference-class.R
\name{resetToInitialOrganism}
\alias{resetToInitialOrganism}
\alias{resetToInitialOrganism,BSRInference-method}
\title{Reset gene names to initial organism providen in first instance}
\usage{
\S4method{resetToInitialOrganism}{BSRInference}(obj, conversion.dict)
}
\arguments{
\item{obj}{BSRInference object}

\item{conversion.dict}{A dictionnary}
}
\value{
An BSRInference object updated for gene names.
The gene names are replaced by the ones from
the organism providen in first instance.
}
\description{
Reset gene names to initial organism providen in first instance
}
\examples{
data(bodyMap.mouse, package = "BulkSignalR")
data(bsrinf.mouse, package = "BulkSignalR")
data(ortholog.dict, package = "BulkSignalR")

#idx <- sample(nrow(bodyMap.mouse), 7500)

#bodyMap.mouse <- bodyMap.mouse[idx,1:3]

#ortholog.dict <- findOrthoGenes(
#    from_organism = "mmusculus",
#    from_values = rownames(bodyMap.mouse)
#)

#matrix.expression.human <- convertToHuman(
#    counts = bodyMap.mouse,
#    dictionary = ortholog.dict
#)

#bsrdm <- BSRDataModel(
#    counts = matrix.expression.human,
#    species = "mmusculus",
#    conversion.dict = ortholog.dict
#)

#bsrdm <- learnParameters(bsrdm,
#    quick = TRUE  
#)

#reactSubset <- getResource(resourceName = "Reactome",
#cache = TRUE)

#subset <- c("REACTOME_BASIGIN_INTERACTIONS",
#"REACTOME_SYNDECAN_INTERACTIONS",
#"REACTOME_ECM_PROTEOGLYCANS",
#"REACTOME_CELL_JUNCTION_ORGANIZATION")

#reactSubset <- reactSubset[
#reactSubset$`Reactome name` \%in\% subset,]

#bsrinf.mouse <- BSRInference(bsrdm,reference="REACTOME")

bsrinf <- resetToInitialOrganism(bsrinf.mouse, 
conversion.dict = ortholog.dict)
}
