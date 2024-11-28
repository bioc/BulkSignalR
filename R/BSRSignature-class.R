#' BulkSignalR ligand-receptor signature Object
#'
#' S4 class to represent gene signatures
#' of inferred ligand-receptor interactions, including
#' their reduced versions.
#'
#' @slot ligands   A list of ligands, one entry per LR interaction.
#' @slot receptors   A list of receptors, one entry per LR interaction.
#' @slot tg.genes  A list of target genes, one entry per LR interaction.
#' @slot pathways  An atomic vector of pathway names, one per interaction.
#' @slot tg.corr  A list of target genes correlation.
#'
#' @export
#' @examples
#' new("BSRSignature")
setClass("BSRSignature",
    slots = c(
        pathways = "character",
        ligands = "list",
        receptors = "list",
        tg.genes = "list",
        tg.corr = "list"
    ),
    prototype = list(
        pathways = "path 1",
        ligands = list("A"),
        receptors = list("B"),
        tg.genes = list(c("a", "b", "c")),
        tg.corr = list(c(0.1, 0.2, 0.3))
    )
)

setValidity(
    "BSRSignature",
    function(object) {
        if (!is.character(object@pathways)) {
            return("pathways is not character")
        }
        if (!is.list(object@ligands)) {
            return("ligands is not a list")
        }
        if (!is.list(object@receptors)) {
            return("receptors is not a list")
        }
        if (!is.list(object@tg.genes)) {
            return("tg.genes is not a list")
        }
        if (!is.list(object@tg.corr)) {
            return("tg.corr is not a list")
        }
        TRUE
    }
)
setMethod("show", "BSRSignature", function(object) {
    print(data.frame(
        L = vapply(
            object@ligands, function(x) paste(x, collapse = ";"),
            character(1)
        ),
        R = vapply(
            object@receptors, function(x) paste(x, collapse = ";"),
            character(1)
        ),
        pathways = object@pathways,
        tgGenes = vapply(
            object@tg.genes, function(x) paste(x, collapse = ";"),
            character(1)
        )
    )[seq_len(min(5, length(object@ligands))), ])
})


# Accessors & setters ========================================================

setGeneric("pathways", signature="x",
    function(x) standardGeneric("pathways")
)
#' pathways accessor
#'
#' @name pathways
#' @aliases pathways,BSRSignature-method
#' @param x BSRSignature
#' @return pathways
#' @examples
#' bsr.sig <- new("BSRSignature")
#' pathways(bsr.sig)
#' @export
setMethod("pathways", "BSRSignature", function(x) x@pathways)

#' ligands accessor
#'
#' @name ligands
#' @aliases ligands,BSRSignature-method
#' @param x BSRSignature
#' @return ligands
#' @examples
#' bsr.sig <- new("BSRSignature")
#' ligands(bsr.sig)
#' @export
setMethod("ligands", "BSRSignature", function(x) x@ligands)

#' receptors accessor
#'
#' @name receptors
#' @aliases receptors,BSRSignature-method
#' @param x BSRSignature
#' @return receptors
#' @examples
#' bsr.sig <- new("BSRSignature")
#' ligands(bsr.sig)
#' @export
setMethod("receptors", "BSRSignature", function(x) x@receptors)

#' Target genes accessor
#'
#' @name tgGenes
#' @aliases tgGenes,BSRSignature-method
#' @param x BSRSignature
#' @return tgGenes
#' @examples
#' bsr.sig <- new("BSRSignature")
#' tgGenes(bsr.sig)
#' @export
setMethod("tgGenes", "BSRSignature", function(x) x@tg.genes)

#' Target gene correlations accessor
#'
#' @name tgCorr
#' @aliases tgCorr,BSRSignature-method
#' @param x BSRSignature
#' @return tgCorr
#' @examples
#' bsr.sig <- new("BSRSignature")
#' tgCorr(bsr.sig)
#' @export
setMethod("tgCorr", "BSRSignature", function(x) x@tg.corr)

