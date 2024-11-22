#' BulkSignalR ligand-receptor signature object for cluster comparisons
#'
#' S4 class to represent gene signatures associated with ligand-receptor
#' interactions that were inferred from the comparison of two clusters
#' of samples. This class inherits from BSRSignature.
#'
#' @slot cmp.name   The name of the comparison.
#' @slot tg.pval    A list of target genes P-values.
#' @slot tg.logFC   A list of target genes logFC.
#' @slot tg.expr   A list of target genes expression
#'
#' @export
#' @examples
#' new("BSRSignatureComp")
#'
setClass("BSRSignatureComp",
    contains = c("BSRSignature"),
    slots = c(
        cmp.name = "character",
        tg.expr = "list",
        tg.pval = "list",
        tg.logFC = "list"
    ),
    prototype = list(
        cmp.name = "myComparison",
        tg.expr = list(c(1, 2, 3)),
        tg.pval = list(c(0.05, 0.1, 0.008)),
        tg.logFC = list(c(-1, 0, 2))
    )
)

setValidity(
    "BSRSignatureComp",
    function(object) {
        if (!is.character(object@cmp.name)) {
            return("cmp.name is not character")
        }
        if (length(object@cmp.name) == 0) {
            return("cmp.name must have a length > 0")
        }
        if (!is.list(object@tg.expr)) {
            return("tg.expr is not a list")
        }
        if (!is.list(object@tg.pval)) {
            return("tg.pval is not a list")
        }
        if (!is.list(object@tg.logFC)) {
            return("tg.logFC is not a list")
        }

        TRUE
    }
)

setMethod("show", "BSRSignatureComp", function(object) {
    callNextMethod()
    cat("Cluster comparison name:", object@cmp.name, "\n")
})


# Accessors & setters ========================================================

#' Target gene expression accessor
#'
#' @name tgExpr
#' @aliases tgExpr,BSRSignatureComp-method
#' @param x BSRSignatureComp object
#' @return tg.expr
#' @examples
#' if(FALSE){
#' }
#' @export
setMethod("tgExpr", "BSRSignatureComp", function(x) x@tg.expr)

#' Target gene P-values accessor
#'
#' @name tgPval
#' @aliases tgPval,BSRSignatureComp-method
#' @param x BSRSignatureComp object
#' @return tg.pval
#' @export
setMethod("tgPval", "BSRSignatureComp", function(x) x@tg.pval)


#' Target gene logFC accessor
#'
#' @name tgLogFC
#' @aliases tgLogFC,BSRSignatureComp-method
#' @param x BSRSignatureComp object
#' @return tg.logFC
#' @export
setMethod("tgLogFC", "BSRSignatureComp", function(x) x@tg.logFC)

#' Comparison name accessor
#'
#' @name comparisonName
#' @aliases comparisonName,BSRSignatureComp-method
#' @param x BSRSignatureComp object
#' @return cmp.name
#' 
#' @export
setMethod("comparisonName", "BSRSignatureComp", function(x) x@cmp.name)
