#' BulkSignalR Cluster Comparison Object
#'
#' An S4 class to represent the comparison of two clusters of samples to
#' infer LR interactions based on the resulting P-values,
#' log-fold-changes (logFC), and expression values.
#'
#' @slot colA   Column indices for the samples in cluster A.
#' @slot colB   Column indices for the samples in cluster B.
#' @slot stats  Comparison statistics A versus B as a data.frame and
#' containing at least two columns named 'pval', 'logFC', and 'expr'.
#'
#' @export
#' @examples
#' if(FALSE){
#' bsrdm <- new("BSRDataModel",
#'     ncounts = matrix(1.5,
#'         nrow = 2, ncol = 4,
#'         dimnames = list(c("A", "B"), c("E", "F", "G", "H"))
#'     ),
#'     log.transformed = TRUE, normalization = "TC"
#' )
#' bsrdm.comp <- as.BSRDataModelComp(bsrdm)
#' colA <- as.integer(1:2)
#' colB <- as.integer(3:4)
#' n <- nrow(ncounts(bsrdm.comp))
#' edger.stats <- data.frame(
#'     pval = runif(n), logFC = rnorm(n, 0, 2),
#'     expr = c(1, 2)
#' )
#' rownames(edger.stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClusterComp(bsrdm.comp, colA, colB, edger.stats)
#' bsrdm.comp <- addClusterComp(bsrdm.comp, bsrcc, "my_comparison")
#' }
setClass("BSRClusterComp",
    slots = c(
        colA = "integer",
        colB = "integer",
        stats = "data.frame"
    ),
    prototype = list(
        colA = as.integer(c(1, 2)),
        colB = as.integer(c(3, 4)),
        stats = data.frame(
            pval = c(0.01, 0.01),
            logFC = c(1, -1), expr = c(1, 2)
        )
    )
)

setValidity(
    "BSRClusterComp",
    function(object) {
        if (!is.integer(object@colA)) {
            return("colA indices are not all integers")
        }
        if (length(object@colA) < 1) {
            return("colA empty")
        }
        if (!is.integer(object@colB)) {
            return("colB indices are not all integers")
        }
        if (length(object@colB) < 1) {
            return("colB empty")
        }
        if (length(intersect(object@colA, object@colB)) > 0) {
            return("colA and colB are not disjoint")
        }
        if (!is.data.frame(object@stats)) {
            return("specified stats are not a data.frame")
        }

        TRUE
    }
)

setMethod(
    "show", "BSRClusterComp",
    function(object) {
        if (length(object@colA) > 5) {
            cat("Cluster A columns:", object@colA[seq_len(5)], "...\n")
        } else {
            cat(
                "Cluster A columns:",
                object@colA[seq_len(length(object@colA))], "\n"
            )
        }
        if (length(object@colB) > 5) {
            cat("Cluster B columns:", object@colB[seq_len(5)], "...\n")
        } else {
            cat(
                "Cluster B columns:",
                object@colB[seq_len(length(object@colB))], "\n"
            )
        }
        print(utils::head(object@stats))
    }
)


# Accessors & setters ========================================================

setGeneric("colA", signature="x",
    function(x) standardGeneric("colA")
)
#' Cluster A columns accessor
#'
#' @name colA
#' @aliases colA,BSRClusterComp-method
#' @param x object BSRClusterComp
#' @return colA
#' @examples
#' bsrdm <- new("BSRDataModel",
#'     ncounts = matrix(1.5,
#'         nrow = 2, ncol = 4,
#'         dimnames = list(c("A", "B"), c("E", "F","G","H"))
#'     ),
#'     log.transformed = TRUE, normalization = "TC"
#' )
#' bsrdm.comp <- as.BSRDataModelComp(bsrdm)
#' colA <- as.integer(1:2)
#' colB <- as.integer(3:4)
#' n <- nrow(ncounts(bsrdm.comp))
#' edger.stats <- data.frame(
#'     pval = runif(n), logFC = rnorm(n, 0, 2),
#'     expr = c(1, 2)
#' )
#' rownames(edger.stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClusterComp(bsrdm.comp, colA, colB, edger.stats)
#' colA(bsrcc)
#' @export
setMethod("colA", "BSRClusterComp", function(x) x@colA)

setGeneric("colA<-", signature=c("x", "value"),
    function(x, value) standardGeneric("colA<-")
)
#' Cluster A columns setter (internal use only)
#'
#' @param x object BSRClusterComp
#' @param value value to be set for BSRClusterComp
#' @return returns \code{NULL}
#' @keywords internal
setMethod("colA<-", "BSRClusterComp", function(x, value) {
    x@colA <- value
    methods::validObject(x)
    x
})

setGeneric("colB", signature="x",
    function(x) standardGeneric("colB")
)
#' Cluster B columns accessor
#'
#' @name colB
#' @aliases colB,BSRClusterComp-method
#' @param x object BSRClusterComp
#' @return colB
#' @examples
#' bsrdm <- new("BSRDataModel",
#'     ncounts = matrix(1.5,
#'         nrow = 2, ncol = 4,
#'         dimnames = list(c("A", "B"), c("E", "F","G","H"))
#'     ),
#'     log.transformed = TRUE, normalization = "TC"
#' )
#' bsrdm.comp <- as.BSRDataModelComp(bsrdm)
#' colA <- as.integer(1:2)
#' colB <- as.integer(3:4)
#' n <- nrow(ncounts(bsrdm.comp))
#' edger.stats <- data.frame(
#'     pval = runif(n), logFC = rnorm(n, 0, 2),
#'     expr = c(1, 2)
#' )
#' rownames(edger.stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClusterComp(bsrdm.comp, colA, colB, edger.stats)
#' colB(bsrcc)
#' @export
setMethod("colB", "BSRClusterComp", function(x) x@colB)

setGeneric("colB<-", signature=c("x", "value"),
    function(x, value) standardGeneric("colB<-")
)
#' Cluster B columns setter (internal use only)
#'
#' @param x object BSRClusterComp
#' @param value value to be set for BSRClusterComp
#' @return returns \code{NULL}
#' @keywords internal
setMethod("colB<-", "BSRClusterComp", function(x, value) {
    x@colB <- value
    methods::validObject(x)
    x
})


setGeneric("stats", signature="x",
    function(x) standardGeneric("stats")
)
#' Cluster comparison statistics accessor
#'
#' @name stats
#' @aliases stats,BSRClusterComp-method
#' @param x BSRClusterComp object
#' @return stats
#' @examples
#' bsrdm <- new("BSRDataModel",
#'     ncounts = matrix(1.5,
#'         nrow = 2, ncol = 4,
#'         dimnames = list(c("A", "B"), c("E", "F","G","H"))
#'     ),
#'     log.transformed = TRUE, normalization = "TC"
#' )
#' bsrdm.comp <- as.BSRDataModelComp(bsrdm)
#' colA <- as.integer(1:2)
#' colB <- as.integer(3:4)
#' n <- nrow(ncounts(bsrdm.comp))
#' edger.stats <- data.frame(
#'     pval = runif(n), logFC = rnorm(n, 0, 2),
#'     expr = c(1, 2)
#' )
#' rownames(edger.stats) <- rownames(ncounts(bsrdm.comp))
#' bsrcc <- defineClusterComp(bsrdm.comp, colA, colB, edger.stats)
#' stats(bsrcc)
#' @export
setMethod("stats", "BSRClusterComp", function(x) x@stats)

setGeneric("stats<-", signature=c("x", "value"),
    function(x, value) standardGeneric("stats<-")
)
#' Cluster comparison statistics setter (internal use only)
#'
#' @param x object BSRClusterComp
#' @param value value to be set for BSRClusterComp
#' @return returns \code{NULL}
#' @keywords internal
setMethod("stats<-", "BSRClusterComp", function(x, value) {
    x@stats <- value
    methods::validObject(x)
    x
})
