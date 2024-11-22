#' BulkSignalR Cluster Comparison Object
#'
#' An S4 class to represent the comparison of two clusters of samples to
#' infer LR interactions based on the resulting P-values,
#' log-fold-changes (logFC), and expression values.
#'
#' @slot col.clusterA   Column indices for the samples in cluster A.
#' @slot col.clusterB   Column indices for the samples in cluster B.
#' @slot differential.stats  Comparison statistics A versus B 
#' as a data.frame and
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
        col.clusterA = "integer",
        col.clusterB = "integer",
        differential.stats = "data.frame"
    ),
    prototype = list(
        col.clusterA = as.integer(c(1, 2)),
        col.clusterB = as.integer(c(3, 4)),
        differential.stats = data.frame(
            pval = c(0.01, 0.01),
            logFC = c(1, -1), expr = c(1, 2)
        )
    )
)

setValidity(
    "BSRClusterComp",
    function(object) {
        if (!is.integer(object@col.clusterA)) {
            return("col.clusterA indices are not all integers")
        }
        if (length(object@col.clusterA) < 1) {
            return("col.clusterA empty")
        }
        if (!is.integer(object@col.clusterB)) {
            return("col.clusterB indices are not all integers")
        }
        if (length(object@col.clusterB) < 1) {
            return("col.clusterB empty")
        }
        if (length(intersect(object@col.clusterA, 
        object@col.clusterB)) > 0) {
            return("col.cluster1 and col.clusterB are not disjoint")
        }
        if (!is.data.frame(object@differential.stats)) {
            return("specified stats are not a data.frame")
        }

        TRUE
    }
)

setMethod(
    "show", "BSRClusterComp",
    function(object) {
        if (length(object@col.clusterA) > 5) {
            cat("Cluster A columns:", 
            object@col.clusterA[seq_len(5)], "...\n")
        } else {
            cat(
                "Cluster A columns:",
                object@col.clusterA[
                    seq_len(length(object@col.clusterA))], "\n"
            )
        }
        if (length(object@col.clusterB) > 5) {
            cat("Cluster B columns:", object@col.clusterB[seq_len(5)], "...\n")
        } else {
            cat(
                "Cluster B columns:",
                object@col.clusterB[
                    seq_len(length(object@col.clusterB))], "\n"
            )
        }
        print(utils::head(object@differential.stats))
    }
)


# Accessors & setters ========================================================

setGeneric("colClusterA", signature="x",
    function(x) standardGeneric("colClusterA")
)
#' Cluster A columns accessor
#'
#' @name colClusterA
#' @aliases colClusterA,BSRClusterComp-method
#' @param x object BSRClusterComp
#' @return col.clusterA
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
#' colClusterA(bsrcc)
#' @export
setMethod("colClusterA", "BSRClusterComp", function(x) x@col.clusterA)

setGeneric("colClusterA<-", signature=c("x", "value"),
    function(x, value) standardGeneric("colClusterA<-")
)
#' Cluster A columns setter (internal use only)
#'
#' @param x object BSRClusterComp
#' @param value value to be set for BSRClusterComp
#' @return returns \code{NULL}
#' @keywords internal
setMethod("colClusterA<-", "BSRClusterComp", function(x, value) {
    x@col.clusterA <- value
    methods::validObject(x)
    x
})

setGeneric("colClusterB", signature="x",
    function(x) standardGeneric("colClusterB")
)
#' Cluster B columns accessor
#'
#' @name colClusterB
#' @aliases colClusterB,BSRClusterComp-method
#' @param x object BSRClusterComp
#' @return col.clusterB
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
#' colClusterB(bsrcc)
#' @export
setMethod("colClusterB", "BSRClusterComp", function(x) x@col.clusterB)

setGeneric("colClusterB<-", signature=c("x", "value"),
    function(x, value) standardGeneric("colClusterB<-")
)
#' Cluster B columns setter (internal use only)
#'
#' @param x object BSRClusterComp
#' @param value value to be set for BSRClusterComp
#' @return returns \code{NULL}
#' @keywords internal
setMethod("colClusterB<-", "BSRClusterComp", function(x, value) {
    x@col.clusterB <- value
    methods::validObject(x)
    x
})


setGeneric("differentialStats", signature="x",
    function(x) standardGeneric("differentialStats")
)
#' Cluster comparison statistics accessor
#'
#' @name differentialStats
#' @aliases differentialStats,BSRClusterComp-method
#' @param x BSRClusterComp object
#' @return diffferential.stats
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
#' differentialStats(bsrcc)
#' @export
setMethod("differentialStats", "BSRClusterComp", function(x) x@differential.stats)

setGeneric("differentialStats<-", signature=c("x", "value"),
    function(x, value) standardGeneric("differentialStats<-")
)
#' Cluster comparison statistics setter (internal use only)
#'
#' @param x object BSRClusterComp
#' @param value value to be set for BSRClusterComp
#' @return returns \code{NULL}
#' @keywords internal
setMethod("differentialStats<-", "BSRClusterComp", function(x, value) {
    x@differential.stats <- value
    methods::validObject(x)
    x
})
