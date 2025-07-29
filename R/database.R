#' Fetch the LR reference database from internet
#'
#' Fetch LR database from remote location.
#'
#' @param onRequest logical TRUE to force
#' downloading again. This will overwrite the
#' pre-existing local database. Default is TRUE.
#' @param verbose Logical
#' @return Returns `NULL`, invisibly. 
#' 
#' @import httr
#' @importFrom curl has_internet
#' @importFrom cli cli_alert_danger cli_alert
#' @export
#' @examples
#' print("Function already called elsewhere by cacheClear()")
#' # createDatabase(onRequest = FALSE)
createDatabase <- function(onRequest = TRUE, verbose = FALSE) {
    # Default directory
    cacheDir <- .SignalR$BulkSignalR_CACHEDIR
    databaseCacheDir <- paste(cacheDir, "database", sep = "/")
    url <- .SignalR$BulkSignalR_DB_URL
    # databaseFilePath <- paste(databaseCacheDir
    #    ,basename(url)
    #    ,sep = "/")
    hasInternet <- tryCatch(expr={curl::has_internet()}, 
        error = FALSE)

    if (!hasInternet & 
        !file.exists(databaseCacheDir)) {
        cli::cli_alert_danger("Your internet connection is off:")
        stop(
        "- Remote database cannot be downloaded."
        )   
    }

    if (!hasInternet & 
        onRequest) {
        cli::cli_alert_danger("Your internet connection is off:")
        stop(
        "- Remote database cannot be downloaded.\n"
        )   
    }

    if (!file.exists(databaseCacheDir) | onRequest) {
        # isDownloaded <- .downloadDatabase(url,databaseFilePath)
        # if(!isDownloaded)
        # stop("Ligand-Receptor database was not downloaded successfully.")
        #cacheVersion()
        if(hasInternet)
            .cacheAdd(
                fpath = url,
                cacheDir = databaseCacheDir,
                resourceName = basename(url),
                verbose = verbose,
                download = TRUE
            )


    }

    cacheVersion(dir="database")

    return(invisible(NULL))
}
