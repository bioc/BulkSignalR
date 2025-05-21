#' Fetch the database from internet.
#'
#' Fetch LR database from remote location.
#'
#' @param onRequest logical True if you force
#' download again. This will overwrite
#' pre-existing database. Default is True.
#' @param verbose Logical TRUE/FALSE
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
        cli::cli_alert_danger("Your internet connection is off :")
        stop(
        "- Remote database can't be downloaded."
        )   
    }

    if (!hasInternet & 
        onRequest) {
        cli::cli_alert_danger("Your internet connection is off :")
        stop(
        "- Remote database can't be downloaded.\n"
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