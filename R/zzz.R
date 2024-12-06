#nameEnv <- ".SignalR-Env"
.SignalR<- new.env(parent = emptyenv()) 

.onLoad <- function(...) {
    
    # handle directory creation over different OS
    cacheDir <- tools::R_user_dir("BulkSignalR", which="cache")
    assign("BulkSignalR_CACHEDIR", cacheDir, envir = .SignalR)

    url <- "https://partage-dev.montp.inserm.fr:9192"
    urlDatabase <- paste0(url,
        "/CBSB/SignalR/database/SignalR.db")
    assign("BulkSignalR_DB_URL", urlDatabase, envir = .SignalR)

    createDatabase(onRequest = FALSE)

    BulkSignalR_LRdb <- getInteractions()

    assign("BulkSignalR_LRdb", BulkSignalR_LRdb, envir = .SignalR)

    ################################
    ##   Resource Cache Files   ###
    ################################
    urlGo <- paste0(url,
        "/CBSB/SignalR/resources/gobp.rds")
    urlReactome <- paste0(url,
        "/CBSB/SignalR/resources/reactome.rds")
    urlNetwork <- paste0(url,
        "/CBSB/SignalR/resources/Network.rds")

    assign("BulkSignalR_GO_URL", urlGo, 
    envir = .SignalR)
    assign("BulkSignalR_Reactome_URL", urlReactome, 
    envir = .SignalR)
    assign("BulkSignalR_Network_URL", urlNetwork, 
    envir = .SignalR)

    createResources(onRequest = FALSE)

    BulkSignalR_Reactome <- getResource(resourceName = "Reactome",
        cache = TRUE)
    BulkSignalR_Gobp <- getResource(resourceName = "GO-BP",
        cache = TRUE)
    BulkSignalR_Network <- getResource(resourceName = "Network",
        cache = TRUE)

    assign("BulkSignalR_Reactome", BulkSignalR_Reactome, 
    envir = .SignalR)
    assign("BulkSignalR_Gobp", BulkSignalR_Gobp, 
    envir = .SignalR)
    assign("BulkSignalR_Network", BulkSignalR_Network, 
    envir = .SignalR)

}
