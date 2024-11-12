

test_that("prepareDataset", {
  data(sdc, package = "BulkSignalR")
  bsrdm <- prepareDataset(sdc)
  expect_s4_class(bsrdm,"BSRDataModel")
  
})


#prepareDataset <- function(
#counts, normalize = TRUE, symbol.col = NULL, min.count = 10,
#prop = 0.1, method = c("UQ", "TC"), 
#log.transformed = FALSE, min.LR.found = 80,
#species = "hsapiens", conversion.dict = NULL,
#UQ.pc = 0.75) 