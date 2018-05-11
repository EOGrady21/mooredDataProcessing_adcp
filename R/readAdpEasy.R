##read.adp.easy

#' @title ADCP PRocessing step 2.0
#'
#' @description Load adp data into R with list that includes all metadata from mooring sheets
#'
#'
#' searches working directory for .000 files, reads into adp objects and adds metadata supplied by list
#' @param metadata list specifying all necessary metadata from log sheet
#'
#'

read.adp.easy <- function(metadata){
  raw <- list.files(path = ".", pattern = "*.000")
  if (raw == 0){
    # stop()
    # warning('No raw ADCP files found')
  }
  adp <- read.adp(raw, latitude = metadata[['latitude']], longitude =metadata[['longitude']] ) #insert lat and lon from mooring logs
  if (missing(metadata)){
    warning('no metadata supplied')
  }
  if (!missing(metadata)) {
    for (m in seq_along(metadata)) {
      adp <- oceSetMetadata(adp, names(metadata)[m], metadata[[m]])
    }
    return(adp)
  }
}
