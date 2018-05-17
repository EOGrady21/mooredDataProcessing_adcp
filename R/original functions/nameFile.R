##creating standard format file names

#' @title File naming format
#'
#' @description pulls metadata from an adp oce object and creates a file name in standard BIO format
#'
#' @param adp an oce class - adp object
#'
#' required metadata : cruise number, mooring number, serial number, sampling interval!! Make sure all metadata is present when experiencing problems
#'
#'
#'
name.file <- function(adp){
  name <- paste('MADCP', adp[['cruiseNumber']], adp[['mooringNumber']], adp[['serialNumber']], adp[['samplingInterval']], sep = '_')
  name
}
