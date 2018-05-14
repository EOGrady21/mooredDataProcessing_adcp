##read.adp.easy

#' @title ADCP PRocessing step 2.0
#'
#' @description Load adp data into R with list that includes all metadata from mooring sheets
#'
#'
#' @param file raw ADCP file (.000 format)
#' @param metadata csv metadata file from template
#'
#'

read.adp.easy <- function(file, metadata){
  if (missing(metadata)){
    warning('no metadata supplied')
  }
  metad <- read.csv(metadata, header = TRUE)

  mn <- as.character(metad$Name)
  mv <- as.character(metad$Value)


  md <- as.list(mv)
  names(md) <- mn

  adp <- read.adp(file, latitude = md[['latitude']], longitude =md[['longitude']] ) #insert lat and lon from mooring logs

  if (!missing(md)) {
    for (m in seq_along(md)) {
      adp <- oceSetMetadata(adp, names(md)[m], md[[m]])
    }
    return(adp)
  }
}

