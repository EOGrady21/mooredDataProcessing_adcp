##read metadata

#' @title ADCP process step 2.1
#'
#'
#' @description Read in all metadata from csv template to adp object, sourced from log sheets
#'
#' @param file csv file name
#' @param obj adp oce object to assign metadata to
#'
#'
#'

read.meta <- function(file, obj){
  md <- read.csv(file, header = TRUE)
  mn <- as.character(names(md[1]))
  mv <- as.character(names(md[2]))
  meta <- as.list(mv)
  names(meta) <- mn

  for (m in seq_along(meta)) {
    obj <- oceSetMetadata(obj, names(meta)[m], meta[[m]])
  }
  return(obj)
}
