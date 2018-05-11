#time cut off
#' @export
#'not working
#'no errors
#'not changing velocity data in adp object
#'
#'
#'


limit_time <- function(t_first, t_final){
  t1 <- as.POSIXct(t_first)
  t1 <- as.numeric(t1)
  t <- adp[['time', "numeric"]]
  t <- as.numeric(t, tz = 'UTC')
  adp[['v']][t < t1] <- NA
  t2 <- as.POSIXct(t_final)
  t2 <- as.numeric(t2)
  adp[['v']][t > t2] <- NA
}
