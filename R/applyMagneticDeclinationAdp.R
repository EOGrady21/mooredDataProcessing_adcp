###applyMAgneticDeclination

#' @title ADCP Processing step 3.8
#'
#' @description apply magnetic declination to ADCP data
#'
#'
#' @param x adp object from oce-class
#' @param lat latitude
#' @param lon longitude
#' @param st start/ deployment time
#' @param et end/ recovery time
#' @param type average or interpolated
#'
#'
#'
#' @return adp object with magnetic declination applied given coordinates
#'
#' If type == average, an average is taken from start and end time declinations and applied uniformly
#'
#' If type == interpolated, the rate of declination is used over time series
#'
#'


applyMagneticDeclinationAdp <- function(x, lat = x[['latitude']], lon = x[['longitude']], st = x[['deploymentTime']], et = x[['recoveryTime']], type = 'average'){
  if (type =='average'){
    #ifs for different time formats? tz argument?
    if (!is.na(lat) & !is.na(lon)){
      if(!is.null(st) & !is.null(et)){
    s <- as.POSIXct(st, tz = 'UTC')
    e <- as.POSIXct(et, tz = 'UTC')
    a <- magneticField(lon, lat, s)
    b <- magneticField(lon, lat, e)
    c <- round(mean(c(a$declination, b$declination)),digits = 2)
    x <- enuToOther(x, heading = c)
    x <- oceSetMetadata(x, 'magneticVariation', c)
    x <- oceSetMetadata(x, 'oceCoordinate', 'enu')
      }
    }
    else {
      warning('Missing required arguments! No processing performed!')
    }
    #x$processingLog <- processingLogAppend(x$processingLog, value = paste0('magnetic variation applied; declination =', c, 'degrees') )
  }
    if (type == 'interpolate'){
      ;
    }
  return(x)
  return(x@metadata$magneticVariation)
  ##still not returning metadata updates
}

