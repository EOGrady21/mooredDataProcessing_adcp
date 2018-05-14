###depth

#' @title Limit Depth ADCP PRocessing Step 3.1
#'
#'@description Functions to limit the bounds of ADCP depth data
#'
#'Can be limited by rmax or time (leaving both for now, one may be deleted in final versions)
#'
#'
#'


#'@title limit depth by rmax
#'
#'@description Use maximum acceptable range values to determine acceptable depth values
#'Uses Teledyne RDI equation, Rmax = Dcos(x)
#'
#'@param x oce object of class adp to be limited
#'@param lat latitude of instrument during sampling
#'


limit_depthbyrmax <- function(x, lat = x[['latitude']]){
  if (!inherits(x, "adp")){
    stop("method is only for objects of class '", "adp", "'")
  }
  d <- x[['depth']]
  deg2rad <- function(deg) {(deg * pi) / (180)}
  rmax <- d *(cos(deg2rad(x@metadata$beamAngle)))


  if(!missing(lat)){
  d <- swDepth(x@data$pressure, latitude = lat, eos = getOption("oceEOS", default = "gsw"))
  d[d < rmax] <- NA
  mdt <- round(mean(d), digits = 2)
  x@metadata$sensorDepth <- mdt
  x@data$depth <- d

  }
  if (missing(lat)){
    if(!is.na(x@metadata$latitude)){
      lat <- x@metadata$latitude
      d <- swDepth(x@data$pressure, latitude = lat, eos = getOption("oceEOS", default = "gsw"))

      d[d < rmax] <- NA
       mdt <- round(mean(x@data$depth, na.rm = TRUE), digits = 2)
       x@metadata$sensorDepth <- mdt
       x@data$depth <- d
    }
    if (is.na(x@metadata$latitude)){
      stop()
    warning('No latitude provided; returning object as is')
    }
    adp@processingLog$time <-processingLogAppend(adp@processingLog, date() )
    adp@processingLog <- processingLogAppend(adp@processingLog, paste0('depth limited by maximum acceptable distance, calulated with Rmax = Dcos(x)'))
    adp@processingLog <- processingLogAppend(adp@processingLog, paste0('Sensor depth and mean depth set to  ', mdt , '  based on trimmed depth values'))

    return(x)
  }
}



#'@title limit depth by time
#'
#'@description Uses deployment and recovery times to limit depth within times that
#'insturment was actively and properly sampliing
#'
#'@param adp oce object of class adp to be limited
#'
#'requires certain meta data features to compute
#'including pressure, latitude, time, deploymentTime, recoveryTime
#'


limit_depthbytime <- function(adp){
  if (!inherits(adp, "adp")){
    stop("method is only for objects of class '", "adp", "'")
  }
  adp[['depth']] <- swDepth(adp[['pressure']], latitude = adp[['latitude']], eos = getOption("oceEOS", default = "gsw"))
  depth <- adp[['depth']]
  depth[as.POSIXct(adp[['time']], tz) <= as.POSIXct(adp[['deploymentTime']], tz) | as.POSIXct(adp[['time']], tz) >= as.POSIXct(adp[['recoveryTime']], tz)] <- NA

 mdt <- round(mean(depth, na.rm = TRUE), digits = 2)
 adp@metadata$sensorDepth <- mdt
 adp@metadata$depthMean <- mdt
  adp@data$depth <- depth
  adp@processingLog$time <-processingLogAppend(adp@processingLog, date() )
  adp@processingLog <- processingLogAppend(adp@processingLog, paste0('depth limited by deployment (', adp[['deploymentTime']], ') and recovery  (', adp[['recoveryTime']], ')  times'))
  adp@processingLog <- processingLogAppend(adp@processingLog, paste0('Sensor depth and mean depth set to  ', mdt , '  based on trimmed depth values'))

  return(adp)
}
