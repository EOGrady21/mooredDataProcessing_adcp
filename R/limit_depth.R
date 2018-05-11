###depth

#' @title Limit Depth ADCP PRocessing Step 3.1
#'
#'
#'


limit_depth <- function(x, lat = x[['latitude']]){
  if (!inherits(x, "adp")){
    stop("method is only for objects of class '", "adp", "'")
}
  deg2rad <- function(deg) {(deg * pi) / (180)}
  rmax <- x@data$depth *(cos(deg2rad(x@metadata$beamAngle)))
  d <- x[['depth']]

  if(!missing(lat)){
  d <- swDepth(x@data$pressure, latitude = lat, eos = getOption("oceEOS", default = "gsw"))
  x@data$depth[d < rmax] <- NA
  mdt <- mean(d)
  x <- oceSetMetadata(x, 'sensorDepth', mdt )

  }
  if (missing(lat)){
    if(!is.na(x@metadata$latitude)){
      lat <- x@metadata$latitude
      d <- swDepth(x@data$pressure, latitude = lat, eos = getOption("oceEOS", default = "gsw"))
       x@data$depth[d < rmax] <- NA
       mdt <- mean(d)
       x <- oceSetMetadata(x, 'sensorDepth', mdt )
    }
    if (is.na(x@metadata$latitude)){
    warning('No latitude provided; returning object as is')
   }
    return(x)
  }
}

