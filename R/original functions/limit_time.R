#time cut off
#'
#' ADCP Processing step 3.4
#'
#'
#'@description Function limits time from before deployment and after recovery of ADCP
#'
#'@param x adp object from oce-class adp
#'@param tz time zone, default is 'UTC'
#'@param dt deployment time of ADCP, default pulls value from metadata
#'@param rt recovery time of ADCP, default pulls value from metadata
#'
#'@return adp object with velocities limited to active ADCP measurement times (after deployment and before recovery)
#'
#'
#'


limit_time <- function(x, tz = 'UTC', dt = x[['deploymentTime']], rt = x[['recoveryTime']]){

  if (!inherits(x, "adp")){
    stop("method is only for objects of class '", "adp", "'")
  }
        #FIX ME : if deployment or recovery time is out of bounds

if(!missing(dt)){
  t1 <- as.POSIXct(dt, tz = tz)
  t <- x[['time', "numeric"]]
  t <- as.POSIXct(t, tz = tz)
  x[['v']][t < t1] <- NA
}
  else if(missing(dt)){
    if (!is.null(x@metadata$deploymentTime)){
    dt <- x@metadata$deploymentTime
    t1 <- as.POSIXct(dt, tz = tz)
    t <- x[['time', "numeric"]]
    t <- as.POSIXct(t, tz = tz)
    x[['v']][t < t1] <- NA
    }
    if (is.null(x@metadata$deploymentTime)){
      warning('No deployment Time provided!')
    }


if(!missing(rt)){
 t2 <- as.POSIXct(rt)
 x[['v']][t > t2] <- NA
}
    else if (missing(rt))
      if (!is.null(x@metadata$recoveryTime)){
        t2 <- as.POSIXct(rt)
        x[['v']][t > t2] <- NA
      }
    if (is.null(x@metadata$recoveryTime)){
      warning('No recovery time provided!')
    }


  return(x)
  }
}
