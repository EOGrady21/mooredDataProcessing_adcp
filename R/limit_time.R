#time cut off
#'
#'Limit Time ADCP Processing step 3.2
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


limit_time <- function(x, tz = 'UTC', dt = adp[['deploymentTime']], rt = adp[['recoveryTime']]){

if (!missing(dt)){
  t1 <- as.POSIXct(dt, tz = tz)
  t <- adp[['time', "numeric"]]
  t <- as.numeric(t, tz = tz)
  adp[['v']][t < t1] <- NA
}
  else if (missing(dt)){
    warning('No deployment time provided ; untrimmed')}
if( !missing(rt)){
  t2 <- as.POSIXct(rt)
  adp[['v']][t > t2] <- NA
}
  else if (missing(rt)){
    warning('No recovery time provided ; untrimmed')
  }


  return(adp)
}
