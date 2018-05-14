####create flag####

#' ADCP processing Step 3.8
#'
#'
#'
#' adpFlag function
#' flag an adp object based on a series of parameters including percent good, and error velocity.
#'
#' This function also defaults to flagging depth values based on the Teledyne RDI standard
#'
#'  Rmax = Dcosx
#'
#'    where Rmax is the maximum acceptable distance range from the ADCP, D is total depth and x is the beam angle of the ADCP.
#' @param adp, an adp object, oce-class
#' @param flagScheme, scheme of flags that will be followed, BODC, MEDS, etc
#' @param pg, The minimum percent good for evaluating beams one + four of the adcp (BIO standard is 25)
#' @param er, Maximum error velocity for evaluating adcp data (BIO standard is 0.46)
#'

#
#'

####adcp process function


adpFlag <- function(adp, flagScheme, pg, er){
  require(oce)
  if (!inherits(adp, "adp")){
    stop("method is only for objects of class '", "adp", "'")
}
  deg2rad <- function(deg) {(deg * pi) / (180)}
  rmax <- adp[['depth']] *(cos(deg2rad(adp[['beamAngle']])))

  #create matrix of maximum acceptable depth
  r <- matrix(rmax, ncol=length(adp[['distance']]), nrow=length(rmax))


  #create matrix of distance (adp to surface)
  dist <- adp[['distance']]
  t <- adp[['time']]
  d <- t(matrix(dist, ncol = length(t), nrow = length(dist)))

  #read in pg per beam
  g <- adp[['g', "numeric"]]

  #combine beam 1 and 4
  lowpg <- g[,,1]+g[,,4]

  #extract error velocities
  ERRV <- adp[['v']][,,4]

  #create logical array of flagged values based on low percent good, high error velocity or surface contamination
  dim = dim(adp[['v']])
  flag <- array(FALSE, dim= dim)
  for (i in 1:3)
    flag[,,i] <- (lowpg<pg) | (abs(ERRV) > er) | r < d

  #initialize and set flag scheme
  adp <- initializeFlags(adp, name = 'v', value = 0)
if(flagScheme == 'BODC'){
  #set adp flags where logical array = TRUE, flag value = 4 (see flag scheme, BODC)
  adp <- setFlags(adp, name = 'v', i= flag, value = 4)
}

if(flagScheme == 'MEDS'){
  adp <- setFlags(adp, name = 'v', i= flag, value = 4)
}

  #adp@processingLog$time <- processingLogAppend(adp@processingLog, date())
  adp@processingLog <- processingLogAppend(paste0(adp@processingLog, 'Quality control flags set based on flag scheme  ', flagScheme))
  return(adp)
}
