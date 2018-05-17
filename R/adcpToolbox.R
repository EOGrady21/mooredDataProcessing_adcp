
require(oce)
require(ncdf4)



##read.adp.easy

#'  ADCP PRocessing step 2.0
#'
#'  Load adp data into R with list that includes all metadata from mooring sheets
#'
#'
#' @param file raw ADCP file (.000 format)
#' @param metadata csv metadata file from template
#'
#' returns an object of class adp (from oce package)
#' uses \code{\link[MASS:oce]{read.adp}}

read.adp.easy <- function(file, metadata){
  if (missing(metadata)){
    warning('no metadata supplied')
  }
  metad <- read.csv(metadata, header = TRUE)

  mn <- as.character(metad[,1])
  mv <- as.character(metad[,2])


  md <- as.list(mv)
  names(md) <- mn

  adp <- read.adp(file, latitude = md[['latitude']], longitude =md[['longitude']] ) #insert lat and lon from mooring logs

  if (!missing(md)) {
    for (m in seq_along(md)) {
      adp <- oceSetMetadata(adp, names(md)[m], md[[m]], note = NULL)
    }
    adp@metadata$latitude <- as.numeric(adp[['latitude']])
    adp@metadata$longitude <- as.numeric(adp[['longitude']])
    adp@processingLog$time <-processingLogAppend(adp@processingLog, date() )
    adp@processingLog <- processingLogAppend(adp@processingLog, 'metadata read in from log sheet')

    return(adp)
  }
}


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
  mn <- as.character(md[,1])
  mv <- as.character(md[,2])
  meta <- as.list(mv)
  names(meta) <- mn

  for (m in seq_along(meta)) {
    obj <- oceSetMetadata(obj, names(meta)[m], meta[[m]], note = NULL)
  }
  return(obj)
}



###applyMAgneticDeclination

#' @title ADCP Processing step 3.2
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


applyMagneticDeclinationAdp <- function(x, lat = x[['latitude']], lon = x[['longitude']], st = x[['deploymentTime']], et = x[['recoveryTime']],tz = 'UTC', type = 'average'){
  if (!inherits(x, "adp")){
    stop("method is only for objects of class '", "adp", "'")
  }
  if (is.na(lat) | is.na(lon)){
    warning('No latitude or longitude values provided!')
  }
  if (is.null(st) | is.null(et)){
    warning('deployment and recovery times not provided!')
  }
  if (type =='average'){
    #ifs for different time formats? tz argument?
    if (!is.na(lat) & !is.na(lon)){
      if(!is.null(st) & !is.null(et)){
        s <- as.POSIXct(st, tz = tz)
        e <- as.POSIXct(et, tz = tz)
        a <- magneticField(lon, lat, s)
        b <- magneticField(lon, lat, e)
        c <- round(mean(c(a$declination, b$declination)),digits = 2)
        coord <- x@metadata$oceCoordinate


        if (coord == 'enu'){
          x <- enuToOther(x, heading = c)
          x@metadata$magneticVariation <- c
          x@metadata$oceCoordinate <- 'enu'
        }
        if (coord != 'enu'){
          warning('Function cannot handle objects in ', coord, 'Object returned as is ; please convert to enu first')

        }
      }
      else {
        warning('Missing required arguments! No processing performed!')
      }
      x@processingLog$time <- processingLogAppend(x@processingLog, date())
      x@processingLog <- processingLogAppend(x@processingLog, value = paste0('magnetic variation applied; declination =', c, 'degrees') )
    }
    if (type == 'interpolate'){
      ;
    }
    return(x)


  }
}



#'@title ADCP Processing step 3.3
#'Limit depth by rmax
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
      warning('No latitude provided; returning object as is')

      stop()
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


limit_depthbytime <- function(adp, tz = 'UTC'){
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

####create flag####

#' ADCP processing Step 3.5
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


adpFlag <- function(adp,  pg, er){
  require(oce)
  if (!inherits(x = adp, "adp")){
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
    flag[,,i] <- (lowpg< pg) | (abs(ERRV) > er) | r < d

  #initialize and set flag scheme
  adp <- initializeFlags(adp, name = 'v', value = 0)


  #set adp flags where logical array = TRUE, flag value = 4 (see flag scheme, BODC)
  adp <- setFlags(adp, name = 'v', i= flag, value = 4)


  #adp@processingLog$time <- processingLogAppend(adp@processingLog, date())
  adp@processingLog <- processingLogAppend(adp@processingLog, 'Quality control flags set based on flag scheme from BODC')
  return(adp)
}




#' @title  ADCP Processing step 4.1
#'
#'
#' @description Exports an adp object to a net cdf using variables and metadata within adp combined with optional additional metatdata
#'
#'
#' @param obj an adp object from the oce class
#' @param name name of the NetCDF file to be produced
#' @param metadata csv file listing metadata names and values to be inserted into global attributes of net CDF








##nc metadata
oceNc_create <- function(adp, name,  metadata){
  if (!inherits(adp, "adp")){
    stop("method is only for objects of class '", "adp", "'")
  }
  if (missing(metadata)){
    warning('no metadata supplied')
  }
  #file name and path
  ncpath <- "./"
  ncname <- name
  ncfname <- paste(ncpath, ncname, ".nc", sep = "")

  ####setting dimensions and definitions####
  #dimension variables from adp object
  time <- adp[['time']]
  dist <- adp[['distance', 'numeric']]
  lon <- adp[['longitude']]
  lat <- adp[['latitude']]

  #create dimensions
  timedim <- ncdim_def("time", "POSIXct", as.double(time))    #time formatting FIX
  depthdim <- ncdim_def("depth", "m", as.double(dist))
  londim <- ncdim_def("lon", "degrees_east" , as.double(lon))
  latdim <- ncdim_def("lat", "degrees_north", as.double(lat))

  #set fill value
  FillValue <- 1e35

  #define variables
  dlname <- "eastward_sea_water_velocity"
  u_def <- ncvar_def("EWCT", "m/sec", list(timedim, depthdim, londim, latdim), FillValue, dlname, prec = "float")

  dlname <- "northward_sea_water_velocity"
  v_def <- ncvar_def("NSCT", "m/sec", list(timedim, depthdim, londim, latdim), FillValue, dlname, prec = "float")

  dlname <- "upward_sea_water_velocity"
  w_def <- ncvar_def("VCSP", "m/sec", list(timedim, depthdim, londim, latdim), FillValue, dlname, prec = "float")

  dlname <- "time_02"
  t_def <- ncvar_def("SYTM", "POSIXct", list(timedim, londim, latdim), FillValue, dlname, prec = "float")

  dlname <- "error_velocity_in_sea_water"
  e_def <- ncvar_def("ERRV", "m/sec", list(timedim, depthdim, londim, latdim), FillValue, dlname, prec = "float")

  dlname <- "ADCP_echo_intensity_beam_1"
  b1_def <- ncvar_def("BEAM_01", "a", list(timedim, depthdim, londim, latdim), FillValue, dlname, prec = "float")

  dlname <- "ADCP_echo_intensity_beam_2"
  b2_def <- ncvar_def("BEAM_02", "a", list(timedim, depthdim, londim, latdim), FillValue, dlname, prec = "float")

  dlname <- "ADCP_echo_intensity_beam_3"
  b3_def <- ncvar_def("BEAM_03", "a", list(timedim, depthdim, londim, latdim), FillValue, dlname, prec = "float")

  dlname <- "ADCP_echo_intensity_beam_4"
  b4_def <- ncvar_def("BEAM_04", "a", list(timedim, depthdim, londim, latdim), FillValue, dlname, prec = "float")

  dlname <- "percent_good_beam_1"
  pg1_def <- ncvar_def("PGDP_01", "counts", list(timedim, depthdim, londim, latdim), FillValue, dlname, prec = "float")

  dlname <- "percent_good_beam_2"
  pg2_def <- ncvar_def("PGDP_02", "counts", list(timedim, depthdim, londim, latdim), FillValue, dlname, prec = "float")

  dlname <- "percent_good_beam_3"
  pg3_def <- ncvar_def("PGDP_03", "counts", list(timedim, depthdim, londim, latdim), FillValue, dlname, prec = "float")

  dlname <- "percent_good_beam_4"
  pg4_def <- ncvar_def("PGDP_04", "counts", list(timedim, depthdim, londim, latdim), FillValue, dlname, prec = "float")

  dlname <- "pitch"
  p_def <- ncvar_def("PTCH", "degrees", list(timedim, londim, latdim), FillValue, dlname, prec = "float")

  dlname <- "roll"
  r_def <- ncvar_def("ROLL", "degrees", list(timedim, londim, latdim), FillValue, dlname, prec = "float")

  dlname <- "height of sea surface"
  hght_def <- ncvar_def("hght", "m", list(timedim, londim, latdim), FillValue, dlname, prec = "float")

  dlname <- "ADCP Transducer Temp."
  Tx_def <- ncvar_def("Tx", "degrees", list(timedim, londim, latdim), FillValue, dlname, prec = "float")

  dlname <- "instrument depth"
  D_def <- ncvar_def("D", "m", list(timedim, londim, latdim), FillValue, dlname, prec = "float")


  dlname <- "quality_flag u"
  qc_u_def <- ncvar_def("QC_flag_u", "", list(timedim, depthdim, londim, latdim), FillValue, dlname, prec = "float")

  dlname <- "quality_flag v"
  qc_v_def <- ncvar_def("QC_flag_v", "", list(timedim, depthdim, londim, latdim), FillValue, dlname, prec = "float")

  dlname <- "quality_flag w"
  qc_w_def <- ncvar_def("QC_flag_w", "", list(timedim, depthdim, londim, latdim), FillValue, dlname, prec = "float")

  ####writing net CDF####
  #write out definitions to new nc file
  ncout <- nc_create(ncfname, list(u_def, v_def, w_def, e_def, t_def, b1_def, b2_def, b3_def, b4_def, pg1_def, pg2_def, pg3_def, pg4_def, p_def, r_def, hght_def, Tx_def, D_def, qc_u_def, qc_v_def, qc_w_def), force_v4 = TRUE)

  #insert variables into nc file
  ncvar_put(ncout, u_def, adp[['v']][,,1])
  ncvar_put(ncout, v_def, adp[['v']][,,2])
  ncvar_put(ncout, w_def, adp[['v']][,,3])
  ncvar_put(ncout, e_def, adp[['v']][,,4])
  ncvar_put(ncout, t_def, adp[['time']])
  ncvar_put(ncout, b1_def, adp[['a', 'numeric']][,,1])
  ncvar_put(ncout, b2_def, adp[['a', 'numeric']][,,2])
  ncvar_put(ncout, b3_def, adp[['a', 'numeric']][,,3])
  ncvar_put(ncout, b4_def, adp[['a', 'numeric']][,,4])
  ncvar_put(ncout, pg1_def, adp[['g', 'numeric']][,,1])
  ncvar_put(ncout, pg2_def, adp[['g', 'numeric']][,,2])
  ncvar_put(ncout, pg3_def, adp[['g', 'numeric']][,,3])
  ncvar_put(ncout, pg4_def, adp[['g', 'numeric']][,,4])
  ncvar_put(ncout, p_def, adp[['pitch']])
  ncvar_put(ncout, r_def, adp[['roll']])
  ncvar_put(ncout, hght_def, adp[['depth']])
  ncvar_put(ncout, Tx_def, adp[['temperature']])
  ncvar_put(ncout, D_def, adp[['depth']])
  ncvar_put(ncout, qc_u_def, adp@metadata$flags$v[,,1])
  ncvar_put(ncout, qc_v_def, adp@metadata$flags$v[,,2])
  ncvar_put(ncout, qc_w_def, adp@metadata$flags$v[,,3])
  ###metadata###

  ####pulled from adp object####
  ncatt_put(ncout, 0, "MOORING", adp[['mooringNumber']])
  ncatt_put(ncout, 0, "Deployment_date", adp[['deploymentTime']])
  ncatt_put(ncout, 0, "Recovery_date", adp[['recoveryTime']])
  ncatt_put(ncout, 0, "firmware_version", adp[['firmwareVersion']])
  ncatt_put(ncout, 0, "frequency", adp[['frequency']])
  ncatt_put(ncout, 0, "beam_pattern", adp[['beamPattern']])
  ncatt_put(ncout, 0, "janus", adp[['numberOfBeams']])
  ncatt_put(ncout, 0, "pings_per_ensemble", adp[['pingsPerEnsemble']])
  ncatt_put(ncout, 0, "valid_correlation_range", adp[['lowCorrThresh']])
  ncatt_put(ncout, 0, "minmax_percent_good", adp[['percentGdMinimum']])
  ncatt_put(ncout, 0,"minmax_percent_good", "100")
  ncatt_put(ncout, 0, "error_velocity_threshold", adp[['errorVelocityMaximum']])
  ncatt_put(ncout, 0, "transmit_pulse_length_cm", adp[['xmitPulseLength']])
  ncatt_put(ncout, 0, "false_target_reject_values", adp[['falseTargetThresh']])
  ncatt_put(ncout, 0, "ADCP_serial_number", adp[['serialNumber']])
  ncatt_put(ncout, 0, "transform", adp[['oceCoordinate']])
  ncatt_put(ncout, 0, "DATA_TYPE", adp[['instrumentType']])
  ncatt_put(ncout, 0, "COORD_SYSTEM", adp[['oceCoordinate']])
  ncatt_put(ncout, 0, "longitude", adp[['longitude']])
  ncatt_put(ncout, 0, "latitude", adp[['latitude']])
  ncatt_put(ncout, 0, "magnetic_variation", adp[['magneticVariation']])
  ncatt_put(ncout, 0, "platform", adp[['platform']])
  ncatt_put(ncout, 0, "sounding", adp[['sounding']])
  ncatt_put(ncout, 0, "chief_scientist", adp[['chiefScientist']])
  ncatt_put(ncout, "depth", "xducer_offset_from_bottom", adp[['depthOffBottom']])
  ncatt_put(ncout, "depth", "bin_size", adp[['cellSize']])
  ncatt_put(ncout, "EWCT", "sensor_type", adp[['instrumentType']])
  ncatt_put(ncout, "EWCT", "sensor_depth", adp[['sensorDepth']])
  ncatt_put(ncout, "EWCT", "serial_number", adp[['serialNumber']])
  ncatt_put(ncout, "NSCT", "sensor_type", adp[['instrumentType']])
  ncatt_put(ncout, "NSCT", "sensor_depth", adp[['sensorDepth']])
  ncatt_put(ncout, "NSCT", "serial_number", adp[['serialNumber']])
  ncatt_put(ncout, "VCSP", "sensor_type", adp[['instrumentType']])
  ncatt_put(ncout, "VCSP", "sensor_depth", adp[['sensorDepth']])
  ncatt_put(ncout, "VCSP", "serial_number", adp[['serialNumber']])
  ncatt_put(ncout, "ERRV", "sensor_type", adp[['instrumentType']])
  ncatt_put(ncout, "ERRV", "sensor_depth", adp[['sensorDepth']])
  ncatt_put(ncout, "ERRV", "serial_number", adp[['serialNumber']])
  ncatt_put(ncout, "BEAM_01", "sensor_type", adp[['instrumentType']])
  ncatt_put(ncout, "BEAM_01", "sensor_depth", adp[['sensorDepth']])
  ncatt_put(ncout, "BEAM_01", "serial_number", adp[['serialNumber']])
  ncatt_put(ncout, "BEAM_02", "sensor_type", adp[['instrumentType']])
  ncatt_put(ncout, "BEAM_02", "sensor_depth", adp[['sensorDepth']])
  ncatt_put(ncout, "BEAM_02", "serial_number", adp[['serialNumber']])
  ncatt_put(ncout, "BEAM_03", "sensor_type", adp[['instrumentType']])
  ncatt_put(ncout, "BEAM_03", "sensor_depth", adp[['sensorDepth']])
  ncatt_put(ncout, "BEAM_03", "serial_number", adp[['serialNumber']])
  ncatt_put(ncout, "BEAM_04", "sensor_type", adp[['instrumentType']])
  ncatt_put(ncout, "BEAM_04", "sensor_depth", adp[['depthMean']])
  ncatt_put(ncout, "BEAM_04", "serial_number", adp[['serialNumber']])
  ncatt_put(ncout, "PGDP_01", "sensor_type", adp[['instrumentType']])
  ncatt_put(ncout, "PGDP_01", "sensor_depth", adp[['sensorDepth']])
  ncatt_put(ncout, "PGDP_01", "serial_number", adp[['serialNumber']])
  ncatt_put(ncout, "PGDP_02", "sensor_type", adp[['instrumentType']])
  ncatt_put(ncout, "PGDP_02", "sensor_depth", adp[['sensorDepth']])
  ncatt_put(ncout, "PGDP_02", "serial_number", adp[['serialNumber']])
  ncatt_put(ncout, "PGDP_03", "sensor_type", adp[['instrumentType']])
  ncatt_put(ncout, "PGDP_03", "sensor_depth", adp[['sensorDepth']])
  ncatt_put(ncout, "PGDP_03", "serial_number", adp[['serialNumber']])
  ncatt_put(ncout, "PGDP_04", "sensor_type", adp[['instrumentType']])
  ncatt_put(ncout, "PGDP_04", "sensor_depth", adp[['sensorDepth']])
  ncatt_put(ncout, "PGDP_04", "serial_number", adp[['serialNumber']])
  ncatt_put(ncout, "EWCT", "generic_name", "u")
  ncatt_put(ncout, "NSCT", "generic_name", "v")
  ncatt_put(ncout, "VCSP", "generic_name", "w")
  ncatt_put(ncout, "ERRV", "generic_name", "w")       #issue in current NC protocol
  ncatt_put(ncout, "BEAM_01", "generic_name", "AGC")
  ncatt_put(ncout, "BEAM_02", "generic_name", "AGC")
  ncatt_put(ncout, "BEAM_03", "generic_name", "AGC")
  ncatt_put(ncout, "BEAM_04", "generic_name", "AGC")
  ncatt_put(ncout, "PGDP_01", "generic_name", "PGd")
  ncatt_put(ncout, "PGDP_02", "generic_name", "PGd")
  ncatt_put(ncout, "PGDP_03", "generic_name", "PGd")
  ncatt_put(ncout, "PGDP_04", "generic_name", "PGd")
  ncatt_put(ncout, "hght", "generic_name", "height")
  ncatt_put(ncout, "hght", "sensor_type", adp[['instrumentType']])
  ncatt_put(ncout, "hght", "sensor_depth", adp[['depthMean']])
  ncatt_put(ncout, "hght", "serial_number", adp[['serialNumber']])
  ncatt_put(ncout, "D", "generic_name", "depth")
  ncatt_put(ncout, "D", "sensor_type", adp[['instrumentType']])
  ncatt_put(ncout, "D", "sensor_depth", adp[['depthMean']])
  ncatt_put(ncout, "D", "serial_number", adp[['serialNumber']])
  ncatt_put(ncout, "Tx", "generic_name", "temp")
  ncatt_put(ncout, "Tx", "sensor_type", adp[['instrumentType']])
  ncatt_put(ncout, "Tx", "sensor_depth", adp[['depthMean']])
  ncatt_put(ncout, "Tx", "serial_number", adp[['serialNumber']])
  ncatt_put(ncout, "QC_flag_u", "comment", "Quality flag resulting from quality control")
  ncatt_put(ncout, "QC_flag_u", "flag_meanings",adp[['flagMeaning']])
  ncatt_put(ncout, "QC_flag_u", "flag_values", adp[['flagValues']])
  ncatt_put(ncout, "QC_flag_u", "References", adp[['flagReferences']])
  ncatt_put(ncout, "QC_flag_v", "comment", "Quality flag resulting from quality control")
  ncatt_put(ncout, "QC_flag_v", "flag_meanings", adp[['flagMeaning']])
  ncatt_put(ncout, "QC_flag_v", "flag_values", adp[['flagValues']])
  ncatt_put(ncout, "QC_flag_v", "References", adp[['flagReferences']])
  ncatt_put(ncout, "QC_flag_w", "comment", "Quality flag resulting from quality control")
  ncatt_put(ncout, "QC_flag_w", "flag_meanings", adp[['flagMeaning']])
  ncatt_put(ncout, "QC_flag_w", "flag_values", adp[['flagValues']])
  ncatt_put(ncout, "QC_flag_w", "References", adp[['flagReferences']])
  ncatt_put(ncout, 0, "DATA_ORIGIN", adp[['institution']])
  ncatt_put(ncout, 0, "WATER_DEPTH", adp[['waterDepth']])
  ncatt_put(ncout, 0, "DELTA_T_sec",adp[['TimeOffset']])
  ncatt_put(ncout, 0, "pred_accuracy", adp[['velocityResolution']])
  #CF conventions

  ncatt_put(ncout, 0, 'Conventions', 'CF-1.7')
  ncatt_put(ncout, 0, "creator_type", "person")
  ncatt_put(ncout, 0, "creator_institution", adp[['institution']])
  ncatt_put(ncout, 0, "program", adp[['program']])
  ncatt_put(ncout, 0, "sea_name", adp[['seaName']])
  ncatt_put(ncout, 0, "time_coverage_start", adp[['deploymentTime']])
  ncatt_put(ncout, 0, "time_coverage_end", adp[['recoveryTime']])
  ncatt_put(ncout, 0, "geospatial_lat_min", adp[['latitude']])
  ncatt_put(ncout, 0, "geospatial_lat_max", adp[['latitude']])
  ncatt_put(ncout, 0, "geosptial_lat_units", "degrees_north")
  ncatt_put(ncout, 0, "geospatial_lon_min", adp[['longitude']])
  ncatt_put(ncout, 0, "geosptial_lon_max", adp[['longitude']])
  ncatt_put(ncout, 0, "geosptial_lon_units", "degrees_east")
  ncatt_put(ncout, 0, "geospatial_vertical_min", min(adp[['depth']]))
  ncatt_put(ncout, 0, "geosptial_vertical_max", max(adp[['depth']]))
  ncatt_put(ncout, 0, "geosptial_vertical_units", "metres")
  ncatt_put(ncout, 0, "geosptial_vertical_positive", adp[['orientation']])     #eg up or down
  ncatt_put(ncout, 0, "institution", adp[['institution']])
  ncatt_put(ncout, 0, "creator_name", adp[['creatorName']])
  ncatt_put(ncout, 0, "creator_url", adp[['creatorUrl']])
  ncatt_put(ncout, 0, "creator_email", adp[['creatorEmail']])
  ncatt_put(ncout, 0, "project", adp[['project']])
  ncatt_put(ncout, 0, "processing_level", adp[['processingLevel']])
  ncatt_put(ncout, 0 , "flag_meanings", adp[['flagMeaning']])
  ncatt_put(ncout, 0 , "flag_values", adp[['flagValues']])
  ncatt_put(ncout, 0, "source", "R code: adcpProcess, github:")
  ncatt_put(ncout, 0, "date_modified", date())
  ncatt_put(ncout,0, "_FillValue", "1e35")
  ncatt_put(ncout, 0, "featureType", "timeSeriesProfile") #link to oce object? ..... if adp == timeSeriesProfile
  ncatt_put

  #BODC P01 names
  ncatt_put(ncout, "EWCT", "sdn_parameter_urn", "SDN:P01::LCEWAP01")
  ncatt_put(ncout, "NSCT", "sdn_parameter_urn", "SDN:P01::LCNSAP01")
  ncatt_put(ncout, "VCSP", "sdn_parameter_urn", "SDN:P01::LRZAAP01")
  ncatt_put(ncout, "ERRV", "sdn_parameter_urn", "SDN:P01::LERRAP01")
  ncatt_put(ncout, "BEAM_01", "sdn_parameter_urn", "SDN:P01::ASAMAP00")
  ncatt_put(ncout, "BEAM_02", "sdn_parameter_urn", "SDN:P01::ASAMAP02")
  ncatt_put(ncout, "BEAM_03", "sdn_parameter_urn", "SDN:P01::ASAMAP03")
  ncatt_put(ncout, "BEAM_04", "sdn_parameter_urn", "SDN:P01::ASAMAP04")
  ncatt_put(ncout, "PGDP_01", "sdn_parameter_urn", "SDN:P01::PCGDAP00")
  ncatt_put(ncout, "PGDP_02", "sdn_parameter_urn", "SDN:P01::PCGDAP02")
  ncatt_put(ncout, "PGDP_03", "sdn_parameter_urn", "SDN:P01::PCGDAP03")
  ncatt_put(ncout, "PGDP_04", "sdn_parameter_urn", "SDN:P01::PCGDAP04")
  #ncatt_put(ncout, "hght", "sdn_parameter_urn", "SDN:P01::")
  ncatt_put(ncout, "D", "sdn_parameter_urn", "SDN:P01::ADEPZZ01")
  ncatt_put(ncout, "Tx", "sdn_parameter_urn", "SDN:P01::TEMPPR01")
  #ncatt_put(ncout, "time_02", "sdn_parameter_urn", "SDN:P01::")
  ncatt_put(ncout, "PTCH", "sdn_parameter_urn", "SDN:P01::PTCHEI01")
  ncatt_put(ncout, "ROLL", "sdn_parameter_urn", "SDN:P01::ROLLFEI01")
  ncatt_put(ncout, "lon", "sdn_parameter_urn", "SDN:P01::ALONZZ01")
  ncatt_put(ncout, "lat", "sdn_parameter_urn", "SDN:P01::ALATZZ01")



  ncatt_put(ncout, "EWCT", "sdn_parameter_name", "Eastward current velocity (Eulerian) in the water body by moored acoustic doppler current profiler (ADCP)")
  ncatt_put(ncout, "NSCT", "sdn_parameter_name", "Northward current velocity (Eulerian) in the water body by moored acoustic doppler current profiler (ADCP)")
  ncatt_put(ncout, "VCSP", "sdn_parameter_name", "Upward current velocity in the water body by moored acoustic doppler current profiler (ADCP)")
  ncatt_put(ncout, "ERRV", "sdn_parameter_name", "Current velocity error in the water body by moored acoustic doppler current profiler (ADCP)")
  ncatt_put(ncout, "BEAM_01", "sdn_parameter_name", "Signal return amplitude from the water body by moored acoustic doppler current profiler (ADCP) beam 1")
  ncatt_put(ncout, "BEAM_02", "sdn_parameter_name", "Signal return amplitude from the water body by moored acoustic doppler current profiler (ADCP) beam 2")
  ncatt_put(ncout, "BEAM_03", "sdn_parameter_name", "Signal return amplitude from the water body by moored acoustic doppler current profiler (ADCP) beam 3")
  ncatt_put(ncout, "BEAM_04", "sdn_parameter_name", "Signal return amplitude from the water body by moored acoustic doppler current profiler (ADCP) beam 4")
  ncatt_put(ncout, "PGDP_01", "sdn_parameter_name", "Acceptable proportion of signal returns by moored acoustic doppler current profiler (ADCP) beam 1")
  ncatt_put(ncout, "PGDP_02", "sdn_parameter_name", "Acceptable proportion of signal returns by moored acoustic doppler current profiler (ADCP) beam 2")
  ncatt_put(ncout, "PGDP_03", "sdn_parameter_name", "Acceptable proportion of signal returns by moored acoustic doppler current profiler (ADCP) beam 3")
  ncatt_put(ncout, "PGDP_04", "sdn_parameter_name", "Acceptable proportion of signal returns by moored acoustic doppler current profiler (ADCP) beam 4")
  ncatt_put(ncout, "D", "sdn_parameter_name", "Depth below surface of the water body")
  ncatt_put(ncout, "Tx", "sdn_parameter_name", "Temperature of the water body")
  ncatt_put(ncout, "PTCH", "sdn_parameter_name", "Orientation (pitch) of measurement platform by inclinometer")
  ncatt_put(ncout, "ROLL", "sdn_parameter_name", "Orientation (roll angle) of measurement platform by inclinometer (second sensor)")
  ncatt_put(ncout, "lon", "sdn_parameter_name", "Longitude east")
  ncatt_put(ncout, "lat", "sdn_parameter_name", "Latitude north")




  ncatt_put(ncout, "EWCT", "sdn_uom_urn", "SDN:P06::UVAA")
  ncatt_put(ncout, "NSCT", "sdn_uom_urn", "SDN:P06::UVAA")
  ncatt_put(ncout, "VCSP", "sdn_uom_urn", "SDN:P06::UVAA")
  ncatt_put(ncout, "ERRV", "sdn_uom_urn", "SDN:P06::UVAA")
  ncatt_put(ncout, "BEAM_01", "sdn_uom_urn", "SDN:P06::UCNT")
  ncatt_put(ncout, "BEAM_02", "sdn_uom_urn", "SDN:P06::UCNT")
  ncatt_put(ncout, "BEAM_03", "sdn_uom_urn", "SDN:P06::UCNT")
  ncatt_put(ncout, "BEAM_04", "sdn_uom_urn", "SDN:P06::UCNT")
  ncatt_put(ncout, "PGDP_01", "sdn_uom_urn", "SDN:P06::UPCT")
  ncatt_put(ncout, "PGDP_02", "sdn_uom_urn", "SDN:P06::UPCT")
  ncatt_put(ncout, "PGDP_03", "sdn_uom_urn", "SDN:P06::UPCT")
  ncatt_put(ncout, "PGDP_04", "sdn_uom_urn", "SDN:P06::UPCT")
  ncatt_put(ncout, "hght", "sdn_uom_urn", "SDN:P06::ULAA")
  ncatt_put(ncout, "D", "sdn_uom_urn", "SDN:P06:ULAA")
  ncatt_put(ncout, "Tx", "sdn_uom_urn", "SDN:P06::UPAA")
  ncatt_put(ncout, "PTCH", "sdn_uom_urn", "SDN:P06:UAAA")
  ncatt_put(ncout, "ROLL", "sdn_uom_urn", "SDN:P06:UAAA")
  ncatt_put(ncout, "lon", "sdn_uom_urn", "SDN:P06::DEGE")
  ncatt_put(ncout, "lat", "sdn_uom_urn", "SDN:P06:DEGN")


  ncatt_put(ncout, "EWCT", "sdn_uom_name", "Metres per second")
  ncatt_put(ncout, "NSCT", "sdn_uom_name", "Metres per second")
  ncatt_put(ncout, "VCSP", "sdn_uom_name", "Metres per second")
  ncatt_put(ncout, "ERRV", "sdn_uom_name", "Metres per second")
  ncatt_put(ncout, "BEAM_01", "sdn_uom_name", "Counts")
  ncatt_put(ncout, "BEAM_02", "sdn_uom_name", "Counts")
  ncatt_put(ncout, "BEAM_03", "sdn_uom_name", "Counts")
  ncatt_put(ncout, "BEAM_04", "sdn_uom_name", "Counts")
  ncatt_put(ncout, "PGDP_01", "sdn_uom_name", "Percent")
  ncatt_put(ncout, "PGDP_02", "sdn_uom_name", "Percent")
  ncatt_put(ncout, "PGDP_03", "sdn_uom_name", "Percent")
  ncatt_put(ncout, "PGDP_04", "sdn_uom_name", "Percent")
  ncatt_put(ncout, "hght", "sdn_uom_name", "Metres")
  ncatt_put(ncout, "D", "sdn_uom_name", "Metres")
  ncatt_put(ncout, "Tx", "sdn_uom_name", "Celsius degree")
  ncatt_put(ncout, "PTCH", "sdn_uom_name", "Degrees")
  ncatt_put(ncout, "ROLL", "sdn_uom_name", "Degrees")
  ncatt_put(ncout, "lon", "sdn_uom_name", "Degrees east")
  ncatt_put(ncout, "lat", "sdn_uom_name", "Degrees north")


  #CF standard names
  ncatt_put(ncout, "EWCT", "standard_name", "eastward_sea_water_velocity")
  ncatt_put(ncout, "NSCT", "standard_name", "northward_sea_water_velocity")
  ncatt_put(ncout, "VCSP", "standard_name", "upward_sea_water_velocity")


  ncatt_put(ncout, "lat", "standard_name", "latitude")
  ncatt_put(ncout, "lon", "standard_name", "longitude")
  ncatt_put(ncout, "D", "standard_name", "depth")
  ncatt_put(ncout, "depth", "poitive", "down")     #direction of depth axis
  ncatt_put(ncout, "depth", "axis", "y")
  ncatt_put(ncout, "time", "axis", "x")


  ####
  ncatt_put(ncout, "EWCT", "data_max", max(adp[['v']][,,1], na.rm = TRUE))
  ncatt_put(ncout, "EWCT", "data_min", min(adp[['v']][,,1], na.rm = TRUE))
  ncatt_put(ncout, "NSCT", "data_max", max(adp[['v']][,,2], na.rm = TRUE))
  ncatt_put(ncout, "NSCT", "data_min", min(adp[['v']][,,2], na.rm = TRUE))
  ncatt_put(ncout, "VCSP", "data_max", max(adp[['v']][,,3], na.rm = TRUE))
  ncatt_put(ncout, "VCSP", "data_min", min(adp[['v']][,,3], na.rm = TRUE))
  ncatt_put(ncout, "ERRV", "data_max", max(adp[['v']][,,4], na.rm = TRUE))
  ncatt_put(ncout, "ERRV", "data_min", min(adp[['v']][,,4], na.rm = TRUE))
  ncatt_put(ncout, "BEAM_01", "data_min", min(adp[['a', 'numeric']][,,1], na.rm= TRUE))
  ncatt_put(ncout, "BEAM_01", "data_max", max(adp[['a', 'numeric']][,,1], na.rm= TRUE))
  ncatt_put(ncout, "BEAM_02", "data_min", min(adp[['a' ,'numeric']][,,2], na.rm= TRUE))
  ncatt_put(ncout, "BEAM_02", "data_max", max(adp[['a', 'numeric']][,,2], na.rm= TRUE))
  ncatt_put(ncout, "BEAM_03", "data_min", min(adp[['a', 'numeric']][,,3], na.rm= TRUE))
  ncatt_put(ncout, "BEAM_03", "data_max", max(adp[['a', 'numeric']][,,3], na.rm= TRUE))
  ncatt_put(ncout, "BEAM_04", "data_min", min(adp[['a', 'numeric']][,,4], na.rm= TRUE))
  ncatt_put(ncout, "BEAM_04", "data_max", max(adp[['a', 'numeric']][,,4], na.rm= TRUE))
  ncatt_put(ncout, "PGDP_01", "data_min", min(adp[['g', 'numeric']][,,1], na.rm= TRUE))
  ncatt_put(ncout, "PGDP_01", "data_max", max(adp[['g', 'numeric']][,,1], na.rm= TRUE))# eg min 25 % good
  ncatt_put(ncout, "PGDP_02", "data_min", min(adp[['g', 'numeric']][,,2], na.rm= TRUE))
  ncatt_put(ncout, "PGDP_02", "data_max", max(adp[['g' ,'numeric']][,,2], na.rm= TRUE))
  ncatt_put(ncout, "PGDP_03", "data_min", min(adp[['g' ,'numeric']][,,3], na.rm= TRUE))
  ncatt_put(ncout, "PGDP_03", "data_max", max(adp[['g', 'numeric']][,,3], na.rm= TRUE))
  ncatt_put(ncout, "PGDP_04", "data_min", min(adp[['g', 'numeric']][,,4], na.rm= TRUE))
  ncatt_put(ncout, "PGDP_04", "data_max", max(adp[['g', 'numeric']][,,4], na.rm= TRUE))
  ncatt_put(ncout, "hght", "data_min", min(adp[['depth', 'numeric']]))
  ncatt_put(ncout, "hght", "data_max", max(adp[['depth', 'numeric']]))
  ncatt_put(ncout, "D", "data_min", min(adp[['depth']]))
  ncatt_put(ncout, "D", "data_max", max(adp[['depth']]))
  ncatt_put(ncout, "Tx", "data_min", min(adp[['temperature']]))
  ncatt_put(ncout, "Tx", "data_max", max(adp[['temperature']]))

  metad <- read.csv(metadata, header = TRUE)

  mn <- as.character(metad[,1])
  mv <- as.character(metad[,2])


  md <- as.list(mv)
  names(md) <- mn

  if (!missing(md)) {
    for (m in seq_along(md)) {
      ncatt_put(ncout, 0, names(md)[m], md[[m]])
    }
    nc_close(ncout)


  }
}


