
require(oce)
require(ncdf4)



##read.adp.easy

#'ADCP Processing step 2.0
#'
#'@family processing
#'
#'
#'Load adp data into R with list that includes all metadata from mooring sheets
#'#' returns an object of class adp (from oce package) uses
#'\code{\link[oce:read.adp]{read.adp}}
#'
#'@param file raw ADCP file (.000 format)
#'@param metadata csv metadata file from template
#'


read.adp.easy <- function(file, metadata){
  if (missing(metadata)){
    warning('no metadata supplied')
  }
  metad <- read.csv(metadata, header = TRUE)

  mn <- as.character(metad[,1])
  mv <- as.character(metad[,2])


  md <- as.list(mv)
  names(md) <- mn

  adp <- read.adp(file, latitude = as.numeric(md[['latitude']]), longitude = as.numeric(md[['longitude']]) ) #insert lat and lon from mooring logs

  if (!missing(md)) {
    for (m in seq_along(md)) {
      adp <- oceSetMetadata(adp, names(md)[m], md[[m]], note = NULL)
    }

    adp@metadata$source <- 'raw'

    adp@metadata$latitude <- as.numeric(adp[['latitude']])
    adp@metadata$longitude <- as.numeric(adp[['longitude']])
    adp@processingLog <- processingLogAppend(adp@processingLog, 'metadata read in from log sheet')

    return(adp)
  }
}


##read metadata

#' ADCP process step 2.1
#'
#'@family processing
#' Read in all metadata from csv template to adp object, sourced from log sheets
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

#' ADCP Processing step 3.2
#'
#'@family processing
#'
#' apply magnetic declination to ADCP data uses
#' \code{\link[oce:magneticField]{magneticField}} to calculate declination and
#' \code{\link[oce:enuToOther]{enuToOTher}} to apply variation to data set
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
#'   If type == average, an average is taken from start and end time
#'   declinations and applied uniformly
#'
#'   If type == interpolated, the rate of declination is used over time series
#'
#'


applyMagneticDeclinationAdp <- function(x, lat = x[['latitude']], lon = x[['longitude']], st = x[['time_coverage_start']], et = x[['time_coverage_end']],tz = 'UTC', type = 'average'){
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
          x@metadata$magnetic_variation <- c
          x@metadata$oceCoordinate <- 'enu'
        }
        if (coord != 'enu'){
          warning('Function cannot handle objects in ', coord, 'Object returned as is ; please convert to enu first')

        }
      }
      else {
        warning('Missing required arguments! No processing performed!')
      }
      x@processingLog <- processingLogAppend(x@processingLog, value = paste0('magnetic variation applied; declination =', c, 'degrees') )
    }
    if (type == 'interpolate'){
      ;
    }
    return(x)


  }
}



#'ADCP Processing step 3.3 Limit depth by rmax
#'
#'@family processing
#'
#'Use maximum acceptable range values to determine acceptable depth values Uses
#'Teledyne RDI equation, Rmax = Dcos(x)
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
    x@metadata$sensor_depth <- mdt
    x@data$depth <- d

  }
  if (missing(lat)){
    if(!is.na(x@metadata$latitude)){
      lat <- x@metadata$latitude
      d <- swDepth(x@data$pressure, latitude = lat, eos = getOption("oceEOS", default = "gsw"))

      d[d < rmax] <- NA
      mdt <- round(mean(x@data$depth, na.rm = TRUE), digits = 2)
      x@metadata$sensor_depth <- mdt
      x@data$depth <- d
    }
    if (is.na(x@metadata$latitude)){
      warning('No latitude provided; returning object as is')

      stop()
    }
    adp@processingLog <- processingLogAppend(adp@processingLog, paste0('depth values adjusted to sea water depth using pressure and latitude'))
    adp@processingLog <- processingLogAppend(adp@processingLog, paste0('depth limited by maximum acceptable distance, calulated with Rmax = Dcos(x)'))
    adp@processingLog <- processingLogAppend(adp@processingLog, paste0('Sensor depth and mean depth set to  ', mdt , '  based on trimmed depth values'))

    return(x)
  }
}



#'@title limit depth by time
#'
#'@family processing
#'
#'@description Uses deployment and recovery times to limit depth within times
#'  that insturment was actively and properly sampliing
#'
#'@param adp oce object of class adp to be limited
#'
#'  requires certain meta data features to compute including pressure, latitude,
#'  time, time_coverage_start, time_coverage_end
#'


limit_depthbytime <- function(adp, tz = 'UTC'){
  if (!inherits(adp, "adp")){
    stop("method is only for objects of class '", "adp", "'")
  }
  adp[['depth']] <- swDepth(adp[['pressure']], latitude = adp[['latitude']], eos = getOption("oceEOS", default = "gsw"))
  depth <- adp[['depth']]
  depth[as.POSIXct(adp[['time']], tz) <= as.POSIXct(adp[['time_coverage_start']], tz) | as.POSIXct(adp[['time']], tz) >= as.POSIXct(adp[['time_coverage_end']], tz)] <- NA

  mdt <- round(mean(depth, na.rm = TRUE), digits = 2)
  adp@metadata$sensor_depth <- mdt

  adp@data$depth <- depth
  adp@processingLog <- processingLogAppend(adp@processingLog, paste0('depth limited by deployment (', adp[['time_coverage_start']], ') and recovery  (', adp[['time_coverage_end']], ')  times'))
  adp@processingLog <- processingLogAppend(adp@processingLog, paste0('Sensor depth and mean depth set to  ', mdt , '  based on trimmed depth values'))

  return(adp)
}

#time cut off
#'
#'ADCP Processing step 3.4
#'
#'@family processing
#'
#'@description Function limits variable (time, salinity, pressure, temperature,
#'  pitch, roll, heading) from before deployment and after recovery of ADCP
#'
#'@param x adp object from oce-class adp
#'@param tz time zone, default is 'UTC'
#'@param dt deployment time of ADCP, default pulls value from metadata (time_coverage_start)
#'@param rt recovery time of ADCP, default pulls value from metadata (time_coverage_end)
#'
#'@return adp object with velocities limited to active ADCP measurement times
#'  (after deployment and before recovery)
#'
#'
#'


limit_time <- function(x, tz = 'UTC', dt = x[['time_coverage_start']], rt = x[['time_coverage_end']]){

  if (!inherits(x, "adp")){
    stop("method is only for objects of class '", "adp", "'")
  }
  #FIX ME : if deployment or recovery time is out of bounds

  if(!missing(dt)){
    t1 <- as.POSIXct(dt, tz = tz)
    t <- x[['time', "numeric"]]
    t <- as.POSIXct(t, tz = tz)


    x[['v']][t < t1] <- NA
    x[['pressure']][t < t1] <- NA
    x[['salinity']][t < t1] <- NA
    x[['temperature']][t < t1] <- NA
    x[['pitch']][t < t1] <- NA
    x[['roll']][t < t1] <- NA
    x[['heading']][t < t1] <- NA
  }
  else if(missing(dt)){
    if (!is.null(x@metadata$time_coverage_start)){
      dt <- x@metadata$time_coverage_start
      t1 <- as.POSIXct(dt, tz = tz)
      t <- x[['time', "numeric"]]
      t <- as.POSIXct(t, tz = tz)
      x[['v']][t < t1] <- NA
      x[['pressure']][t < t1] <- NA
      x[['salinity']][t < t1] <- NA
      x[['temperature']][t < t1] <- NA
      x[['pitch']][t < t1] <- NA
      x[['roll']][t < t1] <- NA
      x[['heading']][t < t1] <- NA
    }
    if (is.null(x@metadata$time_coverage_start)){
      warning('No deployment Time provided!')
    }


    if(!missing(rt)){
      t2 <- as.POSIXct(rt, tz = tz)
      x[['v']][t > t2] <- NA
      x[['pressure']][t > t2] <- NA
      x[['salinity']][t > t2] <- NA
      x[['temperature']][t > t2] <- NA
      x[['pitch']][t > t2] <- NA
      x[['roll']][t > t2] <- NA
      x[['heading']][t > t2] <- NA
    }
    else if (missing(rt))
      if (!is.null(x@metadata$time_coverage_end)){
        t2 <- as.POSIXct(rt, tz = tz)
        x[['v']][t > t2] <- NA
        x[['pressure']][t > t2] <- NA
        x[['salinity']][t > t2] <- NA
        x[['temperature']][t > t2] <- NA
        x[['pitch']][t > t2] <- NA
        x[['roll']][t > t2] <- NA
        x[['heading']][t > t2] <- NA
      }
    if (is.null(x@metadata$time_coverage_end)){
      warning('No recovery time provided!')
    }

    adp@processingLog <- processingLogAppend(adp@processingLog, paste0('Data has been cut off before deployment at  ', dt, '  and after recovery at ', rt))

    return(x)
  }
}

####create flag####

#' ADCP processing Step 3.5
#'
#'@family processing
#'
#' adpFlag function flag an adp object based on a series of parameters including
#' percent good, and error velocity.
#'
#' This function also defaults to flagging depth values based on the Teledyne
#' RDI standard
#'
#' Rmax = Dcosx where Rmax is the maximum acceptable distance range from the
#' ADCP, D is total depth and x is the beam angle of the ADCP.
#'
#' This function also flags data outside the deployment and recovery time
#' bounds.
#'
#' Sets any data not flagged by this processing to 'good' flag value
#'
#' This function uses \code{\link[oce:initializeFlags]{initializeFlags}} to
#' initialize blank flagging scheme for values to be inserted in Then
#' \code{\link[oce: setFlags]{setFlags}} to set flag values based on desired
#' scheme
#'
#' @param adp, an adp object, oce-class
#' @param flagScheme, scheme of flags that will be followed, BODC, MEDS, etc
#' @param pg, The minimum percent good for evaluating beams one + four of the
#'   adcp (BIO standard is 25)
#' @param er, Maximum error velocity for evaluating adcp data (BIO standard is
#'   0.46)
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
    flag[,,i] <- (lowpg< pg) | (abs(ERRV) > er) | r < d | adp[['time']] < adp[['time_coverage_start']] | adp[['time']] > adp[['time_coverage_end']]

  #initialize and set flag scheme
  adp <- initializeFlags(adp, name = 'v', value = 0)


  #set adp flags where logical array = TRUE, flagged based on error, percent good or Rmax, value = 4 (see flag scheme, BODC)
  adp <- setFlags(adp, name = 'v', i= flag, value = 4)
  good <- adp[['flags']][['v']] == 0
  adp <- setFlags(adp, name = 'v', i = good, value = 1)


  adp@processingLog <- processingLogAppend(adp@processingLog, 'Quality control flags set based on flag scheme L20 from BODC')
  return(adp)
}


##creating standard format file names

#' @title File naming format
#'
#'@family general
#'
#' @description pulls metadata from an adp oce object and creates a file name in
#'   standard BIO format
#'
#' @param adp an oce class - adp object
#'
#'   required metadata : cruise number, mooring number, serial number, sampling
#'   interval!! Make sure all metadata is present when experiencing problems
#'
#'
#'
name.file <- function(adp){

  name <- paste('MADCP', adp[['experiment']], adp[['mooring_number']], adp[['serialNumber']], adp[['sampling_interval']], sep = '_')

  name
}

#' ADCP Processing ODF to NetCDF
#'
#'@family odf
#'
#' @description Converting individual odf bins to Net cdf standard format
#'
#' @param files list of odf files
#' @param metadata any extra metadata to be added to net cdf as list form
odf2adp <- function(files, metadata) {
  require(oce)
  require(abind)

  files <- if (is.list(files)) unlist(files) else files

  nd <- length(files)
  ## read the first one to get length of time:
  d <- read.odf(files[1])
  nt <- length(d[['time']])
  vars <- names(d@data)
  vars <- vars[-which(vars == 'time')]
  for (vr in vars) {
    assign(vr, array(NA, dim=c(nt, nd)))
  }
  depth <- NULL
  for (f in 1:length(files)) {
    d <- read.odf(files[f])
    t <- d[['time']]
    depth[f] <- d[['depthMin']]
    for (vr in vars) {
      eval(parse(text=paste0(vr, "[, f] <- d[['", vr, "']]")))
    }
  }

  ## need to sort the depths because of file sorting ...
  o <- order(depth, decreasing = TRUE)
  depth <- depth[o]
  for (vr in vars) {
    eval(parse(text=paste0(vr, "<- ", vr, "[, o]")))
  }

  distance <- max(depth) - depth
  adp <- as.adp(t, distance, v=abind(u, v, w, error, along=3), a=a, q=unknown)
  for (m in names(d@metadata)) {
    if (m != 'units' & m != 'flags' & m != 'dataNamesOriginal') {
      adp <- oceSetMetadata(adp, m, d[[m]])
    }
  }

  ## depthMinMax
  adp <- oceSetMetadata(adp, 'depthMin', min(depth))
  adp <- oceSetMetadata(adp, 'depthMax', max(depth))

  ## add in any extra supplied metadata items
  if (!missing(metadata)) {
    for (m in seq_along(metadata)) {
      adp <- oceSetMetadata(adp, names(metadata)[m], metadata[[m]])
    }
  }
  adp@metadata$source <- 'odf'
  adp@processingLog <- processingLogAppend(adp@processingLog, 'Creation : Data and metadata read into adp object from ODF file')

  return(adp)
}

#'ADCP Processing step 4.1
#'
#'@family processing
#'@family NC
#'
#'
#'@description Exports an adp object to a net cdf using variables and metadata
#'  within adp combined with optional additional metatdata see details in
#'  \code{\link[ncdf4:ncdf4]{ncdf4}} package
#'
#'@param obj an adp object from the oce class
#'@param name name of the NetCDF file to be produced
#'@param metadata csv file listing metadata names and values to be inserted into
#'  global attributes of net CDF




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
  timedim <- ncdim_def("time", "seconds since 1970-01-01T00:00:00Z", as.double(time))    #time formatting FIX
  depthdim <- ncdim_def("depth", "metres", as.double(dist))
  stationdim <- ncdim_def("station", "", as.numeric(adp[['mooring_number']]))
  londim <- ncdim_def("lon", "degrees_east" , as.double(lon))
  latdim <- ncdim_def("lat", "degrees_north", as.double(lat))
  dimnchar <- ncdim_def('nchar', '', 1:23, create_dimvar = FALSE)


  #set fill value
  FillValue <- 1e35

  if (adp@metadata$source == 'raw'){

    #define variables

    dlname <- 'lon'
    lon_def <- ncvar_def(longname= "longitude", units = 'degrees_east', dim = stationdim, name = dlname, prec = 'double')

    dlname <- 'lat'
    lat_def <- ncvar_def( longname = 'latitude', units = 'degrees_north', dim =  stationdim, name = dlname, prec = 'double')

    dlname <- "eastward_sea_water_velocity"
    u_def <- ncvar_def("EWCT", "m/sec", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")

    dlname <- "northward_sea_water_velocity"
    v_def <- ncvar_def("NSCT", "m/sec", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")

    dlname <- "upward_sea_water_velocity"
    w_def <- ncvar_def("VCSP", "m/sec", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")

    dlname <- "time_02"
    t_def <- ncvar_def("ELTMEP01", "seconds since 1970-01-01T00:00:00Z", list( stationdim, timedim), FillValue, dlname, prec = "float")

    dlname <- "error_velocity_in_sea_water"
    e_def <- ncvar_def("ERRV", "m/sec", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")

    dlname <- "ADCP_echo_intensity_beam_1"

    b1_def <- ncvar_def("BEAM_01", "counts", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")

    dlname <- "ADCP_echo_intensity_beam_2"
    b2_def <- ncvar_def("BEAM_02", "counts", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")

    dlname <- "ADCP_echo_intensity_beam_3"
    b3_def <- ncvar_def("BEAM_03", "counts", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")

    dlname <- "ADCP_echo_intensity_beam_4"
    b4_def <- ncvar_def("BEAM_04", "counts", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")


    dlname <- "percent_good_beam_1"
    pg1_def <- ncvar_def("PGDP_01", "counts", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")

    dlname <- "percent_good_beam_2"
    pg2_def <- ncvar_def("PGDP_02", "counts", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")

    dlname <- "percent_good_beam_3"
    pg3_def <- ncvar_def("PGDP_03", "counts", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")

    dlname <- "percent_good_beam_4"
    pg4_def <- ncvar_def("PGDP_04", "counts", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")

    dlname <- "pitch"
    p_def <- ncvar_def("PTCH", "degrees", list( stationdim, timedim), FillValue, dlname, prec = "float")

    dlname <- "roll"
    r_def <- ncvar_def("ROLL", "degrees", list( stationdim, timedim ), FillValue, dlname, prec = "float")

    dlname <- "height of sea surface"
    hght_def <- ncvar_def("hght", "m", list( stationdim, depthdim ), FillValue, dlname, prec = "float")

    dlname <- "ADCP Transducer Temp."
    Tx_def <- ncvar_def("Tx", "degrees", list(stationdim, timedim), FillValue, dlname, prec = "float")

    dlname <- "instrument depth"
    D_def <- ncvar_def("D", "m", list(timedim, stationdim), FillValue, dlname, prec = "float")

    dlname <- "heading"
    head_def <- ncvar_def("HEAD", "degrees", list(timedim, stationdim), FillValue, dlname, prec = "float")

    dlname <- "pressure"
    pres_def <- ncvar_def("PRES", "decibars", list(timedim, stationdim), FillValue, dlname, prec = "float")

    dlname <- "speed of sound"
    svel_def <- ncvar_def("SVEL", "m/s", list(timedim, stationdim), FillValue, dlname, prec = "float")

    dlname <- "time_string"
    ts_def <- ncvar_def("DTUT8601", units = "",dim =  list(timedim, dimnchar), missval = NULL, name =  dlname, prec = "char")


    FillValue <- 0

    dlname <- "quality_flag u"
    qc_u_def <- ncvar_def("QC_flag_u", "", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "integer")

    dlname <- "quality_flag v"
    qc_v_def <- ncvar_def("QC_flag_v", "", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "integer")

    dlname <- "quality_flag w"
    qc_w_def <- ncvar_def("QC_flag_w", "", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "integer")

    ####writing net CDF####
    #write out definitions to new nc file
    ncout <- nc_create(ncfname, list(u_def, v_def, w_def, e_def, t_def, b1_def, b2_def, b3_def, b4_def, pg1_def, pg2_def, pg3_def, pg4_def, p_def, r_def, hght_def, Tx_def, D_def, qc_u_def, qc_v_def, qc_w_def, lon_def, lat_def, head_def, pres_def, svel_def, ts_def), force_v4 = TRUE)
  }

  if (adp@metadata$source == 'odf'){
    #define variables

    dlname <- 'lon'
    lon_def <- ncvar_def(longname= "longitude", units = 'degrees_east', dim = stationdim, name = dlname, prec = 'double')

    dlname <- 'lat'
    lat_def <- ncvar_def( longname = 'latitude', units = 'degrees_north', dim =  stationdim, name = dlname, prec = 'double')

    dlname <- "eastward_sea_water_velocity"
    u_def <- ncvar_def("EWCT", "m/sec", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")

    dlname <- "northward_sea_water_velocity"
    v_def <- ncvar_def("NSCT", "m/sec", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")

    dlname <- "upward_sea_water_velocity"
    w_def <- ncvar_def("VCSP", "m/sec", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")

    dlname <- "time_02"
    t_def <- ncvar_def("ELTMEP01", "seconds since 1970-01-01T00:00:00Z", list( stationdim, timedim), FillValue, dlname, prec = "float")

    dlname <- "error_velocity_in_sea_water"
    e_def <- ncvar_def("ERRV", "m/sec", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")

    dlname <- "ADCP_echo_intensity_beam_1"

    b1_def <- ncvar_def("BEAM_01", "counts", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")


    dlname <- "percent_good_beam_1"
    pg1_def <- ncvar_def("PGDP_01", "counts", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")

    dlname <- "time_string"
    ts_def <- ncvar_def("DTUT8601", units = "",dim =  list(timedim, dimnchar), missval = NULL, name =  dlname, prec = "char")

    ####writing net CDF####
    #write out definitions to new nc file
    ncout <- nc_create(ncfname, list(u_def, v_def, w_def, e_def, t_def, b1_def,  pg1_def, lon_def, lat_def, ts_def), force_v4 = TRUE)


  }


  #insert variables into nc file
  ncvar_put(ncout, u_def, adp[['v']][,,1])
  ncvar_put(ncout, v_def, adp[['v']][,,2])
  ncvar_put(ncout, w_def, adp[['v']][,,3])
  ncvar_put(ncout, e_def, adp[['v']][,,4])
  ncvar_put(ncout, t_def, as.POSIXct(adp[['time']], tz = 'UTC', origin = '1970-01-01 00:00:00'))
  ncvar_put(ncout, lon_def, adp[['longitude']])
  ncvar_put(ncout, lat_def, adp[['latitude']])

  if (adp@metadata$source == 'raw'){
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
    ncvar_put(ncout, hght_def, (adp[['sensor_depth']]- adp[['distance']]))
    ncvar_put(ncout, Tx_def, adp[['temperature']])
    ncvar_put(ncout, D_def, adp[['depth']])
    ncvar_put(ncout, qc_u_def, adp@metadata$flags$v[,,1])
    ncvar_put(ncout, qc_v_def, adp@metadata$flags$v[,,2])
    ncvar_put(ncout, qc_w_def, adp@metadata$flags$v[,,3])
    ncvar_put(ncout, head_def, adp[['heading']])
    ncvar_put(ncout, pres_def, adp[['pressure']])
    ncvar_put(ncout, svel_def, adp[['soundSpeed']])
    ncvar_put(ncout, ts_def, adp[['time']])
  }
  if (adp@metadata$source == 'odf'){
    ncvar_put(ncout, b1_def, adp[['a', 'numeric']])
    ncvar_put(ncout, pg1_def, adp[['q', 'numeric']])
    ncvar_put(ncout, ts_def, adp[['time']])

  }

  ###metadata###

  ncatt_put(ncout, 'station', attname = 'cf_role',attval =  'timeseries_id')
  ncatt_put(ncout, 'time', attname = 'cf_role', attval = 'profile_id')
  ncatt_put(ncout, 'station', 'standard_name', 'platform_name')
  ncatt_put(ncout, 'time' , 'calendar', 'gregorian')
  ncatt_put(ncout, 'time_string', 'note', 'time values as ISO8601 string, YY-MM-DD hh:mm:ss')
  ncatt_put(ncout, 'time_string', 'time_zone', 'UTC')
  nc_att_put(ncout,0, 'processing_level', adp[['processingLog']])


  if (adp@metadata$source == 'raw'){
    ####pulled from adp object####
    ncatt_put(ncout, 0, "mooring_number", adp[['mooring_number']])

    #       deprecated --- Diana Cardoso 06/01/2018
    #ncatt_put(ncout, 0, "deployment_date", adp[['deployment_time']])
    #ncatt_put(ncout, 0, "recovery_date", adp[['recovery_time']])


    ncatt_put(ncout, 0, "firmware_version", adp[['firmwareVersion']])
    ncatt_put(ncout, 0, "frequency", adp[['frequency']])
    ncatt_put(ncout, 0, "beam_pattern", adp[['beamPattern']])
    ncatt_put(ncout, 0, "janus", adp[['numberOfBeams']])
    ncatt_put(ncout, 0, "pings_per_ensemble", adp[['pingsPerEnsemble']])
    ncatt_put(ncout, 0, "valid_correlation_range", adp[['lowCorrThresh']])
    ncatt_put(ncout, 0, "minmax_percent_good", adp[['percentGdMinimum']])
    ncatt_put(ncout, 0,"minmax_percent_good", "100")
    ncatt_put(ncout, 0, "error_velocity_threshold", adp[['errorVelocityMaximum']])
    ncatt_put(ncout, 0, "transmit_pulse_length_cm", adp[['xmitPulseLength']]*100)
    ncatt_put(ncout, 0, "false_target_reject_values", adp[['falseTargetThresh']])
    ncatt_put(ncout, 0, "serial_number", adp[['serialNumber']])

    #     deprecated --- Diana Cardoso 06/01/2018
    #ncatt_put(ncout, 0, "transform", adp[['oceCoordinate']])


    ncatt_put(ncout, 0, "data_type", adp[['instrumentType']])
    ncatt_put(ncout, 0, "data_subtype", adp[['data_subtype']])
    ncatt_put(ncout, 0, "coord_system", adp[['oceCoordinate']])
    ncatt_put(ncout, 0, "longitude", adp[['longitude']])
    ncatt_put(ncout, 0, "latitude", adp[['latitude']])
    ncatt_put(ncout, 0, "magnetic_variation", adp[['magnetic_variation']])
    ncatt_put(ncout, 0, "platform", adp[['platform']])
    ncatt_put(ncout, 0, "sounding", adp[['sounding']])
    ncatt_put(ncout, 0, "chief_scientist", adp[['chief_scientist']])
    ncatt_put(ncout, 0, "data_origin", adp[['institution']])
    ncatt_put(ncout, 0, "water_depth", adp[['sounding']])
    ncatt_put(ncout, 0, "delta_t_sec",adp[['sampling_interval']])
    ncatt_put(ncout, 0, "pred_accuracy", adp[['velocityResolution']]*1000)
    ncatt_put(ncout, "depth", "xducer_offset_from_bottom", as.numeric(adp[['sounding']]) - adp[['sensor_depth']])
    ncatt_put(ncout, "depth", "bin_size", adp[['cellSize']])
    ncatt_put(ncout, "EWCT", "sensor_type", adp[['instrumentType']])
    ncatt_put(ncout, "EWCT", "sensor_depth", adp[['sensor_depth']])
    ncatt_put(ncout, "EWCT", "serial_number", adp[['serialNumber']])
    ncatt_put(ncout, "NSCT", "sensor_type", adp[['instrumentType']])
    ncatt_put(ncout, "NSCT", "sensor_depth", adp[['sensor_depth']])
    ncatt_put(ncout, "NSCT", "serial_number", adp[['serialNumber']])
    ncatt_put(ncout, "VCSP", "sensor_type", adp[['instrumentType']])
    ncatt_put(ncout, "VCSP", "sensor_depth", adp[['sensor_depth']])
    ncatt_put(ncout, "VCSP", "serial_number", adp[['serialNumber']])
    ncatt_put(ncout, "ERRV", "sensor_type", adp[['instrumentType']])
    ncatt_put(ncout, "ERRV", "sensor_depth", adp[['sensor_depth']])
    ncatt_put(ncout, "ERRV", "serial_number", adp[['serialNumber']])
    ncatt_put(ncout, "BEAM_01", "sensor_type", adp[['instrumentType']])
    ncatt_put(ncout, "BEAM_01", "sensor_depth", adp[['sensor_depth']])
    ncatt_put(ncout, "BEAM_01", "serial_number", adp[['serialNumber']])
    ncatt_put(ncout, "BEAM_02", "sensor_type", adp[['instrumentType']])
    ncatt_put(ncout, "BEAM_02", "sensor_depth", adp[['sensor_depth']])
    ncatt_put(ncout, "BEAM_02", "serial_number", adp[['serialNumber']])
    ncatt_put(ncout, "BEAM_03", "sensor_type", adp[['instrumentType']])
    ncatt_put(ncout, "BEAM_03", "sensor_depth", adp[['sensor_depth']])
    ncatt_put(ncout, "BEAM_03", "serial_number", adp[['serialNumber']])
    ncatt_put(ncout, "BEAM_04", "sensor_type", adp[['instrumentType']])
    ncatt_put(ncout, "BEAM_04", "sensor_depth", adp[['sensor_depth']])
    ncatt_put(ncout, "BEAM_04", "serial_number", adp[['serialNumber']])
    ncatt_put(ncout, "PGDP_01", "sensor_type", adp[['instrumentType']])
    ncatt_put(ncout, "PGDP_01", "sensor_depth", adp[['sensor_depth']])
    ncatt_put(ncout, "PGDP_01", "serial_number", adp[['serialNumber']])
    ncatt_put(ncout, "PGDP_02", "sensor_type", adp[['instrumentType']])
    ncatt_put(ncout, "PGDP_02", "sensor_depth", adp[['sensor_depth']])
    ncatt_put(ncout, "PGDP_02", "serial_number", adp[['serialNumber']])
    ncatt_put(ncout, "PGDP_03", "sensor_type", adp[['instrumentType']])
    ncatt_put(ncout, "PGDP_03", "sensor_depth", adp[['sensor_depth']])
    ncatt_put(ncout, "PGDP_03", "serial_number", adp[['serialNumber']])
    ncatt_put(ncout, "PGDP_04", "sensor_type", adp[['instrumentType']])
    ncatt_put(ncout, "PGDP_04", "sensor_depth", adp[['sensor_depth']])
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
    ncatt_put(ncout, "hght", "sensor_depth", adp[['sensor_depth']])
    ncatt_put(ncout, "hght", "serial_number", adp[['serialNumber']])
    ncatt_put(ncout, "D", "generic_name", "depth")
    ncatt_put(ncout, "D", "sensor_type", adp[['instrumentType']])
    ncatt_put(ncout, "D", "sensor_depth", adp[['sensor_depth']])
    ncatt_put(ncout, "D", "serial_number", adp[['serialNumber']])
    ncatt_put(ncout, "Tx", "generic_name", "temp")
    ncatt_put(ncout, "Tx", "sensor_type", adp[['instrumentType']])
    ncatt_put(ncout, "Tx", "sensor_depth", adp[['sensor_depth']])
    ncatt_put(ncout, "Tx", "serial_number", adp[['serialNumber']])
    ncatt_put(ncout, "QC_flag_u", "comment", "Quality flag resulting from quality control")
    ncatt_put(ncout, "QC_flag_u", "flag_meanings",adp[['flag_meaning']])
    ncatt_put(ncout, "QC_flag_u", "flag_values",c(1:9))
    ncatt_put(ncout, "QC_flag_u", "References", adp[['flag_references']])
    ncatt_put(ncout, "QC_flag_v", "comment", "Quality flag resulting from quality control")
    ncatt_put(ncout, "QC_flag_v", "flag_meanings", adp[['flag_meaning']])
    ncatt_put(ncout, "QC_flag_v", "flag_values",c(1:9))
    ncatt_put(ncout, "QC_flag_v", "References", adp[['flag_references']])
    ncatt_put(ncout, "QC_flag_w", "comment", "Quality flag resulting from quality control")
    ncatt_put(ncout, "QC_flag_w", "flag_meanings", adp[['flag_meaning']])
    ncatt_put(ncout, "QC_flag_w", "flag_values",c(1:9))
    ncatt_put(ncout, "QC_flag_w", "References", adp[['flag_references']])

    #CF conventions

    ncatt_put(ncout, 0, 'Conventions', 'CF-1.7')
    ncatt_put(ncout, 0, "creator_type", "person")
    ncatt_put(ncout, 0, "creator_institution", adp[['institution']])
    ncatt_put(ncout, 0, "program", adp[['program']])
    ncatt_put(ncout, 0, "sea_name", adp[['sea_name']])
    ncatt_put(ncout, 0, "time_coverage_start", adp[['time_coverage_start']])
    ncatt_put(ncout, 0, "time_coverage_end", adp[['time_coverage_end']])
    ncatt_put(ncout, 0, "geospatial_lat_min", adp[['latitude']])
    ncatt_put(ncout, 0, "geospatial_lat_max", adp[['latitude']])
    ncatt_put(ncout, 0, "geospatial_lat_units", "degrees_north")
    ncatt_put(ncout, 0, "geospatial_lon_min", adp[['longitude']])
    ncatt_put(ncout, 0, "geospatial_lon_max", adp[['longitude']])
    ncatt_put(ncout, 0, "geospatial_lon_units", "degrees_east")

    if (adp[['orientation']] == 'up'){
      ncatt_put(ncout, 0, "geospatial_vertical_min", adp[['sensor_depth']] + max(adp[['distance']], na.rm = TRUE))
      ncatt_put(ncout, 0, "geospatial_vertical_max", adp[['sensor_depth']] + min(adp[['distance']], na.rm = TRUE))
    }
    if (adp[['orientation']] == 'down'){
      ncatt_put(ncout, 0, "geospatial_vertical_min", adp[['sensor_depth']] + min(adp[['distance']], na.rm = TRUE))
      ncatt_put(ncout, 0, "geospatial_vertical_max", adp[['sensor_depth']] + max(adp[['distance']], na.rm = TRUE))
    }
    ncatt_put(ncout, 0, "geospatial_vertical_units", "metres")
    ncatt_put(ncout, 0, "geospatial_vertical_positive", 'down')
    ncatt_put(ncout, 0, "institution", adp[['institution']])
    ncatt_put(ncout, 0, "creator_name", adp[['creator_name']])
    ncatt_put(ncout, 0, "creator_url", adp[['creator_url']])
    ncatt_put(ncout, 0, "creator_email", adp[['creator_email']])
    ncatt_put(ncout, 0, "project", adp[['project']])
    ncatt_put(ncout, 0, "processing_level", adp[['processing_level']])
    ncatt_put(ncout, 0 , "flag_meanings", adp[['flag_meaning']])
    ncatt_put(ncout, 0 , "flag_values", c(1:9))
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
    ncatt_put(ncout, "BEAM_01", "sdn_parameter_urn", "SDN:P01::TNIHCE01")
    ncatt_put(ncout, "BEAM_02", "sdn_parameter_urn", "SDN:P01::TNIHCE02")
    ncatt_put(ncout, "BEAM_03", "sdn_parameter_urn", "SDN:P01::TNIHCE03")
    ncatt_put(ncout, "BEAM_04", "sdn_parameter_urn", "SDN:P01::TNIHCE04")
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
    ncatt_put(ncout, "HEAD", "sdn_parameter_urn", "SDN:P01::HEADCM01")
    ncatt_put(ncout, "PRES", "sdn_parameter_urn", "SDN:P01::PRESPR01")
    ncatt_put(ncout, "SVEL", "sdn_parameter_urn", "SDN:P01::SVELCV01")
    ncatt_put(ncout, "time_string", "sdn_parameter_urn", "SDN:P01::DTUT8601")



    ncatt_put(ncout, "EWCT", "sdn_parameter_name", "Eastward current velocity (Eulerian) in the water body by moored acoustic doppler current profiler (ADCP)")
    ncatt_put(ncout, "NSCT", "sdn_parameter_name", "Northward current velocity (Eulerian) in the water body by moored acoustic doppler current profiler (ADCP)")
    ncatt_put(ncout, "VCSP", "sdn_parameter_name", "Upward current velocity in the water body by moored acoustic doppler current profiler (ADCP)")
    ncatt_put(ncout, "ERRV", "sdn_parameter_name", "Current velocity error in the water body by moored acoustic doppler current profiler (ADCP)")
    ncatt_put(ncout, "BEAM_01", "sdn_parameter_name", "Echo intensity from the water body by moored acoustic doppler current profiler (ADCP) beam 1")
    ncatt_put(ncout, "BEAM_02", "sdn_parameter_name", "Echo intensity from the water body by moored acoustic doppler current profiler (ADCP) beam 2")
    ncatt_put(ncout, "BEAM_03", "sdn_parameter_name", "Echo intensity from the water body by moored acoustic doppler current profiler (ADCP) beam 3")
    ncatt_put(ncout, "BEAM_04", "sdn_parameter_name", "Echo intensity from the water body by moored acoustic doppler current profiler (ADCP) beam 4")
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
    ncatt_put(ncout, "HEAD", "sdn_parameter_name", "Orientation (horizontal relative to true north) of measurement device {heading}")
    ncatt_put(ncout, "PRES", "sdn_parameter_name", "Pressure (spatial co-ordinate) exerted by the water body by profiling pressure sensor and corrected to read zero at sea level")
    ncatt_put(ncout, "SVEL", "sdn_parameter_name", "Sound velocity in the water body by computation from temperature and salinity by unspecified algorithm")
    ncatt_put(ncout, 'ELTMEP01', "sdn_parameter_name", "Elapsed time (since 1970-01-01T00:00:00Z)")
    ncatt_put(ncout, 'time_string', "sdn_parameter_name", "String corresponding to format 'YYYY-MM-DDThh:mm:ss.sssZ' or other valid ISO8601 string")




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
    ncatt_put(ncout, "HEAD", "sdn_uom_urn", "SDN:P06:UAAA")
    ncatt_put(ncout, "PRES", "sdn_uom_urn", "SDN:P06:UPDB")
    ncatt_put(ncout, "SVEL", "sdn_uom_urn", "SDN:P06:UVAA")
    ncatt_put(ncout, "ELTMEP01", "sdn_uom_urn", "SDN:P06::UTBB")
    ncatt_put(ncout, "time_string", "sdn_uom_urn", "SDN:P06::TISO")


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
    ncatt_put(ncout, "HEAD", "sdn_uom_name", "Degrees")
    ncatt_put(ncout, "PRES", "sdn_uom_name", "Decibars")
    ncatt_put(ncout, "SVEL", "sdn_uom_name", "Metres per second")
    ncatt_put(ncout, "ELTMEP01", "sdn_uom_name", "Seconds")
    ncatt_put(ncout, "time_string", "sdn_uom_name", "ISO8601")


    #CF standard names
    ncatt_put(ncout, "EWCT", "standard_name", "eastward_sea_water_velocity")
    ncatt_put(ncout, "NSCT", "standard_name", "northward_sea_water_velocity")
    ncatt_put(ncout, "VCSP", "standard_name", "upward_sea_water_velocity")
    ncatt_put(ncout, "ELTMEP01", "standard_name", "time")
    ncatt_put(ncout, "lat", "standard_name", "latitude")
    ncatt_put(ncout, "lon", "standard_name", "longitude")
    ncatt_put(ncout, "D", "standard_name", "depth")
    ncatt_put(ncout, "PTCH", "standard_name", "platform_pitch_angle")
    ncatt_put(ncout, "ROLL", "standard_name", "platform_roll_angle")
    ncatt_put(ncout, "PRES", "standard_name", "sea_water_pressure")
    ncatt_put(ncout, "SVEL", "standard_name", "speed_of_sound_in_sea_water")

  }

  if (adp@metadata$source == 'odf'){
    ncatt_put(ncout, 0, "mooring_number", adp[['mooring_number']])

    #       deprecated --- Diana Cardoso 06/01/2018
    #ncatt_put(ncout, 0, "deployment_date", adp[['deployment_time']])
    #ncatt_put(ncout, 0, "recovery_date", adp[['recovery_time']])


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
    ncatt_put(ncout, 0, "serial_number", adp[['serialNumber']])

    #deprecated --- Diana Cardoso 06/01/2018
    #ncatt_put(ncout, 0, "transform", adp[['oceCoordinate']])

    ncatt_put(ncout, 0, "data_type", adp[['instrumentType']])
    ncatt_put(ncout, 0, "data_subtype", adp[['model']])
    ncatt_put(ncout, 0, "coord_system", adp[['oceCoordinate']])
    ncatt_put(ncout, 0, "longitude", adp[['longitude']])
    ncatt_put(ncout, 0, "latitude", adp[['latitude']])
    ncatt_put(ncout, 0, "magnetic_variation", adp[['magneticVariation']])
    ncatt_put(ncout, 0, "platform", adp[['platform']])
    ncatt_put(ncout, 0, "sounding", adp[['sounding']])
    ncatt_put(ncout, 0, "chief_scientist", adp[['chief_scientist']])
    ncatt_put(ncout, 0, "data_origin", adp[['institution']])
    ncatt_put(ncout, 0, "water_depth", adp[['water_depth']])
    ncatt_put(ncout, 0, "delta_t_sec",adp[['sampling_interval']])
    ncatt_put(ncout, 0, "pred_accuracy", adp[['velocityResolution']])

    #FIXME: should be pulled from odf...not in object... issue with oce read.odf
    ncatt_put(ncout, "depth", "xducer_offset_from_bottom", adp[['depth_off_bottom']])

    ncatt_put(ncout, "depth", "bin_size", adp[['cellSize']])
    ncatt_put(ncout, "EWCT", "sensor_type", adp[['instrumentType']])
    ncatt_put(ncout, "EWCT", "sensor_depth", adp[['sensor_depth']])
    ncatt_put(ncout, "EWCT", "serial_number", adp[['serialNumber']])
    ncatt_put(ncout, "NSCT", "sensor_type", adp[['instrumentType']])
    ncatt_put(ncout, "NSCT", "sensor_depth", adp[['sensor_depth']])
    ncatt_put(ncout, "NSCT", "serial_number", adp[['serialNumber']])
    ncatt_put(ncout, "VCSP", "sensor_type", adp[['instrumentType']])
    ncatt_put(ncout, "VCSP", "sensor_depth", adp[['sensor_depth']])
    ncatt_put(ncout, "VCSP", "serial_number", adp[['serialNumber']])
    ncatt_put(ncout, "ERRV", "sensor_type", adp[['instrumentType']])
    ncatt_put(ncout, "ERRV", "sensor_depth", adp[['sensor_depth']])
    ncatt_put(ncout, "ERRV", "serial_number", adp[['serialNumber']])
    ncatt_put(ncout, "BEAM_01", "sensor_type", adp[['instrumentType']])
    ncatt_put(ncout, "BEAM_01", "sensor_depth", adp[['sensor_depth']])
    ncatt_put(ncout, "BEAM_01", "serial_number", adp[['serialNumber']])
    ncatt_put(ncout, "PGDP_01", "sensor_type", adp[['instrumentType']])
    ncatt_put(ncout, "PGDP_01", "sensor_depth", adp[['sensor_depth']])
    ncatt_put(ncout, "PGDP_01", "serial_number", adp[['serialNumber']])
    ncatt_put(ncout, "EWCT", "generic_name", "u")
    ncatt_put(ncout, "NSCT", "generic_name", "v")
    ncatt_put(ncout, "VCSP", "generic_name", "w")
    ncatt_put(ncout, "ERRV", "generic_name", "w")       #issue in current NC protocol
    ncatt_put(ncout, "BEAM_01", "generic_name", "AGC")
    ncatt_put(ncout, "PGDP_01", "generic_name", "PGd")
    #CF

    ncatt_put(ncout, 0, 'Conventions', 'CF-1.7')
    ncatt_put(ncout, 0, "creator_type", "person")
    ncatt_put(ncout, 0, "creator_institution", adp[['institution']])
    ncatt_put(ncout, 0, "program", adp[['program']])
    ncatt_put(ncout, 0, "sea_name", adp[['sea_name']])
    ncatt_put(ncout, 0, "time_coverage_start", adp[['deployment_time']])
    ncatt_put(ncout, 0, "time_coverage_end", adp[['recovery_time']])
    ncatt_put(ncout, 0, "geospatial_lat_min", adp[['latitude']])
    ncatt_put(ncout, 0, "geospatial_lat_max", adp[['latitude']])
    ncatt_put(ncout, 0, "geospatial_lat_units", "degrees_north")
    ncatt_put(ncout, 0, "geospatial_lon_min", adp[['longitude']])
    ncatt_put(ncout, 0, "geospatial_lon_max", adp[['longitude']])
    ncatt_put(ncout, 0, "geospatial_lon_units", "degrees_east")
    if (adp[['orientation']] == 'up'){
      ncatt_put(ncout, 0, "geospatial_vertical_min", adp[['sensor_depth']] + max(adp[['distance']], na.rm = TRUE))
      ncatt_put(ncout, 0, "geospatial_vertical_max", adp[['sensor_depth']] + min(adp[['distance']], na.rm = TRUE))
    }
    if (adp[['orientation']] == 'down'){
      ncatt_put(ncout, 0, "geospatial_vertical_min", adp[['sensor_depth']] + min(adp[['distance']], na.rm = TRUE))
      ncatt_put(ncout, 0, "geospatial_vertical_max", adp[['sensor_depth']] + max(adp[['distance']], na.rm = TRUE))
    }
    ncatt_put(ncout, 0, "geospatial_vertical_units", "metres")
    ncatt_put(ncout, 0, "geospatial_vertical_positive", adp[['orientation']])     #eg up or down
    ncatt_put(ncout, 0, "institution", adp[['institution']])
    ncatt_put(ncout, 0, "creator_name", adp[['creator_name']])
    ncatt_put(ncout, 0, "creator_url", adp[['creator_url']])
    ncatt_put(ncout, 0, "creator_email", adp[['creator_email']])
    ncatt_put(ncout, 0, "project", adp[['project']])
    ncatt_put(ncout, 0, "processing_level", adp[['processingLevel']])
    ncatt_put(ncout, 0 , "flag_meanings", adp[['flag_meaning']])
    ncatt_put(ncout, 0 , "flag_values", c(1:9))
    ncatt_put(ncout, 0, "source", "R code: adcpProcess, github:")
    ncatt_put(ncout, 0, "date_modified", date())
    ncatt_put(ncout,0, "_FillValue", "1e35")
    ncatt_put(ncout, 0, "featureType", "timeSeriesProfile") #link to oce object? ..... if adp == timeSeriesProfile

    #BODC P01 names
    ncatt_put(ncout, "EWCT", "sdn_parameter_urn", "SDN:P01::LCEWAP01")
    ncatt_put(ncout, "NSCT", "sdn_parameter_urn", "SDN:P01::LCNSAP01")
    ncatt_put(ncout, "VCSP", "sdn_parameter_urn", "SDN:P01::LRZAAP01")
    ncatt_put(ncout, "ERRV", "sdn_parameter_urn", "SDN:P01::LERRAP01")
    ncatt_put(ncout, "BEAM_01", "sdn_parameter_urn", "SDN:P01::THINCE01")
    ncatt_put(ncout, "PGDP_01", "sdn_parameter_urn", "SDN:P01::PCGDAP00")
    ncatt_put(ncout, "lon", "sdn_parameter_urn", "SDN:P01::ALONZZ01")
    ncatt_put(ncout, "lat", "sdn_parameter_urn", "SDN:P01::ALATZZ01")
    ncatt_put(ncout, "EWCT", "sdn_parameter_name", "Eastward current velocity (Eulerian) in the water body by moored acoustic doppler current profiler (ADCP)")
    ncatt_put(ncout, "NSCT", "sdn_parameter_name", "Northward current velocity (Eulerian) in the water body by moored acoustic doppler current profiler (ADCP)")
    ncatt_put(ncout, "VCSP", "sdn_parameter_name", "Upward current velocity in the water body by moored acoustic doppler current profiler (ADCP)")
    ncatt_put(ncout, "ERRV", "sdn_parameter_name", "Current velocity error in the water body by moored acoustic doppler current profiler (ADCP)")
    ncatt_put(ncout, "BEAM_01", "sdn_parameter_name", "Echo intensity from the water body by moored acoustic doppler current profiler (ADCP) beam 1")
    ncatt_put(ncout, "PGDP_01", "sdn_parameter_name", "Acceptable proportion of signal returns by moored acoustic doppler current profiler (ADCP) beam 1")
    ncatt_put(ncout, "lon", "sdn_parameter_name", "Longitude east")
    ncatt_put(ncout, "lat", "sdn_parameter_name", "Latitude north")


    #P06
    ncatt_put(ncout, "EWCT", "sdn_uom_urn", "SDN:P06::UVAA")
    ncatt_put(ncout, "NSCT", "sdn_uom_urn", "SDN:P06::UVAA")
    ncatt_put(ncout, "VCSP", "sdn_uom_urn", "SDN:P06::UVAA")
    ncatt_put(ncout, "ERRV", "sdn_uom_urn", "SDN:P06::UVAA")
    ncatt_put(ncout, "BEAM_01", "sdn_uom_urn", "SDN:P06::UCNT")
    ncatt_put(ncout, "PGDP_01", "sdn_uom_urn", "SDN:P06::UPCT")
    ncatt_put(ncout, "lon", "sdn_uom_urn", "SDN:P06::DEGE")
    ncatt_put(ncout, "lat", "sdn_uom_urn", "SDN:P06:DEGN")

    ncatt_put(ncout, "EWCT", "sdn_uom_name", "Metres per second")
    ncatt_put(ncout, "NSCT", "sdn_uom_name", "Metres per second")
    ncatt_put(ncout, "VCSP", "sdn_uom_name", "Metres per second")
    ncatt_put(ncout, "ERRV", "sdn_uom_name", "Metres per second")
    ncatt_put(ncout, "BEAM_01", "sdn_uom_name", "Counts")
    ncatt_put(ncout, "PGDP_01", "sdn_uom_name", "Percent")
    ncatt_put(ncout, "lon", "sdn_uom_name", "Degrees east")
    ncatt_put(ncout, "lat", "sdn_uom_name", "Degrees north")

    #CF standard names
    ncatt_put(ncout, "EWCT", "standard_name", "eastward_sea_water_velocity")
    ncatt_put(ncout, "NSCT", "standard_name", "northward_sea_water_velocity")
    ncatt_put(ncout, "VCSP", "standard_name", "upward_sea_water_velocity")


    ncatt_put(ncout, "lat", "standard_name", "latitude")
    ncatt_put(ncout, "lon", "standard_name", "longitude")

  }
  if(!is.null(adp[['publisher_name']])){
    ncatt_put(ncout, 0, "publisher_name", adp[['publisher_name']])
  }
  if(!is.null(adp[['publisher_url']])){
    ncatt_put(ncout, 0, "publisher_url", adp[['publisher_url']])
  }
  if(!is.null(adp[['publisher_email']])){
    ncatt_put(ncout, 0, "publisher_email", adp[['publisher_email']])
  }

  ####
  ncatt_put(ncout, "EWCT", "data_max", max(adp[['v']][,,1], na.rm = TRUE))
  ncatt_put(ncout, "EWCT", "data_min", min(adp[['v']][,,1], na.rm = TRUE))
  ncatt_put(ncout, "EWCT", "valid_max", 1000)
  ncatt_put(ncout, "EWCT", "valid_min", -1000)

  ncatt_put(ncout, "NSCT", "data_max", max(adp[['v']][,,2], na.rm = TRUE))
  ncatt_put(ncout, "NSCT", "data_min", min(adp[['v']][,,2], na.rm = TRUE))
  ncatt_put(ncout, "NSCT", "valid_max", 1000)
  ncatt_put(ncout, "NSCT", "valid_min", -1000)

  ncatt_put(ncout, "VCSP", "data_max", max(adp[['v']][,,3], na.rm = TRUE))
  ncatt_put(ncout, "VCSP", "data_min", min(adp[['v']][,,3], na.rm = TRUE))
  ncatt_put(ncout, "VCSP", "valid_max", 1000)
  ncatt_put(ncout, "VCSP", "valid_min", -1000)

  ncatt_put(ncout, "ERRV", "data_max", max(adp[['v']][,,4], na.rm = TRUE))
  ncatt_put(ncout, "ERRV", "data_min", min(adp[['v']][,,4], na.rm = TRUE))
  ncatt_put(ncout, "ERRV", "valid_max", 2000)
  ncatt_put(ncout, "ERRV", "valid_min", -2000)

  if(adp@metadata$source == 'raw'){
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
    ncatt_put(ncout, "PTCH", "data_min", min(adp[['pitch']]))
    ncatt_put(ncout, "PTCH", "data_max", max(adp[['pitch']]))
    ncatt_put(ncout, "ROLL", "data_min", min(adp[['roll']]))
    ncatt_put(ncout, "ROLL", "data_max", max(adp[['roll']]))
    ncatt_put(ncout, "HEAD", "data_min", min(adp[['heading']]))
    ncatt_put(ncout, "HEAD", "data_max", max(adp[['heading']]))
    ncatt_put(ncout, "PRES", "data_min", min(adp[['pressure']]))
    ncatt_put(ncout, "PRES", "data_max", max(adp[['pressure']]))
    ncatt_put(ncout, "SVEL", "data_min", min(adp[['soundSpeed']]))
    ncatt_put(ncout, "SVEL", "data_max", max(adp[['soundSpeed']]))

  }
  if( adp@metadata$source == 'odf'){
    ncatt_put(ncout, "BEAM_01", "data_min", min(adp[['a', 'numeric']], na.rm= TRUE))
    ncatt_put(ncout, "BEAM_01", "data_max", max(adp[['a', 'numeric']], na.rm= TRUE))
    ncatt_put(ncout, "PGDP_01", "data_min", min(adp[['q', 'numeric']], na.rm= TRUE))
    ncatt_put(ncout, "PGDP_01", "data_max", max(adp[['q', 'numeric']], na.rm= TRUE))

  }


  if (!missing(metadata)) {
    metad <- read.csv(metadata, header = TRUE)

    mn <- as.character(metad[,1])
    mv <- as.character(metad[,2])


    md <- as.list(mv)
    names(md) <- mn

    for (m in seq_along(md)) {
      ncatt_put(ncout, 0, names(md)[m], md[[m]])
    }
    nc_close(ncout)


  }
}

####adpCombine####

####create netCDF from combined sources####
#from ODF : processed data (u, v, w, werr, beam01, pgdp01, time)
#from archivede netCDF : metadata
#from RAW file: beam2-4, pgdp2-4, ptch, roll, hght, tx, d, heading, pressure, soundspeed
#instrument metadata


#combine all file sources into single adp object
#adp: from ODF list (odf2adp)
#raw file
#archive netCDF

#added metadata to meet standards

#'
#'
#'
#'ADCP Combine
#'
#'@family odf
#'@family NC
#'@family processing
#'
#'@description Combines archived data and metadata from raw, odf and netCDF
#'  sources into adp object which can be exported as netCDF or saved
#'
#'@details \bold{a)	Raw file (.000)}
#'
#'Contains:
#'
#'
#' \bold{Data:} heading, pressure, sound speed data, BEAM_02-_04, PGDP_02-_04, PTCH,
#' ROLL, HGHT, Tx, D
#'
#'
#'\bold{ Metadata: } firmwareVersion, frequency, beamPattern, orientation, beamAngle,
#' numberOfBeams (janus), pingsPerEnsemble, velocityResolution (pred_accuracy),
#' lowCorrThresh (valid_correlation_range), percentGdMinimum
#' (minmax_percent_good), errorVelocityMaximum (error_velocity_threshold),
#' xmitPulseLength (transmit_pulse_length), falseTargetThresh
#' (false_target_reject_values), serialNumber, instrumentType(data_type)
#'
#'
#' b)	\bold{Archived netCDF}
#'
#' Contains:
#'
#'
#'      \bold{ Data:} BEAM_02-_04, PGDP_02-_04, PTCH, ROLL, HGHT, Tx, D
#' ***note this DATA is pulled from the RAW file in case NC file is corrupted or missing
#'
#'
#' \bold{Metadata: }creation_date,  start_time, stop_time, inst_type, history,
#' starting_water_layer, ending_water_layer, depth_note, transform,
#' data_subtype, coord_system, water_mass, pos_const, depth_const, drifter,
#' var_fill, experiment, project, descript, data_cmnt, fill_flag, composite,
#' magnetic_variation, , delta_t_sec, time_between_ping_groups, depth:
#' xducer_offset_from_bottom, depth: bin_size
#'
#'
#'
#' \bold{c)	Archived ODFs (series for each mooring)}
#'
#' Contains:
#'
#'      \bold{Data:} PROCESSED EWCT, NSCT, VCSP, ERRV, BEAM_01 (!which is
#'      actually the average echo intensity), PGDP_04, time, distance
#'
#'
#'     \bold{ Metadata:}  units (v, distance), cellSize, numberOfBeams, orientation,
#'      model, type, serialNumber, ship (platform), scientist
#'      (chief_scientist), institute (data_origin), cruise (cruise_name), station (mooring),
#'      countryInstituteCode, cruiseNumber, startTime, latitude, longitude,
#'      waterDepth, sounding

#'
#'@param adp an adp object sourced from ODF files using
#'  \code{\link[ADCP:odf2adp]{odf2adp}}
#'@param raw a raw ADCP file (.000)
#'@param ncin an archived netCDF file (.nc)
#'
####adpCombine####
adpCombine <- function(adp, raw, ncin = ''){

  a <- read.adp(raw)
  #####pull metadata from RAW####

  firmware_version <- a[['firmwareVersion']]
  frequency <- a[['frequency']]
  beam_pattern <- a[['beamPattern']]
  orientation <- a[['orientation']]
  beam_angle <- a[['beamAngle']]
  janus <- a[['numberOfBeams']]
  pings_per_ensemble <- a[['pingsPerEnsemble']]
  pred_accuracy <- (a[['velocityResolution']]*1000)
  valid_correlation_range <- a[['lowCorrThresh']]
  minmax_percent_good <- a[['percentGdMinimum']]
  error_velocity_threshold <- a[['errorVelocityMaximum']]
  #time_between_ping_groups <- a[['']]
  transmit_pulse_length_cm <-( a[['xmitPulseLength']]*100)
  false_target_reject_values <- a[['falseTargetThresh']]
  serial_number <- a[['serialNumber']]
  data_type <- a[['instrumentType']]

  adp <- oceSetMetadata(adp, 'firmware_version', firmware_version)
  adp <- oceSetMetadata(adp, 'frequency', frequency)
  adp <- oceSetMetadata(adp, 'beam_pattern' , beam_pattern)
  adp <- oceSetMetadata(adp, 'orientation', orientation)
  adp <- oceSetMetadata(adp, 'beam_angle', beam_angle)
  adp <- oceSetMetadata(adp, 'janus', janus)
  adp <- oceSetMetadata(adp, 'pings_per_ensemble', pings_per_ensemble)
  adp <- oceSetMetadata(adp, 'pred_accuracy', pred_accuracy)
  adp <- oceSetMetadata(adp, 'valid_correlation_range', valid_correlation_range)
  adp <- oceSetMetadata(adp, 'minmax_percent_good', minmax_percent_good)
  adp <- oceSetMetadata(adp, 'error_velocity_threshold', error_velocity_threshold)
  #adp <- oceSetMetadata(adp, 'time_between_ping_groups', )
  adp <- oceSetMetadata(adp, 'transmit_pulse_length_cm', transmit_pulse_length_cm)
  adp <- oceSetMetadata(adp, 'false_target_reject_values', false_target_reject_values)
  adp <- oceSetMetadata(adp, 'serial_number', serial_number)
  adp <- oceSetMetadata(adp, 'data_type', data_type)

  #####pull metadata from archive NC####

  if(!missing(ncin)){
    ni <- nc_open(ncin)
    #pull log sheet metadata from incoming netCDF

    creation_date <- ncatt_get(ni, 0, 'CREATION_DATE')
    time_coverage_start <- ncatt_get(ni, 0,   'start_time')
    time_coverage_end <- ncatt_get(ni, 0,  'stop_time')
    inst_type <- ncatt_get(ni, 0, 'INST_TYPE')
    historyadp <- ncatt_get(ni, 0,  'history')
    starting_water_layer <- ncatt_get(ni,  0, 'starting_water_layer')
    ending_water_layer <- ncatt_get(ni, 0,  'ending_water_layer')
    depth_note <- ncatt_get(ni, 0,  'depth_note')

    #     deprecated --- Diana Cardoso 06/01/2018
    #transform <- ncatt_get(ni, 0,  'transform')

    data_subtype <- ncatt_get(ni, 0,  'DATA_SUBTYPE')
    coord_system <- ncatt_get(ni, 0,  'COORD_SYSTEM')
    water_mass <- ncatt_get(ni, 0,  'WATER_MASS')
    pos_const <- ncatt_get(ni, 0,  'POS_CONST')
    depth_const <- ncatt_get(ni, 0,  'DEPTH_CONST')
    drifter <- ncatt_get(ni, 0,  'DRIFTER')
    FillValue <- ncatt_get(ni, 0,  'VAR_FILL')
    experiment <- ncatt_get(ni, 0,  'EXPERIMENT')
    project <- ncatt_get(ni, 0,  'PROJECT')
    description <- ncatt_get(ni, 0,  'DESCRIPT')
    data_comment <- ncatt_get(ni,  0, 'DATA_CMNT')
    fill_flag <- ncatt_get(ni, 0,  'FILL_FLAG')
    composite <- ncatt_get(ni, 0,  'COMPOSITE')
    magnetic_variation <- ncatt_get(ni, 0,  'magnetic_variation')
    delta_t_sec <- ncatt_get(ni, 0, 'DELTA_T_sec')
    ping_interval <- ncatt_get(ni, 0, 'time_between_ping_groups')
    data_origin <- ncatt_get(ni, 0, 'DATA_ORIGIN')
    xducer_offset_from_bottom <- ncatt_get(ni, 'depth', 'xducer_offset_from_bottom')
    bin_size <- ncatt_get(ni, 'depth', 'bin_size')

    nc_close(ni)

    adp <- oceSetMetadata(adp, 'creation_date', creation_date$value)
    adp <- oceSetMetadata(adp, 'data_origin', data_origin$value)
    adp <- oceSetMetadata(adp, 'time_coverage_start', time_coverage_start$value)
    adp <- oceSetMetadata(adp, 'time_coverage_end', time_coverage_end$value)
    adp <- oceSetMetadata(adp, 'inst_type', inst_type$value)
    adp <- oceSetMetadata(adp, 'history', historyadp$value)
    adp <- oceSetMetadata(adp, 'starting_water_layer', starting_water_layer$value)
    adp <- oceSetMetadata(adp, 'ending_water_layer', ending_water_layer$value)
    adp <- oceSetMetadata(adp, 'depth_note', depth_note$value)

    #     deprecated --- Diana Cardoso 06/01/2018
    #adp <- oceSetMetadata(adp, 'transform', transform$value)

    adp <- oceSetMetadata(adp, 'data_subtype', data_subtype$value)
    adp <- oceSetMetadata(adp, 'coord_system', coord_system$value)
    adp <- oceSetMetadata(adp, 'water_mass', water_mass$value)
    adp <- oceSetMetadata(adp, 'pos_const', pos_const$value)
    adp <- oceSetMetadata(adp, 'depth_const', depth_const$value)
    adp <- oceSetMetadata(adp, 'drifter', drifter$value)
    adp <- oceSetMetadata(adp, 'FillValue', FillValue$value)
    adp <- oceSetMetadata(adp, 'experiment', experiment$value)
    adp <- oceSetMetadata(adp, 'project', project$value)
    adp <- oceSetMetadata(adp, 'description', description$value)
    adp <- oceSetMetadata(adp, 'data_comment', data_comment$value)
    adp <- oceSetMetadata(adp, 'fill_flag', fill_flag$value)
    adp <- oceSetMetadata(adp, 'composite', composite$value)
    adp <- oceSetMetadata(adp, 'magnetic_variation', magnetic_variation$value)
    adp <- oceSetMetadata(adp, 'delta_t_sec', delta_t_sec$value)
    adp <- oceSetMetadata(adp, 'xducer_offset_from_bottom', xducer_offset_from_bottom$value)
    adp <- oceSetMetadata(adp, 'bin_size', bin_size$value)
    adp <- oceSetMetadata(adp, 'ping_interval', ping_interval$value)
    adp <- oceSetMetadata(adp, 'sample_interval', pings_per_ensemble * ping_interval$value)



    #set metadata source

    adp <- oceSetMetadata(adp, 'source', 'netCDF, Raw, ODF combined')

  }

  if(missing(ncin)){
    warning('NC file not provided, object is missing metadata')
    warning(
      ' please provide creation_date,  start_time, stop_time, inst_type, history, starting_water_layer, ending_water_layer, depth_note, transform,
      data_subtype, coord_system, water_mass, pos_const, depth_const, drifter,
      var_fill, experiment, project, descript, data_cmnt, fill_flag, composite,
      magnetic_variation, , delta_t_sec, time_between_ping_groups, depth:
      xducer_offset_from_bottom, depth: bin_size'
    )

    adp<- oceSetMetadata(adp, 'source', 'Raw, ODF combined')
  }

  #####pull data from raw file#####
  a <- read.adp(raw)
  BEAM_01 <- a[['a', 'numeric']][,,1]
  BEAM_02 <- a[['a', 'numeric']][,,2]
  BEAM_03 <- a[['a', 'numeric']][,,3]
  BEAM_04 <- a[['a', 'numeric']][,,4]
  PGDP_01 <- a[['g', 'numeric']][,,1]
  PGDP_02 <- a[['g', 'numeric']][,,2]
  PGDP_03 <- a[['g', 'numeric']][,,3]

  PTCH <- a[['pitch']]
  ROLL <- a[['roll']]
  HGHT <- a[['distance']]
  Tx <- a[['temperature']]
  D <- a[['depth']]
  HEAD <- a[['heading']]
  PRES <- a[['pressure']]
  SVEL <- a[['soundSpeed']]


  #####limit dimensions to match odf files####

  ####apply time offset####


  t <-  ( a[['time']] + (adp[['sample_interval']]/2))
  a <- oceSetData(a, 'time', t)

  #limit by time
  limitmat <- matrix(0, nrow = length(a[['time']]), ncol = length(a[['distance']]))
  limitvec <- matrix(0, ncol = length(a[['time']]))



  #create 'flag mask' where 4 = bad vlaue (outside bounds)
  limitmat[as.POSIXct(a[['time']], tz = 'UTC') < as.POSIXct(adp[['time']][[1]], tz = 'UTC') | as.POSIXct(a[['time']], tz = 'UTC') > as.POSIXct(adp[['time']][[length(adp[['time']])]], tz = 'UTC')] <- 4
  limitvec[as.POSIXct(a[['time']], tz = 'UTC') < as.POSIXct(adp[['time']][[1]], tz = 'UTC') | as.POSIXct(a[['time']], tz = 'UTC') > as.POSIXct(adp[['time']][[length(adp[['time']])]], tz = 'UTC')] <- 4


  #limit time variable
  a[['time']][limitvec == 4] <- NA

  #limit other transferable data
  PTCH[limitvec == 4] <- NA
  ROLL[limitvec == 4] <- NA
  Tx[limitvec == 4] <- NA
  D[limitvec == 4] <- NA
  HEAD[limitvec == 4] <- NA
  PRES[limitvec == 4] <- NA
  SVEL[limitvec == 4] <- NA

  BEAM_01[limitmat == 4] <- NA
  BEAM_02[limitmat == 4] <- NA
  BEAM_03[limitmat == 4] <- NA
  BEAM_04[limitmat == 4] <- NA
  PGDP_01[limitmat == 4] <- NA
  PGDP_02[limitmat == 4] <- NA
  PGDP_03[limitmat == 4] <- NA




  #####insert into adp####

  #create an array
  x <- nrow(adp[['a']])
  y <- ncol(adp[['a']])
  z <- 4
  aa <- array(dim = c(x, y, z))


  #combine beams into a single array using dimensions of odf data
  aa[,,1] <- na.omit(BEAM_01[, 1:length(adp[['distance']])])
  aa[,,2] <- na.omit(BEAM_02[, 1:length(adp[['distance']])])
  aa[,,3] <- na.omit(BEAM_03[, 1:length(adp[['distance']])])
  aa[,,4] <- na.omit(BEAM_04[, 1:length(adp[['distance']])])

  #put array into adp object
  adp <- oceSetData(adp, 'a', aa)

  #create a array
  l <- nrow(adp[['q']])
  m <- ncol(adp[['q']])
  n <- 4
  gg <- array(dim = c(l, m, n))

  #combine beams into a single array using dimensions of odf data
  gg[,,1] <- na.omit(PGDP_01[, 1:length(adp[['distance']])])
  gg[,,2] <- na.omit(PGDP_02[, 1:length(adp[['distance']])])
  gg[,,3] <- na.omit(PGDP_03[, 1:length(adp[['distance']])])
  gg[,,4] <- adp[['q', 'numeric']]

  #put array into adp object
  adp <- oceSetData(adp, 'g', gg)

  #insert other data

  adp <- oceSetData(adp, 'pitch', na.omit(PTCH))
  adp <- oceSetData(adp, 'roll', na.omit(ROLL))
  adp <- oceSetData(adp, 'hght', (HGHT[ 1:length(adp[['distance']])]))
  adp <- oceSetData(adp, 'temperature', na.omit(Tx))
  adp <- oceSetData(adp, 'depth', na.omit(D))
  adp <- oceSetData(adp, 'heading', na.omit(HEAD))
  adp <- oceSetData(adp, 'pressure', na.omit(PRES))
  adp <- oceSetData(adp, 'soundSpeed', na.omit(SVEL))

  ####set sensor_depth
  adp <- oceSetMetadata(adp, 'sensor_depth', mean(adp[['depth']]))

  ###fix event qualifier pulled from odf
  adp <- oceSetMetadata(adp, 'eventQualifier', adp[['serialNumber']])
  return(adp)

}


####create netCDF file from combined adp source####
#' NetCDF creation from adp object
#'
#'@family NC
#'
#' @description Creates standardized netCDF file from adp object (produced from
#'   \code{\link[ADCP:adpCombine]{adpCombine}})
#'
#'   Standardized name of file can be created with:
#'   (\code{\link[ADCP:name.file]{name.file}})
#'
#'   product will meet CF compliance, ERDDAP standards, BODC/SDN, and DFO/MEDS
#'   standards
#'
#' @param adp an adp object
#' @param name text string which will name netCDF file

adpNC <- function(adp, name){
  if (!inherits(adp, "adp")){
    stop("method is only for objects of class '", "adp", "'")
  }
  if(missing(name)){
    name <- paste('MADCP', adp[['experiment']], adp[['station']], adp[['serial_number']], adp[['delta_t_sec']], sep = '_')
  }
  #file name and path
  ncpath <- "./"
  ncfname <- paste(ncpath, name, ".nc", sep = "")


  ####setting dimensions and definitions####
  #dimension variables from adp object
  time <- adp[['time']]
  dist <- adp[['distance', 'numeric']]
  lon <- adp[['longitude']]
  lat <- adp[['latitude']]


  #create dimensions
  timedim <- ncdim_def("time", "seconds since 1970-01-01T00:00:00Z", as.double(time))    #time formatting FIX
  depthdim <- ncdim_def("depth", "metres", as.double(dist))
  stationdim <- ncdim_def("station", "counts", as.numeric(adp[['station']]))
  londim <- ncdim_def("lon", "degrees_east" , as.double(lon))
  latdim <- ncdim_def("lat", "degrees_north", as.double(lat))
  dimnchar <- ncdim_def('nchar', '', 1:23, create_dimvar = FALSE)

  #set fill value
  FillValue <- 1e35
  #####define variables####

  dlname <- 'lon'
  lon_def <- ncvar_def(longname= "longitude", units = 'degrees_east', dim = stationdim, name = dlname, prec = 'double')

  dlname <- 'lat'
  lat_def <- ncvar_def( longname = 'latitude', units = 'degrees_north', dim =  stationdim, name = dlname, prec = 'double')

  dlname <- "eastward_sea_water_velocity"
  u_def <- ncvar_def("EWCT", "m/sec", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")

  dlname <- "northward_sea_water_velocity"
  v_def <- ncvar_def("NSCT", "m/sec", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")

  dlname <- "upward_sea_water_velocity"
  w_def <- ncvar_def("VCSP", "m/sec", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")

  dlname <- "time_02"
  t_def <- ncvar_def("ELTMEP01", "seconds since 1970-01-01T00:00:00Z", list( stationdim, timedim), FillValue, dlname, prec = "float")

  dlname <- "error_velocity_in_sea_water"
  e_def <- ncvar_def("ERRV", "m/sec", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")

  dlname <- "ADCP_echo_intensity_beam_1"

  b1_def <- ncvar_def("BEAM_01", "counts", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")

  dlname <- "ADCP_echo_intensity_beam_2"
  b2_def <- ncvar_def("BEAM_02", "counts", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")

  dlname <- "ADCP_echo_intensity_beam_3"
  b3_def <- ncvar_def("BEAM_03", "counts", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")

  dlname <- "ADCP_echo_intensity_beam_4"
  b4_def <- ncvar_def("BEAM_04", "counts", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")


  dlname <- "percent_good_beam_1"
  pg1_def <- ncvar_def("PGDP_01", "counts", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")

  dlname <- "percent_good_beam_2"
  pg2_def <- ncvar_def("PGDP_02", "counts", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")

  dlname <- "percent_good_beam_3"
  pg3_def <- ncvar_def("PGDP_03", "counts", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")

  dlname <- "percent_good_beam_4"
  pg4_def <- ncvar_def("PGDP_04", "counts", list(timedim, depthdim, stationdim), FillValue, dlname, prec = "float")

  dlname <- "pitch"
  p_def <- ncvar_def("PTCH", "degrees", list(  timedim, stationdim), FillValue, dlname, prec = "float")

  dlname <- "roll"
  r_def <- ncvar_def("ROLL", "degrees", list(  timedim , stationdim), FillValue, dlname, prec = "float")

  dlname <- "height of sea surface"
  hght_def <- ncvar_def("hght", "m", list(  depthdim, stationdim ), FillValue, dlname, prec = "float")

  dlname <- "ADCP Transducer Temp."
  Tx_def <- ncvar_def("Tx", "degrees", list( timedim , stationdim), FillValue, dlname, prec = "float")

  dlname <- "instrument depth"
  D_def <- ncvar_def("D", "m", list(timedim, stationdim), FillValue, dlname, prec = "float")

  dlname <- "heading"
  head_def <- ncvar_def("HEAD", "degrees", list(timedim, stationdim), FillValue, dlname, prec = "float")

  dlname <- "pressure"
  pres_def <- ncvar_def("PRES", "decibars", list(timedim, stationdim), FillValue, dlname, prec = "float")

  dlname <- "speed of sound"
  svel_def <- ncvar_def("SVEL", "m/s", list(timedim, stationdim), FillValue, dlname, prec = "float")

  dlname <- "time_string"
  ts_def <- ncvar_def("DTUT8601", units = "",dim =  list(timedim, dimnchar), missval = NULL, name =  dlname, prec = "char")


  #####write out definitions to new nc file####
  ncout <- nc_create(ncfname, list(u_def, v_def, w_def, e_def, t_def, b1_def, b2_def, b3_def, b4_def, pg1_def, pg2_def, pg3_def, pg4_def, p_def, r_def, hght_def, Tx_def, D_def, lon_def, lat_def, head_def, pres_def, svel_def, ts_def), force_v4 = TRUE)
  ncvar_put(ncout, u_def, adp[['v']][,,1])
  ncvar_put(ncout, v_def, adp[['v']][,,2])
  ncvar_put(ncout, w_def, adp[['v']][,,3])
  ncvar_put(ncout, e_def, adp[['v']][,,4])
  ncvar_put(ncout, t_def, as.POSIXct(adp[['time']], tz = 'UTC', origin = '1970-01-01 00:00:00'))
  ncvar_put(ncout, lon_def, adp[['longitude']])
  ncvar_put(ncout, lat_def, adp[['latitude']])
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
  ncvar_put(ncout, hght_def, (adp[['sensor_depth']]- adp[['distance']]))
  ncvar_put(ncout, Tx_def, adp[['temperature']])
  ncvar_put(ncout, D_def, adp[['depth']])
  ncvar_put(ncout, head_def, adp[['heading']])
  ncvar_put(ncout, pres_def, adp[['pressure']])
  ncvar_put(ncout, svel_def, adp[['soundSpeed']])
  ncvar_put(ncout, ts_def, adp[['time']])

  ####metadata####
  ####dimensions####
  ncatt_put(ncout, 'station', attname = 'cf_role',attval =  'timeseries_id')
  ncatt_put(ncout, 'station', 'longitude', adp[['longitude']])
  ncatt_put(ncout, 'station', 'latitiude', adp[['latitude']])
  ncatt_put(ncout, 'time', attname = 'cf_role', attval = 'profile_id')
  ncatt_put(ncout, 'station', 'standard_name', 'platform_name')
  ncatt_put(ncout, 'time' , 'calendar', 'gregorian')
  ncatt_put(ncout, 'time_string', 'note', 'time values as ISO8601 string, YY-MM-DD hh:mm:ss')
  ncatt_put(ncout, 'time_string', 'time_zone', 'UTC')
  ####global####
  ncatt_put(ncout, 0, "mooring_number", adp[['station']])


  #     deprecated --- Diana Cardoso 06/01/2018
  #ncatt_put(ncout, 0, "deployment_date", adp[['deployment_date']])
  #ncatt_put(ncout, 0, "recovery_date", adp[['recovery_date']])


  ncatt_put(ncout, 0, "firmware_version", adp[['firmware_version']])
  ncatt_put(ncout, 0, "frequency", adp[['frequency']])
  ncatt_put(ncout, 0, "beam_pattern", adp[['beam_pattern']])
  ncatt_put(ncout, 0, "janus", adp[['janus']])
  ncatt_put(ncout, 0, "pings_per_ensemble", adp[['pings_per_ensemble']])
  ncatt_put(ncout, 0, "valid_correlation_range", adp[['valid_correlation_range']])
  ncatt_put(ncout, 0, "minmax_percent_good", adp[['minmax_percent_good']])
  ncatt_put(ncout, 0,"minmax_percent_good", "100")
  ncatt_put(ncout, 0, "error_velocity_threshold", adp[['error_velocity_threshold']])
  ncatt_put(ncout, 0, "transmit_pulse_length_cm", adp[['transmit_pulse_length_cm']])
  ncatt_put(ncout, 0, "false_target_reject_values", adp[['false_target_reject_values']])
  ncatt_put(ncout, 0, "serial_number", adp[['serial_number']])

  #     deprecated --- Diana Cardoso 06/01/2018
  #ncatt_put(ncout, 0, "transform", adp[['transform']])

  ncatt_put(ncout, 0, "data_type", adp[['data_type']])
  ncatt_put(ncout, 0, "data_subtype", adp[['data_subtype']])
  ncatt_put(ncout, 0, "coord_system", adp[['coord_system']])
  ncatt_put(ncout, 0, "longitude", adp[['longitude']])
  ncatt_put(ncout, 0, "latitude", adp[['latitude']])
  ncatt_put(ncout, 0, "magnetic_variation", adp[['magnetic_variation']])
  ncatt_put(ncout, 0, "platform", adp[['ship']])
  ncatt_put(ncout, 0, "sounding", adp[['sounding']])
  ncatt_put(ncout, 0, "chief_scientist", adp[['scientist']])

  ncatt_put(ncout, 0, "water_depth", adp[['sounding']])
  ncatt_put(ncout, 0, "delta_t_sec",adp[['delta_t_sec']])
  ncatt_put(ncout, 0, "pred_accuracy", adp[['pred_accuracy']])
  ncatt_put(ncout, 0, "history", adp[['history']])
  ncatt_put(ncout, 0, "starting_water_layer", adp[['starting_water_layer']])
  ncatt_put(ncout, 0, "ending_water_layer", adp[['ending_water_layer']])
  ncatt_put(ncout, 0, "pos_const", adp[['pos_const']])
  ncatt_put(ncout, 0, "depth_const", adp[['depth_const']])
  ncatt_put(ncout, 0, "drifter", adp[['drifter']])
  ncatt_put(ncout, 0, "experiment", adp[['experiment']])
  ncatt_put(ncout, 0, "cruise_name", adp[['cruise']])





  ####variables####

  ncatt_put(ncout, "depth", "xducer_offset_from_bottom", adp[['xducer_offset_from_bottom']])
  ncatt_put(ncout, "depth", "bin_size", adp[['bin_size']])

  ncatt_put(ncout, "EWCT", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "EWCT", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "EWCT", "serial_number", adp[['serial_number']])
  ncatt_put(ncout, "NSCT", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "NSCT", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "NSCT", "serial_number", adp[['serial_number']])
  ncatt_put(ncout, "VCSP", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "VCSP", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "VCSP", "serial_number", adp[['serial_number']])
  ncatt_put(ncout, "ERRV", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "ERRV", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "ERRV", "serial_number", adp[['serial_number']])
  ncatt_put(ncout, "BEAM_01", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "BEAM_01", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "BEAM_01", "serial_number", adp[['serial_number']])
  ncatt_put(ncout, "BEAM_02", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "BEAM_02", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "BEAM_02", "serial_number", adp[['serial_number']])
  ncatt_put(ncout, "BEAM_03", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "BEAM_03", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "BEAM_03", "serial_number", adp[['serial_number']])
  ncatt_put(ncout, "BEAM_04", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "BEAM_04", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "BEAM_04", "serial_number", adp[['serial_number']])
  ncatt_put(ncout, "PGDP_01", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "PGDP_01", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "PGDP_01", "serial_number", adp[['serial_number']])
  ncatt_put(ncout, "PGDP_02", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "PGDP_02", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "PGDP_02", "serial_number", adp[['serial_number']])
  ncatt_put(ncout, "PGDP_03", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "PGDP_03", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "PGDP_03", "serial_number", adp[['serial_number']])
  ncatt_put(ncout, "PGDP_04", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "PGDP_04", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "PGDP_04", "serial_number", adp[['serial_number']])
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
  ncatt_put(ncout, "hght", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "hght", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "hght", "serial_number", adp[['serial_number']])
  ncatt_put(ncout, "D", "generic_name", "depth")
  ncatt_put(ncout, "D", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "D", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "D", "serial_number", adp[['serial_number']])
  ncatt_put(ncout, "Tx", "generic_name", "temp")
  ncatt_put(ncout, "Tx", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "Tx", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "Tx", "serial_number", adp[['serial_number']])
  ncatt_put(ncout, "HEAD", "generic_name", "heading")
  ncatt_put(ncout, "HEAD", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "HEAD", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "HEAD", "serial_number", adp[['serial_number']])
  ncatt_put(ncout, "PRES", "generic_name", "pressure")
  ncatt_put(ncout, "PRES", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "PRES", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "PRES", "serial_number", adp[['serial_number']])
  ncatt_put(ncout, "SVEL", "generic_name", "sound speed")
  ncatt_put(ncout, "SVEL", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "SVEL", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "SVEL", "serial_number", adp[['serial_number']])

  ####CF conventions & BODC standards####
  ncatt_put(ncout, 0, 'Conventions', 'CF-1.7')
  ncatt_put(ncout, 0, "creator_type", "person")

  ncatt_put(ncout, 0, "program", adp[['description']])
  ncatt_put(ncout, 0, "time_coverage_start", adp[['time_coverage_start']])
  ncatt_put(ncout, 0, "time_coverage_end", adp[['time_coverage_end']])
  ncatt_put(ncout, 0, "geospatial_lat_min", adp[['latitude']])
  ncatt_put(ncout, 0, "geospatial_lat_max", adp[['latitude']])
  ncatt_put(ncout, 0, "geospatial_lat_units", "degrees_north")
  ncatt_put(ncout, 0, "geospatial_lon_min", adp[['longitude']])
  ncatt_put(ncout, 0, "geospatial_lon_max", adp[['longitude']])
  ncatt_put(ncout, 0, "geospatial_lon_units", "degrees_east")

  if (adp[['orientation']] == 'up' ){
    ncatt_put(ncout, 0, "geospatial_vertical_min", adp[['sensor_depth']] + max(adp[['distance']], na.rm = TRUE))
    ncatt_put(ncout, 0, "geospatial_vertical_max", adp[['sensor_depth']] + min(adp[['distance']], na.rm = TRUE))
  }
  if (adp[['orientation']] == 'down' ){
    ncatt_put(ncout, 0, "geospatial_vertical_min", adp[['sensor_depth']] + min(adp[['distance']], na.rm = TRUE))
    ncatt_put(ncout, 0, "geospatial_vertical_max", adp[['sensor_depth']] + max(adp[['distance']], na.rm = TRUE))
  }
  ncatt_put(ncout, 0, "geospatial_vertical_units", "metres")
  ncatt_put(ncout, 0, "geospatial_vertical_positive", 'down')

  ncatt_put(ncout, 0, "project", adp[['project']])
  ncatt_put(ncout,0, "_FillValue", "1e35")
  ncatt_put(ncout, 0, "featureType", "timeSeriesProfile")
  ncatt_put(ncout, 0, "date_modified", date())

  #added meta to meet conventions (not found in archive) #to be inserted manually
  #??????
  if(!is.null(adp[['sea_name']])){
    ncatt_put(ncout, 0, "sea_name", adp[['sea_name']])
  }
  if(!is.null(adp[['creator_name']])){
    ncatt_put(ncout, 0, "creator_name", adp[['creator_name']])
  }
  if( !is.null(adp[['creator_url']])){
    ncatt_put(ncout, 0, "creator_url", adp[['creator_url']])
  }
  if(!is.null(adp[['creator_email']])){
    ncatt_put(ncout, 0, "creator_email", adp[['creator_email']])
  }
  if(!is.null(adp[['processing_level']])){
  ncatt_put(ncout, 0, "processing_level", adp[['processing_level']])
  }
  ncatt_put(ncout, 0, "source", "R code: adcpProcess, github:") ##update with link to code
  if(!is.null(adp[['publisher_name']])){
  ncatt_put(ncout, 0, "publisher_name", adp[['publisher_name']])
  }
  if(!is.null(adp[['publisher_url']])){
  ncatt_put(ncout, 0, "publisher_url", adp[['publisher_url']])
  }
  if(!is.null(adp[['publisher_email']])){
  ncatt_put(ncout, 0, "publisher_email", adp[['publisher_email']])
  }

  #redundant metadata to be condensed
  ncatt_put(ncout, 0, "institution", adp[['data_origin']])
  ncatt_put(ncout, 0, "creator_institution", adp[['data_origin']])
  ncatt_put(ncout, 0, "data_origin", adp[['data_origin']])

  ####BODC P01 names####
  ncatt_put(ncout, "EWCT", "sdn_parameter_urn", "SDN:P01::LCEWAP01")
  ncatt_put(ncout, "NSCT", "sdn_parameter_urn", "SDN:P01::LCNSAP01")
  ncatt_put(ncout, "VCSP", "sdn_parameter_urn", "SDN:P01::LRZAAP01")
  ncatt_put(ncout, "ERRV", "sdn_parameter_urn", "SDN:P01::LERRAP01")
  ncatt_put(ncout, "BEAM_01", "sdn_parameter_urn", "SDN:P01::TNIHCE01")
  ncatt_put(ncout, "BEAM_02", "sdn_parameter_urn", "SDN:P01::TNIHCE02")
  ncatt_put(ncout, "BEAM_03", "sdn_parameter_urn", "SDN:P01::TNIHCE03")
  ncatt_put(ncout, "BEAM_04", "sdn_parameter_urn", "SDN:P01::TNIHCE04")
  ncatt_put(ncout, "PGDP_01", "sdn_parameter_urn", "SDN:P01::PCGDAP00")
  ncatt_put(ncout, "PGDP_02", "sdn_parameter_urn", "SDN:P01::PCGDAP02")
  ncatt_put(ncout, "PGDP_03", "sdn_parameter_urn", "SDN:P01::PCGDAP03")
  ncatt_put(ncout, "PGDP_04", "sdn_parameter_urn", "SDN:P01::PCGDAP04")
  #ncatt_put(ncout, "hght", "sdn_parameter_urn", "SDN:P01::")
  ncatt_put(ncout, "D", "sdn_parameter_urn", "SDN:P01::ADEPZZ01")
  ncatt_put(ncout, "Tx", "sdn_parameter_urn", "SDN:P01::TEMPPR01")
  ncatt_put(ncout, "ELTMEP01", "sdn_parameter_urn", "SDN:P01::ELTMEP01")
  ncatt_put(ncout, "PTCH", "sdn_parameter_urn", "SDN:P01::PTCHEI01")
  ncatt_put(ncout, "ROLL", "sdn_parameter_urn", "SDN:P01::ROLLFEI01")
  ncatt_put(ncout, "lon", "sdn_parameter_urn", "SDN:P01::ALONZZ01")
  ncatt_put(ncout, "lat", "sdn_parameter_urn", "SDN:P01::ALATZZ01")
  ncatt_put(ncout, "HEAD", "sdn_parameter_urn", "SDN:P01::HEADCM01")
  ncatt_put(ncout, "PRES", "sdn_parameter_urn", "SDN:P01::PRESPR01")
  ncatt_put(ncout, "SVEL", "sdn_parameter_urn", "SDN:P01::SVELCV01")
  ncatt_put(ncout, "time_string", "sdn_parameter_urn", "SDN:P01::DTUT8601")


  ncatt_put(ncout, "EWCT", "sdn_parameter_name", "Eastward current velocity (Eulerian) in the water body by moored acoustic doppler current profiler (ADCP)")
  ncatt_put(ncout, "NSCT", "sdn_parameter_name", "Northward current velocity (Eulerian) in the water body by moored acoustic doppler current profiler (ADCP)")
  ncatt_put(ncout, "VCSP", "sdn_parameter_name", "Upward current velocity in the water body by moored acoustic doppler current profiler (ADCP)")
  ncatt_put(ncout, "ERRV", "sdn_parameter_name", "Current velocity error in the water body by moored acoustic doppler current profiler (ADCP)")
  ncatt_put(ncout, "BEAM_01", "sdn_parameter_name", "Echo intensity from the water body by moored acoustic doppler current profiler (ADCP) beam 1")
  ncatt_put(ncout, "BEAM_02", "sdn_parameter_name", "Echo intensity from the water body by moored acoustic doppler current profiler (ADCP) beam 2")
  ncatt_put(ncout, "BEAM_03", "sdn_parameter_name", "Echo intensity from the water body by moored acoustic doppler current profiler (ADCP) beam 3")
  ncatt_put(ncout, "BEAM_04", "sdn_parameter_name", "Echo intensity from the water body by moored acoustic doppler current profiler (ADCP) beam 4")
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
  ncatt_put(ncout, "HEAD", "sdn_parameter_name", "Orientation (horizontal relative to true north) of measurement device {heading}")
  ncatt_put(ncout, "PRES", "sdn_parameter_name", "Pressure (spatial co-ordinate) exerted by the water body by profiling pressure sensor and corrected to read zero at sea level")
  ncatt_put(ncout, "SVEL", "sdn_parameter_name", "Sound velocity in the water body by computation from temperature and salinity by unspecified algorithm")
  ncatt_put(ncout, 'ELTMEP01', "sdn_parameter_name", "Elapsed time (since 1970-01-01T00:00:00Z)")
  ncatt_put(ncout, 'time_string', "sdn_parameter_name", "String corresponding to format 'YYYY-MM-DDThh:mm:ss.sssZ' or other valid ISO8601 string")


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
  ncatt_put(ncout, "HEAD", "sdn_uom_urn", "SDN:P06:UAAA")
  ncatt_put(ncout, "PRES", "sdn_uom_urn", "SDN:P06:UPDB")
  ncatt_put(ncout, "SVEL", "sdn_uom_urn", "SDN:P06:UVAA")
  ncatt_put(ncout, "ELTMEP01", "sdn_uom_urn", "SDN:P06::UTBB")
  ncatt_put(ncout, "time_string", "sdn_uom_urn", "SDN:P06::TISO")

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
  ncatt_put(ncout, "HEAD", "sdn_uom_name", "Degrees")
  ncatt_put(ncout, "PRES", "sdn_uom_name", "Decibars")
  ncatt_put(ncout, "SVEL", "sdn_uom_name", "Metres per second")
  ncatt_put(ncout, "ELTMEP01", "sdn_uom_name", "Seconds")
  ncatt_put(ncout, "time_string", "sdn_uom_name", "ISO8601")


  #####CF standard names####
  ncatt_put(ncout, "EWCT", "standard_name", "eastward_sea_water_velocity")
  ncatt_put(ncout, "NSCT", "standard_name", "northward_sea_water_velocity")
  ncatt_put(ncout, "VCSP", "standard_name", "upward_sea_water_velocity")
  ncatt_put(ncout, "ELTMEP01", "standard_name", "time")
  ncatt_put(ncout, "lat", "standard_name", "latitude")
  ncatt_put(ncout, "lon", "standard_name", "longitude")
  ncatt_put(ncout, "D", "standard_name", "depth")
  #ncatt_put(ncout, "Tx", "standard_name", "")
  ncatt_put(ncout, "PTCH", "standard_name", "platform_pitch_angle")
  ncatt_put(ncout, "ROLL", "standard_name", "platform_roll_angle")
  ncatt_put(ncout, "PRES", "standard_name", "sea_water_pressure")
  ncatt_put(ncout, "SVEL", "standard_name", "speed_of_sound_in_sea_water")


  ####data max and min####
  ncatt_put(ncout, "EWCT", "data_max", max(adp[['v']][,,1], na.rm = TRUE))
  ncatt_put(ncout, "EWCT", "data_min", min(adp[['v']][,,1], na.rm = TRUE))
  ncatt_put(ncout, "EWCT", "valid_max", 1000)
  ncatt_put(ncout, "EWCT", "valid_min", -1000)

  ncatt_put(ncout, "NSCT", "data_max", max(adp[['v']][,,2], na.rm = TRUE))
  ncatt_put(ncout, "NSCT", "data_min", min(adp[['v']][,,2], na.rm = TRUE))
  ncatt_put(ncout, "NSCT", "valid_max", 1000)
  ncatt_put(ncout, "NSCT", "valid_min", -1000)

  ncatt_put(ncout, "VCSP", "data_max", max(adp[['v']][,,3], na.rm = TRUE))
  ncatt_put(ncout, "VCSP", "data_min", min(adp[['v']][,,3], na.rm = TRUE))
  ncatt_put(ncout, "VCSP", "valid_max", 1000)
  ncatt_put(ncout, "VCSP", "valid_min", -1000)

  ncatt_put(ncout, "ERRV", "data_max", max(adp[['v']][,,4], na.rm = TRUE))
  ncatt_put(ncout, "ERRV", "data_min", min(adp[['v']][,,4], na.rm = TRUE))
  ncatt_put(ncout, "ERRV", "valid_max", 2000)
  ncatt_put(ncout, "ERRV", "valid_min", -2000)

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

  ncatt_put(ncout, "PTCH", "data_min", min(adp[['pitch']]))
  ncatt_put(ncout, "PTCH", "data_max", max(adp[['pitch']]))

  ncatt_put(ncout, "ROLL", "data_min", min(adp[['roll']]))
  ncatt_put(ncout, "ROLL", "data_max", max(adp[['roll']]))

  ncatt_put(ncout, "HEAD", "data_min", min(adp[['heading']]))
  ncatt_put(ncout, "HEAD", "data_max", max(adp[['heading']]))

  ncatt_put(ncout, "PRES", "data_min", min(adp[['pressure']]))
  ncatt_put(ncout, "PRES", "data_max", max(adp[['pressure']]))

  ncatt_put(ncout, "SVEL", "data_min", min(adp[['soundSpeed']]))
  ncatt_put(ncout, "SVEL", "data_max", max(adp[['soundSpeed']]))


  ####nc close####
  nc_close(ncout)




}





