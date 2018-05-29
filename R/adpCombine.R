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


adpCombine <- function(adp, raw, ncin ){



  #pull data from raw file
  a <- read.adp(raw)
  BEAM_02 <- a[['a', 'numeric']][,,2]
  BEAM_03 <- a[['a', 'numeric']][,,3]
  BEAM_04 <- a[['a', 'numeric']][,,4]
  PGDP_02 <- a[['g', 'numeric']][,,2]
  PGDP_03 <- a[['g', 'numeric']][,,3]
  PGDP_04 <- a[['g', 'numeric']][,,4]
  PTCH <- a[['pitch']]
  ROLL <- a[['roll']]
  HGHT <- a[['distance']]
  Tx <- a[['temperature']]
  D <- a[['depth']]
  HEAD <- a[['heading']]
  PRES <- a[['pressure']]
  SVEL <- a[['soundSpeed']]


  #limit dimensions to match odf files

  #limit by depth bins thrown out
  dim <- (length(adp[['distance']]): length(a[['distance']]))
  BEAM_02[,dim] <- NA
  BEAM_03[,dim] <- NA
  BEAM_04[,dim] <- NA
  PGDP_02[,dim] <- NA
  PGDP_03[,dim] <- NA
  PGDP_04[,dim] <- NA
  HGHT[dim] <- NA



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

  BEAM_02[limitmat == 4] <- NA
  BEAM_03[limitmat == 4] <- NA
  BEAM_04[limitmat == 4] <- NA
  PGDP_02[limitmat == 4] <- NA
  PGDP_03[limitmat == 4] <- NA
  PGDP_04[limitmat == 4] <- NA



  #insert into adp

  #create an array
  x <- nrow(adp[['a']])
  y <- ncol(adp[['a']])
  z <- 4
  aa <- array(dim = c(x, y, z))


  #combine beams into a single array using dimensions of odf data
  aa[,,1] <- adp[['a', 'numeric']]
  aa[,,2] <- na.omit(BEAM_02)
  aa[,,3] <- na.omit(BEAM_03)
  aa[,,4] <- na.omit(BEAM_04)

  #put array into adp object
  adp <- oceSetData(adp, 'a', aa)

  #create a array
  l <- nrow(adp[['q']])
  m <- ncol(adp[['q']])
  n <- 4
  qq <- array(dim = c(l, m, n))

  #combine beams into a single array using dimensions of odf data
  qq[,,1] <- adp[['q', 'numeric']]
  qq[,,2] <- na.omit(PGDP_02)
  qq[,,3] <- na.omit(PGDP_03)
  qq[,,4] <- na.omit(PGDP_04)

  #put array into adp object
  adp <- oceSetData(adp, 'q', qq)


  tz <- 'UTC'
  adp <- oceSetData(adp, 'pitch', na.omit(PTCH))


  adp <- oceSetData(adp, 'roll', na.omit(ROLL))
  adp <- oceSetData(adp, 'hght', na.omit(HGHT))
  adp <- oceSetData(adp, 'temperature', na.omit(Tx))
  adp <- oceSetData(adp, 'depth', na.omit(D))
  adp <- oceSetData(adp, 'heading', na.omit(HEAD))
  adp <- oceSetData(adp, 'pressure', na.omit(PRES))
  adp <- oceSetData(adp, 'soundSpeed', na.omit(SVEL))


  #pull metadata from RAW

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
  ADCP_serial_number <- a[['serialNumber']]

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
  adp <- oceSetMetadata(adp, 'ADCP_serial_number', ADCP_serial_number)


  #pull metadata from archive NC
  ni <- nc_open(ncin)
  #pull log sheet metadata from incoming netCDF

  creation_date <- ncatt_get(ni, 0, 'CREATION_DATE')
  mooring <- ncatt_get(ni, 0, 'MOORING')
  deployment_date <- ncatt_get(ni, 0,   'start_time')
  recovery_date <- ncatt_get(ni, 0,  'stop_time')
  inst_type <- ncatt_get(ni, 0, 'INST_TYPE')
  historyadp <- ncatt_get(ni, 0,  'history')
  starting_water_layer <- ncatt_get(ni,  0, 'starting_water_layer')
  ending_water_layer <- ncatt_get(ni, 0,  'ending_water_layer')
  depth_note <- ncatt_get(ni, 0,  'depth_note')
  transform <- ncatt_get(ni, 0,  'transform')
  data_type <- ncatt_get(ni, 0,  'DATA_TYPE')
  data_subtype <- ncatt_get(ni, 0,  'DATA_SUBTYPE')
  data_origin <- ncatt_get(ni,  0, 'DATA_ORIGIN')
  coord_system <- ncatt_get(ni, 0,  'COORD_SYSTEM')
  water_mass <- ncatt_get(ni, 0,  'WATER_MASS')
  pos_const <- ncatt_get(ni, 0,  'POS_CONST')
  depth_const <- ncatt_get(ni, 0,  'DEPTH_CONST')
  drifter <- ncatt_get(ni, 0,  'DRIFTER')
  FillValue <- ncatt_get(ni, 0,  'VAR_FILL')
  experiment <- ncatt_get(ni, 0,  'EXPERIMENT')
  project <- ncatt_get(ni, 0,  'PROJECT')
  description <- ncatt_get(ni, 0,  'DESCRIPT')
  longitude <- ncatt_get(ni, 0,  'longitude')
  latitude <- ncatt_get(ni, 0,  'latitude')
  data_comment <- ncatt_get(ni,  0, 'DATA_CMNT')
  fill_flag <- ncatt_get(ni, 0,  'FILL_FLAG')
  composite <- ncatt_get(ni, 0,  'COMPOSITE')
  magnetic_variation <- ncatt_get(ni, 0,  'magnetic_variation')
  platform <- ncatt_get(ni, 0,  'Platform')
  sounding<- ncatt_get(ni, 0,  'Sounding')
  chief_scientist <- ncatt_get(ni, 0,  'Chief_Scientist')
  delta_t_sec <- ncatt_get(ni, 0, 'DELTA_T_sec')

  xducer_offset_from_bottom <- ncatt_get(ni, 'depth', 'xducer_offset_from_bottom')
  bin_size <- ncatt_get(ni, 'depth', 'bin_size')

  nc_close(ni)

  adp <- oceSetMetadata(adp, 'creation_date', creation_date$value)
  adp <- oceSetMetadata(adp, 'mooring', mooring$value)
  adp <- oceSetMetadata(adp, 'deployment_date', deployment_date$value)
  adp <- oceSetMetadata(adp, 'recovery_date', recovery_date$value)
  adp <- oceSetMetadata(adp, 'inst_type', inst_type$value)
  adp <- oceSetMetadata(adp, 'history', historyadp$value)
  adp <- oceSetMetadata(adp, 'starting_water_layer', starting_water_layer$value)
  adp <- oceSetMetadata(adp, 'ending_water_layer', ending_water_layer$value)
  adp <- oceSetMetadata(adp, 'depth_note', depth_note$value)
  adp <- oceSetMetadata(adp, 'transform', transform$value)
  adp <- oceSetMetadata(adp, 'data_type', data_type$value)
  adp <- oceSetMetadata(adp, 'data_subtype', data_subtype$value)
  adp <- oceSetMetadata(adp, 'data_origin', data_origin$value)
  adp <- oceSetMetadata(adp, 'coord_system', coord_system$value)
  adp <- oceSetMetadata(adp, 'water_mass', water_mass$value)
  adp <- oceSetMetadata(adp, 'pos_const', pos_const$value)
  adp <- oceSetMetadata(adp, 'depth_const', depth_const$value)
  adp <- oceSetMetadata(adp, 'drifter', drifter$value)
  adp <- oceSetMetadata(adp, 'FillValue', FillValue$value)
  adp <- oceSetMetadata(adp, 'experiment', experiment$value)
  adp <- oceSetMetadata(adp, 'project', project$value)
  adp <- oceSetMetadata(adp, 'description', description$value)
  adp <- oceSetMetadata(adp, 'longitude', longitude$value)
  adp <- oceSetMetadata(adp, 'latitude', latitude$value)
  adp <- oceSetMetadata(adp, 'data_comment', data_comment$value)
  adp <- oceSetMetadata(adp, 'fill_flag', fill_flag$value)
  adp <- oceSetMetadata(adp, 'composite', composite$value)
  adp <- oceSetMetadata(adp, 'magnetic_variation', magnetic_variation$value)
  adp <- oceSetMetadata(adp, 'platform', platform$value)
  adp <- oceSetMetadata(adp, 'sounding', sounding$value)
  adp <- oceSetMetadata(adp, 'chief_scientist', chief_scientist$value)
  adp <- oceSetMetadata(adp, 'delta_t_sec', delta_t_sec$value)
  adp <- oceSetMetadata(adp, 'xducer_offset_from_bottom', xducer_offset_from_bottom$value)
  adp <- oceSetMetadata(adp, 'bin_size', bin_size$value)
  adp <- oceSetMetadata(adp, 'sensor_depth', mean(adp[['depth']]))
  return(adp)

}


####create netCDF file from combined adp source####


adpNC <- function(adp, name){
  if (!inherits(adp, "adp")){
    stop("method is only for objects of class '", "adp", "'")
  }
  if(missing(name)){
    name <- paste('MADCP', adp[['experiment']], adp[['mooring']], adp[['ADCP_serial_number']], adp[['delta_t_sec']], sep = '_')
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
  stationdim <- ncdim_def("station", "counts", as.numeric(adp[['mooring']]))
  londim <- ncdim_def("lon", "degrees_east" , as.double(lon))
  latdim <- ncdim_def("lat", "degrees_north", as.double(lat))
  dimnchar <- ncdim_def('nchar', '', 1:23, create_dimvar = FALSE)

  #set fill value
  FillValue <- 1e35
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
  ts_def <- ncvar_def("DTUT8601", units = "",dim =  list(dimnchar, timedim), missval = NULL, name =  dlname, prec = "char")


  #write out definitions to new nc file
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
  ncvar_put(ncout, pg1_def, adp[['q', 'numeric']][,,1])
  ncvar_put(ncout, pg2_def, adp[['q', 'numeric']][,,2])
  ncvar_put(ncout, pg3_def, adp[['q', 'numeric']][,,3])
  ncvar_put(ncout, pg4_def, adp[['q', 'numeric']][,,4])
  ncvar_put(ncout, p_def, adp[['pitch']])
  ncvar_put(ncout, r_def, adp[['roll']])
  ncvar_put(ncout, hght_def, (adp[['sensor_depth']]- adp[['distance']]))
  ncvar_put(ncout, Tx_def, adp[['temperature']])
  ncvar_put(ncout, D_def, adp[['depth']])
  ncvar_put(ncout, head_def, adp[['heading']])
  ncvar_put(ncout, pres_def, adp[['pressure']])
  ncvar_put(ncout, svel_def, adp[['soundSpeed']])
  ncvar_put(ncout, ts_def, adp[['time']])


  ncatt_put(ncout, 'station', attname = 'cf_role',attval =  'timeseries_id')
  ncatt_put(ncout, 'station', 'longitude', adp[['longitude']])
  ncatt_put(ncout, 'station', 'latitiude', adp[['latitude']])
  ncatt_put(ncout, 'time', attname = 'cf_role', attval = 'profile_id')
  ncatt_put(ncout, 'station', 'standard_name', 'platform_name')
  ncatt_put(ncout, 'time' , 'calendar', 'gregorian')
  ncatt_put(ncout, 'time_string', 'note', 'time values as ISO8601 string')
  ncatt_put(ncout, 'time_string', 'time_zone', 'UTC')
  ncatt_put(ncout, 0, "mooring_number", adp[['mooring']])
  ncatt_put(ncout, 0, "deployment_date", adp[['deployment_date']])
  ncatt_put(ncout, 0, "recovery_date", adp[['recovery_date']])
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
  ncatt_put(ncout, 0, "serial_number", adp[['ADCP_serial_number']])
  ncatt_put(ncout, 0, "transform", adp[['transform']])
  ncatt_put(ncout, 0, "data_type", adp[['data_type']])
  ncatt_put(ncout, 0, "data_subtype", adp[['data_subtype']])
  ncatt_put(ncout, 0, "coord_system", adp[['coord_system']])
  ncatt_put(ncout, 0, "longitude", adp[['longitude']])
  ncatt_put(ncout, 0, "latitude", adp[['latitude']])
  ncatt_put(ncout, 0, "magnetic_variation", adp[['magnetic_variation']])
  ncatt_put(ncout, 0, "platform", adp[['platform']])
  ncatt_put(ncout, 0, "sounding", adp[['sounding']])
  ncatt_put(ncout, 0, "chief_scientist", adp[['chief_scientist']])
  ncatt_put(ncout, 0, "data_origin", adp[['data_origin']])
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





  ncatt_put(ncout, "depth", "xducer_offset_from_bottom", adp[['xducer_offset_from_bottom']])
  ncatt_put(ncout, "depth", "bin_size", adp[['bin_size']])

  ncatt_put(ncout, "EWCT", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "EWCT", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "EWCT", "serial_number", adp[['ADCP_serial_number']])
  ncatt_put(ncout, "NSCT", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "NSCT", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "NSCT", "serial_number", adp[['ADCP_serial_number']])
  ncatt_put(ncout, "VCSP", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "VCSP", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "VCSP", "serial_number", adp[['ADCP_serial_number']])
  ncatt_put(ncout, "ERRV", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "ERRV", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "ERRV", "serial_number", adp[['ADCP_serial_number']])
  ncatt_put(ncout, "BEAM_01", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "BEAM_01", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "BEAM_01", "serial_number", adp[['ADCP_serial_number']])
  ncatt_put(ncout, "BEAM_02", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "BEAM_02", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "BEAM_02", "serial_number", adp[['ADCP_serial_number']])
  ncatt_put(ncout, "BEAM_03", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "BEAM_03", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "BEAM_03", "serial_number", adp[['ADCP_serial_number']])
  ncatt_put(ncout, "BEAM_04", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "BEAM_04", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "BEAM_04", "serial_number", adp[['ADCP_serial_number']])
  ncatt_put(ncout, "PGDP_01", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "PGDP_01", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "PGDP_01", "serial_number", adp[['ADCP_serial_number']])
  ncatt_put(ncout, "PGDP_02", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "PGDP_02", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "PGDP_02", "serial_number", adp[['ADCP_serial_number']])
  ncatt_put(ncout, "PGDP_03", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "PGDP_03", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "PGDP_03", "serial_number", adp[['ADCP_serial_number']])
  ncatt_put(ncout, "PGDP_04", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "PGDP_04", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "PGDP_04", "serial_number", adp[['ADCP_serial_number']])
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
  ncatt_put(ncout, "hght", "serial_number", adp[['ADCP_serial_number']])
  ncatt_put(ncout, "D", "generic_name", "depth")
  ncatt_put(ncout, "D", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "D", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "D", "serial_number", adp[['ADCP_serial_number']])
  ncatt_put(ncout, "Tx", "generic_name", "temp")
  ncatt_put(ncout, "Tx", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "Tx", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "Tx", "serial_number", adp[['ADCP_serial_number']])
  ncatt_put(ncout, "HEAD", "generic_name", "heading")
  ncatt_put(ncout, "HEAD", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "HEAD", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "HEAD", "serial_number", adp[['ADCP_serial_number']])
  ncatt_put(ncout, "PRES", "generic_name", "pressure")
  ncatt_put(ncout, "PRES", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "PRES", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "PRES", "serial_number", adp[['ADCP_serial_number']])
  ncatt_put(ncout, "SVEL", "generic_name", "sound speed")
  ncatt_put(ncout, "SVEL", "sensor_type", adp[['inst_type']])
  ncatt_put(ncout, "SVEL", "sensor_depth", adp[['sensor_depth']])
  ncatt_put(ncout, "SVEL", "serial_number", adp[['ADCP_serial_number']])


  ncatt_put(ncout, 0, 'Conventions', 'CF-1.7')
  ncatt_put(ncout, 0, "creator_type", "person")
  ncatt_put(ncout, 0, "creator_institution", adp[['data_origin']])
  ncatt_put(ncout, 0, "program", adp[['description']])
  ncatt_put(ncout, 0, "time_coverage_start", adp[['deployment_date']])
  ncatt_put(ncout, 0, "time_coverage_end", adp[['recovery_date']])
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
  ncatt_put(ncout, 0, "institution", adp[['data_origin']])
  ncatt_put(ncout, 0, "project", adp[['project']])
  ncatt_put(ncout,0, "_FillValue", "1e35")
  ncatt_put(ncout, 0, "featureType", "timeSeriesProfile")
  ncatt_put(ncout, 0, "date_modified", date())

  #added meta to meet conventions (not found in archive) #to be inserted manually
  #??????
  # ncatt_put(ncout, 0, "sea_name", adp[['sea_name']])
  # ncatt_put(ncout, 0, "creator_name", adp[['creator_name']])
  # ncatt_put(ncout, 0, "creator_url", adp[['creator_url']])
  # ncatt_put(ncout, 0, "creator_email", adp[['creator_email']])
  ncatt_put(ncout, 0, "processing_level", adp[['processing_level']])
  # ncatt_put(ncout, 0, "source", "R code: adcpProcess, github:")

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


  #CF standard names
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

  ncatt_put(ncout, "PGDP_01", "data_min", min(adp[['q', 'numeric']][,,1], na.rm= TRUE))
  ncatt_put(ncout, "PGDP_01", "data_max", max(adp[['q', 'numeric']][,,1], na.rm= TRUE))# eg min 25 % good

  ncatt_put(ncout, "PGDP_02", "data_min", min(adp[['q', 'numeric']][,,2], na.rm= TRUE))
  ncatt_put(ncout, "PGDP_02", "data_max", max(adp[['q' ,'numeric']][,,2], na.rm= TRUE))

  ncatt_put(ncout, "PGDP_03", "data_min", min(adp[['q' ,'numeric']][,,3], na.rm= TRUE))
  ncatt_put(ncout, "PGDP_03", "data_max", max(adp[['q', 'numeric']][,,3], na.rm= TRUE))

  ncatt_put(ncout, "PGDP_04", "data_min", min(adp[['q', 'numeric']][,,4], na.rm= TRUE))
  ncatt_put(ncout, "PGDP_04", "data_max", max(adp[['q', 'numeric']][,,4], na.rm= TRUE))

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

  nc_close(ncout)




}



