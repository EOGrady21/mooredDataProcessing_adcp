####netCDF from ONLY ODF data#####
oceNc_create_odf <- function(adp, name,  metadata){
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

#set fill value
FillValue <- 1e35

#define variables

dlname <- 'lon'
lon_def <- ncvar_def(longname= "longitude", units = 'degrees_east', dim = stationdim, name = dlname, prec = 'double')

dlname <- 'lat'
lat_def <- ncvar_def( longname = 'latitude', units = 'degrees_north', dim =  stationdim, name = dlname, prec = 'double')

dlname <- "eastward_sea_water_velocity"
u_def <- ncvar_def("EWCT", "m/sec", list(stationdim, depthdim, timedim), FillValue, dlname, prec = "float")

dlname <- "northward_sea_water_velocity"
v_def <- ncvar_def("NSCT", "m/sec", list(stationdim, depthdim, timedim), FillValue, dlname, prec = "float")

dlname <- "upward_sea_water_velocity"
w_def <- ncvar_def("VCSP", "m/sec", list(stationdim, depthdim, timedim), FillValue, dlname, prec = "float")

dlname <- "time_02"
t_def <- ncvar_def("SYTM", "seconds since 1970-01-01T00:00:00Z", list( stationdim, timedim), FillValue, dlname, prec = "float")

dlname <- "error_velocity_in_sea_water"
e_def <- ncvar_def("ERRV", "m/sec", list(stationdim, depthdim, timedim), FillValue, dlname, prec = "float")

dlname <- "ADCP_echo_intensity_beam_1"

b1_def <- ncvar_def("BEAM_01", "counts", list(stationdim, depthdim, timedim), FillValue, dlname, prec = "float")


dlname <- "percent_good_beam_1"
pg1_def <- ncvar_def("PGDP_01", "counts", list(stationdim, depthdim, timedim), FillValue, dlname, prec = "float")

####writing net CDF####
#write out definitions to new nc file
ncout <- nc_create(ncfname, list(u_def, v_def, w_def, e_def, t_def, b1_def,  pg1_def, lon_def, lat_def), force_v4 = TRUE)


}


#insert variables into nc file
ncvar_put(ncout, u_def, adp[['v']][,,1])
ncvar_put(ncout, v_def, adp[['v']][,,2])
ncvar_put(ncout, w_def, adp[['v']][,,3])
ncvar_put(ncout, e_def, adp[['v']][,,4])
ncvar_put(ncout, t_def, as.POSIXct(adp[['time']], tz = 'UTC', origin = '1970-01-01 00:00:00'))
ncvar_put(ncout, lon_def, adp[['longitude']])
ncvar_put(ncout, lat_def, adp[['latitude']])
ncvar_put(ncout, b1_def, adp[['a', 'numeric']])
ncvar_put(ncout, pg1_def, adp[['q', 'numeric']])

###metadata###

ncatt_put(ncout, 'station', attname = 'cf_role',attval =  'timeseries_id')
ncatt_put(ncout, 'time', attname = 'cf_role', attval = 'profile_id')
ncatt_put(ncout, 'station', 'standard_name', 'platform_name')
ncatt_put(ncout, 'time' , 'calendar', 'gregorian')

ncatt_put(ncout, 0, "mooring_number", adp[['mooring_number']])
ncatt_put(ncout, 0, "deployment_date", adp[['deployment_time']])
ncatt_put(ncout, 0, "recovery_date", adp[['recovery_time']])
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
ncatt_put(ncout, 0, "transform", adp[['oceCoordinate']])
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
ncatt_put(ncout, 0, "delta_t_sec",adp[['time_offset']])
ncatt_put(ncout, 0, "pred_accuracy", adp[['velocityResolution']])
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

ncatt_put(ncout, "BEAM_01", "data_min", min(adp[['a', 'numeric']], na.rm= TRUE))
ncatt_put(ncout, "BEAM_01", "data_max", max(adp[['a', 'numeric']], na.rm= TRUE))
ncatt_put(ncout, "PGDP_01", "data_min", min(adp[['q', 'numeric']], na.rm= TRUE))
ncatt_put(ncout, "PGDP_01", "data_max", max(adp[['q', 'numeric']], na.rm= TRUE))


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
