####create netCDF from combined sources####
#from ODF : processed data (u, v, w, werr, beam01, pgdp01, time)
#from archivede netCDF : metadata
#from RAW file: beam2-4, pgdp2-4, ptch, roll, hght, tx, d, heading, pressure, soundspeed
  #instrument metadata


#combine all file sources into single adp object
#adp: from ODF list (odf2adp)
#raw file
#archive netCDF



adpCombine <- function(adp, raw, ncin){

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

  #insert into adp
  adp <- oceSetData(adp, a[,,2], BEAM_02)
  adp <- oceSetData(adp, a[,,3], BEAM_03)
  adp <- oceSetData(adp, a[,,4], BEAM_04)
  adp <- oceSetData(adp, g[,,2], PGDP_02)
  adp <- oceSetData(adp, g[,,3], PGDP_03)
  adp <- oceSetData(adp, g[,,4], PGDP_04)
  adp <- oceSetData(adp, 'pitch', PTCH)
  adp <- oceSetData(adp, 'roll', ROLL)
  adp <- oceSetData(adp, 'hght', HGHT)
  adp <- oceSetData(adp, 'temperature', Tx)
  adp <- oceSetData(adp, 'depth', D)
  adp <- oceSetData(adp, 'heading', HEAD)
  adp <- oceSetData(adp, 'pressure', PRES)
  adp <- oceSetData(adp, 'soundSpeed', SVEL)


  #pull metadata from archive NC
  ni <- nc_open(ncin)
  #pull log sheet metadata from incoming netCDF

  creation_date <- ncatt_get(ni, 0, 'CREATION_DATE')
  mooring <- ncatt_get(ni, 0, 'MOORING')
  deployment_date <- ncatt_get(ni, 0,   'Dployment_date')
  recovery_date <- ncatt_get(ni, 0,  'Recovery_date')
  inst_type <- ncatt_get(ni, 0, 'INST_TYPE')
  history <- ncatt_get(ni, 0,  'history')
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

  nc_close(ni)

  adp <- oceSetMetadata(adp, 'creation_date', creation_date)
  adp <- oceSetMetadata(adp, 'mooring', mooring)
  adp <- oceSetMetadata(adp, 'deployment_date', deployment_date)
  adp <- oceSetMetadata(adp, 'recovery_date', recovery_date)
  adp <- oceSetMetadata(adp, 'inst_type', inst_type)
  adp <- oceSetMetadata(adp, 'history', history)
  adp <- oceSetMetadata(adp, 'starting_water_layer', starting_water_layer)
  adp <- oceSetMetadata(adp, 'ending_water_layer', ending_water_layer)
  adp <- oceSetMetadata(adp, 'depth_note', depth_note)
  adp <- oceSetMetadata(adp, 'transform', transform)
  adp <- oceSetMetadata(adp, 'data_type', data_type)
  adp <- oceSetMetadata(adp, 'data_subtype', data_subtype)
  adp <- oceSetMetadata(adp, 'data_origin', data_origin)
  adp <- oceSetMetadata(adp, 'coord_system', coord_system)
  adp <- oceSetMetadata(adp, 'water_mass', water_mass)
  adp <- oceSetMetadata(adp, 'pos_const', pos_const)
  adp <- oceSetMetadata(adp, 'depth_const', depth_const)
  adp <- oceSetMetadata(adp, 'drifter', drifter)
  adp <- oceSetMetadata(adp, 'FillValue', FillValue)
  adp <- oceSetMetadata(adp, 'experiment', experiment)
  adp <- oceSetMetadata(adp, 'project', project)
  adp <- oceSetMetadata(adp, 'description', description)
  adp <- oceSetMetadata(adp, 'longitude', longitude)
  adp <- oceSetMetadata(adp, 'latitude', latitude)
  adp <- oceSetMetadata(adp, 'data_comment', data_comment)
  adp <- oceSetMetadata(adp, 'fill_flag', fill_flag)
  adp <- oceSetMetadata(adp, 'composite', composite)
  adp <- oceSetMetadata(adp, 'magnetic_variation', magnetic_variation)
  adp <- oceSetMetadata(adp, 'platform', platform)
  adp <- oceSetMetadata(adp, 'sounding', sounding)
  adp <- oceSetMetadata(adp, 'chief_scientist', chief_scientist)



}
