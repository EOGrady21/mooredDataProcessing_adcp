###pull meta from netCDF combine with ODF data, output new NetCDF
#reads metadata from netCDF and puts it into adp
#not working yet


#' Copy metadata from netCDF
#'
#' pulls global attributes from netCDF to be put into new format netCDF
#' @param ncin old version netCDF
#' @param ncout new netCDF version

copy.nc.meta <- function(ncin, ncout){
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
  remove(ni)
  #create list of all varibles to be put into new NC

  # var <- as.list(creation_date, mooring, deployment_date, recovery_date, inst_type, history,
  #                starting_water_layer, ending_water_layer, depth_note, transform, data_type,
  #                data_subtype, data_origin, coord_system, water_mass, pos_const, depth_const,
  #                drifter, FillValue, experiment, project, description, longitude, latitude,
  #                data_comment, fill_flag, composite, magnetic_variation, platform, sounding,
  #                chief_scientist )

  #put into outgoing netCDF

  no <- nc_open(ncout)

  for (vr in var){
    if(vr$hasatt == TRUE){
      ncatt_put(no, 0, vr$value)
    }
    else{
      warning( paste('Missing', print(names(var[vr])), '!', sep = "  "))
    }
  }

}
