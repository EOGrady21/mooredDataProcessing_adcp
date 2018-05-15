###example processing procedure

library(oce)
library(ncdf4)

source('adcpToolbox.R')


adp <- read.adp.easy('M1996000.000', 'metadataTemplate.csv')
adp[['depth']] <- swDepth(adp[['pressure']], (adp[['latitude']]))
adp <- applyMagneticDeclinationAdp(adp)
adp <- limit_depthbytime(adp)
adp <- limit_time(adp)
adp <- adpFlag(adp, 25, 0.46)
adp <- handleFlags(adp, flags = 4, actions = list('NA'))
save(adp, file =paste( 'adp', adp[['cruiseNumber']], adp[['mooringNumber']], sep = '_'))
oceNc_create(adp, 'testerb',  'metadataTemplate.csv')
