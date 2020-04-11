## script to get weather covariates for a given sample

# download NEON soil physical properties (dp.10086.001), build site/plot time-series
# clear environment, load packages and set file paths
rm(list=ls())
library(zoo)
library(dplyr)
library(neonUtilities)
library(ggplot2)
library(scales)
library(geoNEON)
source('/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/merge_left.r')

output.path <- "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/dp.10086.soil_phys.rds"

# download data from each site
all_dat <- loadByProduct(dpID="DP1.10086.001", site=c("HARV","OSBS","CPER","DSNY","STER"), package="basic", check.size = F)

# read in tables
soil_pH <- all_dat$sls_soilpH 
soil_moisture <- all_dat$sls_soilMoisture
soil_core <- all_dat$sls_soilCoreCollection

# merge downloaded files
soil_pH_core.merge <- merge_left(soil_core, soil_pH, all.x = T)
soil_phys_merge <- merge_left(soil_pH_core.merge, soil_moisture, all.x = T)

# create new date columns, as both a Date class object and a dateID (YYYY-MM)
soil_phys_merge$asDate <- as.Date(as.yearmon(soil_phys_merge$collectDate))
soil_phys_merge$dateID <- substr(soil_phys_merge$asDate, 1, 7)

# Add precise geolocation data
soil_phys_merge <- geoNEON::def.calc.geo.os(soil_phys_merge, 'sls_soilCoreCollection')
# pull out columns of interest
soil_phys <- soil_phys_merge[,c("sampleID","siteID", "plotID","dateID","asDate", "soilInCaClpH", "litterDepth", "soilTemp", "soilMoisture","soilTemp","sampleBottomDepth", "sampleTopDepth","standingWaterDepth", "adjDecimalLongitude", "adjDecimalLatitude","adjElevation","sampleTiming","geneticSampleID","horizon","nlcdClass","soilInWaterpH")]
soil_phys <- soil_phys[order(soil_phys$dateID, soil_phys$plotID),]

saveRDS(soil_phys, output.path)
