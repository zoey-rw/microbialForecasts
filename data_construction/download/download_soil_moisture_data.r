library(neonUtilities)
outdir <- "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEONSoilMoist_raw"
outpath <- "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEONSoilMoist_raw/rawMoistStacked"
# Download soil temperature data
zipsByProduct(dpID = "DP1.00094.001", site = "all", 
              # enddate = "2015-12", 
              avg = 30, check.size = F, 
              savepath = outdir)
out <- stackByTable(filepath = paste0(outdir, "/filesToStack00094"), savepath = outpath, folder = TRUE, 
                    nCores = 16, saveUnzippedFiles = FALSE)


NEON_moist_allsites <- loadByProduct(dpID = "DP1.00094.001", site = c("CPER","OSBS","HARV","DSNY","STER"), enddate="2018-12", check.size = F, avg = 30)
NEON_moist_30m <- NEON_moist_allsites$SWS_30_minute
NEON_moist_30m$sensorID <- paste0(NEON_moist_30m$siteID, "_", NEON_moist_30m$horizontalPosition, ".", NEON_moist_30m$verticalPosition)
NEON_moist_30m$month <- substr(NEON_moist_30m$startDateTime, 1, 7)
NEON_moist_30m$depth <- substr(NEON_moist_30m$verticalPosition, 3, 3)
NEON_moist_30m$sensor <- substr(NEON_moist_30m$horizontalPosition, 3, 3)
NEON_moist_30m$day <- as.Date(format(as.POSIXct(NEON_moist_30m$startDateTime,format='%Y-%m/%d %H:%M:%S'),format="%Y-%m-%d"))
NEON_moist_30m$date_time <-  substr(NEON_moist_30m$startDateTime, 1, 13)
NEON_moist_30m$date_time <- gsub(" ", "T", NEON_moist_30m$date_time)
saveRDS(NEON_moist_30m, "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEONSoilMoist_raw/allRawMoist_5sites.rds")


# Remove extra rows in sensor data
NEON_moist_allsites <- loadByProduct(dpID = "DP1.00094.001", site ="all", startdate = "2018-12", enddate="2018-12", check.size = F, avg = 30)
sensors <- NEON_moist_allsites$sensor_positions_00094[which(!is.na(soil.raw$sensor_positions_00094$referenceLatitude)),]
sensors$sensorID <- paste(sensors$siteID, sensors$HOR.VER, sep ="_")
sensors$sensorDepth <- abs(sensors$zOffset*100)
saveRDS(sensors, "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEONSoilMoist_raw/sensorLocations.rds")