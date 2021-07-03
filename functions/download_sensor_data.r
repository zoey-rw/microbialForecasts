##' @name download_NEON_soilmoist
##' @description: 
##' Download NEON Soil Water Content and Soil Salinity data by date and site name
##' 
##' @param site four letter NEON site code name(s). If no site is specified, it will download all of them (chr) (e.g "BART" or c("SRER", "KONA", "BART"))
##' @param avg averaging interval (minutes): 1, 30, or both ("all") . default returns both
##' @param var variable of interest: "SWC" (soil water content) or "SIC" (soil ion content) or both ("all") default returns both.
##'     Both variables will be saved in outdir automatically (chr)
##' @param startdate start date as YYYY-mm. If left empty, all data available will be downloaded (chr)
##' @param enddate start date as YYYY-mm. If left empty, all data available will be downloaded (chr) 
##' @param outdir out directory to store the following data:
##'     .rds list files of SWC and SIC data for each site and sensor position, 
##'     sensor positions .csv for each site, 
##'     variable description .csv file,
##'     readme .csv file
##' @return List of specified variable(s) AND prints the path to output folder
##' 
##' @author Juliette Bateman
##' 
##' @example
##' \run{
##' test <- download_NEON_soilmoisture(
##'   site = c("SRER", "BART", "KONA"),
##'   avg = 30,
##'   var = "SWC",
##'   startdate = "2019-01",
##'   enddate = "2020-01",
##'   outdir = getwd())}

## Install NEON libs
#devtools::install_github("NEONScience/NEON-geolocation/geoNEON")
#devtools::install_github("NEONScience/NEON-utilities/neonUtilities", force = TRUE)
#install.packages("BiocManager")
# BiocManager::install("rhdf5")



site = c("CPER") 
site = c("HARV") 
var = "all"
# startdate = "2019-02"
# enddate = "2019-03" 

startdate = "2018-01"
enddate = "2020-06" 
startdate = NA
enddate = NA 
avg = 30
outdir = "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/"

site <- c("OSBS","DSNY","STER","HARV","CPER")
download_NEON_soilmoist_daily <- function(site, avg = "all", var = "all",
                                    startdate = NA, enddate = NA,
                                    outdir) {
  
  sitenames <- paste(site, collapse="_")
  require(dplyr)
  #################### Data Download from NEON #################### 
  soil.raw = neonUtilities::loadByProduct(dpID = "DP1.00094.001", site = site, avg = avg, startdate = startdate, enddate = enddate, check.size = FALSE)
  
  # Export into new folder in outdir
  dir = paste0(outdir, "/NEONSoilMoist", "_daily")
  dir.create(dir)
  
  #################### Clean-up Data Observations ####################
  # Only select data from list and remove flagged observations 
  if (avg == 30) {
    data.raw = soil.raw$SWS_30_minute %>% stats::na.omit()
  } else if (avg == 1) {
    data.raw = soil.raw$SWS_1_minute %>% stats::na.omit()
  } else {
    data.raw = list(soil.raw$SWS_1_minute, soil.raw$SWS_30_minute) %>% stats::na.omit()
  }
  
  # Reformat day variable
  data.raw$day <- format(as.POSIXct(data.raw$startDateTime,format='%Y-%m/%d %H:%M:%S'),format="%Y-%m-%d")

  data.raw$sensorID <- paste0(data.raw$siteID, "_", data.raw$horizontalPosition, ".", data.raw$verticalPosition)
  
  # Separate variables, omit flagged data obs
  data.raw.SWC = (split(data.raw, data.raw$VSWCFinalQF))$'0' %>%
    dplyr::select(c("domainID", "siteID", "sensorID", "day", "horizontalPosition", "verticalPosition", "startDateTime", "endDateTime", "VSWCMean", "VSWCMinimum", "VSWCMaximum", "VSWCVariance", "VSWCNumPts", "VSWCExpUncert", "VSWCStdErMean"))
  data.raw.SIC = (split(data.raw, data.raw$VSICFinalQF))$'0' %>%
    dplyr::select(c("domainID", "siteID", "sensorID", "day", "horizontalPosition", "verticalPosition", "startDateTime", "endDateTime","VSICMean", "VSICMinimum", "VSICMaximum", "VSICVariance", "VSICNumPts", "VSICExpUncert", "VSICStdErMean"))
  
  # Summarize by day
  daily.SWC <- data.raw.SWC %>% group_by(siteID, day, sensorID) %>% mutate_at(.vars = c("VSWCMean", "VSWCMinimum", "VSWCMaximum", "VSWCNumPts"), .funs = mean) %>% 
    distinct(siteID, sensorID, .keep_all=TRUE)
  daily.SIC <- data.raw.SIC %>% group_by(siteID, day, sensorID) %>% mutate_at(.vars = c("VSICMean", "VSICMinimum", "VSICMaximum", "VSICNumPts"), .funs = mean) %>% 
    distinct(siteID, day, sensorID, .keep_all=TRUE)
  
  # Remove extra rows in sensor data
  sensors <- soil.raw$sensor_positions_00094[which(!is.na(soil.raw$sensor_positions_00094$referenceLatitude)),]
  sensors$sensorID <- paste(sensors$siteID, sensors$HOR.VER, sep ="_")
  sensors$sensorDepth <- abs(sensors$zOffset*100)
  
  #################### Save data into files ####################
  # all sensor data in one file
  utils::write.csv(sensors, file = paste0(dir, "/sensor_positions.csv"))
  utils::write.csv(soil.raw$variables_00094, file = paste0(dir, "/variable_description.csv"))
  saveRDS(daily.SIC, file = paste0(dir, "/", sitenames, startdate, "_", enddate, "_SIC_data.rds"))
  saveRDS(daily.SWC, file = paste0(dir, "/", sitenames, startdate, "_", enddate, "_SWC_data.rds"))
 
  # Return file path to data and print lists of 
  cat("Done! NEON soil moisture data has been downloaded and stored in ", paste0(dir), ".")
  
}



# saveRDS(data.raw.SWC, file = paste0(dir, "/raw_CPER_SWC_data.rds"))
# saveRDS(data.raw.SWC, file = paste0(dir, "/raw_HARV_SWC_data.rds"))
# 
# utils::write.csv(sensors, file = paste0(dir, "/HARV_sensor_positions.csv"))
# saveRDS(daily.SIC, file = paste0(dir, "/", startdate, "_", enddate, "_HARV_SIC_data.rds"))
# saveRDS(daily.SWC, file = paste0(dir, "/", startdate, "_", enddate, "_HARV_SWC_data.rds"))
