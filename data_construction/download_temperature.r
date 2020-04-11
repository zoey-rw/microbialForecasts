# download NEON temperature and combine it with the pre-downloaded ERA5 estimates.

rm(list=ls())
library(neonUtilities)
library(zoo)
library(rjags)
library(coda)

#### 1. Download NEON data 

# create output dataframe
NEON_temp_all.sites <- setNames(data.frame(matrix(ncol = 5, nrow = 1)), 
                                c("siteID", "month", "monthly_temp", "monthly_temp_sd", "source"))
sites <- c("DSNY", "HARV", "OSBS", "CPER", "STER")

# download NEON data for all 5 fies  
NEON_temp_expanded <- loadByProduct(dpID="DP1.00002.001", site=sites, 
                                    package="expanded", check.size=F, avg=30)
NEON_temp_30m <- NEON_temp_expanded$SAAT_30min
# save
saveRDS(NEON_temp_30m, "NEON_temp_data_raw.RDS")


# clean up data
NEON_temp_30m$startDateTime <- as.POSIXct(NEON_temp_30m$startDateTime, 
                                          format="%Y-%m-%d T %H:%M:%S Z", tz="GMT")
NEON_temp_30m$month <- substr(NEON_temp_30m$startDateTime, 1, 7)
NEON_temp_30m <- NEON_temp_30m[NEON_temp_30m$finalQF!=1,]
NEON_temp_30m <- NEON_temp_30m[!is.na(NEON_temp_30m$tempSingleMean) & 
                                 !is.na(NEON_temp_30m$tempSingleExpUncert) & 
                                 !is.na(NEON_temp_30m$tempSingleStdErMean),]

# Read in ERA5 ensembles of historical data
clim_ensembles <- readRDS("climate_ensembles_13-19.rds")

# save as model input data
save(clim_ensembles, NEON_temp_30m, file = "temp_time_series_input.Rdata")





output.path <- "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/dp.00002.SAAT_30min.rds"

# Read in data
temp30_sites <- read.csv('NEON-pheno-temp-timeseries/temp/SAAT_30min.csv', stringsAsFactors = FALSE)

# create new dataframe without NAs
temp30_noNA <- temp30 %>%
  drop_na(tempSingleMean)  # tidyr function

sum(is.na(temp30_noNA$tempSingleMean))


temp30_noNA$startDateTime <- as.POSIXct(temp30_noNA$startDateTime,
                                        format = "%Y-%m-%dT%H:%M:%SZ", tz = "GMT")
# check that conversion worked
str(temp30_noNA$startDateTime)


res <- googleway::google_timezone(c(lat, long), time1, key = NULL)




# create output dataframe
NEON_temp_all.sites <- setNames(data.frame(matrix(ncol = 5, nrow = 1)), 
                                c("siteID", "month", "monthly_temp", "monthly_temp_sd", "source"))
sites <- c("DSNY", "HARV", "OSBS", "CPER", "STER")

site <- "DSNY"
NEON_temp_basic <- loadByProduct(dpID="DP1.00002.001", site=site, nCores=3, startdate = "2017-01",enddate="2017-06",
                                    package="basic", check.size=F, avg=30)
NEON_temp_30m <- NEON_temp_expanded$SAAT_30min
# download NEON data for all 5 fies  
NEON_temp_expanded <- loadByProduct(dpID="DP1.00002.001", site=site, 
                                    package="expanded", check.size=F, avg=30)
NEON_temp_30m <- NEON_temp_expanded$SAAT_30min
saveRDS(NEON_temp_30m, "NEON_temp_data_raw.RDS")
# clean up data
NEON_temp_30m$startDateTime <- as.POSIXct(NEON_temp_30m$startDateTime, 
                                          format="%Y-%m-%d T %H:%M:%S Z", tz="GMT")
NEON_temp_30m$month <- substr(NEON_temp_30m$startDateTime, 1, 7)
NEON_temp_30m <- NEON_temp_30m[NEON_temp_30m$finalQF!=1,]
NEON_temp_30m <- NEON_temp_30m[!is.na(NEON_temp_30m$tempSingleMean) & 
                                 !is.na(NEON_temp_30m$tempSingleExpUncert) & 
                                 !is.na(NEON_temp_30m$tempSingleStdErMean),]

# save data
save(NEON_temp_30m, output.path)

