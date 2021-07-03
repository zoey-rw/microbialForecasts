
# 16 million rows. oy vey.
NEON_temp_30m <- read.csv("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEONSoilTemp_daily/dailyTempStacked.rds/stackedFiles/ST_30_minute.csv")
NEON_temp_30m <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEONSoilTemp_raw_allsites.rds")
NEON_temp_30m$sensorID <- paste0(NEON_temp_30m$siteID, "_", NEON_temp_30m$horizontalPosition, ".", NEON_temp_30m$verticalPosition)

NEON_temp_30m$month <- substr(NEON_temp_30m$startDateTime, 1, 7)
NEON_temp_30m$depth <- substr(NEON_temp_30m$verticalPosition, 3, 3)
NEON_temp_30m$sensor <- substr(NEON_temp_30m$horizontalPosition, 3, 3)
#NEON_temp_30m$day <- as.Date(format(as.POSIXct(NEON_temp_30m$startDateTime,format='%Y-%m/%dT%H:%M:%S'),format="%Y-%m-%d"))
NEON_temp_30m$day <- as.Date(substr(as.character(NEON_temp_30m$startDateTime), 1, 10))
NEON_temp_30m$date_time <-  substr(NEON_temp_30m$startDateTime, 1, 13)
NEON_temp_30m$date_time <- gsub(" ", "T", NEON_temp_30m$date_time)

# Decide which sensor/depth is most correlated with observations
dat_soil <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/soilData.rds")
dat_soil$day <- as.Date(dat_soil$collectDate)
dat_soil$date_time <- substr(dat_soil$collectDate, 1, 13)
dat_soil <- dat_soil[which(!is.na(dat_soil$soilTemp)),]

sens.subset <- NEON_temp_30m[NEON_temp_30m$day %in% dat_soil$day,]
sens.subset <- sens.subset[sens.subset$date_time %in% dat_soil$date_time,]

merged <- merge(sens.subset, dat_soil, by = c("date_time","siteID"), all.x=T)
merged <- merged[which(!is.na(merged$soilTemp)),]

ggplot(merged) + geom_point(aes(y = soilTempMean, x = soilTemp, color = depth)) + facet_grid(rows = vars(sensor), cols=vars(plotID)) + geom_abline()

ggplot(merged) + geom_point(aes(y = soilTempMean, x = soilTemp, color = as.factor(depth))) + facet_grid(rows = vars(sensor)) + geom_abline()

ggplot(merged) + geom_point(aes(y = soilTempMean, x = soilTemp))+ geom_abline()

ggplot(merged[merged$depth %in%  c(1:5),]) + geom_point(aes(y = soilTempMean, x = soilTemp, color = depth)) + facet_grid(rows = vars(sensor), cols=vars(plotID)) + geom_abline()

# Use mean of all sensors for the site
sensor_by_site <- merged[merged$depth %in% c(1:5),] %>% group_by(siteID, depth, day.y) %>% mutate(by_site = mean(soilTempMean))
ggplot(sensor_by_site) + geom_point(aes(y = by_site, x = soilTemp, color = depth)) + facet_grid(~siteID) + geom_abline()


sensor_by_site <- merged[merged$depth %in% c(1:5),] %>% group_by(siteID, date_time) %>% mutate(by_site = mean(soilTempMean))
ggplot(sensor_by_site) + geom_point(aes(y = by_site, x = soilTemp)) + facet_grid(~siteID) + geom_abline()

sensor_by_site <- merged[merged$depth %in% c(3),] %>% group_by(siteID, date_time) %>% mutate(by_site = mean(soilTempMean))
ggplot(sensor_by_site) + geom_point(aes(y = by_site, x = soilTemp)) + #facet_grid(~siteID) + 
	geom_abline()
fit <- lm(by_site ~ soilTemp, sensor_by_site)
cor_plot <- ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
	geom_point() +
	stat_smooth(method = "lm", col = "red") +
	labs(title = "Calibration fits", subtitle = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 4))) + 
	xlab("Microbial soil samples") + ylab("Mean of soil sensor array") + 
	theme_minimal() + theme (text = element_text(size = 12)) 
cor_plot

# Use mean for specific sensors
temp_by_sensor <- merged[merged$depth %in% c(1:5),] %>% group_by(siteID, horizontalPosition, day.y) %>% mutate(by_sensor = mean(soilTempMean))
ggplot(temp_by_sensor) + geom_point(aes(y = by_sensor, x = soilTemp, color = as.factor(horizontalPosition))) + geom_abline() + facet_grid(~siteID)

shallow <- sensor_by_site[sensor_by_site$depth %in% c(1)]

# 
# 
# # Remove extra rows in sensor data
# sensors <- NEON_temp_allsites$sensor_positions_00041[which(!is.na(soil.raw$sensor_positions_00041$referenceLatitude)),]
# sensors$sensorID <- paste(sensors$siteID, sensors$HOR.VER, sep ="_")
# sensors$sensorDepth <- abs(sensors$zOffset*100)
# saveRDS(sensors, "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEONSoilTemp_daily/sensorLocations.rds")
# 
# 



# Trying to get sensor location data from API

library(httr)
library(jsonlite)
library(dplyr, quietly=T)
library(downloader)


site <- "CPER"
# Request all location data from a specific site
all.locs.site <- GET(paste0("http://data.neonscience.org/api/v0/locations/", site))
all.locs.text <- fromJSON(content(all.locs.site, as="text"), simplifyDataFrame=T, flatten=T)
site.urls <- all.locs.text$data$locationChildrenUrls

# get data availability from location/date of interest
soil.ar <- GET(site.urls[grep("SOILAR", site.urls)])
soil.ar.files <- fromJSON(content(soil.ar, as="text"))
sens.urls <- soil.ar.files$data$locationChildrenUrls

# Loop through for each sensor
for (s in 1:5){
  # Specific soil sensor
  sens <- GET(sens.urls[[s]])
  sens.files <- fromJSON(content(sens, as="text"), simplifyDataFrame=T, flatten=T)
  temp.urls <- sens.files$data$locationChildrenUrls[grep("SOILTP", sens.files$data$locationChildrenUrls)]
  
  # Loop through for each depth
  for (d in 1:9){
    # Specific vertical soil sensor depth
    depth <- GET(temp.urls[[d]])
    depth.files <- fromJSON(content(depth, as="text"))
    depth.url <- depth.files$data$locationChildrenUrls
    
    # Finally, lat/lon for that sensor
    loc <- GET(depth.url)
    loc.file <- fromJSON(content(loc, as="text"))
    
    # Specifically, soil temperature sensor
    tmp2 <- GET(urls2[[2]])
    tmp.files2 <- fromJSON(content(tmp2, as="text"))
    urls3 <- unlist(tmp.files2$data$locationChildrenUrls)
    
    # Location data for that one sensor
    tmp3 <- GET(urls3[[1]])
    tmp.files3 <- fromJSON(content(tmp3, as="text"))
  }
  
}

