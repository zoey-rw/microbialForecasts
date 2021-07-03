# Workflow to compare NEON soil sensor data on volumetric water content with sampled moisture data
# Used to decide which soil sensors most accurately reflect the conditions experienced by soil microbes
# Decided on 2nd sensor depth, averaged by NEON site

library(neonUtilities)
library(dplyr)
library(fuzzyjoin)
library(cutpointr)
library(tidyverse)
library(gridExtra)

source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/functions/helperFunctions.r")

# Code outline: assumes you have already downloaded raw sensor data and sample data
## 1. Download bulk density/particle size data ##
## 2. Calculate bulk density for each plot/horizon ##
## 3. Convert gravimetric to volumetric soil moisture ##
## 4. Decide which horizon each sensor belongs to ##
## 5. Match up sensor data with sample data, by hour ##
## 6. Visualize data relationships ##


##### 1. Download bulk density/particle size data #####

# # Get initial characterization for each site. This includes bulk density for many (but not all) plots.
# init_phys <- loadByProduct(dpID = "DP1.10047.001", site = "all", check.size = F)
# ps <- init_phys$spc_particlesize
# ps$coarseFragPercent <- (ps$coarseFrag2To5 + ps$coarseFrag5To20) * .001
# bd <- init_phys$spc_bulkdensity %>% select(-c(domainID, nrcsDescriptionID, bulkDensProcessedDate, bulkDensMethodPub, bulkDensMethod, bulkDensIDnrcs, bulkDensSampleType, bulkDensCenterDepth, laboratoryName, remarks, publicationDate, namedLocation, uid))

# # Combine particle size and bulk density tables
# bd_ps <- merge(bd, ps[,c("horizonID","coarseFragPercent")], by = "horizonID", all.x = T)
# bd_ps$fineFragPercent <- 1 - bd_ps$coarseFragPercent
# saveRDS(bd_ps, "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/bulkDensity_allsites.rds")

##### 2. Calculate bulk density for each plot/horizon #####

bd_ps <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/bulkDensity_allsites.rds")

# label horizons by O or M
bd_ps$horizon <- ifelse(grepl("O", bd_ps$horizonName), "O", "M")

#sites defined elsewhere (TO DO: delete later)
bd_ps <- bd_ps[bd_ps$siteID %in% sites,]

# Subset to shallower bulk density measurements
bd_shallow <- bd_ps[bd_ps$bulkDensTopDepth < 25 & bd_ps$bulkDensBottomDepth < 35,] 

# Use third bar measurements whenever possible (user guide recommends that)
bd_shallow$bulkDens <- ifelse(!is.na(bd_shallow$bulkDensThirdBar), 
                              bd_shallow$bulkDensThirdBar, bd_shallow$bulkDensFieldMoist)

# Multiply by <2mm fraction (for the same horizon) to get more accurate value
bd_shallow$bulkDensFineFrag <- bd_shallow$bulkDens * bd_shallow$fineFragPercent

# Summarize by site or plot
site_plot_bd <- bd_shallow %>% group_by(siteID, horizon) %>%  #depth) %>% 
  mutate(mean_bd_FineFrag_site = mean(bulkDensFineFrag, na.rm=T), 
         mean_topDepth_site = mean(bulkDensTopDepth, na.rm=T),  
         mean_bottomDepth_site = mean(bulkDensBottomDepth, na.rm=T)) %>% 
  ungroup() %>% 
  group_by(plotID, horizon) %>%  #depth) %>% 
  mutate(mean_bd_FineFrag_plot = mean(bulkDensFineFrag, na.rm=T)) %>% 
  distinct(plotID, horizon, .keep_all=TRUE) %>% as.data.frame()


##### 3. Convert gravimetric to volumetric soil moisture #####

# Read in soil samples
dat_soil <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/soilData.rds")
dat_soil$day <- as.Date(dat_soil$collectDate)
dat_soil$date_time <- substr(dat_soil$collectDate, 1, 13)

# Merge with plot bulk density valies
scc_bd <- merge(dat_soil, site_plot_bd, all.x  = T, by = c("siteID", "plotID","horizon"))
# Use site-level bulk density.
scc_bd$vol_SoilMoisture_site_FineFrag <- scc_bd$soilMoisture * scc_bd$mean_bd_FineFrag_site
# Use plot-level bulk density.
scc_bd$vol_SoilMoisture_FineFrag <- scc_bd$soilMoisture * scc_bd$mean_bd_FineFrag_plot

# here someone else does what I started doing to geolocate sensors:
# https://github.com/cran/Z10/blob/41948458ae0094b1d0810dd68de9fa096d17c668/R/daily_soil_temp_mean.R

##### 4. Decide which horizon each sensor belongs to #####

# craete match-up table for soil horizon depths
all_horizon_depths <- bd_ps %>% dplyr::select(c(siteID,plotID,bulkDensTopDepth,bulkDensBottomDepth,horizon))

# # View variety of fine-scale horizon classifications. To justify taking a median cutpoint for the site.
# ggplot(bd_ps[bd_ps$siteID=="OSBS",]) + geom_errorbar(mapping=aes(x=horizon, ymin=bulkDensTopDepth, ymax=bulkDensBottomDepth, color=horizonName), width=0.2, size=3) + facet_wrap(~plotID, nrow=3)

# Identify sites that even have organic & mineral layers
has_organic <- unique(all_horizon_depths[all_horizon_depths$horizon == "O",]$siteID)
has_mineral <- unique(all_horizon_depths[all_horizon_depths$horizon == "M",]$siteID)
has_both <- intersect(has_organic, has_mineral)

# Keep only the shallowest mineral layer, which we'll use to determine a good cutpoint between organic and mineral
horizon_depths <- all_horizon_depths[all_horizon_depths$bulkDensTopDepth < 30,] %>% 
  group_by(plotID, horizon) %>% 
  dplyr::slice(which.min(bulkDensTopDepth)) %>% dplyr::filter(siteID %in% has_both)

# View variety of organic vs. mineral horizon classifications. To justify taking a median cutpoint for the site.
# ggplot(horizon_depths[horizon_depths$siteID=="HARV",]) + geom_errorbar(mapping=aes(x=horizon, ymin=bulkDensTopDepth, ymax=bulkDensBottomDepth, color=horizon), width=0.2, size=3) + facet_wrap(~plotID, nrow=2)

# Expand data frame to long form (TO DO: don't understand the error but seems to work)
expanded <- horizon_depths %>%
  nest(bulkDensTopDepth, bulkDensBottomDepth) %>%
  mutate(data = map(data, ~seq(unique(.x$bulkDensTopDepth), unique(.x$bulkDensBottomDepth), 1))) %>%
  unnest(data) %>% dplyr::rename(depth = "data") %>% ungroup()

# Decide organic-mineral transition depth point for each site
cutpoints <- cutpointr(expanded, depth, horizon, siteID,
                       pos_class = "M", direction = ">=",
                       method = maximize_boot_metric,
                       metric = accuracy)
key1 <- cutpoints %>% dplyr::rename(siteID = "subgroup", bulkDensBottomDepth = "optimal_cutpoint") %>% mutate(bulkDensTopDepth = 0, horizon = "O") %>% dplyr::select(siteID,bulkDensTopDepth, bulkDensBottomDepth, horizon)
key2 <- cutpoints %>% dplyr::rename(siteID = "subgroup", bulkDensTopDepth = "optimal_cutpoint") %>% mutate(bulkDensBottomDepth = 30, horizon = "M") %>% dplyr::select(siteID,bulkDensTopDepth, bulkDensBottomDepth, horizon)

# Keep only the shallowest mineral layer, which we'll use to determine a good cutpoint between organic and mineral
horizon_depth_key <- all_horizon_depths[all_horizon_depths$bulkDensTopDepth < 30,] %>% 
  group_by(siteID, horizon) %>% 
	dplyr::slice(which.min(bulkDensTopDepth)) %>% dplyr::filter(!siteID %in% has_both) %>% mutate(bulkDensBottomDepth = 30)
horizon_depth_key <- bind_rows(key1, key2, horizon_depth_key)

# Read in soil moisture sensor locations
sensors <- read.csv("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEONSoilMoist_daily/rawMoistStacked/stackedFiles/sensor_positions_00094.csv")
sensors$sensorID <- paste(sensors$siteID, sensors$HOR.VER, sep ="_")
sensors$sensorDepth <- abs(sensors$zOffset*100)

# According to NEON's ATBD, sensors integrate 5cm and 5cm below, so we'll cut sensors off at 25cm.
sensors <- sensors %>% dplyr::select(name, referenceName, siteID, sensorID, sensorDepth) %>% filter(sensorDepth < 25)
# Match to horizons
sensors.to.keep <- fuzzy_left_join(
  sensors, horizon_depth_key,
  by = c(
    "siteID" = "siteID",
    "sensorDepth" = "bulkDensTopDepth",
    "sensorDepth" = "bulkDensBottomDepth"
  ),
  match_fun = list(`==`, `>=`, `<=`)
) 

##### 5. Match up sensor data with sample data, by hour #####

# Read in soil moisture, and subset to shallow sensors
daily.SWC <- data.table::fread("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEONSoilMoist_daily/rawMoistStacked.rds/stackedFiles/SWS_30_minute.csv")
daily.SWC <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEONSoilMoist_raw_allsites.rds")

daily.SWC$sensorID <- paste0(daily.SWC$siteID, "_", daily.SWC$horizontalPosition, ".", daily.SWC$verticalPosition) # can delete later
daily.SWC$date_time <-  substr(daily.SWC$startDateTime, 1, 13)
daily.SWC$date_time <- gsub(" ", "T", daily.SWC$date_time) # get date-time to the hour
daily.SWC$depth <- substr(daily.SWC$verticalPosition, 3, 3)

#Merge w/ sensor depth data
#daily.SWC.ID <- daily.SWC
daily.SWC.ID <- merge(daily.SWC, sensors.to.keep[,c("sensorID", "sensorDepth","horizon")], by = "sensorID", all=F)

# Prep sensor data for merging
daily.SWC.merge <- daily.SWC.ID %>%  filter(date_time %in% scc_bd$date_time) %>% select(sensorID,siteID,horizontalPosition,VSWCMean,VSWCFinalQF,date_time,sensorDepth,horizon,depth)
#daily.SWC.merge <- daily.SWC.ID %>%  filter(date_time %in% scc_bd$date_time) %>% select(sensorID,siteID,horizontalPosition,VSWCMean,VSWCFinalQF,date_time,depth)

# Remove data that didn't pass QF <- removes too much (for now)
#daily.SWC.merge <- daily.SWC.merge %>% filter(VSWCFinalQF == 0)

# Merge sensor and sample data.
merged <- merge(scc_bd, daily.SWC.merge, by = c("date_time","siteID","horizon"))
merged <- merged[which(!is.na(merged$vol_SoilMoisture_site_FineFrag) & !is.na(merged$VSWCMean)),]


##### 6. Visualize data relationships #####

# Compare all sensor measurements against all sampled data
ggplot(merged) + geom_point(aes(y = vol_SoilMoisture_FineFrag, x = VSWCMean, col = as.factor(horizontalPosition))) + geom_abline()

# Color by sensor depth
ggplot(merged) + geom_point(aes(y = VSWCMean, x = vol_SoilMoisture_FineFrag, color = sensorDepth)) + facet_grid(rows = vars(horizontalPosition), cols=vars(siteID), scales="free") + geom_abline()

# Only shallowest 2 sensors
ggplot(merged[merged$depth <= 2,]) + geom_point(aes(y = VSWCMean, x = vol_SoilMoisture_FineFrag, color = horizon)) + facet_grid(rows = vars(horizontalPosition), cols=vars(siteID), scales="free") + geom_abline()

# Only keep mineral horizons
ggplot(merged[merged$depth <=2 & merged$horizon == "M",]) + geom_point(aes(y = VSWCMean, x = vol_SoilMoisture_FineFrag, color = as.factor(depth))) + facet_grid(rows = vars(horizontalPosition), cols=vars(siteID), scales="free") + geom_abline()

# Get mean of all depths per sensors 
mois_by_sensor <- merged[merged$horizon == "M",] %>% group_by(siteID, horizontalPosition, date_time) %>% mutate(by_sensor = mean(VSWCMean))
ggplot(mois_by_sensor) + geom_point(aes(y = by_sensor, x = vol_SoilMoisture_FineFrag, color = as.factor(horizontalPosition))) + geom_abline() + facet_grid(~siteID)

# Mean of sensors by depth
mois_by_depth <- merged %>% group_by(siteID, depth, date_time) %>% mutate(by_depth = mean(VSWCMean))
ggplot(mois_by_depth) + geom_point(aes(y = by_depth, x = vol_SoilMoisture_FineFrag, color = as.factor(depth))) + geom_abline() + facet_wrap(~siteID, nrow=2)

# Mean of 2nd shallowest sensor layer per site
mois_by_site <- merged[merged$depth == 2,] %>% group_by(siteID, date_time) %>% mutate(by_sensor = mean(VSWCMean))
ggplot(mois_by_site) + geom_point(aes(y = by_sensor, x = vol_SoilMoisture_FineFrag)) + geom_abline() + facet_wrap(~siteID, nrow=2)