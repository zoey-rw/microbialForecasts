library(tidyr)
library(dplyr)

# Read in and prep the temperature data
temp_raw <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/temp_time_series_output.rds")

# subset to second half of 2013, for testing
temp_raw <- temp_raw[grepl("2013", temp_raw$month),]
temp_raw <- temp_raw[which(as.numeric(substr(temp_raw$month,6,7)) > 5),]
temp <- temp_raw %>% mutate(dateID = stringr::str_replace_all(month, "-", "")) %>% 
  #mutate(siteID = as.numeric(factor(siteID))) %>% 
  dplyr::select(-c(source,month))
temp <- temp[,c(1,4,2,3)]

precip_raw <- readRDS("/usr3/graduate/zrwerbin/temporal_forecasts/data/precip_time_series_output.rds")
precip_raw <- precip_raw[grepl("2013", precip_raw$month),]
precip_raw <- precip_raw[which(as.numeric(substr(precip_raw$month,6,7)) > 5),]

precip <- precip_raw %>% mutate(dateID = stringr::str_replace_all(month, "-", "")) %>% 
  #mutate(siteID = as.numeric(factor(siteID))) %>% 
  dplyr::select(-c(source,month))
precip <- precip[,c(1,4,2,3)]

weather <- merge(temp, precip)
weather <- weather %>% dplyr::select(-c(monthly_temp_sd, monthly_precip_sd))

phys_raw <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/dp.10086.soil_phys.rds")
plots2013 <- as.character(unique(phys_raw[grepl("2013", phys_raw$dateID),]$plotID))
plot_pH <-  phys_raw %>% 
  group_by(siteID, plotID) %>% summarize(pH = mean(soilInCaClpH, na.rm=T))
plot_pH <- plot_pH[which(plot_pH$plotID %in% plots2013),]

plot_pH[plot_pH$siteID == "CPER" & is.na(plot_pH$pH),]$pH <- 6.37  
plot_pH[plot_pH$siteID == "STER" & is.na(plot_pH$pH),]$pH <- 6.39  
plot_pH[plot_pH$siteID == "DSNY" & is.na(plot_pH$pH),]$pH <- 3.50  
plot_pH[plot_pH$siteID == "OSBS" & is.na(plot_pH$pH),]$pH <- 3.95
plot_pH[plot_pH$siteID == "HARV" & is.na(plot_pH$pH),]$pH <- 3.46

saveRDS(list(weather, plot_pH), "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/EE509_model_covariates.rds")
