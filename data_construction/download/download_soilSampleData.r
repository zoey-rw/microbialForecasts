# Download all soils data using function from Clara Qin's microbial R package
library(geoNEON)
library(dplyr)
# Source relevant scripts
# source("/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/R/utils.R")
# dat_soil <- downloadRawSoilData(startYrMo = "2013-06", endYrMo = "2020-09", 
#                     outDir = "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/")
dat_soil <- data.table::fread("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/raw//sls_soilData_2020-09-23.csv")


# remove analytical replicates
dat_soil <- dat_soil[which(!dat_soil$analyticalRepNumber %in% c(2:4)),]

# Fix unmerged %C and %N columns in soil moisture data
dup <- dat_soil[duplicated(dat_soil$cnSampleID) & !is.na(dat_soil$cnSampleID),]
dup <- dat_soil[which(dat_soil$cnSampleID %in% dup$cnSampleID),] # all duplicated rows, including first appearances
c.only.ind <- which(is.na(dup$nitrogenPercent) & !is.na(dup$organicCPercent))
n.only.ind <- which(!is.na(dup$nitrogenPercent) & is.na(dup$organicCPercent))
c_n_merged <- full_join(dplyr::select(dup[c.only.ind,], -nitrogenPercent),
                        dplyr::select(dup[n.only.ind,], sampleID, nitrogenPercent))
c_n_merged$CNratio <- round(c_n_merged$organicCPercent/c_n_merged$nitrogenPercent,1)
c_n_merged$CNvals_merged <- TRUE # so we know for later that these values are different

# merge back in with rest of C/N data
dat_soil$CNvals_merged <- FALSE
dat_soil <- dat_soil[!(dat_soil$sampleID %in% c_n_merged$sampleID),]
dat_soil <- rbind(dat_soil, c_n_merged)

# Add geolocation data to soilCoreCollection
dat_soil <- getLocTOS(dat_soil, "sls_soilCoreCollection")

# Remove duplicates and select output columns 
out <- dat_soil[!duplicated(dat_soil),]
out <- out %>% dplyr::select(domainID, siteID, plotID, namedLocation, plotType, nlcdClass, adjDecimalLatitude, adjDecimalLongitude, 
                      adjElevation, collectDate, sampleTiming, standingWaterDepth, sampleID, horizon, soilTemp, soilMoisture,
                      litterDepth, sampleTopDepth, sampleBottomDepth, geneticSampleID, nitrogenPercent, organicCPercent, CNratio, CNvals_merged, soilInWaterpH, soilInCaClpH
)
out$dateID <- substr(out$collectDate,1,7)
out$year <- substr(out$collectDate,1,4)

saveRDS(out, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/soilSample_data_allsites.rds")
