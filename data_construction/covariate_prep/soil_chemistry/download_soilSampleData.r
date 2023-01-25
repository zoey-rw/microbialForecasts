# Download all soils data using function from Clara Qin's microbial R package
library(geoNEON)
library(dplyr)
# Source relevant scripts
# source("https://raw.githubusercontent.com/claraqin/neonMicrobe/master/R/get_neon_data.R")
# dat_soil <- downloadRawSoilData(startYrMo = "2013-06", endYrMo = "2020-09",
#                     outDir = "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/")
dat_soil <- data.table::fread("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/sls_soilData_2021-06-16.csv")
dat_soil <- data.table::fread("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/sls_soilData_2022-11-07.csv")

# # Bug in geoNEON package, no longer calculating lat/long

# plot.all <- geoNEON::getLocByName(dat_soil, locCol="namedLocation", locOnly=T, token="eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJ6cndlcmJpbkBidS5lZHUiLCJzY29wZSI6InJhdGU6cHVibGljIiwiaXNzIjoiaHR0cHM6Ly9kYXRhLm5lb25zY2llbmNlLm9yZy8iLCJleHAiOjE3NTc4ODg2NzAsImlhdCI6MTYwMDIwODY3MCwiZW1haWwiOiJ6cndlcmJpbkBidS5lZHUifQ.8eW8vxUOiton-kQ_Xyvva0QSHD_BDd2E5IGeNKW3WHib-m7UpTnEhGFAUlAHGdsUyz-dKE1jMOAGS5A_NRYXGg")
# # Remove duplicates and select output columns
# out <- dat_soil[!duplicated(dat_soil),]
# out <- out %>% dplyr::select(domainID, siteID, plotID, namedLocation, plotType, nlcdClass, adjDecimalLatitude, adjDecimalLongitude,
# 														 adjElevation, collectDate, sampleTiming, standingWaterDepth, sampleID, horizon, soilTemp, soilMoisture,
# 														 litterDepth, sampleTopDepth, sampleBottomDepth, geneticSampleID, nitrogenPercent, organicCPercent, CNratio, CNvals_merged, soilInWaterpH, soilInCaClpH
# )

# Remove duplicates and select output columns
out <- dat_soil[!duplicated(dat_soil),]
out <- out %>% dplyr::select(domainID, siteID, plotID, namedLocation, plotType, nlcdClass, decimalLatitude, decimalLongitude,
                      elevation, collectDate, sampleTiming, standingWaterDepth, sampleID, horizon, soilTemp, soilMoisture,
                      litterDepth, sampleTopDepth, sampleBottomDepth, geneticSampleID, nitrogenPercent, organicCPercent, soilInWaterpH, soilInCaClpH
)
out$dateID <- substr(out$collectDate,1,7)
out$year <- substr(out$collectDate,1,4)

saveRDS(out, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/soilSample_data_allsites.rds")
