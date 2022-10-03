# library(neonMicrobe)
#
# soil_Harv <- downloadSoilData(
#   startYrMo = "2014-07",
#   sites = c("HARV"), outDir = "/projectnb2/talbot-lab-data/zrwerbin/interactions/data/neonMicrobe_metadata/"
# )
soil_Harv <- read.csv("/projectnb2/talbot-lab-data/zrwerbin/interactions/data/neonMicrobe_metadata/sls_soilData_2022-07-27.csv")

soil_Harv$horizon <- substr(soil_Harv$sampleID, 10, 10)
soil_Harv$dateID <- substr(soil_Harv$collectDate, 1, 7)
select_columns <- c("siteID", "plotID", "dateID","horizon", "nlcdClass", "collectDate",
                    "sampleTiming", "standingWaterDepth", "sampleID",
                    "soilTemp", "litterDepth", "sampleTopDepth", "sampleBottomDepth", "soilMoisture",
                    "geneticSampleID", "soilInCaClpH","soilInWaterpH","CNratio", "nitrogenPercent", "organicCPercent")
soils <- soil_Harv[, select_columns]

soils$day <- substr(soils$collectDate, 1, 10)



# Nitrogen rate data

library(neonUtilities)
library(neonNTrans)
soilData <- loadByProduct(site = "HARV", dpID = "DP1.10086.001", startdate = "2014-01-01", package = "basic", check.size = F, token="eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJ6cndlcmJpbkBidS5lZHUiLCJzY29wZSI6InJhdGU6cHVibGljIiwiaXNzIjoiaHR0cHM6Ly9kYXRhLm5lb25zY2llbmNlLm9yZy8iLCJleHAiOjE3NTc4ODg2NzAsImlhdCI6MTYwMDIwODY3MCwiZW1haWwiOiJ6cndlcmJpbkBidS5lZHUifQ.8eW8vxUOiton-kQ_Xyvva0QSHD_BDd2E5IGeNKW3WHib-m7UpTnEhGFAUlAHGdsUyz-dKE1jMOAGS5A_NRYXGg")

#soilData <- readRDS("/projectnb/talbot-lab-data/zrwerbin/interactions/data/NEON_soilData.rds")
out <- def.calc.ntrans(kclInt = soilData$ntr_internalLab,
                       kclIntBlank = soilData$ntr_internalLabBlanks,
                       kclExt = soilData$ntr_externalLab,
                       soilMoist = soilData$sls_soilMoisture,
                       dropAmmoniumFlags = "blanks exceed sample value",
                       dropNitrateFlags = "blanks exceed sample value",
                       dropConditions = c("deprecatedMethod", "other"))
summary <- out$data_summary
summary$siteID <- substr(summary$sampleID, 1, 4)
summary$plotID <- substr(summary$sampleID, 1, 8)
summary$horizon <- substr(summary$sampleID, 10, 10)
summary$dateID <- substr(summary$collectDate, 1, 7)
summary$collectDate <- NULL



all_N_data <- out$all_data %>% select(siteID,plotID,sampleID,incubationPairID,nTransBoutType,soilMoisture,soilAmmoniumNugPerGram,soilNitrateNitriteNugPerGram,soilInorganicNugPerGram,netNminugPerGramPerDay,netNitugPerGramPerDay)
#
#
# N_wide = all_N_data %>% pivot_wider(id_cols = c(siteID, plotID, incubationPairID),
# 	values_from =
# 																			c(soilMoisture, soilAmmoniumNugPerGram, soilNitrateNitriteNugPerGram, soilInorganicNugPerGram),
# 																		names_from = nTransBoutType, #names_prefix = nTransBoutType,
# 																		names_expand=T, values_fn = list)

initial <- all_N_data[all_N_data$nTransBoutType=="tInitial",]  %>%
	rename(soilMoisture_initial = soilMoisture,
				 soilAmmoniumNugPerGram_initial = soilAmmoniumNugPerGram,
				 soilNitrateNitriteNugPerGram_initial = soilNitrateNitriteNugPerGram,
				 soilInorganicNugPerGram_initial=soilInorganicNugPerGram) %>%
	select(-c(netNminugPerGramPerDay,netNitugPerGramPerDay,nTransBoutType))

final <- all_N_data[all_N_data$nTransBoutType=="tFinal",] %>%
	rename(soilMoisture_final = soilMoisture,
				 soilAmmoniumNugPerGram_final = soilAmmoniumNugPerGram,
				 soilNitrateNitriteNugPerGram_final = soilNitrateNitriteNugPerGram,
				 soilInorganicNugPerGram_final=soilInorganicNugPerGram,
				 sampleID_final = sampleID) %>% select(-nTransBoutType)

N_rates = merge(initial, final, by=c("siteID","plotID","incubationPairID")) %>% filter(!is.na(incubationPairID))


#merged <- merge(summary, soils, all.x=T)
#soils <- merge(summary, soils, all.y=T)


# Add in sensor values
daily.SWC <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEONSoilMoist_raw_allsites.rds")
daily.harv_SWC <- daily.SWC %>%
  #filter(siteID=="HARV" & month %in% c("2018-04","2018-05", "2018-06", "2018-07", "2018-08","2018-09","2018-10")) %>%
  filter(siteID=="HARV") %>%
  group_by(siteID, day) %>%
  summarize(sensor_moisture=mean(VSWCMean, na.rm=T))

daily.harv_SWC$day <- substr(daily.harv_SWC$day, 1, 10)

soils <- merge(soils, daily.harv_SWC, by.x = c("siteID","day"), by.y = c("siteID","day"), all.x = T)



merged_all <- merge(N_rates, soils, all=T)

#N_rates <- merge(N_rates, daily.harv_SWC, by.x = c("siteID","day"), by.y = c("siteID","day"), all.x = T)
#soils <- merge(soils, harv_SWC, by.x = c("siteID","dateID"), by.y = c("siteID","month"), all.x = T)

#merged <- merge(merged, harv_SWC, by.x = c("siteID","dateID"), by.y = c("siteID","month"), all.x = T)
#soils <- merge(soils, harv_SWC, by.x = c("siteID","dateID"), by.y = c("siteID","month"), all.x = T)

saveRDS(N_rates, "/projectnb2/talbot-lab-data/zrwerbin/interactions/data/HARV_DAMM_val_N_rates.rds")
saveRDS(merged_all, "/projectnb2/talbot-lab-data/zrwerbin/interactions/data/HARV_DAMM_all_data.rds")

