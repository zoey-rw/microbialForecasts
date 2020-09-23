# Combine soil moisture data products and calculate uncertainties for each site-source combo.
library(dplyr)
library(tidyr)
library(lme4)
library(merTools)
library(ggplot2)

##### 1. Read in data sources, which are already aggregated by month #####

## Get list of NEON fieldsites/locations ##
fieldsites <- read.csv("https://www.neonscience.org/science-design/field-sites/export")

## Get SMV data (includes SCAN and SMAP) ##  
smv_month <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/DAAC_SMV/monthly_SMV_allsites.rds")
# Only need surface zone from SCAN and SMAP, not GRACE or rootzones
smv_month <- smv_month %>% dplyr::select(-c(mean_GRACE_s, min_GRACE_s, max_GRACE_s, min_SMAP_r, max_SMAP_r, mean_SMAP_r, min_SCAN_r, mean_SCAN_r, max_SCAN_r ))

## Get NEON moisture data ##
neon <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/NEONSoilMois_monthly_allsites.rds")
neon <- neon %>% mutate(neon_mean = NEON_moist_mean*100,
                        neon_sd = NEON_moist_sd * 100) # units for everything else are 0-100 instead of 0-1

## Get SMOS data ##
smos <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/SMOS/monthly_SMOS_allsites.rds")

#### 2. Create dataframe with all of the sites/dates that we need data for #####
# Get list of site/dates where we have microbial data
microb.avail = Z10::dp.avail("DP1.10108.001")
avail.site.months <- microb.avail %>% 
  unnest(months) %>% unnest(site) %>% rename(siteID = site, month = months) %>% 
  mutate(date = as.Date(paste0(month, "-01"))) 

sites <- unique(avail.site.months$siteID)
thisyear <- cbind.data.frame(siteID = unique(sites), month = "2020-09", date = "2020-09-01") # gap fill until present
avail.site.months <- rbind.data.frame(avail.site.months, thisyear) %>% 
  group_by(siteID) %>% padr::pad() %>% # ignore error, it's probably within the padr package
  mutate(month = substr(date, 1, 7))

##### 3. FIT CALIBRATION MODELS, CONVERT DATA, ESTIMATE UNCERTAINTIES FOR each site x source #####
# same thing for three sources, 
# except the SCAN model has random intercepts and not slopes due to unidentifiability (only 4 sites())
# Using bootstrapping method modified from the example here: https://www.rdocumentation.org/packages/lme4/versions/1.1-23/topics/predict.merMod

# Merge sources together into a main data frame
SMV_neon_merged <- merge(neon, smv_month, all=T)
SMV_neon_merged <- merge(avail.site.months, SMV_neon_merged, all.x=T)
SMV_neon_merged$month <-  gsub("-","", SMV_neon_merged$month)
all_sources_merged <- merge(SMV_neon_merged, smos, by = c("siteID", "month"), all.x=T)

  
  ### FOR SMAP DATA ###
  # remove unnecessary columns, and na.omit due to bootstrapping function
  pred_data <- all_sources_merged[which(!is.na(all_sources_merged$mean_SMAP_s)),] %>% 
    dplyr::select(siteID, month, mean_SMAP_s, neon_mean) %>% 
      mutate(siteID = as.factor(siteID)) %>% 
      na.omit()
  ## random slopes linear model to get predictions
  lmm <- lmer(data = pred_data, neon_mean ~ mean_SMAP_s + (mean_SMAP_s|siteID))
  # model hasn't converged, but fits better than *only random slopes* or *only random intercepts*
  all_sources_merged$SMAP_est <- predict(lmm, newdata = all_sources_merged)
  
  ## Get site-level uncertainties 
  predict.fun <- function(my.lmm) predict(my.lmm, newdata = pred_data)   # predict.merMod 
  # Make predictions in 500 bootstraps, to get standard deviation of each point.
  lmm.boots <- bootMer(lmm, predict.fun, nsim = 500)
  pred_data$SMAP_est_sd <- apply(lmm.boots$t, 2, sd, na.rm=T)
  # Get mean SD for each site.
  SMAP.pred.out <- pred_data %>% group_by(siteID) %>% 
    dplyr::summarize(SMAP_site_sd = mean(SMAP_est_sd, na.rm=T))
  # Add into main dataframe.
  all_sources_merged <- merge(all_sources_merged, SMAP.pred.out, all.x = T)
  
  
  #### FOR SCAN DATA ####
  # remove unnecessary columns, and na.omit due to bootstrapping function
  pred_data <- all_sources_merged[which(!is.na(all_sources_merged$mean_SCAN_s)),] %>% 
    dplyr::select(siteID, month, mean_SCAN_s, neon_mean) %>% 
    mutate(siteID = as.factor(siteID)) %>% 
    na.omit()
  
  SCAN.site.ind <- which(all_sources_merged$siteID %in% unique(pred_data$siteID))
  ## random intercept linear model to get predictions (slopes was unidentifiable)
  lmm <- lmer(data = pred_data, neon_mean ~ mean_SCAN_s + (1|siteID))
  all_sources_merged$SCAN_est <- NA
  all_sources_merged[SCAN.site.ind,]$SCAN_est <- predict(lmm, newdata = all_sources_merged[SCAN.site.ind,])
  
  ## Get site-level uncertainties 
  predict.fun <- function(my.lmm) predict(my.lmm, newdata = pred_data)   # predict.merMod 
  # Make predictions in 500 bootstraps, to get standard deviation of each point.
  lmm.boots <- bootMer(lmm, predict.fun, nsim = 500)
  pred_data$SCAN_est_sd <- apply(lmm.boots$t, 2, sd, na.rm=T)
  # Get mean SD for each site.
  SCAN.pred.out <- pred_data %>% group_by(siteID) %>% 
    dplyr::summarize(SCAN_site_sd = mean(SCAN_est_sd, na.rm=T))
  # Add into main dataframe.
  all_sources_merged <- merge(all_sources_merged, SCAN.pred.out, all.x = T)

  
  ### FOR SMOS DATA ###
  # remove unnecessary columns, and na.omit due to bootstrapping function
  pred_data <- all_sources_merged[which(!is.na(all_sources_merged$SM)),] %>% 
    dplyr::select(siteID, month, SM, neon_mean) %>% 
    mutate(siteID = as.factor(siteID)) %>% 
    na.omit()
 
   ## random slopes linear model to get predictions
  SMOS.site.ind <- which(all_sources_merged$siteID %in% unique(pred_data$siteID))
  ## random intercept linear model to get predictions (slopes was unidentifiable)
  lmm <- lmer(data = pred_data, neon_mean ~ SM + (SM|siteID))
  all_sources_merged$SMOS_est <- NA
  all_sources_merged[SMOS.site.ind,]$SMOS_est <- predict(lmm, newdata = all_sources_merged[SMOS.site.ind,])
  
  ## Get site-level uncertainties  
  predict.fun <- function(my.lmm) predict(my.lmm, newdata = pred_data)   # predict.merMod 
  # Make predictions in 100 bootstraps, to get standard deviation of each point.
  lmm.boots <- bootMer(lmm, predict.fun, nsim = 100)
  pred_data$SMOS_est_sd <- apply(lmm.boots$t, 2, sd, na.rm=T)
  # Get mean SD for each site.
  SMOS.pred.out <- pred_data %>% group_by(siteID) %>% 
    dplyr::summarize(SMOS_site_sd = mean(SMOS_est_sd, na.rm=T))
  # Add into main dataframe.
  all_sources_merged <- merge(all_sources_merged, SMOS.pred.out, all.x = T)
  
  
  #### EVALUATE OUTPUT
  # how well does SCAN do now? # R2 = .5094
  summary(lm(neon_mean ~ SCAN_est, all_sources_merged))
  # how well does SMAP do now? # .57
  summary(lm(neon_mean ~ SMAP_est, all_sources_merged))
  # how well does SMOS do now? # .528
  summary(lm(neon_mean ~ SMOS_est, all_sources_merged))

  # Now select the means and uncertainties to use for final dataset, NEON > SCAN > SMAP > SMOS
  soil.moisture.out <- all_sources_merged %>% mutate(moisture = ifelse(!is.na(neon_mean), neon_mean,
                                                                    ifelse(!is.na(SCAN_est), SCAN_est,
                                                                    ifelse(!is.na(SMAP_est), SMAP_est,
                                                                           SMOS_est))),
                                                     source = ifelse(!is.na(neon_mean), "NEON",
                                                                     ifelse(!is.na(SCAN_est), "SCAN",
                                                                            ifelse(!is.na(SMAP_est), "SMAP",
                                                                                   "SMOS"))),
                                                     moisture_sd = ifelse(!is.na(neon_mean), neon_sd,
                                                                          ifelse(!is.na(SCAN_est), SCAN_site_sd,
                                                                                 ifelse(!is.na(SMAP_est), SMAP_site_sd,
                                                                                        SMOS_site_sd))))
  # Error bars = 1 standard deviation
  soil.moisture.out$hi <- soil.moisture.out$moisture + soil.moisture.out$moisture_sd
  soil.moisture.out$low <- soil.moisture.out$moisture - soil.moisture.out$moisture_sd
  p <- ggplot(soil.moisture.out, aes(x = date, y = moisture, color = source)) + 
    geom_point(show.legend = F) + facet_wrap(~siteID, nrow = 8) + 
    geom_errorbar(aes(ymin=low,ymax=hi),alpha=0.3) + theme_minimal() + theme (legend.position = c(0.7, 0.04)) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 5, alpha = 1))) + ylab("Soil moisture")
  p

  ggsave(p, filename = "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/figures/soil_moisture_calibration.png", device = "png", width = 15, height = 12, units = "in")
  
# Plots look alright: 
# Uncertainty increases for sites with less NEON data.
# As expected, worse data sources have larger uncertainties.
  
  saveRDS(soil.moisture.out, "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/monthly_soil_moisture.rds")
