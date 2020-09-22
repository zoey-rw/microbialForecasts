# Combine soil moisture data products and calculate uncertainties for each site-source combo.
library(dplyr)
library(tidyr)
library(lme4)
library(merTools)


# Get list of NEON fieldsites/locations
fieldsites <- read.csv("https://www.neonscience.org/science-design/field-sites/export")

#### Get SMV data ####
smv_month <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/DAAC_SMV/monthly_SMV_allsites.rds")
# Only need surface zone from SCAN and SMAP
smv_month <- smv_month %>% dplyr::select(-c(mean_GRACE_s, min_GRACE_s, max_GRACE_s, min_SMAP_r, max_SMAP_r, mean_SMAP_r, min_SCAN_r, mean_SCAN_r, max_SCAN_r ))

#### 2. GET NEON MOISTURE DATA ####
neon <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/NEONSoilMois_monthly_allsites.rds")
neon <- neon %>% mutate(neon_mean = simpleMean*100)

#### Get SMOS data ####
smos <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/SMOS/monthly_SMOS_allsites.rds")

#### 3. GET NEON SITES WHERE WE HAVE MICROBES ####
# get list of site/dates where we have microbial data
microb.avail = Z10::dp.avail("DP1.10108.001")
avail.site.months <- microb.avail %>% 
  unnest(months) %>% unnest(site) %>% rename(siteID = site, month = months) %>% 
  mutate(date = as.Date(paste0(month, "-01"))) 

sites <- unique(avail.site.months$siteID)
thisyear <- cbind.data.frame(siteID = unique(sites), month = "2020-09", date = "2020-09-01")
avail.site.months <- rbind.data.frame(avail.site.months, thisyear) %>% 
  group_by(siteID) %>% padr::pad() %>% 
  mutate(month = substr(date, 1, 7))


# Merge things together
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
  all_sources_merged$SMAP_est <- predict(lmm, newdata = all_sources_merged)
  
  ## Get site-level uncertainties 
  predict.fun <- function(my.lmm) predict(my.lmm, newdata = pred_data)   # predict.merMod 
  # Make predictions in 100 bootstraps, to get standard deviation of each point.
  lmm.boots <- bootMer(lmm, predict.fun, nsim = 100)
  pred_data$SMAP_est_sd <- apply(lmm.boots$t, 2, sd, na.rm=T)
  # Get mean SD for each site.
  SMAP.pred.out <- pred_data %>% group_by(siteID) %>% 
    dplyr::summarize(SMAP_site_sd = mean(SMAP_est_sd, na.rm=T))
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
  # Make predictions in 100 bootstraps, to get standard deviation of each point.
  lmm.boots <- bootMer(lmm, predict.fun, nsim = 100)
  pred_data$SCAN_est_sd <- apply(lmm.boots$t, 2, sd, na.rm=T)
  # Get mean SD for each site.
  SCAN.pred.out <- pred_data %>% group_by(siteID) %>% 
    dplyr::summarize(SCAN_site_sd = mean(SCAN_est_sd, na.rm=T))
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
  all_sources_merged <- merge(all_sources_merged, SMOS.pred.out, all.x = T)
  
  
  #### EVALUATE OUTPUT
  # how well does SCAN do now? # R2 = .5094
  summary(lm(neon_mean ~ SCAN_est, all_sources_merged))
  # how well does SMAP do now? # .57
  summary(lm(neon_mean ~ SMAP_est, all_sources_merged))
  # how well does SMOS do now? # .528
  summary(lm(neon_mean ~ SMOS_est, all_sources_merged))

  # Now select the means and uncertainties to use for final dataset
  all_sources_merged$NEON_moist_mean <- all_sources_merged$NEON_moist_mean * 100
  all_sources_merged$NEON_moist_sd <- all_sources_merged$NEON_moist_sd * 100
  soil.moisture.out <- all_sources_merged %>% mutate(moisture = ifelse(!is.na(NEON_moist_mean), NEON_moist_mean,
                                                                    ifelse(!is.na(SCAN_est), SCAN_est,
                                                                    ifelse(!is.na(SMAP_est), SMAP_est,
                                                                           SMOS_est))),
                                                     source = ifelse(!is.na(NEON_moist_mean), "NEON",
                                                                     ifelse(!is.na(SCAN_est), "SCAN",
                                                                            ifelse(!is.na(SMAP_est), "SMAP",
                                                                                   "SMOS"))),
                                                     moisture_sd = ifelse(!is.na(NEON_moist_mean), NEON_moist_sd,
                                                                          ifelse(!is.na(SCAN_est), SCAN_site_sd,
                                                                                 ifelse(!is.na(SMAP_est), SMAP_site_sd,
                                                                                        SMOS_site_sd))))
  # Error bars = 1 standard deviation
  soil.moisture.out$hi <- soil.moisture.out$moisture + soil.moisture.out$moisture_sd
  soil.moisture.out$low <- soil.moisture.out$moisture - soil.moisture.out$moisture_sd
  p <- ggplot(soil.moisture.out, aes(x = date, y = moisture, color = source)) + 
    geom_point(show.legend = F) + facet_wrap(~siteID, nrow = 8) + 
    geom_errorbar(aes(ymin=low,ymax=hi),alpha=0.3) + theme_minimal() + theme (legend.position = c(0.7, 0.04)) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 5, alpha = 1)))
  p

