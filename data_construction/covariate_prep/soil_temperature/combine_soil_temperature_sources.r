
library(dplyr)
library(tidyr)
library(lme4)
library(merTools)
library(ggplot2)

# Read in NEON soil temp
NEON_temp_monthly <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/NEONSoilTemp_monthly_allsites.rds")
#NEON_temp_monthly <- NEON_temp_monthly %>% filter(!grepl("2020", dateID))

# Read in daymet daily/weekly/monthly air temp 
daymet_monthly <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/daymet_monthly.rds")
daymet_monthly <- daymet_monthly %>% filter(!dateID %in% c(201301, 201302, 201303, 201304, 201305)) %>% dplyr::rename(month = dateID)


NEON_temp_monthly$month <- gsub("-","",NEON_temp_monthly$month)
monthly <- merge(NEON_temp_monthly, daymet_monthly, by = c("siteID", "month"), all=T)

# NEON_temp_wide <- NEON_temp_monthly %>% pivot_wider(id_cols = dateID, 
#                                                     names_from = siteID, 
#                                                     values_from = NEON_temp_mean, 
#                                                     names_glue = "{siteID}_temp_mean") %>% 
#   arrange(dateID)

# Get list of site/dates where we have microbial data
microb.avail = Z10::dp.avail("DP1.10108.001")
avail.site.months <- microb.avail %>% 
  unnest(months) %>% unnest(site) %>% dplyr::rename(siteID = site, month = months) %>% 
  mutate(date = as.Date(paste0(month, "-01"))) 

avail.site.months <- rbind(avail.site.months, c("OSBS","2013-06","2013-08-01"))
avail.site.months <- rbind(avail.site.months, c("OSBS","2013-06","2013-07-01"))

sites <- unique(avail.site.months$siteID)
thisyear <- cbind.data.frame(siteID = unique(sites), month = "2020-09", date = "2020-09-01") # gap fill until present
avail.site.months <- rbind.data.frame(avail.site.months, thisyear) %>% 
  group_by(siteID) %>% padr::pad() %>% # ignore error, it's probably within the padr package
  mutate(month = gsub("-","",substr(date, 1, 7)))


##### 3. FIT CALIBRATION MODELS, CONVERT DATA, ESTIMATE UNCERTAINTIES FOR each site x source #####
# Using bootstrapping method modified from the example here: https://www.rdocumentation.org/packages/lme4/versions/1.1-23/topics/predict.merMod

# Merge sources together into a main data frame
sources_merged <- merge(avail.site.months, monthly, all.x=T)



### FOR DAYMET DATA ###
# remove unnecessary columns, and na.omit due to bootstrapping function
pred_data <- sources_merged[which(!is.na(sources_merged$maxtemp_mean)),] %>% 
  dplyr::select(siteID, month, mintemp_mean, maxtemp_mean, precip_mean,  NEON_temp_mean) %>% 
  mutate(siteID = as.factor(siteID)) %>% 
  na.omit()
## random slopes linear model to get predictions
lmm <- lmer(data = pred_data, NEON_temp_mean ~ mintemp_mean + (mintemp_mean|siteID))
# model hasn't converged, but fits better than *only random slopes* or *only random intercepts*
sources_merged$daymet_est <- predict(lmm, newdata = sources_merged)

## Get site-level uncertainties 
predict.fun <- function(my.lmm) predict(my.lmm, newdata = pred_data)   # predict.merMod 
# Make predictions in 500 bootstraps, to get standard deviation of each point.
lmm.boots <- bootMer(lmm, predict.fun, nsim = 500)
pred_data$daymet_est_sd <- apply(lmm.boots$t, 2, sd, na.rm=T)
# Get mean SD for each site.
daymet.pred.out <- pred_data %>% group_by(siteID) %>% 
  dplyr::summarize(daymet_site_sd = mean(daymet_est_sd, na.rm=T))
# Add into main dataframe.
sources_merged <- merge(sources_merged, daymet.pred.out, all.x = T)


#### EVALUATE OUTPUT
# how well does DAYMET do now? # .528
summary(lm(NEON_temp_mean ~ daymet_est, sources_merged))

# Now select the means and uncertainties to use for final dataset, NEON > SCAN > SMAP > SMOS
soil.temperature.out <- sources_merged %>% mutate(temperature_out = ifelse(!is.na(NEON_temp_mean), NEON_temp_mean, daymet_est),
                                                   source = ifelse(!is.na(NEON_temp_mean), "NEON", "DAYMET"),
                                                   temperature_sd_out = ifelse(!is.na(NEON_temp_mean), NEON_temp_sd, daymet_site_sd))
# Error bars = 1 standard deviation
soil.temperature.out$hi <- soil.temperature.out$temperature_out + soil.temperature.out$temperature_sd_out
soil.temperature.out$low <- soil.temperature.out$temperature_out - soil.temperature.out$temperature_sd_out

# Mean-center and scale values for model
temp.scale <- scale(soil.temperature.out$temperature_out)
scale_factor <- attr(temp.scale, "scaled:scale")
soil.temperature.out$temperature <- temp.scale
# SD transformed by same scale factor 
soil.temperature.out$temperature_sd <- soil.temperature.out$temperature_sd_out/scale_factor



saveRDS(soil.temperature.out, "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/monthly_soil_temperature.rds")


# VISUALIIIIZE
p <- ggplot(soil.temperature.out, aes(x = date, y = temperature_out, color = source)) + 
  geom_point(show.legend = F) + facet_wrap(~siteID) +#, nrow = 8) + 
  geom_errorbar(aes(ymin=low,ymax=hi),alpha=0.3) + theme_minimal() + theme (legend.position = c(0.83, 0.04), text = element_text(size = 16)) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 5, alpha = 1))) +
  ylab("Soil temperature") + xlab(NULL)
p
#ggsave(p, filename = "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/figures/soil_temperature_calibration.png", device = "png", width = 15, height = 12, units = "in")

# Plots look alright: 
# Uncertainty is low for all sites. Daymet trends look incredibly similar to NEON trends.




fit <- lm(daymet_est ~ NEON_temp_mean, soil.temperature.out)
p2 <- ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
		geom_point() +
		stat_smooth(method = "lm", col = "red") +
		labs(title = "Calibration fits", subtitle = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 4))) + 
	xlab("Observed values from NEON") + ylab("Predicted values from Daymet") + 
	theme_minimal() + theme (text = element_text(size = 12)) 
p2

library(ggpubr)
ggarrange(p, p2, widths = c(2,1), labels = c("A","B"))

