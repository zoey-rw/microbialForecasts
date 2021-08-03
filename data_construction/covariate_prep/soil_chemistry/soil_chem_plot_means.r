#Core level aggregation.
#Getting core, plot, site and global means and sd's for all observations made at the core level.
#This is currently soil pH, soil %C and soil C:N.
#clear environment, source paths.
library(dplyr)

# Read in microbial data
dat_sequenceMetadata <- read.csv("/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/data/NEON/sequence_metadata/mmg_soilMetadata_all_2020-09-26.csv") %>% select(siteID, plotID) %>% distinct()
# # # get rid of duped dnaSampleIDs.
# dat_sequenceMetadata$geneticSampleID <- gsub("-DNA[1234]","",dat_sequenceMetadata$dnaSampleID)
# dat_sequenceMetadata <- dat_sequenceMetadata[!duplicated(dat_sequenceMetadata$geneticSampleID),]

# Read in soil sample data
dat_soil <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/soilSample_data_allsites.rds")

div_dat <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_ITS.rds")[[1]] %>% select(siteID, plotID, horizon) %>% distinct()

# plot(dat_soil$pH, dat_soil$CNratio)
# plot(dat_soil$pH, dat_soil$organicCPercent)
# Correlation between pH and C:N ratio is too strong, so we can't keep both. 
# Using percent carbon instead.

#finalize columns for core.level.----
core.level <- dat_soil[,c('geneticSampleID','year','sampleID','collectDate','siteID','plotID','dateID','horizon','soilMoisture','nitrogenPercent','soilInCaClpH','organicCPercent')] %>% 
  dplyr::rename(pH = soilInCaClpH,
                pC = organicCPercent)

# Add in plots from microbial data
core.level <- merge(core.level, div_dat, all=T)
core.level <- merge(core.level, dat_sequenceMetadata, all=T)

# Scale by mean and standard deviation
core.level$pH <- scale(core.level$pH)
core.level$pC <- scale(core.level$pC)

# Get summaries at a bunch of different spatial and temporal scales
core.level.plot <- core.level %>% filter(!grepl("LOGISTICAL", sampleID)) %>% 
  group_by(siteID,horizon) %>%  dplyr::mutate(site_anytime_pH = mean(pH, na.rm=T),
                                                     site_anytime_pH_sd = sd(pH, na.rm=T),
                                                     site_anytime_pC = mean(pC, na.rm=T),
                                                     site_anytime_pC_sd = sd(pC, na.rm=T)) %>% ungroup() %>%
  group_by(siteID,year,plotID,horizon) %>% dplyr::mutate(plot_year_pH = mean(pH, na.rm=T),
                                                           plot_year_pH_sd = sd(pH, na.rm=T)) %>% ungroup() %>%
  group_by(siteID,plotID,horizon) %>% dplyr::mutate(plot_anytime_pH = mean(pH, na.rm=T), 
                                                         plot_anytime_pH_sd = sd(pH, na.rm=T),
                                                         plot_anytime_pC = mean(pC, na.rm=T), 
                                                         plot_anytime_pC_sd = sd(pC, na.rm=T)) %>% ungroup() %>% 
  group_by(siteID,dateID,plotID,horizon) %>% dplyr::mutate(plot_pH = mean(pH, na.rm=T), 
                                     plot_pH_sd = sd(pH, na.rm=T)) %>% ungroup() %>% 
  group_by(siteID,plotID) %>% dplyr::mutate(plot_anytime_anyhorizon_pH = mean(pH, na.rm=T), 
                                                    plot_anytime_anyhorizon_pH_sd = sd(pH, na.rm=T)) %>% 
  group_by(siteID) %>% dplyr::mutate(site_anytime_anyhorizon_pC = mean(pC, na.rm=T), 
                                            site_anytime_anyhorizon_pC_sd = sd(pC, na.rm=T)) %>% as.data.frame()

# We're not choosing the best resolution possible: 
# We're using pH at the plot level using all available measurements, so that it is static over time
# Same for pC, which is filled in with a lot of site-level data, and missing for 6 sites
core.level.plot <- core.level.plot %>% mutate(final_pH = ifelse(#!is.na(pH), pH, 
                                                                #ifelse(!is.na(plot_year_pH), plot_year_pH, 
                                                                       #ifelse(
                                                                         !is.na(plot_anytime_pH), plot_anytime_pH, 
                                                                              ifelse(!is.na(plot_anytime_anyhorizon_pH), plot_anytime_anyhorizon_pH, 
                                                                                     ifelse(!is.na(site_anytime_pH), site_anytime_pH, NA))),
                                              final_pH_sd = ifelse(#!is.na(pH), .01, 
                                                                   #ifelse(!is.na(plot_year_pH), plot_year_pH_sd, 
                                                                          #ifelse(
                                                                          	!is.na(plot_anytime_pH), plot_anytime_pH_sd, 
                                                                                 ifelse(!is.na(plot_anytime_anyhorizon_pH), plot_anytime_anyhorizon_pH_sd, 
                                                                                        ifelse(!is.na(site_anytime_pH), site_anytime_pH_sd, NA))),
                                              final_pC = ifelse(!is.na(plot_anytime_pC), plot_anytime_pC, 
                                                                ifelse(!is.na(site_anytime_pC), site_anytime_pC, 
                                                                       ifelse(!is.na(site_anytime_anyhorizon_pC), site_anytime_anyhorizon_pC, NA))),
                                              final_pC_sd = ifelse(!is.na(plot_anytime_pC), plot_anytime_pC_sd, 
                                                                ifelse(!is.na(site_anytime_pC), site_anytime_pC_sd, 
                                                                       ifelse(!is.na(site_anytime_anyhorizon_pC), site_anytime_anyhorizon_pC_sd, NA)))) %>% as.data.frame()
  
# Get the horizon for each site with the most site-date combos
dom_horizons <- dat_soil %>% distinct(siteID, dateID, horizon) %>% 
  group_by(siteID, horizon) %>% dplyr::count() %>% ungroup() %>% group_by(siteID) %>% filter(n == max(n)) 

# If there's a tie, use the mineral layer, which usually has more core observations
duped <- dom_horizons[duplicated(dom_horizons$siteID),]$siteID
dom_duped <- dom_horizons[dom_horizons$siteID %in% duped & dom_horizons$horizon == "M",]
dom_horizons_out <- rbind(dom_horizons[!dom_horizons$siteID %in% duped,], dom_duped)
saveRDS(dom_horizons_out, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/dominantHorizonsSite.rds")
# Subset to the most dominant horizon for each site 
soil.chem.out <- core.level.plot %>% filter(paste(siteID, horizon) %in% paste(dom_horizons_out$siteID, dom_horizons_out$horizon))

#Some plots in core-level not in plot-level.
to_add <- as.character(unique(dat_sequenceMetadata[!(dat_sequenceMetadata$plotID %in% soil.chem.out$plotID),]$plotID))
to_add <- data.frame(to_add)
colnames(to_add) <- 'plotID'
soil.chem.out <- plyr::rbind.fill(soil.chem.out,to_add)
soil.chem.out$siteID <- substring(soil.chem.out$plotID, 1, 4)
soil.chem.out <- soil.chem.out[order(soil.chem.out$plotID),]


# remove plots in plot-level that aren't in microbial data.
# soil.chem.out <- soil.chem.out[soil.chem.out$plotID %in% dat_sequenceMetadata$plotID,]

# now keep only the distinct values for each plot
soil.chem.out.plot <- soil.chem.out %>% distinct(siteID, plotID, final_pH, final_pH_sd, final_pC, final_pC_sd) %>% dplyr::rename(pH = final_pH,
																																																									pH_sd = final_pH_sd,
																																																									pC = final_pC,
																																																									pC_sd = final_pC_sd)
# # there's one missing plot, let's fill it in (clunky method)
# osbs_pH <- mean(soil.chem.out.plot[soil.chem.out.plot$siteID == "OSBS",]$pH)
# osbs_pH_sd <- mean(soil.chem.out.plot[soil.chem.out.plot$siteID == "OSBS",]$pH_sd)
# osbs_pC <- mean(soil.chem.out.plot[soil.chem.out.plot$siteID == "OSBS",]$pC)
# osbs_pC_sd <- mean(soil.chem.out.plot[soil.chem.out.plot$siteID == "OSBS",]$pC_sd)
# 
# soil.chem.out.plot <- rbind(soil.chem.out.plot, c("OSBS", "OSBS_006", osbs_pH, osbs_pH_sd, osbs_pC, osbs_pC_sd))

#save output.----
saveRDS(soil.chem.out.plot, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/soilChemPlot.rds")









# Plot the sample pH and the plot-date-pH
p1 <- ggplot(soil.chem.out.plot[soil.chem.out.plot$siteID %in%  c("BART","NIWO","TALL"),]) + 
	geom_boxplot(aes(x = plotID, y = pH)) +
	geom_point(aes(x = plotID, y = pH)) + facet_grid(~siteID, scales = "free_x") + 
	ggtitle("pH at 3 representative NEON sites", "Plot-level means were treated as constant over time.") +
	ylab("pH in CaCl, soil core") + xlab("Plot ID") + theme(axis.text.x = element_text(angle = 310))

# Plot the sample pC and corresponding plot pC
p2 <- ggplot(soil.chem.out.plot[soil.chem.out.plot$siteID %in% c("BART","NIWO","TALL"),]) + 
	geom_boxplot(aes(x = plotID, y = pC)) +
	geom_point(aes(x = plotID, y = pC)) + facet_grid(~siteID, scales = "free_x", drop = T) + 
	ggtitle("Percent carbon at 3 representative NEON sites", "Plot-level means were treated as constant over time.") +
	ylab("%C, soil core") + xlab("Plot ID") + theme(axis.text.x = element_text(angle = 310))

library(ggpubr)
ggarrange(p1, p2, nrow=2, labels = c("A","B"))



chem_df <- dat_soil %>% select(siteID, plotID, soilInCaClpH, organicCPercent, year)
ggplot(chem_df[chem_df$siteID %in% c("CPER"),]) + 
	geom_point(aes(x = plotID, y = organicCPercent)) + 
	#facet_grid(~siteID, scales = "free_x", drop = T) + 
 theme(axis.text.x = element_text(angle = 310))
