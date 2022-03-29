# New plant script because others got...messy
# Prep plant covariates for forecasting models: tree basal area (averaged per plot, gap-filled by site), and LAI (monthly by site)

# Load R packages, install if necessary.
#devtools::install_github("admahood/neondiversity") 
library(neondiversity) 
library(tidyverse)
library(zoo)
library(corrplot) 

# Output from clean_relEM.r
basal_relEM <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEON_relEM_plot.level.rds")
# Output from clean_plantDiv_data.r
div <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/annual_plant_data.rds")

poss_site_plots <- div %>% tidyr::expand(nesting(siteID, plotID))
basal_relEM <- merge(poss_site_plots, basal_relEM, all=T)

# Fill missing values with site means
relEM_fill_site <- basal_relEM %>% group_by(siteID) %>% mutate(site_mean_relEM = mean(relEM, na.rm=T)) %>% 
	mutate(gap_filled = ifelse(is.na(relEM), TRUE, FALSE),
				 relEM = ifelse(is.na(relEM), site_mean_relEM, relEM)) %>% ungroup()

# Look into why some plots have no trees. What is the land cover type?
missing <- unique(relEM_fill_site[is.na(relEM_fill_site$site_mean_relEM),]$siteID)
site_info <- read.csv("https://www.neonscience.org/sites/default/files/NEON_Field_Site_Metadata_20210226_0.csv")
site_info <- site_info %>% select(field_site_id, field_dominant_nlcd_classes)
site_info[site_info$field_site_id %in% missing,]
# Makes sense that these sites should have no trees
# field_site_id                       field_dominant_nlcd_classes
# 4           BARR                      Emergent Herbaceous Wetlands
# 15          CPER                              Grassland/Herbaceous
# 18          DCFS                              Grassland/Herbaceous
# 30          JORN                                       Shrub/Scrub
# 32          KONA                                  Cultivated Crops
# 47          OAES                  Grassland/Herbaceous|Shrub/Scrub
# 65          STER                                  Cultivated Crops
# 73          TOOL                           Dwarf Scrub|Shrub/Scrub
# 79          WOOD Emergent Herbaceous Wetlands|Grassland/Herbaceous

# Fill missing site-level means with zeros 
relEM_fill_site <- relEM_fill_site %>% 
	mutate(gap_filled = ifelse(is.na(relEM), TRUE, FALSE),
				 relEM = ifelse(is.na(relEM), 0, relEM)) %>% ungroup()

# Save original values
relEM_fill_site$relEM_orig <- relEM_fill_site$relEM

# Mean center and scale data
relEM_fill_site$relEM <- scale(relEM_fill_site$relEM_orig)

saveRDS(relEM_fill_site, "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/mean_relEM_data.rds")




library(tidyverse)
library(data.table)


lai_dat_orig <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/modis/LAI_allsites.rds")

lai_dat <- lai_dat_orig %>% mutate(dateID = substr(calendar_date, 1, 7), 
																	 dateID = gsub("-","",dateID, fixed = T))

# Mean by site/month
lai_site_monthly <- lai_dat %>% group_by(siteID, dateID) %>% 
	summarize(mean_lai = mean(data, na.rm=T)) %>% 
	mutate(missing = ifelse(is.na(mean_lai), TRUE, FALSE))

# Mean center and scale data
lai_site_monthly$LAI <- scale(lai_site_monthly$mean_lai)[,1]
lai_site_monthly$LAI_orig <- lai_site_monthly$mean_lai

saveRDS(lai_site_monthly, "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/mean_LAI_data.rds")


ggplot(lai_site_monthly) + geom_point(aes(x = dateID, y = mean_lai)) + facet_grid(~siteID)
