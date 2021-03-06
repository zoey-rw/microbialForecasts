# Prep predictor data for models. A subsequent function will refine this dataset to only include specific plots/sites.
require(padr)
require(tibble)
require(purrr)
require(dplyr)	
require(tidyr)	
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

# # Read in covariate data chem_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/soilChemPlot.rds")
moisture_in <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/monthly_soil_moisture.rds")
temp_in <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/monthly_soil_temperature.rds")
plant_in <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/annual_plant_data.rds")
dom_soil_horizons_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/dominantHorizonsSite.rds")
chem_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/soilChemPlot.rds")

moisture = moisture_in;
temp = temp_in;
chem = chem_in; 
plant = plant_in;

dom_soil_horizons = dom_soil_horizons_in



# reformat moisture to be site x date
mois_site <- moisture %>% dplyr::select(siteID, month, moisture) %>% 
	arrange(month) %>% unique() %>% group_by(siteID) %>% 
	tidyr::fill(moisture, .direction = "updown") %>% ungroup() %>% 
	pivot_wider(names_from = month, values_from = moisture,
							values_fill = NA) %>% 
	arrange(siteID) %>% 
	 column_to_rownames(var = "siteID")
# mois_site <- mois_site %>% data.matrix() 

# reformat moisture uncertainty to be site x date
mois_site_sd <- moisture %>% dplyr::select(siteID, month, moisture_sd) %>% 
	arrange(month) %>% unique() %>% group_by(siteID) %>% 
	tidyr::fill(moisture_sd, .direction = "updown") %>% ungroup() %>% 
	pivot_wider(names_from = month, values_from = moisture_sd,
							values_fill = NA) %>% 
	arrange(siteID) %>% 
	column_to_rownames(var = "siteID")
#mois_site_sd <- mois_site_sd %>% data.matrix() 

# reformat temperature to be site x date
temp_site <- temp %>% dplyr::select(siteID, month, temperature) %>% 
	arrange(month) %>% unique() %>% group_by(siteID) %>% 
	tidyr::fill(temperature, .direction = "updown") %>% ungroup() %>% 
	pivot_wider(names_from = month, values_from = temperature,
							values_fill = NA) %>% 
	arrange(siteID) %>% 
	column_to_rownames(var = "siteID")
# temp_site <- temp_site %>% data.matrix() 

# reformat temperature uncertainty to be site x date
temp_site_sd <- temp %>% dplyr::select(siteID, month, temperature_sd) %>% 
	arrange(month) %>% unique() %>% group_by(siteID) %>% 
	tidyr::fill(temperature_sd, .direction = "updown") %>% ungroup() %>% 
	pivot_wider(names_from = month, values_from = temperature_sd,
							values_fill = NA) %>% 
	arrange(siteID) %>% 
	column_to_rownames(var = "siteID")
# temp_site_sd <- temp_site_sd %>% data.matrix() 



# reformat pH to be plot x date (even though dates aren't real) - use dates from temperature
pH_plot <- chem  %>% #merge(unique(dat[,c("siteID", "plotID")]), all.y=T) %>% #filter(plotID %in% dat$plotID) %>% 
	dplyr::select(plotID, siteID, pH) %>% merge(unique(temp[,c("siteID", "date","month")]))  %>% 
	arrange(date, plotID) %>% unique() %>% 
	group_by(siteID) %>% tidyr::fill(pH, .direction = "updown") %>% ungroup() %>% 
	dplyr::select(-siteID, -date) %>% 
	pivot_wider(names_from = month, values_from = pH, values_fn = list)  %>% 
	arrange(plotID)  %>% 
	column_to_rownames(var = "plotID") 
pH_plot[pH_plot=="NULL"] <- NA
# pH_plot <- pH_plot %>% data.matrix() 

# reformat pH uncertainty to be plot x date (even though dates aren't real) - use dates from temperature
pH_plot_sd <- chem  %>% #merge(unique(dat[,c("siteID", "plotID")]), all.y=T) %>% #filter(plotID %in% dat$plotID) %>% 
	#mutate(pH = scale(as.numeric(pH_sd))) %>% 
	dplyr::select(plotID, siteID, pH_sd) %>% merge(unique(temp[,c("siteID", "date","month")]))  %>% 
	arrange(date, plotID) %>% unique() %>% 
	group_by(siteID) %>% tidyr::fill(pH_sd, .direction = "updown") %>% ungroup() %>% 
	dplyr::select(-siteID, -date) %>% 
	pivot_wider(names_from = month, values_from = pH_sd, values_fn = list)  %>% 
	arrange(plotID)  %>% 
	column_to_rownames(var = "plotID") 
pH_plot_sd[pH_plot_sd=="NULL"] <- NA
# pH_plot_sd <- pH_plot_sd %>% data.matrix() 


# reformat pC to be plot x date (even though dates aren't real) - use dates from temperature
pC_plot <- chem  %>% #merge(unique(dat[,c("siteID", "plotID")]), all.y=T) %>% #filter(plotID %in% dat$plotID) %>% 
	dplyr::select(plotID, siteID, pC) %>% merge(unique(temp[,c("siteID", "date","month")]))  %>% 
	arrange(date, plotID) %>% unique() %>% 
	group_by(siteID) %>% tidyr::fill(pC, .direction = "updown") %>% ungroup() %>% 
	dplyr::select(-siteID, -date) %>% 
	pivot_wider(names_from = month, values_from = pC, values_fn = list)  %>% 
	arrange(plotID)  %>% 
	column_to_rownames(var = "plotID") 
pC_plot[pC_plot=="NULL"] <- NA
# pC_plot <- pC_plot %>% data.matrix() 

# reformat pC uncertainty to be plot x date (even though dates aren't real) - use dates from temperature
pC_plot_sd <- chem  %>% #merge(unique(dat[,c("siteID", "plotID")]), all.y=T) %>% #filter(plotID %in% dat$plotID) %>% 
	#mutate(pH = scale(as.numeric(pH_sd))) %>% 
	dplyr::select(plotID, siteID, pC_sd) %>% merge(unique(temp[,c("siteID", "date","month")]))  %>% 
	arrange(date, plotID) %>% unique() %>% 
	group_by(siteID) %>% tidyr::fill(pC_sd, .direction = "updown") %>% ungroup() %>% 
	dplyr::select(-siteID, -date) %>% 
	pivot_wider(names_from = month, values_from = pC_sd, values_fn = list)  %>% 
	arrange(plotID)  %>% 
	column_to_rownames(var = "plotID") 
pC_plot_sd[pC_plot_sd=="NULL"] <- NA
# pC_plot_sd <- pC_plot_sd %>% data.matrix() 

# reformat plant data to be plot x date (even though dates are only years)
nspp_plot <- plant %>% #merge(unique(dat[,c("siteID", "plotID")]), all.y=T) %>% 
	mutate(nspp_total = scale(nspp_total, scale=T)) %>% 
	dplyr::select(plotID, nspp_total, year_date) %>% group_by(plotID) %>% 
	pad(interval = "month", start_val = as.Date("2013-06-01"), end_val = as.Date("2020-01-01")) %>% 
	#  	tidyr::fill(nspp_total) %>% 
	group_by(plotID) %>% tidyr::fill(nspp_total, .direction = "updown") %>% ungroup() %>% #dplyr::select(-siteID) %>% 
	mutate(dateID = as.numeric(as.character(stringr::str_replace_all(substr(year_date, 1, 7), "-", "")))) %>% 
	dplyr::select(-c(year_date)) %>% arrange(dateID) %>% 
	pivot_wider(names_from = dateID, values_from = nspp_total, values_fn = mean) %>% arrange(plotID)  %>% 
	column_to_rownames(var = "plotID") 
nspp_plot[nspp_plot=="NULL"] <- NA
# nspp_plot <- nspp_plot %>% data.matrix()

rc_grass_plot <- plant %>% #merge(unique(dat[,c("siteID", "plotID")]), all.y=T) %>% 
	mutate(rc_Poaceae = scale(rc_Poaceae, scale=T)) %>% 
	dplyr::select(plotID, rc_Poaceae, year_date) %>% group_by(plotID) %>% 
	pad(interval = "month", start_val = as.Date("2013-06-01"), end_val = as.Date("2020-01-01")) %>% 
	group_by(plotID) %>% tidyr::fill(rc_Poaceae, .direction = "updown") %>% ungroup() %>% #dplyr::select(-siteID) %>% 
	mutate(dateID = as.numeric(as.character(stringr::str_replace_all(substr(year_date, 1, 7), "-", "")))) %>% 
	dplyr::select(-c(year_date)) %>% arrange(dateID) %>% 
	pivot_wider(names_from = dateID, values_from = rc_Poaceae, values_fn = mean) %>% arrange(plotID)  %>% 
	column_to_rownames(var = "plotID") 
rc_grass_plot[rc_grass_plot=="NULL"] <- NA
# rc_grass_plot <- rc_grass_plot  %>% data.matrix() 


rc_exotic_plot <- plant %>%# merge(unique(chem_in[,c("siteID", "plotID")]), all.y=T) %>% 
	mutate(rel_cover_exotic = scale(rel_cover_exotic, scale=T)) %>% 
	dplyr::select(plotID, rel_cover_exotic, year_date) %>% group_by(plotID) %>% 
	pad(interval = "month", start_val = as.Date("2013-06-01"), end_val = as.Date("2020-01-01")) %>% 
	group_by(plotID) %>% tidyr::fill(rel_cover_exotic, .direction = "updown") %>% ungroup() %>% #dplyr::select(-siteID) %>% 
	mutate(dateID = as.numeric(as.character(stringr::str_replace_all(substr(year_date, 1, 7), "-", "")))) %>% 
	dplyr::select(-c(year_date)) %>% arrange(dateID) %>% 
	pivot_wider(names_from = dateID, values_from = rel_cover_exotic, values_fn = mean) %>% arrange(plotID)  %>% 
	column_to_rownames(var = "plotID") 
rc_exotic_plot[rc_exotic_plot=="NULL"] <- NA
# rc_exotic_plot <- rc_exotic_plot  %>% data.matrix() 

all_predictors <- list("mois" = mois_site,
											 "mois_sd" = mois_site_sd,
											 "temp" = temp_site,
											 "temp_sd" = temp_site_sd,
											 "pH" = pH_plot,
											 "pH_sd" = pH_plot_sd,
											 "pC" = pC_plot,
											 "pC_sd" = pC_plot_sd, 
											 "nspp" = nspp_plot,
											 "rc_grass" = rc_grass_plot,
											 "rc_exotic" = rc_exotic_plot
											 )


saveRDS(all_predictors, "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/all_predictor_data.rds")
