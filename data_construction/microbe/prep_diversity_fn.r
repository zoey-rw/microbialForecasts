## NOTE: removed crib_fun step because my rowSums are ending up >1 (need to fix probably)

# k <- 1
# j <- 1
# 
# # Read in covariate data chem_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/soilChemPlot.rds")
moisture_in <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/monthly_soil_moisture.rds")
temp_in <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/monthly_soil_temperature.rds")
plant_in <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/annual_plant_data.rds")
dom_soil_horizons_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/dominantHorizonsSite.rds")
chem_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/soilChemPlot.rds")
# subset_by_prevalence = T
# min.prev = 4
# fcast_timesteps = 10
moisture = moisture_in;
temp = temp_in;
chem = chem_in; 
plant = plant_in;
#													subset_by_prevalence = TRUE, 
min.prev = 5;
fcast_timesteps = 5;
dom_soil_horizons = dom_soil_horizons_in

prepDivData <- function(rank.df, 
													moisture = moisture_in,
													temp = temp_in,
													chem = chem_in, 
													plant = plant_in,
													#													subset_by_prevalence = TRUE, 
													min.prev = 5,
													fcast_timesteps = 5,
													dom_soil_horizons = dom_soil_horizons_in
){
	require(padr)
	require(tibble)
	require(purrr)
	require(phyloseq)
	require(dplyr)	
	require(tidyr)	
	source("/projectnb2/talbot-lab-data/zrwerbin/NEFI_microbe/NEFI_functions/crib_fun.r")
	
	
	rank.df.orig <- rank.df
	dates <- as.Date(rank.df$dates, format = "%Y%m%d")
	rank.df.orig$dates <- NULL
	dat <- rank.df.orig
	
	
	
	
	# Remove horizons that aren't the dominant one for that site
	keep_hor <- paste0(dom_soil_horizons_in$siteID, dom_soil_horizons_in$horizon)
	dat <- dat %>% 
		mutate(horizon = ifelse(grepl("-M-", sampleID), "M", "O")) %>%  
		mutate(site_hor = paste0(siteID, horizon)) %>% 
		dplyr::filter(site_hor %in% keep_hor) %>% 
		select(-c(horizon, site_hor))
	
	
	# Add coreIDs
	with_coreIDs <- dat %>% tibble::rownames_to_column() %>%  arrange(dateID,siteID,plotID) %>% 
		mutate(number = 1) %>% group_by(plot_date) %>% 
		dplyr::mutate(core = cumsum(number)) %>% ungroup() %>% as.data.frame()
	
	# remove cores 4/5 for any plots
	with_coreIDs <- with_coreIDs[with_coreIDs$core < 4,]
	
	# only keep plots with at least 4 time points
	dates_per_plot <- table(with_coreIDs[!duplicated(with_coreIDs$plot_date),]$plotID)
	keep_plots <- names(dates_per_plot[dates_per_plot > min.prev])
	dat_subset <- with_coreIDs[which(with_coreIDs$plotID %in% keep_plots),]
	
	
	
	start_date <- paste0(substr(min(dates, na.rm = T), 1, 7), "-01")
	poss_dates <- seq.Date(as.Date(start_date), as.Date("2019-09-01"), by = "month")
	poss_dateID <- as.numeric(as.character(stringr::str_replace_all(substr(poss_dates, 1, 7), "-", "")))
	
	# Expand data frame to include all possible plot-date combinations
	all_poss_date_combos <- tidyr::expand(dat_subset, nesting(siteID, plotID, core), poss_dateID) 
	all_poss_date_combos$dateID <- all_poss_date_combos$poss_dateID
	#all_poss_date_combos <- all_poss_date_combos[all_poss_date_combos$core==1,] # actually we only need 1 core per plot.
	# Merge back with actual df
	expanded_dat <- merge(dat_subset, all_poss_date_combos, all = T)
	expanded_dat$plot_date <- paste0(expanded_dat$plotID, "_", expanded_dat$dateID)
	
	
	
	# assign start-indexes for each site, with the earliest (1) being 2013-06
	site_start_index <- list()
	for (s in unique(expanded_dat$siteID)){
		s_dat <- expanded_dat %>% filter(siteID == s & !is.na(expanded_dat[,c(6)]))
		s_init <- min(s_dat$dateID)
		s_index <- which(poss_dateID == s_init)
		site_start_index[[s]] = s_index 
	}
	site_start_index <- stack(site_start_index) %>% 
		select(siteID = ind, start = values) %>% 
		mutate(index = start + 1)
	
	
	# assign start-indexes for each site, with the earliest (1) being 2013-06
	plot_start_index <- list()
	for (p in unique(expanded_dat$plotID)){
		p_dat <- expanded_dat %>% filter(plotID == p & !is.na(expanded_dat[,c(6)]))
		p_init <- min(p_dat$dateID)
		p_index <- which(poss_dateID == p_init)
		plot_start_index[[p]] = p_index 
	}
	plot_start_index <- stack(plot_start_index) %>% 
		select(plotID = ind, start = values) %>% 
		mutate(index = start + 1)
	
	
	
	
	# reformat moisture to be plot x date
	mois_plot <- moisture %>% select(siteID, month, moisture) %>% 
		mutate(moisture = scale(moisture)) %>% 
		merge(unique(dat[,c("siteID", "plotID")])) %>% 
		arrange(month, plotID) %>% unique() %>% group_by(siteID) %>% 
		tidyr::fill(moisture, .direction = "updown") %>% ungroup() %>% select(-siteID) %>% 
		pivot_wider(names_from = month, values_from = moisture, values_fn = list, values_fill = NA) %>% 
		arrange(plotID) %>% 
		column_to_rownames(var = "plotID")
	mois_plot[mois_plot=="NULL"] <- NA
	mois_plot <- mois_plot %>% data.matrix() 
	
	
	mois_plot_uncert <- moisture %>% select(siteID, month, moisture_sd) %>% 
		merge(unique(dat[,c("siteID", "plotID")])) %>% 
		arrange(month, plotID) %>% unique() %>% group_by(siteID) %>% 
		tidyr::fill(moisture_sd, .direction = "updown") %>% ungroup() %>% select(-siteID) %>% 
		pivot_wider(names_from = month, values_from = moisture_sd, values_fn = list) %>% arrange(plotID) %>% 
		column_to_rownames(var = "plotID") 
	mois_plot_uncert[mois_plot_uncert=="NULL"] <- NA
	mois_plot_uncert <- mois_plot_uncert %>% data.matrix()
	# convert NULLs to NAs
	
	# reformat temperature to be plot x date
	temp_plot <- temp %>% select(siteID, month, temperature) %>% merge(unique(dat[,c("siteID", "plotID")])) %>% 
		mutate(temperature = scale(temperature)) %>% 
		arrange(month, plotID) %>% unique() %>% group_by(siteID) %>% 
		tidyr::fill(temperature, .direction = "updown") %>% ungroup() %>% select(-siteID) %>% 
		pivot_wider(names_from = month, values_from = temperature, values_fn = list)  %>% arrange(plotID) %>% 
		column_to_rownames(var = "plotID") 
	temp_plot[temp_plot=="NULL"] <- NA
	temp_plot <- temp_plot %>% data.matrix() 
	
	temp_plot_uncert <- temp %>% select(siteID, month, temperature_sd) %>% merge(unique(dat[,c("siteID", "plotID")])) %>% 
		arrange(month, plotID) %>% unique() %>% group_by(siteID) %>% 
		tidyr::fill(temperature_sd, .direction = "updown") %>% ungroup() %>% select(-siteID) %>% 
		pivot_wider(names_from = month, values_from = temperature_sd, values_fn = list)  %>% arrange(plotID) %>% 
		column_to_rownames(var = "plotID") 
	temp_plot_uncert[temp_plot_uncert=="NULL"] <- NA
	temp_plot_uncert <- temp_plot_uncert %>% data.matrix() 
	
	
	# reformat chemistry to be plot x date (even though dates aren't real) - use dates from temperature
	pH_plot <- chem  %>% merge(unique(dat[,c("siteID", "plotID")]), all.y=T) %>% #filter(plotID %in% dat$plotID) %>% 
		mutate(pH = scale(as.numeric(pH))) %>% 
		select(plotID, siteID, pH) %>% merge(unique(temp[,c("siteID", "date")]))  %>% 
		arrange(date, plotID) %>% unique() %>% 
		group_by(siteID) %>% tidyr::fill(pH, .direction = "updown") %>% ungroup() %>% 
		select(-siteID) %>% 
		pivot_wider(names_from = date, values_from = pH, values_fn = list)  %>% 
		arrange(plotID)  %>% 
		column_to_rownames(var = "plotID") 
	pH_plot[pH_plot=="NULL"] <- NA
	pH_plot <- pH_plot %>% data.matrix() 
	
	
	# reformat chemistry to be plot x date (even though dates aren't real) - use dates from temperature
	pH_plot_sd <- chem  %>% merge(unique(dat[,c("siteID", "plotID")]), all.y=T) %>% #filter(plotID %in% dat$plotID) %>% 
		mutate(pH = scale(as.numeric(pH_sd))) %>% 
		select(plotID, siteID, pH_sd) %>% merge(unique(temp[,c("siteID", "date")]))  %>% 
		arrange(date, plotID) %>% unique() %>% 
		group_by(siteID) %>% tidyr::fill(pH_sd, .direction = "updown") %>% ungroup() %>% 
		select(-siteID) %>% 
		pivot_wider(names_from = date, values_from = pH_sd, values_fn = list)  %>% 
		arrange(plotID)  %>% 
		column_to_rownames(var = "plotID") 
	pH_plot_sd[pH_plot_sd=="NULL"] <- NA
	pH_plot_sd <- pH_plot_sd %>% data.matrix() 
	
	pC_plot <- chem %>% merge(unique(dat[,c("siteID", "plotID")]), all.y=T) %>% 
		mutate(pC = scale(as.numeric(pC))) %>% 
		select(plotID, siteID, pC) %>% merge(unique(temp[,c("siteID", "date")]))  %>% 
		arrange(date, plotID) %>% unique() %>% 
		group_by(siteID) %>% tidyr::fill(pC, .direction = "updown") %>% ungroup() %>% select(-siteID) %>% 
		pivot_wider(names_from = date, values_from = pC, values_fn = list)  %>% arrange(plotID)  %>% 
		column_to_rownames(var = "plotID")
	pC_plot[pC_plot=="NULL"] <- NA
	pC_plot <- pC_plot %>% data.matrix() 
	
	# reformat plant data to be plot x date (even though dates are only years)
	nspp_plot <- plant %>% merge(unique(dat[,c("siteID", "plotID")]), all.y=T) %>% 
		mutate(nspp_total = scale(nspp_total)) %>% 
		select(plotID, nspp_total, year_date) %>% group_by(plotID) %>% 
		pad(interval = "month", start_val = as.Date("2013-06-01"), end_val = as.Date("2020-09-01")) %>% 
		#  	tidyr::fill(nspp_total) %>% 
		group_by(plotID) %>% tidyr::fill(nspp_total, .direction = "updown") %>% ungroup() %>% #select(-siteID) %>% 
		mutate(dateID = as.numeric(as.character(stringr::str_replace_all(substr(year_date, 1, 7), "-", "")))) %>% 
		select(-c(year_date)) %>% arrange(dateID) %>% 
		pivot_wider(names_from = dateID, values_from = nspp_total, values_fn = mean) %>% arrange(plotID)  %>% 
		column_to_rownames(var = "plotID") 
	nspp_plot[nspp_plot=="NULL"] <- NA
	nspp_plot <- nspp_plot %>% data.matrix()
	
	rc_grass_plot <- plant %>% merge(unique(dat[,c("siteID", "plotID")]), all.y=T) %>% 
		mutate(rc_Poaceae = scale(rc_Poaceae)) %>% 
		select(plotID, rc_Poaceae, year_date) %>% group_by(plotID) %>% 
		pad(interval = "month", start_val = as.Date("2013-06-01"), end_val = as.Date("2020-09-01")) %>% 
		group_by(plotID) %>% tidyr::fill(rc_Poaceae, .direction = "updown") %>% ungroup() %>% #select(-siteID) %>% 
		mutate(dateID = as.numeric(as.character(stringr::str_replace_all(substr(year_date, 1, 7), "-", "")))) %>% 
		select(-c(year_date)) %>% arrange(dateID) %>% 
		pivot_wider(names_from = dateID, values_from = rc_Poaceae, values_fn = mean) %>% arrange(plotID)  %>% 
		column_to_rownames(var = "plotID") 
	rc_grass_plot[rc_grass_plot=="NULL"] <- NA
	rc_grass_plot <- rc_grass_plot  %>% data.matrix() 
	
	
	rc_exotic_plot <- plant %>% merge(unique(dat[,c("siteID", "plotID")]), all.y=T) %>% 
		mutate(rel_cover_exotic = scale(rel_cover_exotic)) %>% 
		select(plotID, rel_cover_exotic, year_date) %>% group_by(plotID) %>% 
		pad(interval = "month", start_val = as.Date("2013-06-01"), end_val = as.Date("2020-01-01")) %>% 
		group_by(plotID) %>% tidyr::fill(rel_cover_exotic, .direction = "updown") %>% ungroup() %>% #select(-siteID) %>% 
		mutate(dateID = as.numeric(as.character(stringr::str_replace_all(substr(year_date, 1, 7), "-", "")))) %>% 
		select(-c(year_date)) %>% arrange(dateID) %>% 
		pivot_wider(names_from = dateID, values_from = rel_cover_exotic, values_fn = mean) %>% arrange(plotID)  %>% 
		column_to_rownames(var = "plotID") %>% as.matrix() 
	rc_exotic_plot[rc_exotic_plot=="NULL"] <- NA
	rc_exotic_plot <- rc_exotic_plot  %>% data.matrix() 
	
	
	
	
	
	
	
	
	# subset to plots that have been observed for multiple (min.prev) dates 
	keep_plotIDs <- keep_plots
	
	# now remove those plots from covariates
	mois_plot 				<- mois_plot[rownames(mois_plot) %in% keep_plotIDs,]
	mois_plot_uncert 	<- mois_plot_uncert[rownames(mois_plot_uncert) %in% keep_plotIDs,]
	temp_plot 	<- temp_plot[rownames(temp_plot) %in% keep_plotIDs,]
	temp_plot_uncert 	<- temp_plot_uncert[rownames(temp_plot_uncert) %in% keep_plotIDs,]
	pH_plot 	<- pH_plot[rownames(pH_plot) %in% keep_plotIDs,]
	pH_plot_sd 	<- pH_plot_sd[rownames(pH_plot_sd) %in% keep_plotIDs,]
	pC_plot 	<- pC_plot[rownames(pC_plot) %in% keep_plotIDs,]
	nspp_plot 	<- nspp_plot[rownames(nspp_plot) %in% keep_plotIDs,]
	rc_grass_plot 	<- rc_grass_plot[rownames(rc_grass_plot) %in% keep_plotIDs,]
	rc_exotic_plot 	<- rc_exotic_plot[rownames(rc_exotic_plot) %in% keep_plotIDs,]
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	# reformat abundances into separate matrices
	# also applying cribari-neto transformation to prevent any hard zeroes
	date_dats <- expanded_dat %>% #tibble::rownames_to_column() %>%  
		arrange(dateID,siteID,plotID) %>% 
		mutate(number = 1) %>% group_by(plot_date) %>% 
		dplyr::mutate(core = cumsum(number)) %>% ungroup() %>%  group_split(dateID) %>% 
		map(.f = ~ .x %>% mutate(coreID = paste0(plotID, "-", core)))
	y_list <- date_dats %>% map(.f = ~ .x %>% #dplyr::select(colnames(out_top10)) %>% 
																dplyr::select(!c(siteID, plotID, dateID, core, rowname, sample, asDate, site_date,without_horizon,geneticSampleID,
																								 sampleID, plot_date, number, coreID,poss_dateID)) %>% 
																as.matrix() #%>% interval_transform()
	) 
	names(y_list) <- sort(unique(expanded_dat$dateID))
	y <- array(unlist(y_list), 
						 dim=c(nrow(y_list[[1]]), ncol(y_list[[1]]), length(y_list)),
						 dimnames = list(date_dats[[1]]$plotID, colnames(y_list[[1]]),  names(y_list)))
	
	## CREATE NA VALUES FOR FORECASTING
	# if(!is.null(fcast_timesteps)){
	# 	names <- substr(seq.Date(max(dates, na.rm=T), length.out = fcast_timesteps+1, by = "month"), 1, 7)
	# 	names <- gsub("-","",names)
	# 	names <- names[-1]
	# 	fcast_NA <- matrix(nrow = 3, ncol = fcast_timesteps) 
	# 	colnames(fcast_NA) <- names
	# 	plot_dats <- lapply(plot_dats, function(x) cbind(x, fcast_NA))
	# }
	# 
	#core_plot, plot_site factors for model indexing.
	#expanded_dat$plotID <- factor(paste0(expanded_dat$siteID,'_',expanded_dat$plotID))
	expanded_dat$plotID <- factor(expanded_dat$plotID)
	by_date <- expanded_dat %>% group_split(dateID)
	siteID = as.factor(by_date[[1]]$siteID)
	core_plot = as.factor(by_date[[1]]$plotID)
	plot_site <- as.factor(unique(by_date[[1]][,c('siteID','plotID')])$siteID)
	n.core.by.date <- expanded_dat %>%  group_split(dateID) %>% 
		map(.f = ~ .x %>% nrow())
	
	# Create truth outputs
	truth.site <- expanded_dat %>% group_by(siteID) %>% select(-c(core,number,poss_dateID)) %>% 
		summarize(.cols = across(where(is.numeric), ~ mean(.x, na.rm=T))) %>% 
		select(-1) %>% as.matrix()
	truth.plot <- expanded_dat %>% group_by(plot_date) %>% 
		select(-c(number, rowname, core, poss_dateID)) %>% 
		summarize(across(where(is.numeric), ~ mean(.x, na.rm=T))) %>% ungroup() %>% 
		as.matrix()
	
	# reorganize truth data
	truth.plot.long <- truth.plot %>% as.data.frame() %>% 
		separate(plot_date, sep="_", into=c("siteID","plotID","dateID")) %>% 
		mutate(plotID = paste0(siteID, "_", plotID),
					 plot_num = as.numeric(as.factor(plotID)),
					 date_num = as.numeric(as.factor(dateID))) %>% 
		pivot_longer(cols = 4,names_to = "shannon", values_to = "truth")
	
	
	return(list(y = y, 
							siteID = siteID, 
							plotID = core_plot, 
							plot.truth = truth.plot.long,
							site_start = site_start_index[site_start_index$siteID %in% siteID,],
							plot_start = plot_start_index[plot_start_index$plotID %in% core_plot,],
							n.core.by.date = n.core.by.date,
							mois = mois_plot, 				
							mois_sd = mois_plot_uncert,
							temp = temp_plot,
							temp_sd = temp_plot_uncert,
							pH = pH_plot,
							pH_sd = pH_plot_sd,
							pC = pC_plot,
							nspp = nspp_plot,
							rc_grass = rc_grass_plot,
							rc_exotic = rc_exotic_plot))
}

