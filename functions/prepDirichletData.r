# This function subsets the microbial data by prevalence, then only keeps the corresponding predictor data. 

## NOTE: removed crib_fun step because my rowSums are ending up >1 (need to fix probably)

min.prev = 3;
max.date = "20170101"
dom_soil_horizons <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/dominantHorizonsSite.rds")
predictor_data <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/all_predictor_data.rds")

prepModelData <- function(rank.df, 
													min.prev = 3,
													max.date = "20170101",
													predictor_data = NULL,
													dom_soil_horizons = NULL
){
	require(padr)
	require(tibble)
	require(purrr)
	require(dplyr)	
	require(tidyr)	
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
	
	if (is.null(predictor_data)){
		predictor_data <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/all_predictor_data.rds")
	}
	if (is.null(dom_soil_horizons)){
		dom_soil_horizons <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/dominantHorizonsSite.rds")
	}
	
# 	# for testing
# 	rank.df <- d1$phylum_bac
#   rank.df.orig <- rank.df
#   dat <- rank.df.orig

  dat <- rank.df

  # Subset to dates before "max.date"
  max.date <- as.Date(max.date, format = "%Y%m%d")
  min.date <- as.Date("20130601", format = "%Y%m%d")
  dates <- as.Date(dat$dates, format = "%Y%m%d")
  dat <- dat[which(dates < max.date),]
  dat <- dat[which(dates > min.date),]
  dat$dates <- NULL
  
  
  # Remove horizons that aren't the dominant one for that site
  keep_hor <- paste0(dom_soil_horizons$siteID, dom_soil_horizons$horizon)
  dat <- dat %>% 
  	mutate(horizon = ifelse(grepl("-M-", sampleID), "M", "O")) %>%  
  	mutate(site_hor = paste0(siteID, horizon)) %>% 
  	dplyr::filter(site_hor %in% keep_hor) %>% 
  	select(-c(horizon, site_hor))
  
 
  # Add coreIDs
  with_coreIDs <- dat %>% tibble::rownames_to_column() %>%  #arrange(dateID,siteID,plotID) %>% 
  	mutate(number = 1) %>% group_by(plot_date) %>% 
  	dplyr::mutate(core = cumsum(number)) %>% ungroup() %>% as.data.frame()
  
  # # remove cores 4/5 for any plots
  # with_coreIDs <- with_coreIDs[with_coreIDs$core < 4,]
  
  # Count of dates per plot
  dates_per_plot <- table(with_coreIDs[!duplicated(with_coreIDs$plot_date),]$plotID)
  
  # Only keep plots with a minimum number of time points
  keep_plots <- names(dates_per_plot[dates_per_plot >= min.prev])
  if (length(keep_plots) == 0){
  	errorCondition("No plots meet the minimum-date requirements. Lower 'min.prev' or include more data.")
  	}
  dat_subset <- with_coreIDs[which(with_coreIDs$plotID %in% keep_plots),]

  
  # Fill in any empty dates
  start_date <- paste0(substr(min(dates, na.rm = T), 1, 7), "-01")
  start_date <- "2013-06-01"
  poss_dates <- seq.Date(as.Date(start_date), as.Date(max(dates, na.rm = T)), by = "month")
  poss_dateID <- as.numeric(as.character(stringr::str_replace_all(substr(poss_dates, 1, 7), "-", "")))
  
  # Add dates for forecasts
  #names <- substr(seq.Date(max(poss_dateID, na.rm=T), length.out = fcast_timesteps+1, by = "month"), 1, 7)
  #   	names <- gsub("-","",names)
  #   	names <- names[-1]
  
  # Expand data frame to include all possible plot-date combinations
  all_poss_date_combos <- tidyr::expand(dat_subset, nesting(siteID, plotID, core), poss_dateID) 
  all_poss_date_combos$dateID <- all_poss_date_combos$poss_dateID
  all_poss_date_combos <- all_poss_date_combos[all_poss_date_combos$core==1,] # actually we only need 1 core per plot.
  # Merge back with actual df
  expanded_dat <- merge(dat_subset, all_poss_date_combos, all = T)
  expanded_dat$plot_date <- paste0(expanded_dat$plotID, "_", expanded_dat$dateID)
  expanded_dat <- expanded_dat %>% arrange(siteID, plotID, dateID) %>% filter(plotID %in% keep_plots)
  
  
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
  	p_init <- min(p_dat$dateID, na.rm = T)
  	p_index <- which(poss_dateID == p_init)
  	plot_start_index[[p]] = p_index 
  }
  plot_start_index <- stack(plot_start_index) %>% 
  	select(plotID = ind, start = values) %>% 
  	mutate(index = start + 1)
 
  
  # siteID = as.factor(expanded_dat$siteID)
  # plotID = as.factor(expanded_dat$plotID)
  # plot_site <- as.factor(unique(expanded_dat[,c('siteID','plotID')])$siteID)
  
  y_long <- dat_subset %>% select(-c(rowname, siteID, plotID,dateID,sampleID,plot_date,number,core))
  y_long_trans <- y_long %>% as.matrix() %>% interval_transform()
  y <- y_long_trans
  siteID = dat_subset$siteID
  plotID = dat_subset$plotID
#  plot_site <- as.factor(unique(dat_subset[,c('siteID','plotID')])$siteID)
  
  # Create output timepoints
  expanded_dat$timepoint <- as.numeric(as.factor(expanded_dat$dateID))
  timepoint <- expanded_dat[match(dat_subset$dateID, expanded_dat$dateID),]$timepoint
  names(timepoint) <- expanded_dat[match(dat_subset$dateID, expanded_dat$dateID),]$dateID
  
  # subset covariates to plots/sites that have been observed for multiple (min.prev) dates, and before the max date
  	keep_sites <- unique(substr(keep_plots, 1, 4))
  	# By site
  	mois 				<- predictor_data$mois %>% filter(rownames(predictor_data$mois) %in% keep_sites) %>% data.matrix() 
  	mois_sd 				<- predictor_data$mois_sd %>% filter(rownames(predictor_data$mois_sd) %in% keep_sites) %>% data.matrix() 
  	temp 				<- predictor_data$temp %>% filter(rownames(predictor_data$temp) %in% keep_sites) %>% data.matrix() 
  	temp_sd 				<- predictor_data$temp_sd %>% filter(rownames(predictor_data$temp_sd) %in% keep_sites) %>% data.matrix() 
  	# By plot
  	pH 				<- predictor_data$pH %>% filter(rownames(predictor_data$pH) %in% keep_plots) %>% select(88) %>% data.matrix() 
  	pH_sd 				<- predictor_data$pH_sd %>% filter(rownames(predictor_data$pH_sd) %in% keep_plots) %>%  #data.matrix() %>% 
  		select(88)  %>% data.matrix() 
  	pC 				<- predictor_data$pC %>% filter(rownames(predictor_data$pC) %in% keep_plots) %>% #data.matrix() %>% 
  		select(88) %>% data.matrix() 
  	pC_sd 				<- predictor_data$pC_sd %>% filter(rownames(predictor_data$pC_sd) %in% keep_plots) %>% #data.matrix() %>% 
  		select(88) %>% data.matrix() 
  	nspp 				<- predictor_data$nspp %>% filter(rownames(predictor_data$nspp) %in% keep_plots)  %>% data.matrix() 
  	rc_grass 				<- predictor_data$rc_grass %>% filter(rownames(predictor_data$rc_grass) %in% keep_plots) %>% data.matrix() 
  	rc_exotic 				<- predictor_data$rc_exotic %>% filter(rownames(predictor_data$rc_exotic) %in% keep_plots) %>% data.matrix() 
  	
  	site_start_temp <- site_start_index[site_start_index$siteID %in% keep_sites,]
  	plot_start_temp <- plot_start_index[plot_start_index$plotID %in% keep_plots,]
  	
  	site_start <- site_start_temp$start
  	plot_start <- plot_start_temp$start 
  	plot_index <- plot_start_temp$index 
  	names(plot_index) <- plot_start_temp$plotID
  	names(plot_start) <- plot_start_temp$plotID
  	names(site_start) <- site_start_temp$siteID
  	#identical(rownames(pH), unique(rownames(y[,,1])))
  	
  	plot_site <- substr(names(plot_start), 1, 4)
  	
  	#core_plot, plot_site factors for model indexing.
  	plot_num <- match(plotID, names(plot_start))
  	plot_site_num <- match(plot_site, names(site_start))
  	#identical(rownames(pH), unique(rownames(y[,,1])))
  	
  	
  	# Create truth outputs
  	# truth.site <- expanded_dat %>% group_by(siteID) %>% select(-c(core,number,poss_dateID)) %>%
  	# 	summarize(.cols = across(where(is.numeric), ~ mean(.x, na.rm=T))) %>%
  	# 	select(-1) %>% as.matrix()
  	truth.plot <- expanded_dat %>% group_by(plot_date) %>%
  		select(-c(number, rowname, core, poss_dateID)) %>%
  		summarize(across(where(is.numeric), ~ mean(.x, na.rm=T))) %>% ungroup() %>%
  		as.matrix()
  	
  	# reorganize truth data
  	truth.plot.long <- truth.plot %>% as.data.frame() %>%
  		separate(plot_date, sep="_", into=c("siteID","plotID","dateID")) %>%
  		mutate(plotID = paste0(siteID, "_", plotID),
  					 date_num = as.numeric(as.factor(dateID)),
  					 plot_num = match(plotID, names(plot_start)),
  					 site_num = match(siteID, names(site_start))) %>%
  		relocate(plot_num, date_num, site_num, timepoint, .before=1) %>%
  		pivot_longer(cols = 8:last_col(),names_to = "species", values_to = "truth")
  	
  
  	return(list(y = y, 
  							siteID = siteID, 
  							plotID = plotID, 
  							plot_site = plot_site,
  							site_start = site_start,
  							plot_start = plot_start,
  							plot_index = plot_index,
  							plot_num = plot_num,
  							truth.plot.long = truth.plot.long,
  							plot_site_num = plot_site_num,
  							N.date = length(poss_dates),
  							timepoint = timepoint,
  							mois = mois, 				
  							mois_sd = mois_sd,
  							temp = temp,
  							temp_sd = temp_sd,
  							pH = pH,
  							pH_sd = pH_sd,
  							pC = pC,
  							pC_sd = pC_sd,
  							nspp = nspp,
  							rc_grass = rc_grass,
  							rc_exotic = rc_exotic))
}

  
  
# 
#   # reformat abundances into separate matrices
#   # also applying cribari-neto transformation to prevent any hard zeroes
#   date_dats <- expanded_dat %>% #tibble::rownames_to_column() %>%
#   	arrange(dateID,siteID,plotID) %>%
#   	mutate(number = 1) %>% group_by(plot_date) %>%
#   	dplyr::mutate(core = cumsum(number)) %>% ungroup() %>%  group_split(dateID) %>%
#   	map(.f = ~ .x %>% mutate(coreID = paste0(plotID, "-", core)))
#   y_list <- date_dats %>% map(.f = ~ .x %>% #dplyr::select(colnames(out_top10)) %>%
#   															dplyr::select(!c(siteID, plotID, dateID, core, rowname,
#   																							 sampleID, plot_date, number, coreID,poss_dateID)) %>%
#   															as.matrix() %>% interval_transform()
#   														)
#   names(y_list) <- sort(unique(expanded_dat$dateID))
#   y <- array(unlist(y_list),
#   					 dim=c(nrow(y_list[[1]]), ncol(y_list[[1]]), length(y_list)),
#   					 dimnames = list(date_dats[[1]]$plotID, colnames(y_list[[1]]),  names(y_list)))
  
  ## CREATE NA VALUES FOR FORECASTING
  # if(!is.null(fcast_timesteps)){
#   	names <- substr(seq.Date(max(dates, na.rm=T), length.out = fcast_timesteps+1, by = "month"), 1, 7)
#   	names <- gsub("-","",names)
#   	names <- names[-1]
# fcast_NA <- array(dim=c(nrow(y_list[[1]]), ncol(y_list[[1]]), fcast_timesteps))
#   	colnames(fcast_NA) <- names
#   	plot_dats <- lapply(plot_dats, function(x) cbind(x, fcast_NA))
  # }
  # 
  # #core_plot, plot_site factors for model indexing.
  # expanded_dat$plotID <- factor(expanded_dat$plotID)
  # by_date <- expanded_dat %>% group_split(dateID)
  # siteID = as.factor(by_date[[1]]$siteID)
  # core_plot = as.factor(by_date[[1]]$plotID)
  # plot_site <- as.factor(unique(by_date[[1]][,c('siteID','plotID')])$siteID)
  # n.core.by.date <- expanded_dat %>%  group_split(dateID) %>% 
  # 	map(.f = ~ .x %>% nrow())
  # 

  
#  
#   return(list(y = y, 
#   						siteID = siteID, 
#   						plotID = core_plot, 
#   						plot.truth = truth.plot.long,
#   						site_start = site_start_index[site_start_index$siteID %in% siteID,],
#   						plot_start = plot_start_index[plot_start_index$plotID %in% core_plot,],
#   						n.core.by.date = n.core.by.date,
#   						mois = mois_plot, 				
#   						mois_sd = mois_plot_uncert,
#   						temp = temp_plot,
#   						temp_sd = temp_plot_uncert,
#   						pH = pH_plot,
#   						pH_sd = pH_plot_sd,
#   						pC = pC_plot,
#   						nspp = nspp_plot,
#   						rc_grass = rc_grass_plot,
#   						rc_exotic = rc_exotic_plot))
# }

