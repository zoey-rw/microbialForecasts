## NOTE: removed crib_fun step because my rowSums are ending up >1 (need to fix probably)
#
# k <- 1
# j <- 1
#
# Read in covariate data
# chem_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/soilChemPlot.rds")
# min.prev = 3;
# max.date = "20200101"
# dom_soil_horizons <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/dominantHorizonsSite.rds")
# predictor_data <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/all_predictor_data.rds")


#' @title 			prepFunctionalData
#' @description prepFunctionalData
#' @export
prepFunctionalData <- function(rank.df,
															 min.date = "20130601",
															 max.date = "20170101",
															 predictor_data = NULL,
												min.prev = 5,
												dom_soil_horizons = NULL,
												full_timeseries = F
){
	require(padr)
	require(lubridate)

	if (is.null(predictor_data)){
		predictor_data <- readRDS(here("data","/clean/all_predictor_data.rds"))
	}
	if (is.null(dom_soil_horizons)){
		dom_soil_horizons <- readRDS(here("data","/clean/dominantHorizonsSite.rds"))
	}


	# For testing
	rank.df.orig <- rank.df
	dat <- rank.df.orig
	#Remove sites missing key covariates
	dat <- dat[which(!dat$siteID %in% c("ABBY","LAJA")),]


	# Reorder & remove missing rows
	dat <- dat %>% arrange(siteID, plotID, dateID) %>%
		mutate(plot_date = paste0(plotID, "_", dateID))
	dat <- dat[!is.na(dat$plot_date),]

	# Subset to filtered dates
	min.date <- as.Date(min.date, format = "%Y%m%d"); message("Min: ", (min.date))
	max.date <- as.Date(max.date, format = "%Y%m%d"); message("Max: ", (max.date))
	dates <- dat$dates
	dat <- dat[which(dates <= max.date & dates >= min.date),]
	dat$dates <- NULL

	# Remove horizons that aren't the dominant one for that site
	keep_hor <- paste0(dom_soil_horizons$siteID, dom_soil_horizons$horizon)
	dat <- dat %>%
		mutate(horizon = ifelse(grepl("-M-", sampleID), "M", "O"),
					 site_hor = paste0(siteID, horizon)) %>%
		dplyr::filter(site_hor %in% keep_hor) %>%
		select(-c(horizon, site_hor))


	# Add coreIDs
	with_coreIDs <- dat %>% tibble::rownames_to_column() %>%
		mutate(number = 1) %>% group_by(plot_date) %>%
		dplyr::mutate(core = cumsum(number)) %>% ungroup() %>%
		mutate(rowname=NULL, number = NULL) %>%
		as.data.frame()

	# Count of dates per plot, only keep plots with a minimum number of time points
	dates_per_plot <- with_coreIDs[!duplicated(with_coreIDs$plot_date),]$plotID %>% table()
	keep_plots <- dates_per_plot[dates_per_plot >= min.prev] %>% names()
	keep_sites <- keep_plots %>% substr(1, 4) %>% unique()


	if (length(keep_plots) == 0){
		errorCondition("No plots meet the minimum-date requirements. Lower 'min.prev' or include more data.")
	} else {
		dat_subset <- with_coreIDs[which(with_coreIDs$plotID %in% keep_plots),]
	}


	# Expand data frame to include all possible plot-date combinations
	poss_dateID <- seq.Date(min.date, max.date, by = "month") %>%
		substr(1, 7) %>% str_replace_all("-", "") %>%
		as.character() %>% as.numeric()
	all_poss_date_combos <- tidyr::expand(dat_subset,
																				nesting(siteID, plotID, core),
																				poss_dateID)  %>% rename(dateID = poss_dateID) %>%
		filter(core==1) %>% distinct() %>% mutate(plot_date = paste0(plotID, "_", dateID))
	# Merge back with actual df
	expanded_dat <- merge(dat_subset, all_poss_date_combos, all = T) %>% arrange(siteID, plotID, dateID)

	# Assign start dates for each site, and indices for looping through
	not_na <- dat_subset %>% filter(!is.na(sampleID))
	site_start <- not_na %>% group_by(siteID) %>% summarise(dateID = min(dateID)) %>%
		mutate(start = match(dateID, poss_dateID)) %>% select(siteID, start) %>% with(., split(start, siteID)) %>% unlist()
	plot_start <- not_na %>% group_by(plotID) %>% summarise(dateID = min(dateID)) %>%
		mutate(start = match(dateID, poss_dateID)) %>% select(plotID, start) %>% with(., split(start, plotID)) %>% unlist()
	plot_index <- plot_start + 1




	# Interval transform for model
	y <- dat_subset %>%
		select(-c(core, siteID, plotID,dateID, sampleID, plot_date)) %>%
		as.matrix() %>% interval_transform()


	#observations, core_plot, plot_site factors for model indexing.
	siteID = dat_subset$siteID
	plotID = dat_subset$plotID
	plot_num			<- match(plotID, names(plot_start))
	plot_site 		<- substr(names(plot_start), 1, 4)
	plot_site_num <- match(plot_site, names(site_start))


	# Create output timepoints
	expanded_dat$timepoint <- as.numeric(as.factor(expanded_dat$dateID))
	timepoint <- expanded_dat[match(dat_subset$dateID, expanded_dat$dateID),]$timepoint
	names(timepoint) <- expanded_dat[match(dat_subset$dateID, expanded_dat$dateID),]$dateID

	# Only return entire (mostly empty) timeseries if using for forecasting
	if (full_timeseries){
		timepoint <- expanded_dat$timepoint
		names(timepoint) <- expanded_dat$dateID
		max.predictor.date = "2022-01-01"
	} else {
		timepoint <-
			expanded_dat[match(dat_subset$dateID, expanded_dat$dateID),]$timepoint
		names(timepoint) <-
			expanded_dat[match(dat_subset$dateID, expanded_dat$dateID),]$dateID
		max.predictor.date = NULL
	}

	# # Create truth outputs
	truth.plot <- expanded_dat %>% group_by(plot_date) %>%
		select(-c(core)) %>%
		summarize(across(where(is.numeric), ~ mean(.x, na.rm=T))) %>% ungroup() %>%
		as.matrix()
	# reorganize truth data
	truth.plot.long <- truth.plot %>% as.data.frame() %>%
		separate(plot_date, sep="_", into=c("siteID","plotID","dateID")) %>%
		mutate(plotID = paste0(siteID, "_", plotID),
					 date_num = as.numeric(as.factor(dateID)),
					 plot_num = match(plotID, names(plot_start)),
					 site_num = match(siteID, names(site_start)),
					 timepoint = as.numeric(timepoint)) %>%
		relocate(plot_num, date_num, site_num, timepoint, .before=1) %>%
		pivot_longer(cols = 8:last_col(),names_to = "species", values_to = "truth")



	# Create output list with everything so far
	out.list <- list(y = y,
									 siteID = siteID,
									 plotID = plotID,
									 plot_site = plot_site,
									 site_start = site_start,
									 plot_start = plot_start,
									 plot_index = plot_index,
									 plot_num = plot_num,
									 plot_site_num = plot_site_num,
									 truth.plot.long = truth.plot.long,
									 N.date = length(poss_dateID),
									 timepoint = timepoint,
									 dates_per_plot = dates_per_plot,
									 N.plot =  length(unique(plotID)),
									 N.spp = ncol(y),
									 N.core = nrow(y),
									 N.site = length(unique(siteID))
	)

	# subset covariates to plots/sites that have been observed for multiple (min.prev) dates, and before the max date
	filt_predictor_data <- lapply(predictor_data, filter_date_site, keep_sites = keep_sites,
																keep_plots = keep_plots, min.date = min.date,
																max.date = max.date, max.predictor.date)
	names(filt_predictor_data) <- recode(names(filt_predictor_data), relEM_plot = "relEM")

	# Add sine/cosine
	sin_cos_month <- get_sin_cos(colnames(filt_predictor_data$mois))
	filt_predictor_data$sin_mo = sin_cos_month$sin
	filt_predictor_data$cos_mo = sin_cos_month$cos

	out.list <- c(out.list, filt_predictor_data)

	return(out.list)
}

