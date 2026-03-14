## NOTE: removed crib_fun step because my rowSums are ending up >1 (need to fix probably)
#
# k <- 1
# j <- 1
#
# # Read in covariate data
# chem_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/soilChemPlot.rds")
# dom_soil_horizons <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/dominantHorizonsSite.rds")
# predictor_data <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/all_predictor_data.rds")
# #
# min.prev = 3;
# min.date = "20151101"
# min.date = "20130601"
# max.date = "20180101"
# dom_soil_horizons <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/dominantHorizonsSite.rds")
# predictor_data <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/all_predictor_data.rds")
# full_timeseries = F


#' @title 			prepBetaRegData
#' @description prepBetaRegData
#' @export
prepBetaRegData <- function(rank.df,
															min.date = "20130601",
															max.date = "20170101",
															predictor_data = NULL, #readRDS("./data/clean/all_predictor_data.rds"),
															min.prev = 3,
															dom_soil_horizons = NULL, #readRDS("./data/clean/dominantHorizonsSite.rds"),
															full_timeseries = FALSE,
															keep_vec = NULL
){
	require(padr)
	require(lubridate)

	if (is.null(predictor_data)){
		predictor_data <- readRDS(here("data", "clean/all_predictor_data.rds"))
	}
	if (is.null(dom_soil_horizons)){
		dom_soil_horizons <- readRDS(here("data", "clean/dominantHorizonsSite.rds"))
	}

	# For testing
	rank.df.orig <- rank.df
	dat <- rank.df.orig

	# CRITICAL FIX: Reconstruct missing dateID values from sampleID or dates
	if (!"dateID" %in% colnames(dat)) {
		stop("Required column 'dateID' missing from rank.df input")
	}
	
	# Check for NA dateID and reconstruct if needed
	if (any(is.na(dat$dateID))) {
		na_dateid_count <- sum(is.na(dat$dateID))
		message("Found ", na_dateid_count, " rows with NA dateID - attempting to reconstruct from sampleID")
		
		na_dateid_mask <- is.na(dat$dateID)
		
		# Try to reconstruct from sampleID using parseNEONsampleIDs
		if ("sampleID" %in% colnames(dat)) {
			sampleIDs_to_parse <- dat$sampleID[na_dateid_mask]
			valid_sampleIDs <- sampleIDs_to_parse[!is.na(sampleIDs_to_parse)]
			if (length(valid_sampleIDs) > 0) {
				tryCatch({
					parsed <- parseNEONsampleIDs(valid_sampleIDs)
					if ("dateID" %in% colnames(parsed)) {
						# Match parsed results back to original dat rows by rownames (parseNEONsampleIDs preserves original sampleID in rownames)
						parsed_rownames <- rownames(parsed)
						for (i in which(na_dateid_mask)) {
							sampleID_val <- dat$sampleID[i]
							# Find matching row in parsed (rownames should match original sampleID)
							matched_idx <- which(parsed_rownames == sampleID_val)
							if (length(matched_idx) > 0 && matched_idx[1] <= nrow(parsed)) {
								parsed_dateID <- parsed$dateID[matched_idx[1]]
								if (!is.na(parsed_dateID)) {
									# Convert dateID to numeric format (YYYYMM)
									dat$dateID[i] <- suppressWarnings(as.numeric(as.character(parsed_dateID)))
								}
							}
						}
						reconstructed_count <- sum(!is.na(dat$dateID) & na_dateid_mask)
						if (reconstructed_count > 0) {
							message("Reconstructed ", reconstructed_count, " dateID values from sampleID")
						}
					}
				}, error = function(e) {
					message("Could not reconstruct dateID from sampleID: ", conditionMessage(e))
				})
			}
		}
		
		# If still NA, try to reconstruct from dates column
		if (any(is.na(dat$dateID))) {
			still_na_mask <- is.na(dat$dateID)
			if ("dates" %in% colnames(dat)) {
				dates_vals <- dat$dates[still_na_mask]
				if (inherits(dates_vals, "Date") || (is.character(dates_vals) && any(!is.na(dates_vals)))) {
					tryCatch({
						date_parsed <- as.Date(dates_vals)
						reconstructed_dateID <- as.numeric(format(date_parsed, "%Y%m"))
						valid_reconstruction <- !is.na(reconstructed_dateID)
						if (any(valid_reconstruction)) {
							dat$dateID[still_na_mask][valid_reconstruction] <- reconstructed_dateID[valid_reconstruction]
							message("Reconstructed ", sum(valid_reconstruction), " dateID values from dates column")
						}
					}, error = function(e) {
						message("Could not reconstruct dateID from dates: ", conditionMessage(e))
					})
				}
			}
		}
		
		# After reconstruction attempts, filter out any remaining NA dateID
		if (any(is.na(dat$dateID))) {
			still_na_count <- sum(is.na(dat$dateID))
			warning("Could not reconstruct dateID for ", still_na_count, " rows - filtering out these rows")
			dat <- dat[!is.na(dat$dateID), ]
		}
	}

	if (!is.null(keep_vec)){
		dat <- dat[,colnames(rank.df) %in% keep_vec]
	}
	#	colnames(dat) <- lapply(strsplit(colnames(dat), "\\."), "[[", 1)


	# If more than one species + other, replace species with other column
	if (ncol(dat) > 8) {
		tots <- rowSums(dat[,7:ncol(dat)])
		dat$other <- 1-tots
		# Remove samples where more than 99% of reads are "Other"
		dat <- dat %>% filter(tots > .01)
	}

	#Remove sites missing key covariates
	dat <- dat[which(!dat$siteID %in% c("ABBY","LAJA")),]


	# Reorder & remove missing rows
	# Since we've validated dateID is not NA above, plot_date should be valid
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
																				tidyr::nesting(siteID, plotID, core),
																				poss_dateID)  %>% dplyr::rename(dateID = poss_dateID) %>%
		filter(core==1) %>% distinct() %>% mutate(plot_date = paste0(plotID, "_", dateID))
	# Merge back with actual df
	expanded_dat <- merge(dat_subset, all_poss_date_combos, all = T) %>% arrange(siteID, plotID, dateID)
	
	# CRITICAL FIX: Ensure dates column is properly mapped from dateID
	# Create approximate dates from dateID for all timepoints
	expanded_dat <- expanded_dat %>%
		mutate(dates = paste0(substr(dateID, 1, 4), "-", substr(dateID, 5, 6), "-01"))

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
	expanded_dat$plot_num			<- match(expanded_dat$plotID, names(plot_start))
	
	# CRITICAL FIX: Proper plot-to-site indexing for site-level effects and covariates
	# Each plot belongs to a site, and we need to map plot index p to site index k
	plot_site 		<- substr(names(plot_start), 1, 4)  # Extract site ID from plot ID
	unique_sites <- unique(plot_site)  # Get unique site IDs in order they appear
	site_indices <- 1:length(unique_sites)  # Create sequential site indices
	names(site_indices) <- unique_sites  # Name them with site IDs
	
	# Now map each plot to its site index: plot_site_num[p] gives site index for plot p
	plot_site_num <- site_indices[plot_site]


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
		summarize(across(where(is.numeric), ~ mean(.x, na.rm=T)),
						  dates = as.character(dates[1])) %>% ungroup() # Use first date in each group (all dates in month are equivalent)
	# reorganize truth data
	truth.plot.long <- truth.plot %>%
		separate(plot_date, sep="[-_]", into=c("siteID","plotID","dateID")) %>%
		mutate(plotID = paste0(siteID, "_", plotID),
					 date_num = as.numeric(as.factor(dateID)),
					 plot_num = match(plotID, names(plot_start)),
					 site_num = match(siteID, names(site_start)),
					 timepoint = as.numeric(timepoint)) %>%
		relocate(plot_num, date_num, site_num, timepoint, .before=1) %>%
		pivot_longer(cols = -c(plot_num, date_num, site_num, timepoint, siteID, plotID, dateID, dates), names_to = "species", values_to = "truth")


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
									 N.site = length(unique(siteID)),
									 sample_values = expanded_dat
	)

	# subset covariates to plots/sites that have been observed for multiple (min.prev) dates, and before the max date
	# CRITICAL FIX: When full_timeseries=TRUE, use the full environmental data range instead of taxonomic data range
	if (full_timeseries) {
		# For full time series, use the full environmental data range
		env_min_date <- min(colnames(predictor_data$temp))
		env_max_date <- max(colnames(predictor_data$temp))
		# CRITICAL FIX: Override max.predictor.date to use environmental data range
		env_max_predictor_date <- paste0(substr(env_max_date, 1, 4), "-", substr(env_max_date, 5, 6), "-01")
		filt_predictor_data <- lapply(predictor_data, filter_date_site, keep_sites = keep_sites,
																	keep_plots = keep_plots, min.date = env_min_date,
																	max.date = env_max_date, max.predictor.date = env_max_predictor_date)
	} else {
		# For regular processing, use the taxonomic data range
		filt_predictor_data <- lapply(predictor_data, filter_date_site, keep_sites = keep_sites,
																	keep_plots = keep_plots, min.date = min.date,
																	max.date = max.date, max.predictor.date)
	}
	names(filt_predictor_data) <- recode(names(filt_predictor_data), relEM_plot = "relEM")

	# CRITICAL FIX: Handle LAI data structure - convert data.frame to matrix
	if ("LAI" %in% names(filt_predictor_data)) {
		if (is.data.frame(filt_predictor_data$LAI)) {
			filt_predictor_data$LAI <- as.matrix(filt_predictor_data$LAI)
		}
	}

	# Add sine/cosine
	sin_cos_month <- get_sin_cos(colnames(filt_predictor_data$mois))
	filt_predictor_data$sin_mo = sin_cos_month$sin
	filt_predictor_data$cos_mo = sin_cos_month$cos

	# CRITICAL FIX: Update N.date for full time series
	if (full_timeseries) {
		out.list$N.date <- length(colnames(filt_predictor_data$temp))
	}

	out.list <- c(out.list, filt_predictor_data)

	return(out.list)
}

