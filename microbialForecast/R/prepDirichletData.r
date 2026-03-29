#' @title 			prepDirichletData
#' @description Prepare data for Dirichlet regression models with compositional data
#' @export
prepDirichletData <- function(rank.df,
															min.date = "20130601",
															max.date = "20170101",
															predictor_data = NULL,
															min.prev = 3,
															dom_soil_horizons = NULL,
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

	if (!is.null(keep_vec)){
		dat <- dat[,colnames(rank.df) %in% keep_vec]
	}

	# Convert dates to proper format
	dat$dates <- as.Date(dat$dates, format = "%Y%m%d")
	dat$plot_date <- as.Date(dat$plot_date, format = "%Y%m%d")

	# Filter by date range
	dat <- dat[dat$dates >= as.Date(min.date, format = "%Y%m%d") & 
						 dat$dates <= as.Date(max.date, format = "%Y%m%d"), ]

	# Get unique identifiers
	siteID <- unique(dat$siteID)
	plotID <- unique(dat$plotID)
	dateID <- unique(dat$dateID)

	# For Dirichlet models, we need at least 3 taxa (composition requires multiple components)
	metadata_cols <- c("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date")
	taxa_cols <- setdiff(colnames(dat), metadata_cols)

	# Exclude pre-existing "other" column from prevalence ranking — it's not a real taxon
	real_taxa_cols <- setdiff(taxa_cols, "other")
	taxa_prevalence <- sapply(real_taxa_cols, function(taxa) sum(dat[[taxa]] > 0, na.rm = TRUE))
	taxa_prevalence <- taxa_prevalence[order(taxa_prevalence, decreasing = TRUE)]

	if (length(taxa_prevalence) < 3) {
		stop("Dirichlet model requires at least 3 taxa; this rank has ", length(taxa_prevalence), ".")
	}

	# Keep top 3 real taxa; sum everything else (including pre-existing "other") into "other"
	top_taxa <- names(taxa_prevalence)[1:3]
	other_taxa <- setdiff(taxa_cols, top_taxa)  # includes remaining real taxa + pre-existing "other"
	if (length(other_taxa) > 0) {
		dat$other <- rowSums(dat[, other_taxa, drop = FALSE], na.rm = TRUE)
		keep_taxa <- c(top_taxa, "other")
	} else {
		keep_taxa <- top_taxa
	}

	remove_taxa <- setdiff(taxa_cols, keep_taxa)
	if (length(remove_taxa) > 0) {
		dat <- dat[, !colnames(dat) %in% remove_taxa]
	}
	cat("Keeping", length(keep_taxa), "taxa for Dirichlet modeling:", paste(keep_taxa, collapse=", "), "\n")

	# Note: y matrix will be created after all NA filtering is complete

	# --- CRITICAL FIX: Date Reconstruction ---
	if (!"dateID" %in% colnames(dat)) stop("Required column 'dateID' missing")
	if (any(is.na(dat$dateID))) {
		na_mask <- is.na(dat$dateID)
		if ("sampleID" %in% colnames(dat)) {
			try({
				parsed <- parseNEONsampleIDs(dat$sampleID[na_mask])
				if ("dateID" %in% colnames(parsed)) {
					dat$dateID[na_mask] <- suppressWarnings(as.numeric(as.character(parsed$dateID)))
				}
			}, silent = TRUE)
		}
		if (any(is.na(dat$dateID)) && "dates" %in% colnames(dat)) {
			dat$dateID[is.na(dat$dateID)] <- as.numeric(format(as.Date(dat$dates[is.na(dat$dateID)]), "%Y%m"))
		}
		dat <- dat[!is.na(dat$dateID), ]
	}

	# --- CRITICAL FIX: Continuous Timeline for AR(1) ---
	poss_dateID <- as.numeric(format(seq.Date(as.Date(min.date, format = "%Y%m%d"), as.Date(max.date, format = "%Y%m%d"), by = "month"), "%Y%m"))

	date_to_timepoint <- setNames(1:length(poss_dateID), poss_dateID)
	timepoint <- date_to_timepoint[as.character(dat$dateID)]

	dat <- dat[!is.na(timepoint), ]
	timepoint <- timepoint[!is.na(timepoint)]

	# Create plot and site mappings
	unique_plots <- sort(unique(dat$plotID))
	plot_to_num <- setNames(1:length(unique_plots), unique_plots)
	plot_num <- plot_to_num[as.character(dat$plotID)]

	unique_sites <- sort(unique(dat$siteID))
	site_to_num <- setNames(1:length(unique_sites), unique_sites)
	site_num <- site_to_num[as.character(dat$siteID)]

	# Create plot_site_num mapping
	plot_site_num <- sapply(unique_plots, function(p) {
		site_to_num[as.character(dat$siteID[dat$plotID == p][1])]
	})

	# Create plot_start relative to continuous timeline
	plot_start <- sapply(unique_plots, function(p) {
		min(date_to_timepoint[as.character(dat$dateID[dat$plotID == p])], na.rm = TRUE)
	})

	plot_index <- 1:length(unique_plots)

	# CRITICAL FIX 1: Calculate actual site_start instead of rep(1)
	# This ensures the sanitizer ignores valid NAs that exist before the site's first observation
	site_start <- sapply(unique_sites, function(s) {
		min(date_to_timepoint[as.character(dat$dateID[dat$siteID == s])], na.rm = TRUE)
	})

	truth.plot.long <- dat[, c("plotID", "dateID", "sampleID", "siteID")]
	truth.plot.long$site_num <- as.integer(site_num)

	dates_per_plot <- sapply(unique_plots, function(p) length(unique(dat$dateID[dat$plotID == p])))

	# --- CRITICAL FIX 2: Use filter_date_site to ensure clean environmental matrices ---
	keep_plots <- as.character(unique_plots)
	keep_sites <- as.character(unique_sites)

	# Inject rownames into uncertainty matrices so filter_date_site can actually filter them
	for (var in c("temp_sd", "mois_sd", "LAI")) {
		if (var %in% names(predictor_data) && is.null(rownames(predictor_data[[var]])) && !is.null(rownames(predictor_data$temp))) {
			rownames(predictor_data[[var]]) <- rownames(predictor_data$temp)
		}
	}
	for (var in c("pH_sd", "pC_sd", "relEM_plot")) {
		if (var %in% names(predictor_data) && is.null(rownames(predictor_data[[var]])) && !is.null(rownames(predictor_data$pH))) {
			rownames(predictor_data[[var]]) <- rownames(predictor_data$pH)
		}
	}

	# Only run filter_date_site on env vars the Dirichlet model needs (exclude site_skip, nspp, etc. that lack filterable rownames)
	env_vars_needed <- c("temp", "mois", "LAI", "temp_sd", "mois_sd", "pH", "pC", "relEM_plot", "pH_sd", "pC_sd")
	predictor_data_sub <- predictor_data[names(predictor_data) %in% env_vars_needed]

	if (full_timeseries) {
		# For full time series, use the full environmental data range
		env_min_date <- min(colnames(predictor_data$temp))
		env_max_date <- max(colnames(predictor_data$temp))
		env_max_predictor_date <- paste0(substr(env_max_date, 1, 4), "-", substr(env_max_date, 5, 6), "-01")
		env_data <- lapply(predictor_data_sub, filter_date_site, keep_sites = keep_sites,
																keep_plots = keep_plots, min.date = env_min_date,
																max.date = env_max_date, max.predictor.date = env_max_predictor_date)
	} else {
		# For regular processing, use the taxonomic data range
		env_data <- lapply(predictor_data_sub, filter_date_site, keep_sites = keep_sites,
																keep_plots = keep_plots, min.date = min.date,
																max.date = max.date, max.predictor.date = NULL)
	}
	names(env_data) <- dplyr::recode(names(env_data), relEM_plot = "relEM")

	# CRITICAL FIX: Handle LAI data structure - convert data.frame to matrix
	if ("LAI" %in% names(env_data)) {
		if (is.data.frame(env_data$LAI)) {
			env_data$LAI <- as.matrix(env_data$LAI)
		}
	}

	# Add sine/cosine matching the environmental data columns
	sin_cos_month <- get_sin_cos(colnames(env_data$mois))
	env_data$sin_mo = sin_cos_month$sin
	env_data$cos_mo = sin_cos_month$cos

	# Update N.date for full time series if necessary
	if (full_timeseries) {
		poss_dateID <- as.numeric(gsub("-", "", substr(colnames(env_data$temp), 1, 7)))
	}
	
	# Create the response matrix for Dirichlet (compositional data) AFTER all filtering
	y <- as.matrix(dat[, keep_taxa, drop = FALSE])
	
	# Ensure all values are between 0 and 1 (proportions)
	y <- pmin(pmax(y, 0), 1)

	# Add pseudocount to prevent exact zeros (Dirichlet density is -Inf at zero)
	y[y == 0] <- 1e-6

	# Normalize to sum to 1 for each row (Dirichlet requirement)
	y_sums <- rowSums(y)

	# FIX: Prevent NaN if an observation row is all zeroes
	y[y_sums == 0, ] <- 1 / ncol(y)
	y_sums[y_sums == 0] <- 1

	y <- y / y_sums
	
	# Create the final output list
	output <- list(
		y = y,
		siteID = dat$siteID,
		plotID = dat$plotID,
		plot_site = dat$plotID,  # For compatibility
		site_start = site_start,
		plot_start = plot_start,
		plot_index = plot_index,
		plot_num = plot_num,
		plot_site_num = plot_site_num,
		truth.plot.long = truth.plot.long,
		N.date = length(poss_dateID),
		timepoint = timepoint,
		dates_per_plot = dates_per_plot,
		N.plot = length(unique_plots),
		N.spp = length(keep_taxa),
		N.core = nrow(y),
		N.site = length(unique_sites),
		sample_values = y,  # For compatibility
		keep_taxa = keep_taxa  # Save the taxa names being modeled
	)
	
	# Add environmental data
	output <- c(output, env_data)
	
	return(output)
}
