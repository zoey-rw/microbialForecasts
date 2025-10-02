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

	# For Dirichlet models, we need to identify which taxa to keep
	# Get all non-metadata columns (these are the taxa)
	metadata_cols <- c("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date")
	taxa_cols <- setdiff(colnames(dat), metadata_cols)
	
	# Filter to top 3 taxa plus "other" category (ultra-reduced complexity for Dirichlet)
	taxa_prevalence <- sapply(taxa_cols, function(taxa) sum(dat[[taxa]] > 0, na.rm = TRUE))
	taxa_prevalence <- taxa_prevalence[order(taxa_prevalence, decreasing = TRUE)]
	
	# Keep top 3 taxa
	top_taxa <- names(taxa_prevalence)[1:min(3, length(taxa_prevalence))]
	
	# Create "other" category by summing remaining taxa
	other_taxa <- names(taxa_prevalence)[4:length(taxa_prevalence)]
	if (length(other_taxa) > 0) {
		dat$other <- rowSums(dat[, other_taxa, drop = FALSE], na.rm = TRUE)
		keep_taxa <- c(top_taxa, "other")
	} else {
		# If we have 3 or fewer taxa, just use them all
		keep_taxa <- top_taxa
	}
	
	# Remove the original taxa columns that are not in keep_taxa
	# This prevents duplication issues
	remove_taxa <- setdiff(taxa_cols, keep_taxa)
	if (length(remove_taxa) > 0) {
		dat <- dat[, !colnames(dat) %in% remove_taxa]
	}
	
	if (length(keep_taxa) == 0) {
		stop("No taxa available for Dirichlet modeling")
	}
	cat("Keeping", length(keep_taxa, "\n"), "taxa (top 3 + other) for Dirichlet modeling\n")
	
	# Note: y matrix will be created after all NA filtering is complete
	
	# Create timepoint mapping
	unique_dates <- sort(unique(dat$dateID))
	# Remove any NA values from unique_dates
	unique_dates <- unique_dates[!is.na(unique_dates)]
	date_to_timepoint <- setNames(1:length(unique_dates), unique_dates)
	timepoint <- date_to_timepoint[as.character(dat$dateID)]
	
	# Handle any NA values in timepoint by filtering out those rows
	if (any(is.na(timepoint))) {
		cat("WARNING: Found", sum(is.na(timepoint, "\n")), "NA values in timepoint, filtering out these rows\n")
		# Filter out rows with NA timepoint
		valid_rows <- !is.na(timepoint)
		dat <- dat[valid_rows, ]
		timepoint <- timepoint[valid_rows]
		cat("  Filtered to", nrow(dat, "\n"), "rows after removing NA timepoint values\n")
	}
	
	# Create plot mapping
	unique_plots <- sort(unique(dat$plotID))
	# Remove any NA values from unique_plots
	unique_plots <- unique_plots[!is.na(unique_plots)]
	plot_to_num <- setNames(1:length(unique_plots), unique_plots)
	plot_num <- plot_to_num[as.character(dat$plotID)]
	
	# Handle any NA values in plot_num by filtering out those rows
	if (any(is.na(plot_num))) {
		cat("WARNING: Found", sum(is.na(plot_num, "\n")), "NA values in plot_num, filtering out these rows\n")
		# Filter out rows with NA plot_num
		valid_rows <- !is.na(plot_num)
		dat <- dat[valid_rows, ]
		plot_num <- plot_num[valid_rows]
		timepoint <- timepoint[valid_rows]  # Also update timepoint to match
		cat("  Filtered to", nrow(dat, "\n"), "rows after removing NA plot_num values\n")
	}
	
	# Create site mapping
	unique_sites <- sort(unique(dat$siteID))
	# Remove any NA values from unique_sites
	unique_sites <- unique_sites[!is.na(unique_sites)]
	site_to_num <- setNames(1:length(unique_sites), unique_sites)
	site_num <- site_to_num[as.character(dat$siteID)]
	
	# Handle any NA values in site_num by filtering out those rows
	if (any(is.na(site_num))) {
		cat("WARNING: Found", sum(is.na(site_num, "\n")), "NA values in site_num, filtering out these rows\n")
		# Filter out rows with NA site_num
		valid_rows <- !is.na(site_num)
		dat <- dat[valid_rows, ]
		site_num <- site_num[valid_rows]
		plot_num <- plot_num[valid_rows]  # Also update plot_num to match
		timepoint <- timepoint[valid_rows]  # Also update timepoint to match
		cat("  Filtered to", nrow(dat, "\n"), "rows after removing NA site_num values\n")
	}
	
	# Create plot_site_num mapping
	plot_site_num <- rep(NA, length(unique_plots))
	for (i in 1:length(unique_plots)) {
		plot_id <- unique_plots[i]
		plot_data <- dat[dat$plotID == plot_id, ]
		if (nrow(plot_data) > 0) {
			site_id <- unique(plot_data$siteID)[1]
			plot_site_num[i] <- site_to_num[as.character(site_id)]
		}
	}
	
	# Create plot_start mapping (first timepoint for each plot)
	plot_start <- rep(1, length(unique_plots))
	for (i in 1:length(unique_plots)) {
		plot_id <- unique_plots[i]
		plot_data <- dat[dat$plotID == plot_id, ]
		if (nrow(plot_data) > 0) {
			plot_timepoints <- date_to_timepoint[as.character(plot_data$dateID)]
			plot_start[i] <- min(plot_timepoints, na.rm = TRUE)
		}
	}
	
	# Create plot_index (for compatibility with existing code)
	plot_index <- 1:length(unique_plots)
	
	# Create site_start (for compatibility with existing code)
	site_start <- rep(1, length(unique_sites))
	
	# Create truth.plot.long (for compatibility with existing code)
	truth.plot.long <- dat[, c("plotID", "dateID", "sampleID")]
	
	# Create dates_per_plot (for compatibility with existing code)
	dates_per_plot <- rep(1, length(unique_plots))
	for (i in 1:length(unique_plots)) {
		plot_id <- unique_plots[i]
		plot_data <- dat[dat$plotID == plot_id, ]
		dates_per_plot[i] <- length(unique(plot_data$dateID))
	}
	
	# Get environmental predictors using the same approach as prepBetaRegData
	env_predictors <- c("temp", "mois", "pH", "pC", "relEM", "LAI")
	env_data <- list()
	
	cat("Available predictors in predictor_data:", paste(names(predictor_data, "\n"), collapse = ", "), "\n")
	
	for (pred in env_predictors) {
		# Handle special case where relEM is stored as relEM_plot
		pred_name <- ifelse(pred == "relEM" && "relEM_plot" %in% names(predictor_data), "relEM_plot", pred)
		
		if (pred_name %in% names(predictor_data)) {
			env_data[[pred]] <- predictor_data[[pred_name]]
			cat("Added", pred, "predictor (from", pred_name, ")\n")
		} else {
			cat("WARNING:", pred, "predictor not found in predictor_data\n", "\n")
		}
	}
	
	# Add seasonal predictors
	months <- 1:length(unique_dates)
	env_data$sin_mo <- sin(2 * pi * months / 12)
	env_data$cos_mo <- cos(2 * pi * months / 12)
	
	# Create the response matrix for Dirichlet (compositional data) AFTER all filtering
	y <- as.matrix(dat[, keep_taxa, drop = FALSE])
	
	# Ensure all values are between 0 and 1 (proportions)
	y <- pmin(pmax(y, 0), 1)
	
	# Normalize to sum to 1 for each row (Dirichlet requirement)
	y_sums <- rowSums(y)
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
		N.date = length(unique_dates),
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
