#'  @title calculate_plot_summary_from_samples
#'  @description Calculate plot summary statistics from samples2 data
#'  @param samples2 MCMC samples containing plot estimates
#'  @return List with statistics and quantiles data frames
calculate_plot_summary_from_samples <- function(samples2) {
	# Combine all samples if there are multiple chains
	if (is.mcmc.list(samples2)) {
		# Multiple chains - combine them
		combined_samples <- do.call(rbind, lapply(samples2, as.matrix))
	} else if (is.matrix(samples2)) {
		# Single matrix
		combined_samples <- samples2
	} else {
		# Single samples object
		combined_samples <- as.matrix(samples2)
	}
	
	# Calculate summary statistics
	means <- apply(combined_samples, 2, mean, na.rm = TRUE)
	quantiles <- apply(combined_samples, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
	
	# Create data frames
	means_df <- data.frame(Mean = means, rowname = names(means))
	quantiles_df <- data.frame(t(quantiles), rowname = names(means))
	colnames(quantiles_df) <- c("2.5%", "25%", "50%", "75%", "97.5%", "rowname")
	
	return(list(means_df, quantiles_df))
}

#'  @title summarize_beta_regression
#'  @description Summarize NIMBLE state-space beta regression models for microbial taxa and functional groups
#'  Assumes input RDS files contain a list of:
#'  MCMC samples, parameter summaries, latent state summaries, and model-fitting metadata
#'	@param overwrite want to save new summary files even if there's an existing, recent summary file
#' @export
summarize_beta_model <- function(file_path, save_summary = NULL, overwrite=NULL, drop_other = TRUE){
	message("=== SUMMARIZE_BETA_MODEL FUNCTION CALLED ===")
	require(stringr)
	require(tidyr)
	# Load helper functions
	if (file.exists("microbialForecast/R/helperFunctions.r")) {
		source("microbialForecast/R/helperFunctions.r")
	} else if (file.exists("R/helperFunctions.r")) {
		source("R/helperFunctions.r")
	}
	# Load summarizeModels functions
	if (file.exists("microbialForecast/R/summarizeModels.r")) {
		source("microbialForecast/R/summarizeModels.r")
	} else if (file.exists("R/summarizeModels.r")) {
		source("R/summarizeModels.r")
	}
	if(summary_exists(file_path)) { # checks that a summary is needed (samples files are new)
		if (is.null(overwrite)) {
			return("Summary file already exists")
		}
	}
	# Read in file, assign named contents to global environment
	read_in <- readRDS(file_path)
	message("File loaded successfully")
	#list2env(read_in,globalenv())

	# Read in samples
	samples <- read_in$samples
	param_summary <- read_in$param_summary
	plot_summary <- read_in$plot_summary
	
	# Debug: Check initial param_summary
	message("  DEBUG: Initial param_summary names: ", paste(names(param_summary), collapse = ", "))
	
	# Check if we have samples2 (plot estimates) and use it for plot_summary if available
	samples2 <- read_in$samples2
	
	# Check if existing plot_summary is malformed (contains parameter data instead of plot data)
	plot_summary_valid <- FALSE
	if (!is.null(plot_summary) && length(plot_summary) > 0) {
		# Check if plot_summary contains actual plot data (plot_mu or Ex parameters)
		if (is.list(plot_summary) && length(plot_summary) >= 1) {
			plot_summary_rownames <- rownames(plot_summary[[1]])
			plot_summary_valid <- any(grepl("plot_mu\\[", plot_summary_rownames)) || 
								 any(grepl("plot_rel\\[", plot_summary_rownames)) ||
								 any(grepl("Ex\\[", plot_summary_rownames))
		}
	}
	
	# Regenerate plot_summary from samples2 if needed
	if (!is.null(samples2) && length(samples2) > 0) {
		# Check if samples2 contains plot_mu parameters and has reasonable dimensions
		if (is.mcmc.list(samples2)) {
			sample_cols <- colnames(samples2[[1]])
			n_cols <- ncol(samples2[[1]])
		} else if (is.matrix(samples2)) {
			sample_cols <- colnames(samples2)
			n_cols <- ncol(samples2)
		} else {
			sample_cols <- NULL
			n_cols <- 0
		}
		
		# Process samples2 if it has plot parameters (plot_mu or Ex) and plot_summary is invalid
		has_plot_params <- any(grepl("plot_mu", sample_cols)) || any(grepl("Ex\\[", sample_cols))
		if (!is.null(sample_cols) && has_plot_params && !plot_summary_valid) {
			message("  Regenerating plot_summary from samples2 (", n_cols, " parameters)")
			# Calculate plot_summary from samples2
			plot_summary <- calculate_plot_summary_from_samples(samples2)
		} else if (plot_summary_valid) {
			message("  Using existing valid plot_summary")
		}
	}
	
	# If plot_summary is still null or empty, use samples for plot estimates
	if (is.null(plot_summary) || length(plot_summary) == 0) {
		samples_for_plot <- samples
	} else {
		samples_for_plot <- samples
	}
	
	# Always recalculate param_summary from samples to ensure correct format
	message("  Recalculating param_summary from samples to ensure correct format...")
	
	# Check if samples exist and are valid
	if (is.null(samples) || length(samples) == 0) {
		message("  ERROR: No samples found, cannot calculate param_summary")
		# Create empty param_summary structure
		means_df <- data.frame(Mean = numeric(0), rowname = character(0))
		quantiles_df <- data.frame(rowname = character(0))
		param_summary <- list(means_df, quantiles_df)
		names(param_summary) <- c("means", "quantiles")
	} else {
		# Combine all samples if there are multiple chains
		if (is.list(samples) && length(samples) > 1) {
			# Multiple chains - combine them
			combined_samples <- do.call(rbind, lapply(samples, as.matrix))
		} else if (is.list(samples) && length(samples) == 1) {
			# Single chain in list
			combined_samples <- as.matrix(samples[[1]])
		} else {
			# Single samples object
			combined_samples <- as.matrix(samples)
		}
		
		# Check if combined_samples is valid
		if (is.null(combined_samples) || nrow(combined_samples) == 0 || ncol(combined_samples) == 0) {
			message("  ERROR: Invalid samples data, cannot calculate param_summary")
			# Create empty param_summary structure
			means_df <- data.frame(Mean = numeric(0), rowname = character(0))
			quantiles_df <- data.frame(rowname = character(0))
			param_summary <- list(means_df, quantiles_df)
			names(param_summary) <- c("means", "quantiles")
		} else {
			# Calculate parameter summaries
			means <- apply(combined_samples, 2, mean, na.rm = TRUE)
			quantiles <- apply(combined_samples, 2, quantile, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
			
			# Create param_summary structure with parameters as rows (not columns)
			# This matches the expected format for extract_summary_row function
			means_df <- data.frame(Mean = means)
			quantiles_df <- as.data.frame(t(quantiles))
			
			# Add rowname column for extract_summary_row function
			means_df$rowname <- rownames(means_df)
			quantiles_df$rowname <- rownames(quantiles_df)
			
			param_summary <- list(means_df, quantiles_df)
			names(param_summary) <- c("means", "quantiles")
			
			message("  Successfully recalculated param_summary from samples")
			message("  DEBUG: param_summary names: ", paste(names(param_summary), collapse = ", "))
			message("  DEBUG: quantiles colnames: ", paste(colnames(param_summary[[2]]), collapse = ", "))
		}
	}
	
	# Handle different data structures
	# Newer models have truth.plot.long nested under model_data
	# Older models have the data directly in model_data
	if("truth.plot.long" %in% names(read_in$metadata$model_data)) {
		truth.plot.long <- read_in$metadata$model_data$truth.plot.long
	} else {
		truth.plot.long <- read_in$metadata$model_data
	}
	
	# Ensure truth.plot.long is a data frame
	if (!is.data.frame(truth.plot.long)) {
		stop("truth.plot.long is not a data frame. Class: ", class(truth.plot.long))
	}


	# Extract run information
	info <- basename(file_path) %>% str_split("_") %>% unlist()
	model_id <- basename(file_path) %>%  str_replace("samples_", "") %>%  str_replace(".rds", "")

	parsed_id = parse_model_id(model_id)
	rank.name.eval <- parsed_id[[1]]
	model_name <- parsed_id[[6]]
	summary_type <- parsed_id[[8]]  # Fixed: should be index 8, not 3
	group  <- parsed_id[[5]]
	time_period <- parsed_id[[2]]

	if (length(info) > 0 && !is.na(tail(info,1)) && tail(info,1) == "summary.rds") {
			info <- info %>% head(-1)
	}

	species <- parsed_id[[4]]
	rank_only <- parsed_id[[3]]

# Add columns based on
	if (summary_type=="functional") {
		rank.name <- rank.name.eval
		rank_only <- summary_type
	fg_cat <- assign_fg_categories(species)
	group <- assign_fg_kingdoms(fg_cat)
	} else {

		taxa_key = stack(microbialForecast:::rank_spec_names) %>%
			select(species = values, rank.name = ind)

		rank.name <- as.character(taxa_key[match(species, taxa_key$species),]$rank.name)
		rank_only <-  rank.name  %>% str_split("_") %>% unlist() %>% head(1)
		group <-  rank.name  %>% str_split("_") %>% unlist() %>% tail(1)
		}
	taxon.name = species

	message("Summarizing ", species, ", ", rank.name, ", ", time_period, ", ", model_name)


	cov_key <- switch(model_name,
										"all_covariates" = microbialForecast:::all_covariates_key,
										"env_cov" = microbialForecast:::all_covariates_key,
										"env_cycl" = microbialForecast:::all_covariates_key,
										"cycl_only" = microbialForecast:::cycl_only_key)

	taxon_key <- unique(truth.plot.long$species)
	names(taxon_key) <- seq(1, length(taxon_key))

	sites <- truth.plot.long %>% select(site_num, siteID) %>% unique()
	site_key <- sites[["siteID"]]
	names(site_key) <- sites[["site_num"]]

	# Add some info to observational data for merging
	# Check if dateID column exists, if not use a default date
	if ("dateID" %in% colnames(truth.plot.long)) {
		truth.plot.long <- truth.plot.long %>%
			mutate(dates = fixDate(dateID),
						 truth = as.numeric(truth),
						 model_name = !!model_name,
						 taxon = species,
						 rank_name = rank.name,
						 group = !!group,
						 rank_only = !!rank_only,
						 time_period = !!time_period,
						 fcast_type = !!summary_type,
						 pretty_group = ifelse(!is.na(group) & group %in% c("16S","bac"), "Bacteria", "Fungi"),
						 model_id = !!model_id) %>%
			mutate(time_period =
						 	recode(as.character(!!time_period), !!!microbialForecast:::date_recode))
	} else {
		# If no dateID column, add minimal required columns
		truth.plot.long <- truth.plot.long %>%
			mutate(truth = as.numeric(truth),
						 model_name = !!model_name,
						 taxon = species,
						 rank_name = rank.name,
						 group = !!group,
						 rank_only = !!rank_only,
						 time_period = !!time_period,
						 fcast_type = !!summary_type,
						 pretty_group = ifelse(!is.na(group) & group %in% c("16S","bac"), "Bacteria", "Fungi"),
						 model_id = !!model_id) %>%
			mutate(time_period =
						 	recode(as.character(!!time_period), !!!microbialForecast:::date_recode))
	}

	if (summary_type=="functional") {
		truth.plot.long <- truth.plot.long %>%
			mutate(fg_cat = !!fg_cat,
						 fcast_type = "Functional")
	} else {
		truth.plot.long <- truth.plot.long %>%
			mutate(fg_cat = NA,
						 fcast_type = "Taxonomic")
	}

	if (drop_other) {
		truth.plot.long <- truth.plot.long %>% filter(species != "other")
	}

	if (!is.na(species) && nchar(species) > 0) {
		taxon_key[1] = species
	}

	# Calculate plot median and quantiles
	pred.quantiles <- plot_summary[[2]] %>% parse_plot_mu_vars() %>%
		merge(truth.plot.long, by = c("plot_num", "timepoint"), all = T)


	# For scoring the predictions, need mean and SD
	pred.means <- plot_summary[[1]] %>% parse_plot_mu_vars() %>%
		merge(truth.plot.long, by = c("plot_num", "timepoint"), all = T)

	pred.quantiles$Mean <- pred.means$Mean
	pred.quantiles$SD <- pred.means$SD

	# Get mean values for parameters
	means <- param_summary[[1]]
	quantiles <- param_summary[[2]]
	
	# Debug: Check param_summary structure
	message("  DEBUG: param_summary names: ", paste(names(param_summary), collapse = ", "))
	message("  DEBUG: means colnames: ", paste(colnames(means), collapse = ", "))
	message("  DEBUG: quantiles colnames: ", paste(colnames(quantiles), collapse = ", "))
	
	# Check if param_summary has any data
	if (nrow(means) == 0 || nrow(quantiles) == 0) {
		message("  WARNING: Empty param_summary, creating empty parameter lists")
		eff_list <- data.frame(taxon = species, stringsAsFactors = FALSE)
		eff_list2 <- data.frame(taxon = species, stringsAsFactors = FALSE)
	} else {
		# Use correct parameter names for beta regression models
		eff_list <- lapply(c("precision", "intercept", "rho", "legacy_effect", "site_effect_sd"),
											 function(x) extract_summary_row(param_summary[[1]], var = x)) %>%
			plyr::rbind.fill() %>%
			mutate(taxon = !!species)
		eff_list2 <- lapply(c("precision", "intercept", "rho", "legacy_effect", "site_effect_sd"),
												function(x) extract_summary_row(param_summary[[2]], var = x)) %>%
			plyr::rbind.fill() %>%
			mutate(taxon = !!species)
	}
	
	# Only assign Median if both lists have the same number of rows
	if (nrow(eff_list) > 0 && nrow(eff_list2) > 0 && nrow(eff_list) == nrow(eff_list2)) {
		# Debug: Check if 50% column exists
		if ("50%" %in% colnames(eff_list2)) {
			eff_list$Median = eff_list2[,"50%"]
		} else {
			message("  WARNING: eff_list2 missing 50% column, colnames: ", paste(colnames(eff_list2), collapse = ", "))
			eff_list$Median = NA
		}
	} else if (nrow(eff_list) > 0) {
		eff_list$Median = NA
	}

	# Get site effect sizes
	if (nrow(means) == 0 || nrow(quantiles) == 0) {
		site_eff_out <- data.frame(taxon = species, siteID = "UNKNOWN", stringsAsFactors = FALSE)
		site_eff_out2 <- data.frame(taxon = species, siteID = "UNKNOWN", `50%` = NA, stringsAsFactors = FALSE)
	} else {
		site_eff_out <- tryCatch({
			extract_summary_row(param_summary[[1]], var = "site_effect") %>%
				extract_bracketed_vals(varname1 = "site_num")  %>%
				mutate(taxon = !!species,
							 siteID = recode(site_num, !!!site_key))
		}, error = function(e) {
			cat("  WARNING: Error extracting site effects:", e$message, "\n")
			data.frame(taxon = species, siteID = "UNKNOWN", stringsAsFactors = FALSE)
		})
		
		site_eff_out2 <- tryCatch({
			extract_summary_row(param_summary[[2]], var = "site_effect") %>%
				extract_bracketed_vals(varname1 = "site_num")  %>%
				mutate(taxon = !!species,
							 siteID = recode(site_num, !!!site_key))
		}, error = function(e) {
			cat("  WARNING: Error extracting site effect quantiles:", e$message, "\n")
			data.frame(taxon = species, siteID = "UNKNOWN", `50%` = NA, stringsAsFactors = FALSE)
		})
	}
	
	# Only add Median if the column exists
	if ("50%" %in% colnames(site_eff_out2)) {
		site_eff_out$Median = site_eff_out2[,"50%"]
	} else {
		site_eff_out$Median = NA
	}

	# Get beta sizes per rank
	if (nrow(means) == 0 || nrow(quantiles) == 0) {
		beta_out <- data.frame(taxon = species, stringsAsFactors = FALSE)
		beta_ci <- data.frame(taxon = species, stringsAsFactors = FALSE)
	} else {
		beta_out <- extract_summary_row(param_summary[[1]], var = "beta") %>%
			extract_bracketed_vals(varname1 = "beta_num")  %>%
			mutate(beta = recode(beta_num, !!!cov_key),
						 taxon = !!species)
						# taxon = recode(taxon_num, !!!taxon_key))

		# Use quantiles to assign significance to beta parameters.
		beta_ci <- extract_summary_row(param_summary[[2]], var = "beta") %>%
			extract_bracketed_vals(varname1 = "beta_num")  %>%
			mutate(beta = recode(beta_num, !!!cov_key),
						 taxon = !!species)
						# taxon = recode(taxon_num, !!!taxon_key))
	}
	
	# Only calculate significance if we have valid data
	if (nrow(beta_out) > 0 && nrow(beta_ci) > 0 && all(!is.na(beta_ci$`2.5%`)) && all(!is.na(beta_ci$`97.5%`))) {
		beta_out$significant <- microbialForecast:::is_significant(beta_ci$`2.5%`, beta_ci$`97.5%`)
		beta_out$effSize <- abs(beta_out$Mean)
	} else {
		beta_out$significant <- NA
		beta_out$effSize <- NA
	}

	# Combine parameter estimates into summary
	#if (length(unique(beta_out$taxon)) > 2){
	summary_df <-
		plyr::rbind.fill(beta_out, eff_list, site_eff_out) %>% mutate(rank = rank.name)
	
	# Add metadata columns directly instead of joining
	summary_df <- summary_df %>%
		mutate(model_name = !!model_name,
					 group = !!group,
					 rank_only = !!rank_only,
					 time_period = !!time_period,
					 fcast_type = !!summary_type,
					 pretty_group = ifelse(!is.na(group) & group %in% c("16S","bac"), "Bacteria", "Fungi"),
					 model_id = !!model_id)
	# } else {
	# 	summary_df <-
	# 		plyr::rbind.fill(beta_out, eff_list, site_eff_out)  %>% mutate(rank = rank.name) %>%
	# 		left_join(truth.plot.long[1, 11:18])
	# }


	## Calculate gelman diagnostics to assess convergence
	gd <- tryCatch({
		add_gelman(read_in, rank.name) %>% mutate(rank = rank.name,  taxon = !!species) %>%
			mutate(model_name = !!model_name,
						 group = !!group,
						 rank_only = !!rank_only,
						 time_period = !!time_period,
						 fcast_type = !!summary_type,
						 pretty_group = ifelse(!is.na(group) & group %in% c("16S","bac"), "Bacteria", "Fungi"),
						 model_id = !!model_id)
	}, error = function(e) {
		cat("  WARNING: Error calculating Gelman diagnostics:", e$message, "\n")
		# Create a minimal data frame with required columns
		data.frame(
			parameter = "UNKNOWN",
			Point.est = NA,
			Upper.C.I. = NA,
			effSize = NA,
			rank.name = rank.name,
			taxon = species,
			model_name = model_name,
			group = group,
			rank_only = rank_only,
			time_period = time_period,
			fcast_type = summary_type,
			pretty_group = ifelse(!is.na(group) & group %in% c("16S","bac"), "Bacteria", "Fungi"),
			model_id = model_id,
			stringsAsFactors = FALSE
		)
	})

	if (drop_other) {
		summary_df <- summary_df %>% filter(taxon != "other")
		#pred.means <- pred.means %>% filter(taxon != "other")
		pred.quantiles <- pred.quantiles %>% filter(taxon != "other")
	}

	out <- list(summary_df, pred.means, pred.quantiles, gd)
	if (!is.null(save_summary)) {
		savePath <- gsub("samples","summary",file_path)
		saveRDS(out, savePath)
		message("Saved summary to ", savePath)
		return(TRUE)
	} else {
		return(out)
	}
}