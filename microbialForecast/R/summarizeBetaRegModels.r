#'  @title calculate_plot_summary_from_samples
#'  @description Calculate plot summary statistics from samples2 data
#'  @param samples2 MCMC samples containing plot estimates
#'  @return List with statistics and quantiles data frames
calculate_plot_summary_from_samples <- function(samples2) {
  # normalize names early (mu[...] -> plot_mu[...])
  .norm_names <- function(nm) sub("^mu\\[", "plot_mu[", nm)

  to_mat <- function(x) {
    if (inherits(x, "mcmc")) return(as.matrix(x))
    x <- as.matrix(x)
    storage.mode(x) <- "double"
    x
  }

  if (coda::is.mcmc.list(samples2)) {
    base_cols <- .norm_names(colnames(samples2[[1]]))
    mats <- lapply(samples2, function(ch) {
      m <- to_mat(ch)
      colnames(m) <- .norm_names(colnames(m))
      keep <- intersect(base_cols, colnames(m))
      m[, keep, drop = FALSE]
    })
    common <- Reduce(intersect, lapply(mats, colnames))
    # preserve ordering of base_cols
    common <- base_cols[base_cols %in% common]
    mats <- lapply(mats, function(m) m[, common, drop = FALSE])
    combined_samples <- do.call(rbind, mats)

  } else if (is.matrix(samples2)) {
    combined_samples <- to_mat(samples2)
    colnames(combined_samples) <- .norm_names(colnames(combined_samples))

  } else if (is.list(samples2)) {
    # list-of-matrices or list-of-mcmc
    mats <- lapply(samples2, to_mat)
    # use first element as base
    base_cols <- .norm_names(colnames(mats[[1]]))
    mats <- lapply(mats, function(m) {
      colnames(m) <- .norm_names(colnames(m))
      keep <- intersect(base_cols, colnames(m))
      m[, keep, drop = FALSE]
    })
    common <- Reduce(intersect, lapply(mats, colnames))
    common <- base_cols[base_cols %in% common]
    combined_samples <- do.call(rbind, lapply(mats, function(m) m[, common, drop = FALSE]))

  } else {
    combined_samples <- to_mat(samples2)
    colnames(combined_samples) <- .norm_names(colnames(combined_samples))
  }

  if (is.null(colnames(combined_samples)))
    stop("samples2 has no column names; expected names like 'plot_mu[plot,time]'.")

  means  <- apply(combined_samples, 2, mean, na.rm = TRUE)
  sds    <- apply(combined_samples, 2, sd,   na.rm = TRUE)
  quants <- apply(combined_samples, 2, quantile,
                  probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)

  means_df <- data.frame(Mean = means, SD = sds, check.names = FALSE)
  means_df$param <- names(means)
  rownames(means_df) <- names(means)

  quantiles_df <- as.data.frame(t(quants), check.names = FALSE)
  colnames(quantiles_df) <- c("2.5%", "25%", "50%", "75%", "97.5%")
  quantiles_df$param <- rownames(quantiles_df)

  list(means = means_df, quantiles = quantiles_df)
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
		if (is.null(overwrite) || overwrite == FALSE) {
			return("Summary file already exists")
		}
	}
	# Read in file, assign named contents to global environment
	read_in <- readRDS(file_path)
	message("File loaded successfully")
	message("  DEBUG: File loaded, checking structure...")
	#list2env(read_in,globalenv())

	# Read in samples
	samples <- read_in$samples
	param_summary <- read_in$param_summary
	plot_summary <- read_in$plot_summary
	
	# Debug: Check initial param_summary
	message("  DEBUG: Initial param_summary names: ", paste(names(param_summary), collapse = ", "))
	
	# Handle different param_summary structures for driver uncertainty models
	# Driver uncertainty models have deeply nested structures, so we'll skip the nested handling
	# and let the function recalculate param_summary from samples
	message("  DEBUG: Driver uncertainty model detected, will recalculate param_summary from samples")
	
	message("  DEBUG: Final param_summary names: ", paste(names(param_summary), collapse = ", "))
	message("  DEBUG: Data extracted successfully, proceeding to model ID parsing...")
	
	# Convert summary.mcmc to list format if needed
	if (inherits(plot_summary, "summary.mcmc")) {
		message("  Converting summary.mcmc plot_summary to list format")
		plot_summary <- list(
			plot_summary$statistics,
			plot_summary$quantiles
		)
		names(plot_summary) <- c("statistics", "quantiles")
	}
	
	# Validate plot_summary has plot data
	if (!is.null(plot_summary) && is.list(plot_summary) && length(plot_summary) >= 2) {
		rownames <- rownames(plot_summary[[2]])
		has_plot_data <- any(grepl("^(plot_mu|plot_rel|Ex)\\[", rownames))
		if (!has_plot_data) {
			message("  WARNING: plot_summary does not contain plot data")
			plot_summary <- NULL  # Will create empty estimates later
		}
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
			sds   <- apply(combined_samples, 2, sd,   na.rm = TRUE)
			quants <- apply(combined_samples, 2, quantile, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
			
			# Create param_summary structure with rownames set to parameter names
			means_df <- data.frame(Mean = means, SD = sds, check.names = FALSE)
			rownames(means_df) <- colnames(combined_samples)
			means_df$param <- rownames(means_df)
			
			quantiles_df <- as.data.frame(t(quants), check.names = FALSE)
			colnames(quantiles_df) <- c("2.5%", "25%", "50%", "75%", "97.5%")
			rownames(quantiles_df) <- colnames(combined_samples)
			quantiles_df$param <- rownames(quantiles_df)
			
			param_summary <- list(means = means_df, quantiles = quantiles_df)
			
			message("  Successfully recalculated param_summary from samples")
			message("  DEBUG: param_summary names: ", paste(names(param_summary), collapse = ", "))
			message("  DEBUG: quantiles colnames: ", paste(colnames(param_summary[[2]]), collapse = ", "))
		}
	}
	
	# Handle different data structures
	# Newer models have truth.plot.long nested under model_data
	# Older models have the data directly in model_data
	message("  DEBUG: Extracting truth.plot.long...")
	if("truth.plot.long" %in% names(read_in$metadata$model_data)) {
		truth.plot.long <- read_in$metadata$model_data$truth.plot.long
		message("  DEBUG: Using nested truth.plot.long")
	} else {
		truth.plot.long <- read_in$metadata$model_data
		message("  DEBUG: Using direct model_data")
	}
	message("  DEBUG: truth.plot.long extracted, class:", class(truth.plot.long))
	
	# Ensure truth.plot.long is a data frame
	if (!is.data.frame(truth.plot.long)) {
		stop("truth.plot.long is not a data frame. Class: ", class(truth.plot.long))
	}


	# Extract run information
	info <- basename(file_path) %>% str_split("_") %>% unlist()
	model_id <- basename(file_path) %>%  str_replace("samples_", "") %>%  str_replace(".rds", "")

	message("  DEBUG: About to parse model_id:", model_id)
	tryCatch({
		parsed_id = parse_model_id(model_id)
		message("  DEBUG: Model ID parsed successfully")
	}, error = function(e) {
		message("  ERROR in parse_model_id:", e$message)
		stop(e)
	})
	rank.name.eval <- parsed_id[[1]]
	model_name <- parsed_id[[6]]
	summary_type <- parsed_id[[8]]  # Fixed: should be index 8, not 3
	group  <- parsed_id[[5]]
	time_period <- parsed_id[[2]]
	has_driver_uncertainty <- parsed_id[[9]]  # New driver uncertainty flag

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
	message("  DEBUG: Starting parameter extraction...")


	cov_key <- switch(model_name,
										"all_covariates" = microbialForecast:::all_covariates_key,
										"env_cov" = microbialForecast:::all_covariates_key,
										"env_cycl" = microbialForecast:::all_covariates_key,
										"cycl_only" = microbialForecast:::cycl_only_key)
	message("  DEBUG: cov_key assigned, length:", length(cov_key))

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
						 pretty_group = ifelse(!is.na(group) & !is.na(group %in% c("16S","bac")) & group %in% c("16S","bac"), "Bacteria", "Fungi"),
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
						 pretty_group = ifelse(!is.na(group) & !is.na(group %in% c("16S","bac")) & group %in% c("16S","bac"), "Bacteria", "Fungi"),
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

	# Initialize taxon_key before using it
	taxon_key <- unique(truth.plot.long$species)
	names(taxon_key) <- seq(1, length(taxon_key))
	
	if (!is.na(species) && nchar(species) > 0) {
		taxon_key[1] = species
	}

	# Calculate plot median and quantiles
	# Handle different plot_summary structures
	if (is.list(plot_summary) && length(plot_summary) >= 2) {
		# For standard plot_summary with plot_mu parameters
		# Handle mu[...] parameters by renaming to plot_mu[...] if needed
		rownms <- rownames(plot_summary[[2]])
		if (any(grepl("^mu\\[", rownms)) && !any(grepl("^plot_mu\\[", rownms))) {
			message("  Renaming mu[...] parameters to plot_mu[...] for compatibility")
			rownames(plot_summary[[1]]) <- sub("^mu\\[", "plot_mu[", rownames(plot_summary[[1]]))
			rownames(plot_summary[[2]]) <- sub("^mu\\[", "plot_mu[", rownames(plot_summary[[2]]))
		}
		
		# CRITICAL FIX: Ensure truth values are preserved correctly during merge
		plot_quantiles <- plot_summary[[2]] %>% parse_plot_mu_vars()
		plot_means <- plot_summary[[1]] %>% parse_plot_mu_vars()
		
		# Merge with truth data, being careful about column name conflicts
		pred.quantiles <- plot_quantiles %>%
			merge(truth.plot.long, by = c("plot_num", "timepoint"), all = T, suffixes = c("", "_truth"))
		
		pred.means <- plot_means %>%
			merge(truth.plot.long, by = c("plot_num", "timepoint"), all = T, suffixes = c("", "_truth"))

		# CRITICAL FIX: Ensure we use the correct truth column (not dateID values)
		# If there are duplicate columns, prefer the one from truth.plot.long
		if ("truth_truth" %in% names(pred.quantiles)) {
			pred.quantiles$truth <- pred.quantiles$truth_truth
			pred.quantiles$truth_truth <- NULL
		}
		if ("truth_truth" %in% names(pred.means)) {
			pred.means$truth <- pred.means$truth_truth
			pred.means$truth_truth <- NULL
		}
		
		# CRITICAL FIX: Validate that truth values are in [0,1] range (not dateIDs)
		if (any(pred.quantiles$truth > 1, na.rm = TRUE)) {
			message("  WARNING: Detected corrupted truth values (dateIDs) - fixing...")
			# Find rows where truth values look like dateIDs (> 200000)
			corrupted_rows <- pred.quantiles$truth > 200000
			if (any(corrupted_rows, na.rm = TRUE)) {
				message("  Found ", sum(corrupted_rows, na.rm = TRUE), " corrupted truth values")
				# Set corrupted values to NA
				pred.quantiles$truth[corrupted_rows] <- NA_real_
				pred.means$truth[corrupted_rows] <- NA_real_
			}
		}

		# Safety check before copying Mean/SD
		stopifnot(identical(pred.means$plot_num, pred.quantiles$plot_num) &&
		          identical(pred.means$timepoint, pred.quantiles$timepoint))
		pred.quantiles$Mean <- pred.means$Mean
		pred.quantiles$SD   <- pred.means$SD
	} else {
		# Fallback for other structures
		message("  Unknown plot_summary structure - creating empty plot estimates")
		pred.quantiles <- data.frame(
			plot_num = integer(0),
			timepoint = integer(0),
			Mean = numeric(0),
			SD = numeric(0),
			`2.5%` = numeric(0),
			`25%` = numeric(0),
			`50%` = numeric(0),
			`75%` = numeric(0),
			`97.5%` = numeric(0),
			taxon = character(0)
		)
		pred.means <- pred.quantiles
	}

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
					 pretty_group = ifelse(!is.na(group) & !is.na(group %in% c("16S","bac")) & group %in% c("16S","bac"), "Bacteria", "Fungi"),
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
						 pretty_group = ifelse(!is.na(group) & !is.na(group %in% c("16S","bac")) & group %in% c("16S","bac"), "Bacteria", "Fungi"),
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
			pretty_group = ifelse(!is.na(group) & !is.na(group %in% c("16S","bac")) & group %in% c("16S","bac"), "Bacteria", "Fungi"),
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
	message("  DEBUG: Created output list with", length(out), "elements")
	if (!is.null(save_summary) && save_summary == TRUE) {
		savePath <- gsub("samples","summary",file_path)
		message("  DEBUG: About to save to", savePath)
		saveRDS(out, savePath)
		message("Saved summary to ", savePath)
		return(TRUE)
	} else {
		message("  DEBUG: Not saving, returning output list")
		return(out)
	}
}