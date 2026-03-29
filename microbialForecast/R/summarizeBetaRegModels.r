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

  # CRITICAL FIX: Check for plot_mu columns before calculating
  plot_mu_cols <- grep("^plot_mu\\[|^mu\\[", colnames(combined_samples), value = TRUE)
  if (length(plot_mu_cols) == 0) {
    stop("CRITICAL: No plot_mu or mu columns found in samples2. Cannot calculate plot summary. ",
         "This indicates the model output does not contain plot-level predictions. ",
         "Available columns (first 20): ", paste(head(colnames(combined_samples), 20), collapse=", "))
  }
  
  calc_time <- Sys.time()
  message("  TIMING: Starting quantile calculations on ", ncol(combined_samples), " columns, ", nrow(combined_samples), " rows")
  
  means  <- colMeans(combined_samples, na.rm = TRUE)
  means_time <- Sys.time()
  message("  TIMING: Means calculated (", round(as.numeric(means_time - calc_time, units="secs"), 2), "s)")
  
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    sds <- matrixStats::colSds(combined_samples, na.rm = TRUE)
  } else {
    # Fallback: optimized calculation using colMeans and matrix operations
    # SD = sqrt(mean((x - mean(x))^2))
    means_mat <- matrix(means, nrow = nrow(combined_samples), ncol = ncol(combined_samples), byrow = TRUE)
    sds <- sqrt(colMeans((combined_samples - means_mat)^2, na.rm = TRUE))
  }
  sds_time <- Sys.time()
  message("  TIMING: SDs calculated (", round(as.numeric(sds_time - means_time, units="secs"), 2), "s)")
  
  probs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  n_cols <- ncol(combined_samples)
  
  # For large matrices, process quantiles in chunks to avoid memory issues
  chunk_size <- 1000
  if (n_cols > chunk_size) {
    message("  TIMING: Processing quantiles in chunks of ", chunk_size, " columns")
    n_chunks <- ceiling(n_cols / chunk_size)
    quants_list <- vector("list", n_chunks)
    
    for (chunk_idx in 1:n_chunks) {
      start_col <- (chunk_idx - 1) * chunk_size + 1
      end_col <- min(chunk_idx * chunk_size, n_cols)
      chunk_cols <- start_col:end_col
      
      quants_list[[chunk_idx]] <- apply(combined_samples[, chunk_cols, drop = FALSE], 2, quantile,
                                        probs = probs, na.rm = TRUE)
    }
    quants <- do.call(cbind, quants_list)
    colnames(quants) <- colnames(combined_samples)
  } else {
    # For smaller matrices, use standard apply
    quants <- apply(combined_samples, 2, quantile, probs = probs, na.rm = TRUE)
  }
  
  quants_time <- Sys.time()
  message("  TIMING: Quantiles calculated (", round(as.numeric(quants_time - sds_time, units="secs"), 2), "s)")
  
  # Validate that quantiles are not all NA
  all_na_cols <- sapply(1:ncol(quants), function(i) all(is.na(quants[, i])))
  if (any(all_na_cols)) {
    na_col_names <- colnames(combined_samples)[all_na_cols]
    warning("CRITICAL: Found ", sum(all_na_cols), " columns with all-NA quantiles: ", 
            paste(head(na_col_names, 5), collapse=", "),
            ". This suggests the MCMC samples for these parameters are invalid.")
  }

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
	start_time <- Sys.time()
	message("=== SUMMARIZE_BETA_MODEL FUNCTION CALLED ===")
	require(stringr)
	require(tidyr)
	# Helper functions (helperFunctions.r, summarizeModels.r) are loaded via package namespace
	if(summary_exists(file_path)) { # checks that a summary is needed (samples files are new)
		if (is.null(overwrite) || overwrite == FALSE) {
			return("Summary file already exists")
		}
	}
	# Read in file, assign named contents to global environment
	read_time <- Sys.time()
	read_in <- readRDS(file_path)
	message("File loaded successfully (", round(as.numeric(Sys.time() - read_time, units="secs"), 2), "s)")
	message("  DEBUG: File loaded, checking structure...")
	#list2env(read_in,globalenv())

	# Read in samples
	samples <- read_in$samples
	samples2 <- read_in$samples2  # Raw plot estimates from step 02
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
	
	# Convert summary.mcmc to list format if needed (plot_summary created by step 02)
	if (inherits(plot_summary, "summary.mcmc")) {
		message("  Converting summary.mcmc plot_summary to list format")
		plot_summary <- list(
			plot_summary$statistics,
			plot_summary$quantiles
		)
		names(plot_summary) <- c("statistics", "quantiles")
	}
	
	# Check if plot_summary exists but is invalid (all NA quantiles or wrong structure)
	plot_summary_valid <- FALSE
	if (!is.null(plot_summary) && is.list(plot_summary) && length(plot_summary) >= 2) {
		# Check if plot_summary has valid quantile data
		if (is.data.frame(plot_summary[[2]]) && nrow(plot_summary[[2]]) > 0) {
			# Check if quantile columns exist and have non-NA values
			quant_cols <- c("2.5%", "25%", "50%", "75%", "97.5%")
			has_quant_cols <- any(quant_cols %in% names(plot_summary[[2]]))
			if (has_quant_cols) {
				# Check if any quantiles are non-NA
				med_col <- if("50%" %in% names(plot_summary[[2]])) plot_summary[[2]][["50%"]] else NULL
				if (!is.null(med_col) && sum(is.finite(med_col)) > 0) {
					plot_summary_valid <- TRUE
					message("  plot_summary is valid with", sum(is.finite(med_col)), "finite median values")
				} else {
					message("  WARNING: plot_summary exists but all quantiles are NA - will try to regenerate from samples2")
				}
			} else {
				message("  WARNING: plot_summary exists but missing quantile columns - will try to regenerate from samples2")
			}
		}
	}
	
	# If plot_summary is invalid or missing, try to create it from samples2
	if (!plot_summary_valid && !is.null(samples2)) {
		message("  Attempting to regenerate plot_summary from samples2...")
		tryCatch({
			# Check if samples2 has plot_mu columns
			samples2_matrix <- NULL
			if (is.matrix(samples2)) {
				samples2_matrix <- samples2
			} else if (is.list(samples2) && length(samples2) > 0) {
				# Try to extract matrix from list (could be mcmc.list or list of matrices)
				if (inherits(samples2, "mcmc.list")) {
					samples2_matrix <- as.matrix(do.call(rbind, samples2))
				} else if (is.matrix(samples2[[1]])) {
					samples2_matrix <- do.call(rbind, lapply(samples2, as.matrix))
				}
			}
			
			if (!is.null(samples2_matrix) && ncol(samples2_matrix) > 0) {
				col_names <- colnames(samples2_matrix)
				if (!is.null(col_names)) {
					# Check for plot_mu or mu columns
					has_plot_mu <- any(grepl("^plot_mu\\[", col_names)) || any(grepl("^mu\\[", col_names))
					if (has_plot_mu) {
						message("  Found plot_mu/mu columns in samples2 - generating plot_summary")
						plot_summary <- calculate_plot_summary_from_samples(samples2_matrix)
						plot_summary_valid <- TRUE
						message("  Successfully regenerated plot_summary from samples2")
					} else {
						message("  WARNING: samples2 does not contain plot_mu or mu columns")
						message("    samples2 columns (first 20):", paste(head(col_names, 20), collapse=", "))
						message("    Cannot regenerate plot_summary - will use empty structure")
					}
				}
			}
		}, error = function(e) {
			message("  ERROR regenerating plot_summary from samples2:", e$message)
		})
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
										"all_covariates" = microbialForecast:::env_cycl_covariates_key,
										"env_cov" = microbialForecast:::env_cov_covariates_key,
										"env_cycl" = microbialForecast:::env_cycl_covariates_key,
										"cycl_only" = microbialForecast:::cycl_only_key)
	message("  DEBUG: cov_key assigned, length:", length(cov_key))

	sites <- truth.plot.long %>% select(site_num, siteID) %>% unique()
	site_key <- sites[["siteID"]]
	names(site_key) <- sites[["site_num"]]
	
	# Create plot mapping from plot_num to plotID (truth.plot.long should always have both)
	# CRITICAL FIX: Ensure truth.plot.long has the required columns before creating plot_key
	if (!"plot_num" %in% names(truth.plot.long) || !"plotID" %in% names(truth.plot.long)) {
		stop("truth.plot.long is missing required columns (plot_num or plotID) - cannot create plot mapping")
	}
	
	plot_mapping <- truth.plot.long %>% 
		select(plot_num, plotID) %>% 
		unique() %>%
		filter(!is.na(plot_num) & !is.na(plotID))
	
	# CRITICAL FIX: Check if plot_mapping is empty
	if (nrow(plot_mapping) == 0) {
		stop("plot_mapping is empty - truth.plot.long has no valid plot_num/plotID pairs")
	}
	
	plot_key <- plot_mapping[["plotID"]]
	names(plot_key) <- as.character(plot_mapping[["plot_num"]])

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
	# CRITICAL FIX: Validate plot_summary structure before processing
	if (is.null(plot_summary)) {
		stop("CRITICAL: plot_summary is NULL - cannot create model summary. ",
		     "This indicates the samples file is missing plot_summary or samples2 data.")
	}
	
	if (is.list(plot_summary) && length(plot_summary) >= 2) {
		# CRITICAL FIX: Check if plot_summary is empty
		if (nrow(plot_summary[[1]]) == 0 || nrow(plot_summary[[2]]) == 0) {
			stop("CRITICAL: plot_summary is empty (0 rows). ",
			     "This indicates the model has no plot-level predictions. ",
			     "Cannot create valid model summary without prediction data.")
		}
		
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
		
		# CRITICAL FIX: Validate plot_num consistency before merge
		# plot_num in plot_quantiles comes from plot_mu[plot_num, timepoint] parameter names
		# plot_num in truth.plot.long comes from match(plotID, names(plot_start))
		# They should match because both use the same plot_start vector from model fitting
		if (nrow(plot_quantiles) > 0 && "plot_num" %in% names(plot_quantiles)) {
			plot_quantiles_plot_nums <- unique(plot_quantiles$plot_num[!is.na(plot_quantiles$plot_num)])
			truth_plot_nums <- unique(truth.plot.long$plot_num[!is.na(truth.plot.long$plot_num)])
			
			# Check if all plot_quantiles plot_nums exist in truth.plot.long
			missing_in_truth <- setdiff(plot_quantiles_plot_nums, truth_plot_nums)
			if (length(missing_in_truth) > 0) {
				# Some plot_num values from model output don't exist in truth data
				# This indicates a data integrity issue - the model was fit with different plots
				warning("CRITICAL: Found ", length(missing_in_truth), " plot_num values in plot_quantiles (", 
				        paste(head(missing_in_truth, 5), collapse=", "), 
				        ") that don't exist in truth.plot.long. ",
				        "This suggests plot_num mismatch between model output and truth data. ",
				        "These rows will have NA for all truth columns after merge.")
			}
			
			# Check if truth has plot_nums not in plot_quantiles (less critical, but worth noting)
			missing_in_quantiles <- setdiff(truth_plot_nums, plot_quantiles_plot_nums)
			if (length(missing_in_quantiles) > 0 && length(missing_in_quantiles) > length(truth_plot_nums) * 0.1) {
				# More than 10% of truth plots missing from model output
				warning("Found ", length(missing_in_quantiles), " plot_num values in truth.plot.long (", 
				        paste(head(missing_in_quantiles, 5), collapse=", "), 
				        ") that don't exist in plot_quantiles. ",
				        "This suggests the model output is missing predictions for some plots.")
			}
		}
		
		# Merge with truth data, being careful about column name conflicts
		merge_time <- Sys.time()
		message("  TIMING: Starting merge - plot_quantiles: ", nrow(plot_quantiles), " rows, truth.plot.long: ", nrow(truth.plot.long), " rows")
		pred.quantiles <- plot_quantiles %>%
			merge(truth.plot.long, by = c("plot_num", "timepoint"), all = T, suffixes = c("", "_truth"))
		merge1_time <- Sys.time()
		message("  TIMING: pred.quantiles merge complete (", round(as.numeric(merge1_time - merge_time, units="secs"), 2), "s)")
		
		pred.means <- plot_means %>%
			merge(truth.plot.long, by = c("plot_num", "timepoint"), all = T, suffixes = c("", "_truth"))
		merge2_time <- Sys.time()
		message("  TIMING: pred.means merge complete (", round(as.numeric(merge2_time - merge1_time, units="secs"), 2), "s)")

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

		# CRITICAL FIX: Ensure dateID is properly handled after merge
		# Handle dateID_truth suffix if present
		if ("dateID_truth" %in% names(pred.quantiles)) {
			# Prefer dateID_truth (from truth.plot.long) over dateID (from plot_quantiles, if any)
			if (!"dateID" %in% names(pred.quantiles) || all(is.na(pred.quantiles$dateID))) {
				pred.quantiles$dateID <- pred.quantiles$dateID_truth
			} else {
				# Both exist - prefer truth, but fill NAs from plot_quantiles
				na_idx <- is.na(pred.quantiles$dateID)
				pred.quantiles$dateID[na_idx] <- pred.quantiles$dateID_truth[na_idx]
			}
			pred.quantiles$dateID_truth <- NULL
		}
		
		# CRITICAL FIX: Validate that timepoint maps uniquely to dateID
		# timepoint is GLOBAL (from prepBetaRegData) and should map 1:1 to dateID
		# If we have multiple dateID values for the same timepoint, this indicates a data integrity issue
		# We should STOP and investigate, not silently fix it
		validation_time <- Sys.time()
		if ("timepoint" %in% names(pred.quantiles) && "dateID" %in% names(pred.quantiles)) {
			message("  TIMING: Starting dateID validation on ", nrow(pred.quantiles), " rows")
			timepoint_dateid_check <- pred.quantiles %>%
				filter(!is.na(timepoint) & !is.na(dateID)) %>%
				group_by(timepoint) %>%
				summarise(n_dateids = n_distinct(dateID), .groups = "drop") %>%
				filter(n_dateids > 1)
			validation_done <- Sys.time()
			message("  TIMING: dateID validation complete (", round(as.numeric(validation_done - validation_time, units="secs"), 2), "s)")
			
			if (nrow(timepoint_dateid_check) > 0) {
				# Multiple dateIDs for same timepoint - this is a critical data integrity issue
				# Show which timepoints are affected with detailed diagnostics
				problematic_timepoints <- timepoint_dateid_check$timepoint
				dup_details <- pred.quantiles %>%
					filter(timepoint %in% problematic_timepoints & !is.na(dateID)) %>%
					select(timepoint, dateID, plot_num) %>%
					distinct() %>%
					arrange(timepoint, dateID)
				
				# Check if this came from plot_quantiles or truth.plot.long
				plot_quantiles_has_dateid <- "dateID" %in% names(plot_quantiles)
				truth_has_dateid <- "dateID" %in% names(truth.plot.long)
				
				# Print detailed diagnostic information
				message("CRITICAL DATA INTEGRITY ERROR: Found ", nrow(timepoint_dateid_check), 
				        " timepoint(s) with multiple dateID values.")
				message("  timepoint should map 1:1 to dateID (timepoint is GLOBAL from prepBetaRegData).")
				message("  This indicates the merge between plot_quantiles and truth.plot.long created incorrect mappings.")
				message("  Problematic timepoints: ", paste(head(problematic_timepoints, 10), collapse=", "))
				message("  Detailed mappings (timepoint -> dateID):")
				print(head(dup_details, 20))
				message("  plot_quantiles has dateID column: ", plot_quantiles_has_dateid)
				message("  truth.plot.long has dateID column: ", truth_has_dateid)
				message("  Number of rows in plot_quantiles: ", nrow(plot_quantiles))
				message("  Number of rows in truth.plot.long: ", nrow(truth.plot.long))
				message("  Number of unique (plot_num, timepoint) in plot_quantiles: ", 
				        ifelse("plot_num" %in% names(plot_quantiles) && "timepoint" %in% names(plot_quantiles),
				               nrow(plot_quantiles %>% select(plot_num, timepoint) %>% distinct()),
				               "N/A"))
				message("  Number of unique (plot_num, timepoint) in truth.plot.long: ",
				        ifelse("plot_num" %in% names(truth.plot.long) && "timepoint" %in% names(truth.plot.long),
				               nrow(truth.plot.long %>% select(plot_num, timepoint) %>% distinct()),
				               "N/A"))
				
				stop("CRITICAL DATA INTEGRITY ERROR: timepoint maps to multiple dateID values. ",
				     "This must be investigated and fixed upstream. ",
				     "Do not silently fix this - it indicates a serious problem with the merge or data structure. ",
				     "See diagnostic output above for details.")
			}
		}
		
		# Safety check before copying Mean/SD
		stopifnot(identical(pred.means$plot_num, pred.quantiles$plot_num) &&
		          identical(pred.means$timepoint, pred.quantiles$timepoint))
		pred.quantiles$Mean <- pred.means$Mean
		pred.quantiles$SD   <- pred.means$SD
		
		# Apply same dateID fix to pred.means
		if ("dateID_truth" %in% names(pred.means)) {
			if (!"dateID" %in% names(pred.means) || all(is.na(pred.means$dateID))) {
				pred.means$dateID <- pred.means$dateID_truth
			} else {
				na_idx <- is.na(pred.means$dateID)
				pred.means$dateID[na_idx] <- pred.means$dateID_truth[na_idx]
			}
			pred.means$dateID_truth <- NULL
		}
		
		# CRITICAL FIX: Validate pred.quantiles is not empty after merge
		if (nrow(pred.quantiles) == 0) {
			stop("CRITICAL: pred.quantiles is empty after merge. ",
			     "This indicates no valid plot-timepoint combinations exist. ",
			     "Cannot create valid model summary without prediction data.")
		}
		
		# CRITICAL FIX: Check if quantiles are all NA (indicates calculation failure)
		quantile_cols <- c("2.5%", "25%", "50%", "75%", "97.5%")
		quantile_cols_present <- intersect(quantile_cols, names(pred.quantiles))
		if (length(quantile_cols_present) > 0) {
			all_na_rows <- apply(pred.quantiles[, quantile_cols_present, drop = FALSE], 1, function(x) all(is.na(x)))
			if (all(all_na_rows)) {
				stop("CRITICAL: All quantile values are NA in pred.quantiles. ",
				     "This indicates the quantile calculation failed completely. ",
				     "Cannot create valid model summary with all-NA quantiles.")
			} else if (sum(all_na_rows) > nrow(pred.quantiles) * 0.5) {
				warning("CRITICAL: More than 50% of rows have all-NA quantiles (", 
				        sum(all_na_rows), " out of ", nrow(pred.quantiles), "). ",
				        "This suggests a serious issue with the quantile calculation.")
			}
		}
		
		# CRITICAL FIX: Ensure plotID is always present (truth.plot.long should have it)
		# The merge uses all=T (outer join) with suffixes = c("", "_truth")
		# Since plot_quantiles doesn't have plotID, plotID should come directly from truth.plot.long
		# However, when plot_num doesn't match, merge creates rows with NA for all truth columns
		# We need to ensure plotID is always present, using merge result first, then reconstruction
		
		# Step 1: Check if plotID exists from merge (should be present if plot_num matched)
		if (!"plotID" %in% names(pred.quantiles)) {
			# plotID not in pred.quantiles - check if it came from merge as plotID_truth
			if ("plotID_truth" %in% names(pred.quantiles)) {
				pred.quantiles$plotID <- pred.quantiles$plotID_truth
				pred.quantiles$plotID_truth <- NULL
			}
		}
		
		# Step 2: If plotID exists but has NAs, or doesn't exist, reconstruct from plot_num
		if (!"plotID" %in% names(pred.quantiles) || all(is.na(pred.quantiles$plotID))) {
			if (nrow(pred.quantiles) > 0 && length(plot_key) > 0 && "plot_num" %in% names(pred.quantiles)) {
				# Reconstruct plotID from plot_num using plot_key
				pred_plot_nums <- as.character(pred.quantiles$plot_num)
				matched_plotIDs <- plot_key[pred_plot_nums]
				
				# Check how many matched
				matched_count <- sum(!is.na(matched_plotIDs))
				unmatched_count <- sum(is.na(matched_plotIDs))
				
				if (unmatched_count > 0 && matched_count == 0) {
					# CRITICAL: No plot_num values match - this is a data integrity failure
					unmatched_plot_nums <- unique(pred_plot_nums[is.na(matched_plotIDs)])
					stop("CRITICAL DATA INTEGRITY ERROR: All plot_num values in plot_quantiles (", 
					     paste(head(unmatched_plot_nums, 5), collapse=", "), 
					     ") do not match plot_key from truth.plot.long. ",
					     "This indicates plot_num was assigned inconsistently between model output and truth data. ",
					     "Cannot proceed without valid plotID mapping.")
				} else if (unmatched_count > 0) {
					# Some plot_num values don't match - warn but use what we can
					unmatched_plot_nums <- unique(pred_plot_nums[is.na(matched_plotIDs)])
					warning("Found ", unmatched_count, " rows with plot_num values not in plot_key: ", 
					        paste(head(unmatched_plot_nums, 5), collapse=", "),
					        ". This indicates plot_num mismatch between model output and truth data.")
				}
				
				# Use matched plotIDs, filling NAs from merge if available
				if (!"plotID" %in% names(pred.quantiles)) {
					pred.quantiles$plotID <- matched_plotIDs
				} else {
					# plotID exists but is all NA - replace with matched values
					pred.quantiles$plotID <- matched_plotIDs
				}
				
				# Fill any remaining NAs from merge if available (shouldn't happen if merge worked)
				if (any(is.na(pred.quantiles$plotID)) && "plotID_truth" %in% names(pred.quantiles)) {
					na_idx <- is.na(pred.quantiles$plotID)
					pred.quantiles$plotID[na_idx] <- pred.quantiles$plotID_truth[na_idx]
					pred.quantiles$plotID_truth <- NULL
				}
			} else {
				# Cannot reconstruct - this is a critical error
				if (nrow(pred.quantiles) > 0) {
					stop("CRITICAL: Cannot reconstruct plotID - plot_key is empty or plot_num is missing. ",
					     "pred.quantiles has ", nrow(pred.quantiles), " rows but no way to assign plotID.")
				}
			}
		} else {
			# plotID exists but might have some NAs - fill them from plot_key if possible
			if (any(is.na(pred.quantiles$plotID)) && length(plot_key) > 0 && "plot_num" %in% names(pred.quantiles)) {
				na_idx <- is.na(pred.quantiles$plotID)
				pred_plot_nums <- as.character(pred.quantiles$plot_num[na_idx])
				matched_plotIDs <- plot_key[pred_plot_nums]
				if (sum(!is.na(matched_plotIDs)) > 0) {
					pred.quantiles$plotID[na_idx] <- matched_plotIDs
				}
				# Also check plotID_truth from merge
				if (any(is.na(pred.quantiles$plotID[na_idx])) && "plotID_truth" %in% names(pred.quantiles)) {
					still_na_idx <- na_idx & is.na(pred.quantiles$plotID)
					pred.quantiles$plotID[still_na_idx] <- pred.quantiles$plotID_truth[still_na_idx]
				}
			}
		}
		
		# Final validation: plotID must be present and non-empty
		if (!"plotID" %in% names(pred.quantiles) || (nrow(pred.quantiles) > 0 && all(is.na(pred.quantiles$plotID)))) {
			stop("CRITICAL: plotID is missing or all NA in pred.quantiles after all reconstruction attempts. ",
			     "Cannot create valid model summary without plotID.")
		}
		# Apply same plotID fix to pred.means
		if (!"plotID" %in% names(pred.means)) {
			if ("plotID_truth" %in% names(pred.means)) {
				pred.means$plotID <- pred.means$plotID_truth
				pred.means$plotID_truth <- NULL
			}
		}
		
		if (!"plotID" %in% names(pred.means) || all(is.na(pred.means$plotID))) {
			if (nrow(pred.means) > 0 && length(plot_key) > 0 && "plot_num" %in% names(pred.means)) {
				pred_plot_nums <- as.character(pred.means$plot_num)
				matched_plotIDs <- plot_key[pred_plot_nums]
				matched_count <- sum(!is.na(matched_plotIDs))
				unmatched_count <- sum(is.na(matched_plotIDs))
				
				if (unmatched_count > 0 && matched_count == 0) {
					stop("CRITICAL DATA INTEGRITY ERROR: All plot_num values in pred.means do not match plot_key.")
				}
				
				if (!"plotID" %in% names(pred.means)) {
					pred.means$plotID <- matched_plotIDs
				} else {
					pred.means$plotID <- matched_plotIDs
				}
				
				if (any(is.na(pred.means$plotID)) && "plotID_truth" %in% names(pred.means)) {
					na_idx <- is.na(pred.means$plotID)
					pred.means$plotID[na_idx] <- pred.means$plotID_truth[na_idx]
					pred.means$plotID_truth <- NULL
				}
			} else if (nrow(pred.means) > 0) {
				stop("CRITICAL: Cannot reconstruct plotID for pred.means.")
			}
		}
		
		if (!"plotID" %in% names(pred.means) || (nrow(pred.means) > 0 && all(is.na(pred.means$plotID)))) {
			stop("CRITICAL: plotID is missing or all NA in pred.means after all reconstruction attempts.")
		}
		
		# CRITICAL FIX: Ensure siteID is always present (truth.plot.long should have it)
		# If missing after merge, reconstruct from site_num using site_key, or extract from plotID
		if (!"siteID" %in% names(pred.quantiles)) {
			if (nrow(pred.quantiles) > 0) {
				# First try to reconstruct from site_num using site_key
				if ("site_num" %in% names(pred.quantiles) && length(site_key) > 0) {
					pred.quantiles$siteID <- site_key[as.character(pred.quantiles$site_num)]
					# Check if reconstruction worked (some site_nums might not be in site_key)
					if (any(is.na(pred.quantiles$siteID)) && "plotID" %in% names(pred.quantiles)) {
						# Fill missing values by extracting from plotID
						missing_idx <- is.na(pred.quantiles$siteID)
						pred.quantiles$siteID[missing_idx] <- substr(pred.quantiles$plotID[missing_idx], 1, 4)
					}
				} else if ("plotID" %in% names(pred.quantiles)) {
					# Fallback: extract siteID from plotID (first 4 chars)
					pred.quantiles$siteID <- substr(pred.quantiles$plotID, 1, 4)
				} else {
					pred.quantiles$siteID <- character(0)
				}
			} else {
				pred.quantiles$siteID <- character(0)
			}
		}
		if (!"siteID" %in% names(pred.means)) {
			if (nrow(pred.means) > 0) {
				# First try to reconstruct from site_num using site_key
				if ("site_num" %in% names(pred.means) && length(site_key) > 0) {
					pred.means$siteID <- site_key[as.character(pred.means$site_num)]
					# Check if reconstruction worked (some site_nums might not be in site_key)
					if (any(is.na(pred.means$siteID)) && "plotID" %in% names(pred.means)) {
						# Fill missing values by extracting from plotID
						missing_idx <- is.na(pred.means$siteID)
						pred.means$siteID[missing_idx] <- substr(pred.means$plotID[missing_idx], 1, 4)
					}
				} else if ("plotID" %in% names(pred.means)) {
					# Fallback: extract siteID from plotID (first 4 chars)
					pred.means$siteID <- substr(pred.means$plotID, 1, 4)
				} else {
					pred.means$siteID <- character(0)
				}
			} else {
				pred.means$siteID <- character(0)
			}
		}
	} else {
		# Fallback for other structures
		message("  Unknown plot_summary structure - creating empty plot estimates")
		pred.quantiles <- data.frame(
			plot_num = integer(0),
			timepoint = integer(0),
			plotID = character(0),  # CRITICAL: Include plotID column even when empty
			siteID = character(0),  # CRITICAL: Include siteID column even when empty
			Mean = numeric(0),
			SD = numeric(0),
			`2.5%` = numeric(0),
			`25%` = numeric(0),
			`50%` = numeric(0),
			`75%` = numeric(0),
			`97.5%` = numeric(0),
			taxon = character(0),
			stringsAsFactors = FALSE
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
		# Use parameter names that cover both beta regression and truncated normal models:
		# "precision" = cloglog beta precision; "core_sd" = truncnorm observation noise;
		# "sigma$" = truncnorm process noise (anchored to avoid matching sigma_proc etc.)
		eff_list <- lapply(c("precision", "core_sd", "sigma$", "intercept", "rho",
		                     "legacy_effect", "site_effect_sd"),
											 function(x) extract_summary_row(param_summary[[1]], var = x)) %>%
			plyr::rbind.fill() %>%
			mutate(taxon = !!species)
		eff_list2 <- lapply(c("precision", "core_sd", "sigma$", "intercept", "rho",
		                      "legacy_effect", "site_effect_sd"),
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
		save_time <- Sys.time()
		message("  DEBUG: About to save to", savePath)
		saveRDS(out, savePath)
		save_done <- Sys.time()
		message("Saved summary to ", savePath, " (", round(as.numeric(save_done - save_time, units="secs"), 2), "s)")
		total_time <- Sys.time()
		message("  TIMING: Total function time: ", round(as.numeric(total_time - start_time, units="secs"), 2), "s")
		return(TRUE)
	} else {
		message("  DEBUG: Not saving, returning output list")
		return(out)
	}
}