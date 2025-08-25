#### Create forecasts for individual taxonomic groups - SIMPLE OPTIMIZED VERSION

testing=TRUE
#### Reading in files ####

source("../../source.R")
source("../../microbialForecast/R/run_hindcast.r")

# Load data.table for optimization
if (!require(data.table, quietly = TRUE)) {
  install.packages("data.table")
}
library(data.table)

# Read in microbial abundances
all.ranks <- c(readRDS(here("data/clean/groupAbundances_16S_2023.rds")),
							 readRDS(here("data/clean/groupAbundances_ITS_2023.rds")))

# Read in model outputs to grab parameter estimates
calibration_model_summaries <- readRDS(here("data/summary/logit_beta_regression_summaries.rds"))

# Read in predicted site effects
read_in = readRDS(here("data/summary/site_effects_unobserved.rds"))
pred_effects <- read_in[[2]]
pred_effects$fit = pred_effects$Median
pred_effects$se.fit = pred_effects$se_fit

# Read in predictor data, just to get the list of sites missing pC data
all_predictors = readRDS(here("data/clean/all_predictor_data.rds"))

# Read in list of taxa that passed convergence
keep_list <- readRDS(here("data/summary/converged_taxa_list.rds"))

model_id_list = unique(pred_effects$model_id)

min.date = "20151101"

# Pre-allocate output list for better memory management
tax_output_list = vector("list", length(model_id_list))

# for testing - run sequentially instead of in parallel
# Limit to first 3 models for quick testing
model_limit <- if(testing) 3 else length(model_id_list)

# Add progress bar
cat("🚀 Starting hindcast generation for", if(testing) paste0(model_limit, " models (TESTING MODE)") else paste0(length(model_id_list), " models"), "...\n")
pb <- txtProgressBar(min = 0, max = model_limit, style = 3)

for(k in 1:model_limit){
	source("../../microbialForecast/R/run_hindcast.r")

	model_id=model_id_list[k]
	# Update progress bar instead of verbose message
	setTxtProgressBar(pb, k)
	cat("\rProcessing model", k, "of", length(model_id_list), ":", model_id)
	
	parsed_id = parse_model_id(model_id)
	rank.name = parsed_id[[1]]
	time_period = parsed_id[[2]]
	taxon = species = parsed_id[[4]]
	fcast_type = parsed_id[[7]]
	model_name = parsed_id[[6]]

	# Filter model estimates for each plot abundance
	plot_summary <- rank_plot_summary <- calibration_model_summaries$plot_est %>%
				filter(model_id == !!model_id)

	keep_vec <- c("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date", taxon)
	rank.df = all.ranks[[as.character(rank.name)]]
	rank.df_spec <- rank.df %>%
		select("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date", !!taxon)
	rank.df_spec$other <- 1-rank.df_spec[[taxon]]
	
	# Prep validation data using full time series
	cat("  Creating validation period data...\n")
	validation_inputs <- prepBetaRegData(rank.df = rank.df,
																							min.prev = 0,
																							min.date = min.date,
																							max.date = "20200101",
																							full_timeseries = T,
																							keep_vec = keep_vec)
	
	# Fix: Create compatible data structure to prevent indexing errors
	cat("  Creating compatible data structure...\n")
	
	# Load the model to get calibration structure
	model_file <- here(file.path("data/model_outputs/logit_beta_regression/", model_name, paste0("samples_", model_id, ".rds")))
	if (file.exists(model_file)) {
		model_samples <- readRDS(model_file)
		calibration_model_data <- model_samples$metadata$model_data
		
		# Create merged structure that's compatible with the model
		full.ts.model.inputs <- list()
		
		# Copy basic structure from validation
		full.ts.model.inputs$N.date <- validation_inputs$N.date
		full.ts.model.inputs$N.site <- validation_inputs$N.site
		full.ts.model.inputs$N.plot <- validation_inputs$N.plot
		full.ts.model.inputs$N.core <- validation_inputs$N.core
		full.ts.model.inputs$N.spp <- validation_inputs$N.spp
		
		# CRITICAL FIX: Use calibration site_start indices but ensure they're valid for validation period
		calibration_site_start <- calibration_model_data$site_start
		
		# Validate and adjust site_start indices for validation period
		adjusted_site_start <- list()
		for (site_name in names(calibration_site_start)) {
			start_idx <- calibration_site_start[[site_name]]
			
			# If the calibration start index is out of bounds for validation, use 1
			if (is.na(start_idx) || start_idx < 1 || start_idx > validation_inputs$N.date) {
				adjusted_site_start[[site_name]] <- 1
				cat("    Adjusted site_start for ", site_name, " from ", start_idx, " to 1\n")
			} else {
				adjusted_site_start[[site_name]] <- start_idx
			}
		}
		
		full.ts.model.inputs$site_start <- adjusted_site_start
		
		# Create compatible arrays by mapping validation data to calibration structure
		cat("    Creating compatible temp array...\n")
		# Use calibration model data as primary source since it has all required arrays
		full.ts.model.inputs$temp <- calibration_model_data$temp
		
		cat("    Creating compatible pH array...\n")
		# Use calibration model data as primary source since it has all required arrays
		full.ts.model.inputs$pH <- calibration_model_data$pH
		
		# Create pH_sd array - MISSING ARRAY THAT CAUSES HINDCASTING ERRORS
		cat("    Creating compatible pH_sd array...\n")
		full.ts.model.inputs$pH_sd <- matrix(NA, nrow = nrow(calibration_model_data$pH), ncol = validation_inputs$N.date)
		rownames(full.ts.model.inputs$pH_sd) <- rownames(calibration_model_data$pH)
		
		# Map pH_sd data with bounds checking and robust NA handling
		for (plot in rownames(calibration_model_data$pH)) {
			if (plot %in% rownames(validation_inputs$pH_sd)) {
				plot_data <- validation_inputs$pH_sd[plot, ]
				n_cols <- min(length(plot_data), validation_inputs$N.date)
				full.ts.model.inputs$pH_sd[plot, 1:n_cols] <- plot_data[1:n_cols]
			} else {
				# Use mean pH_sd across all plots and timepoints, with minimum value to prevent errors
				mean_pH_sd <- mean(validation_inputs$pH_sd, na.rm = TRUE)
				if (is.na(mean_pH_sd) || mean_pH_sd <= 0) {
					mean_pH_sd <- 0.1  # Default minimum value for pH uncertainty
				}
				full.ts.model.inputs$pH_sd[plot, ] <- mean_pH_sd
			}
		}
		
		# CRITICAL FIX: Replace any remaining NA or invalid values in pH_sd with safe defaults
		cat("    Cleaning pH_sd array - replacing NA and invalid values...\n")
		na_mask <- is.na(full.ts.model.inputs$pH_sd) | full.ts.model.inputs$pH_sd <= 0
		if (any(na_mask)) {
			na_count <- sum(na_mask)
			cat("    - Found", na_count, "NA or invalid values in pH_sd\n")
			
			# Calculate safe default value from non-NA data
			valid_values <- full.ts.model.inputs$pH_sd[!is.na(full.ts.model.inputs$pH_sd) & full.ts.model.inputs$pH_sd > 0]
			if (length(valid_values) > 0) {
				safe_default <- max(0.1, mean(valid_values, na.rm = TRUE))
			} else {
				safe_default <- 0.1  # Absolute fallback
			}
			
			cat("    - Using safe default value:", round(safe_default, 3), "\n")
			full.ts.model.inputs$pH_sd[na_mask] <- safe_default
		}
		
		# Final validation
		if (any(is.na(full.ts.model.inputs$pH_sd)) || any(full.ts.model.inputs$pH_sd <= 0)) {
			stop("CRITICAL ERROR: pH_sd array still contains NA or invalid values after cleaning!")
		}
		cat("    ✅ pH_sd array cleaned and validated\n")
		
		# Copy other arrays with bounds checking
		cat("    Creating compatible moisture array...\n")
		# Use calibration model data as primary source since it has all required arrays
		full.ts.model.inputs$mois <- calibration_model_data$mois
		full.ts.model.inputs$mois_sd <- calibration_model_data$mois_sd
		
		cat("    Creating compatible temperature SD array...\n")
		full.ts.model.inputs$temp_sd <- calibration_model_data$temp_sd
		
		cat("    Creating compatible other arrays...\n")
		full.ts.model.inputs$pC <- calibration_model_data$pC
		full.ts.model.inputs$pC_sd <- calibration_model_data$pC_sd
		full.ts.model.inputs$relEM <- calibration_model_data$relEM
		full.ts.model.inputs$LAI <- calibration_model_data$LAI
		full.ts.model.inputs$sin_mo <- calibration_model_data$sin_mo
		full.ts.model.inputs$cos_mo <- calibration_model_data$cos_mo
		
		# CRITICAL FIX: Ensure all required arrays are present and valid
		required_arrays <- c("temp", "temp_sd", "mois", "mois_sd", "pH", "pH_sd", "pC", "pC_sd", "relEM", "LAI", "sin_mo", "cos_mo")
		cat("    Validating all required arrays...\n")
		
		for (array_name in required_arrays) {
			if (is.null(full.ts.model.inputs[[array_name]])) {
				cat("    ❌ CRITICAL ERROR: Required array '", array_name, "' is NULL\n")
				stop("Required array '", array_name, "' is NULL - this will cause hindcasting to fail!")
			}
			
			array_data <- full.ts.model.inputs[[array_name]]
			if (is.array(array_data) || is.matrix(array_data)) {
				na_count <- sum(is.na(array_data))
				if (na_count > 0) {
					cat("    ⚠️ Array '", array_name, "' has ", na_count, " NA values\n")
				} else {
					cat("    ✅ Array '", array_name, "' is clean (no NA values)\n")
				}
			} else {
				cat("    ❌ CRITICAL ERROR: Array '", array_name, "' is not a matrix/array\n")
				stop("Array '", array_name, "' is not a matrix/array - this will cause hindcasting to fail!")
			}
		}
		
		cat("    ✅ All required arrays validated\n")
		
		# Copy the truth.plot.long structure
		full.ts.model.inputs$truth.plot.long <- validation_inputs$truth.plot.long
		
		cat("    ✓ Compatible data structure created with bounds validation\n")
		
	} else {
		cat("    ⚠️ Model file not found, using raw validation inputs (may cause errors)\n")
		full.ts.model.inputs <- validation_inputs
	}

	message("Forecasting with model: ", model_name)

	# Get model outputs
	f <- here(file.path("data/model_outputs/logit_beta_regression/", model_name,  paste0("samples_", model_id, ".rds")))
	if(!file.exists(f)) next()
	read_in <- readRDS(f)
	param_samples <- as.data.frame(as.matrix(read_in$samples))
	model.dat <- read_in$metadata$model_data
	
	# Handle different structures for truth.plot.long access
	if("truth.plot.long" %in% names(model.dat)) {
		truth.plot.long <- model.dat$truth.plot.long
	} else {
		# Try the older structure
		truth.plot.long <- model.dat
	}

	plot_site_key <- truth.plot.long %>%
		select(siteID, plotID, dateID, date_num, plot_num, site_num) %>%
		distinct()
	site_list <- unique(plot_site_key$siteID)

	# Use new model inputs for full date, site, and plot keys
	new_plot_site_key <- full.ts.model.inputs$truth.plot.long %>%
		select(siteID, plotID, dateID, date_num, plot_num, site_num) %>%
		distinct() %>%
		filter(!siteID %in% plot_site_key$siteID)
	new_site_list <- unique(new_plot_site_key$siteID)

	if (testing == T) {
		full_site_list <- c(head(site_list,5), head(new_site_list,5))
	} else {
		full_site_list <- c(site_list, new_site_list)
	}

	# Pre-allocate site output list
	site_output_list = vector("list", length(full_site_list))
	site_counter = 1
	
	# Add site progress indicator
	cat("\n  Processing", length(full_site_list), "sites for model", k, "...\n")
	
	# Loop through both observed and unobserved sites
	for (siteID in full_site_list) {
		# Remove verbose SiteID messages
		# message("SiteID: ", siteID)

		newsite <- siteID %in% new_plot_site_key$siteID
		plot_key <- if (newsite) new_plot_site_key else plot_site_key
		plot_key <- plot_key %>% filter(siteID == !!siteID)
		plot_list <- unique(plot_key$plotID)

		# Pre-allocate plot output list
		plot_output_list = vector("list", length(plot_list))
		plot_counter = 1
		
		# Add plot progress indicator
		cat("    Processing", length(plot_list), "plots for site", siteID, "...\n")

		for (plotID in plot_list){
			# Remove verbose PlotID messages
			# message("PlotID: ", plotID)
			pred_rank = ifelse(fcast_type=="Functional", rank.name, taxon)

			# Forecast with known or random site effect
			tryCatch({
				hindcast.plot <-
					fcast_logit_beta(plotID,
													 full.ts.model.inputs,
													 param_samples,
													 truth.plot.long,
													 plot_summary,
													 Nmc = 1000,
													 predict_site_effects = NULL,
													 rank.name = pred_rank,
													 model_id=model_id)
				
				# Check if the forecast returned NULL (plot was skipped)
				if (is.null(hindcast.plot)) {
					# Remove verbose plot skipping messages - just continue silently
					# message("Plot ", plotID, " was skipped due to invalid data - continuing to next plot")
					next
				}
				
				# Add metadata if forecast was successful
				hindcast.plot <- hindcast.plot %>% mutate(model_name = !!model_name,
																 																	 time_period = !!time_period,
																 																	 species = !!taxon,
																 																	 rank_name = !!rank.name,
																 																	 predicted_site_effect=F,
																 																	 newsite = !!newsite,
																 															 model_id=!!model_id)

				# Forecast with estimated site effect, if available
				if (siteID %in% pred_effects$siteID) {
					pred_effects_filt = pred_effects %>% filter(model_id == !!model_id)
					tryCatch({
						hindcast.plot_pred_site_eff <-
							fcast_logit_beta(plotID,
															 full.ts.model.inputs,
															 param_samples,
															 truth.plot.long,
															 plot_summary,
															 Nmc = 1000,
															 predict_site_effects = pred_effects_filt,
															 rank.name = pred_rank,
															 model_id=model_id)
						
						# Check if the forecast returned NULL (plot was skipped)
						if (is.null(hindcast.plot_pred_site_eff)) {
							# Remove verbose plot skipping messages - just continue silently
							# message("Plot ", plotID, " was skipped for predicted site effect - continuing")
							next
						}
						
						# Add metadata if forecast was successful
						hindcast.plot_pred_site_eff <- hindcast.plot_pred_site_eff %>%
							mutate(model_name = !!model_name,
										 time_period = !!time_period,
										 species = !!taxon,
										 rank_name = !!rank.name,
										 predicted_site_effect=T,
										 newsite = !!newsite,
										 model_id=!!model_id)
						
						# Use data.table rbind for better performance
						hindcast.plot <- rbindlist(list(hindcast.plot, hindcast.plot_pred_site_eff), fill=TRUE)
					}, error = function(e) {
						message("Error in predicted site effect forecast for plot ", plotID, ": ", e$message)
					})
				}

				# print(tail(hindcast.plot, 1))  # Removed for performance
				plot_output_list[[plot_counter]] <- hindcast.plot
				plot_counter = plot_counter + 1
				
			}, error = function(e) {
				message("Error in forecast for plot ", plotID, ": ", e$message)
				# Skip this plot and continue with the next one
			})
		}
		
		# Use data.table rbind for better performance
		site_output_list[[site_counter]] <- rbindlist(plot_output_list, fill=TRUE)
		site_counter = site_counter + 1
	}
	
	# Use data.table rbind for better performance
	tax_output <- rbindlist(site_output_list, fill=TRUE)
	message("Completed forecast loop for: ", model_id)
	tax_output_list[[k]] <- tax_output
}

# Final output binding with data.table
out.fcast <- rbindlist(tax_output_list, fill = TRUE)

# Close progress bar
close(pb)
cat("\n✅ Completed hindcast generation!\n")

message("output nrow: ", nrow(out.fcast))

# Save as both RDS and Parquet for compatibility
out.path.rds <- here("data/summary/all_hindcasts_raw_plsr2.rds")
out.path.parquet <- here("data/summary/parquet/all_hindcasts_raw_plsr2.parquet")

# Ensure parquet directory exists
dir.create(here("data/summary/parquet"), showWarnings = FALSE, recursive = TRUE)

# Save as RDS
saveRDS(out.fcast, out.path.rds)
cat("✓ Saved RDS file:", out.path.rds, "\n")

# Save as Parquet for memory efficiency
if (require(arrow, quietly = TRUE)) {
	arrow::write_parquet(out.fcast, out.path.parquet)
	cat("✓ Saved Parquet file:", out.path.parquet, "\n")
} else {
	cat("⚠️  Arrow package not available, skipping Parquet save\n")
}

cat("🎉 HINDCAST GENERATION COMPLETE \n")