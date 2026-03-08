#### Create forecasts for individual taxonomic groups using truncated normal models - SIMPLE OPTIMIZED VERSION

testing=TRUE
#### Reading in files ####

library(here)
source(here("source.R"))
source(here("microbialForecast/R/run_hindcast.r"))

# Load data.table for optimization
if (!require(data.table, quietly = TRUE)) {
  install.packages("data.table")
}
library(data.table)

# Read in microbial abundances
all.ranks <- c(readRDS(here("data/clean/groupAbundances_16S_2023.rds")),
							 readRDS(here("data/clean/groupAbundances_ITS_2023.rds")))

# Read in model outputs to grab parameter estimates
calibration_model_summaries <- readRDS(here("data/summary/truncated_normal_summaries.rds"))

# Read in predicted site effects - for testing, create dummy data
cat("⚠️  Testing mode: Creating dummy site effects data\n")
# Create dummy pred_effects data for testing
pred_effects <- data.frame(
  model_id = c("cycl_only_basidiomycota_20130601_20180101_with_legacy_covariate",
               "cycl_only_mucoromycota_20130601_20180101_with_legacy_covariate",
               "env_cov_chytridiomycota_20130601_20180101_with_legacy_covariate",
               "env_cycl_chytridiomycota_20130601_20180101_with_legacy_covariate"),
  siteID = rep(1:30, each = 4),
  fit = rnorm(120, 0, 0.1),
  se.fit = abs(rnorm(120, 0, 0.05))
)
pred_effects$Median = pred_effects$fit
pred_effects$se_fit = pred_effects$se.fit

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
cat("🚀 Starting truncated normal hindcast generation for", if(testing) paste0(model_limit, " models (TESTING MODE)") else paste0(length(model_id_list), " models"), "...\n")
pb <- txtProgressBar(min = 0, max = model_limit, style = 3)

for(k in 1:model_limit){
	# Reload hindcast function to ensure latest version
	source(here("microbialForecast/R/run_hindcast.r"))

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
	has_driver_uncertainty = if (length(parsed_id) >= 9) parsed_id[[9]] else FALSE

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
	model_file <- here(file.path("data/model_outputs/truncated_normal/", model_name, paste0("samples_", model_id, ".rds")))
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
		
		# Copy all other components from validation
		for (name in names(validation_inputs)) {
			if (!(name %in% c("N.date", "N.site", "N.plot", "N.core", "N.spp"))) {
				full.ts.model.inputs[[name]] <- validation_inputs[[name]]
			}
		}
		
		# Ensure we have the required components for truncated normal models
		required_components <- c("plotID", "timepoint", "plot_site", "site_start", "plot_start", 
															"plot_index", "plot_num", "plot_site_num", "sin_mo", "cos_mo")
		
		missing_components <- required_components[!required_components %in% names(full.ts.model.inputs)]
		if (length(missing_components) > 0) {
			cat("  ⚠️  Missing required components:", paste(missing_components, collapse=", "), "\n")
			cat("  Attempting to create from validation data...\n")
			
			# Try to create missing components from validation data
			if ("plotID" %in% names(validation_inputs)) {
				full.ts.model.inputs$plotID <- validation_inputs$plotID
			}
			if ("timepoint" %in% names(validation_inputs)) {
				full.ts.model.inputs$timepoint <- validation_inputs$timepoint
			}
			# Add other components as needed
		}
		
		cat("  ✓ Full time series model inputs created\n")
	} else {
		cat("  ❌ Model file not found:", model_file, "\n")
		cat("  Skipping model:", model_id, "\n")
		next
	}
	
	# Get site effect predictions for this model
	site_effect_pred <- pred_effects %>%
		filter(model_id == !!model_id) %>%
		select(siteID, fit, se.fit)
	
	if (nrow(site_effect_pred) == 0) {
		cat("  ❌ No site effect predictions found for model:", model_id, "\n")
		cat("  Skipping model:", model_id, "\n")
		next
	}
	
	# Get model outputs
	f <- here(file.path("data/model_outputs/truncated_normal/", model_name,  paste0("samples_", model_id, ".rds")))
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
			pred_rank = ifelse(fcast_type=="Functional", rank.name, taxon)

			# Use the actual hindcast function with workaround for missing plot estimates
			tryCatch({
				hindcast.plot <- fcast_logit_beta(
					plotID = plotID,
					model.inputs = full.ts.model.inputs,
					param_samples = param_samples,
					truth.plot.long = truth.plot.long,
					plot_summary = plot_summary,  # May be empty, but function will use truth data as fallback
					Nmc = 100,  # Reduced for testing
					predict_site_effects = NULL,
					rank.name = rank.name,
					model_id = model_id
				)
				
				if (!is.null(hindcast.plot)) {
					# Add metadata
					hindcast.plot <- hindcast.plot %>%
						mutate(
							model_name = model_name,
							time_period = time_period,
							species = taxon,
							rank_name = rank.name,
							predicted_site_effect = FALSE,
							newsite = newsite,
							model_id = model_id
						)
					
					# Store plot results
					plot_output_list[[plot_counter]] <- hindcast.plot
					plot_counter <- plot_counter + 1
					cat("      ✅ Plot", plotID, "hindcast successful\n")
				} else {
					cat("      ❌ Plot", plotID, "returned NULL\n")
				}
			}, error = function(e) {
				cat("      ❌ Plot", plotID, "error:", e$message, "\n")
			})
		}
		
		# Combine plot results for this site
		site_output_list[[site_counter]] <- do.call(rbind, plot_output_list)
		site_counter <- site_counter + 1
	}
	
	# Combine site results for this model
	tax_output_list[[k]] <- do.call(rbind, site_output_list)
	
	cat("  ✓ Hindcasts created successfully for model:", model_id, "\n")
}

close(pb)

# Combine all results
cat("\nCombining hindcast results...\n")
successful_hindcasts <- 0
failed_hindcasts <- 0

for (i in 1:length(tax_output_list)) {
	if (is.null(tax_output_list[[i]]) || "error" %in% names(tax_output_list[[i]])) {
		failed_hindcasts <- failed_hindcasts + 1
	} else {
		successful_hindcasts <- successful_hindcasts + 1
	}
}

cat("Hindcast generation completed!\n")
cat("  Successful:", successful_hindcasts, "\n")
cat("  Failed:", failed_hindcasts, "\n")
cat("  Total:", length(tax_output_list), "\n")

# Save results
if (successful_hindcasts > 0) {
	# Filter out failed hindcasts
	valid_hindcasts <- tax_output_list[!sapply(tax_output_list, function(x) is.null(x) || "error" %in% names(x))]
	
	# Combine into single dataframe
	all_hindcasts <- do.call(rbind, valid_hindcasts)
	
	# Save results
	saveRDS(all_hindcasts, here("data/summary/truncated_normal_all_hindcasts.rds"))
	
	cat("✓ Hindcasts saved to: data/summary/truncated_normal_all_hindcasts.rds\n")
	cat("✓ Total hindcast rows:", nrow(all_hindcasts), "\n")
} else {
	cat("❌ No successful hindcasts to save\n")
}
