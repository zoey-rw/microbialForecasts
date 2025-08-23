#### Create forecasts for individual taxonomic groups - SIMPLE OPTIMIZED VERSION

testing=TRUE
#### Reading in files ####

source("source.R")
source("microbialForecast/R/run_hindcast.r")

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
	source("microbialForecast/R/run_hindcast.r")

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
	full.ts.model.inputs <- prepBetaRegData(rank.df = rank.df,
																							min.prev = 0,
																							min.date = min.date,
																							max.date = "20200101",
																							full_timeseries = T,
																							keep_vec = keep_vec)

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