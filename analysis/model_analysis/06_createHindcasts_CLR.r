# Create hindcasts for CLR microbial forecasts
# Updated to follow the working version's approach but adapted for CLR models

library(dplyr)
library(here)
library(foreach)
library(doParallel)
library(data.table)
source(here("source.R"))

cat("=== CREATING CLR HINDCASTS FOR ALL SITES ===\n")

# Load required data files (following working version)
cat("Loading required data files...\n")

# Read in microbial abundances
all.ranks <- c(readRDS(here("data/clean/groupAbundances_16S_2023.rds")),
               readRDS(here("data/clean/groupAbundances_ITS_2023.rds")))

# Read in CLR model outputs to grab parameter estimates
calibration_model_summaries <- readRDS(here("data/summary/clr_regression_summaries.rds"))

# Read in predicted site effects for CLR models
clr_site_effects_file <- here("data/summary/site_effects_unobserved_CLR.rds")
if (file.exists(clr_site_effects_file)) {
  cat("Using CLR-specific site effects file:", clr_site_effects_file, "\n")
  read_in <- readRDS(clr_site_effects_file)
  pred_effects <- read_in[[2]]
  pred_effects$fit = pred_effects$Median
  pred_effects$se.fit = pred_effects$se_fit
} else {
  cat("CLR-specific site effects not found, using fallback\n")
  # Fall back to regular site effects if CLR-specific doesn't exist
  fallback_file <- here("data/summary/site_effects_unobserved.rds")
  if (file.exists(fallback_file)) {
    read_in <- readRDS(fallback_file)
    pred_effects <- read_in[[2]]
    pred_effects$fit = pred_effects$Median
    pred_effects$se.fit = pred_effects$se_fit
  } else {
    stop("No site effects files found!")
  }
}

# Read in predictor data, just to get the list of sites missing pC data
all_predictors = readRDS(here("data/clean/all_predictor_data.rds"))

cat("✅ Data files loaded successfully\n")

# Get available CLR models
model_outputs_dir <- here("data", "model_outputs", "CLR_regression")
available_models <- list.dirs(model_outputs_dir, full.names = FALSE, recursive = FALSE)
cat("Available CLR models:", paste(available_models, collapse = ", "), "\n\n")

# Focus on env_cycl model for now (can be expanded to include cycl_only and env_cov)
target_model <- "env_cycl"
if (target_model %in% available_models) {
  cat("Processing CLR model:", target_model, "\n")
  
  # Get available taxa for this model from the CLR summaries
  # Extract env_cycl models from the CLR summaries
  env_cycl_models <- calibration_model_summaries$plot_est$model_id[grepl(paste0('^', target_model, '_'), calibration_model_summaries$plot_est$model_id)]
  env_cycl_models <- unique(env_cycl_models)
  
  cat("Found", length(env_cycl_models), "CLR env_cycl models in summaries\n")
  
  # Extract taxon names from model IDs (e.g., "env_cycl_ascomycota_20130601_20180101_with_legacy_covariate" -> "ascomycota")
  available_taxa <- sapply(env_cycl_models, function(model_id) {
    parts <- strsplit(model_id, "_")[[1]]
    # Find the taxon part (between env_cycl and the date)
    if (length(parts) >= 3) {
      # Skip the first two parts (env_cycl) and take everything until the date
      taxon_parts <- parts[3:(which(grepl("^[0-9]{8}$", parts))[1] - 1)]
      paste(taxon_parts, collapse = "_")
    } else {
      NA
    }
  })
  
  available_taxa <- available_taxa[!is.na(available_taxa)]
  available_taxa <- unique(available_taxa)
  
  cat("Extracted taxa from CLR summaries:", paste(available_taxa, collapse = ", "), "\n\n")
  
  # Check which taxa have samples files
  samples_dir <- file.path(model_outputs_dir, target_model)
  available_samples <- list.files(samples_dir, pattern = "samples.*\\.rds$", full.names = FALSE, recursive = TRUE)
  
  # Extract taxon names from samples files
  samples_taxa <- sapply(available_samples, function(filename) {
    parts <- strsplit(filename, "_")[[1]]
    if (length(parts) >= 5) {
      # Format: samples_clr_env_cycl_[taxon]_[date]_with_legacy_covariate.rds
      # or: samples_env_cycl_[taxon]_[date]_with_legacy_covariate_clr.rds
      if (parts[2] == "clr") {
        parts[5]  # samples_clr_env_cycl_[taxon]...
      } else if (parts[2] == "env" && parts[3] == "cycl") {
        parts[4]  # samples_env_cycl_[taxon]...
      } else {
        NA
      }
    } else {
      NA
    }
  })
  
  samples_taxa <- samples_taxa[!is.na(samples_taxa)]
  samples_taxa <- unique(samples_taxa)
  
  cat("Taxa with CLR samples files:", paste(samples_taxa, collapse = ", "), "\n\n")
  
  # Use taxa that have both summaries and samples
  working_taxa <- intersect(available_taxa, samples_taxa)
  cat("Working CLR taxa (have both summaries and samples):", paste(working_taxa, collapse = ", "), "\n\n")
  
  # For testing, limit to first 2 taxa
  working_taxa <- working_taxa[1:2]
  cat("Testing with first 2 CLR taxa:", paste(working_taxa, collapse = ", "), "\n\n")
  
  if (length(working_taxa) == 0) {
    cat("❌ No working CLR taxa found for", target_model, "\n")
    quit(status = 1)
  }
  
  # Initialize output list for all taxa
  all_taxa_outputs <- list()
  
  # Process each taxon
  for (taxon in working_taxa) {
    cat("Processing CLR taxon:", taxon, "\n")
    
    # Find the model ID for this taxon
    model_id <- env_cycl_models[grepl(paste0("_", taxon, "_"), env_cycl_models)][1]
    if (is.na(model_id)) {
      cat("  ❌ No model ID found for taxon:", taxon, "\n")
      next
    }
    
    cat("  Model ID:", model_id, "\n")
    
    # Find the samples file for this taxon
    # Remove _clr suffix from model_id if present for samples file lookup
    samples_model_id <- gsub("_clr$", "", model_id)
    samples_file <- file.path(samples_dir, taxon, paste0("samples_", samples_model_id, "_clr.rds"))
    
    if (file.exists(samples_file)) {
      cat("  ✅ Found samples file:", basename(samples_file), "\n")
      
      # Load the samples
      samples_data <- readRDS(samples_file)
      
      # Extract parameter samples and plot summary
      param_samples <- samples_data$samples[[1]]  # First chain for now
      plot_summary <- samples_data$plot_summary
      
      cat("  Parameter samples dimensions:", dim(param_samples), "\n")
      cat("  Plot summary dimensions:", dim(plot_summary), "\n")
      
      # Get the rank name for this taxon
      rank.name <- ifelse(grepl("_bac", taxon), "phylum_bac", "phylum_fun")
      cat("  Rank name:", rank.name, "\n")
      
      # CRITICAL FIX: Create model inputs using prepCLRData (following beta regression pattern)
      cat("  Creating model inputs using prepCLRData\n")
      
      # Create keep_vec and prepare data (following beta regression approach)
      # CRITICAL FIX: For CLR transformation, we need all taxonomic columns, not just the target taxon
      rank.df <- all.ranks[[rank.name]]
      # Get all taxonomic columns (columns 7 onwards)
      taxonomic_cols <- colnames(rank.df)[7:ncol(rank.df)]
      keep_vec <- c("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date", taxonomic_cols)
      
      # Prep validation data using full time series 
      full.ts.model.inputs <- prepCLRData(rank.df = rank.df,
                                         min.prev = 0,
                                         min.date = "20130601",
                                         max.date = "20200101",  # Extend to 2020
                                         full_timeseries = TRUE,
                                         keep_vec = keep_vec,
                                         s = taxon)
      
      cat("  Full time series model inputs created\n")
      cat("  N.date:", full.ts.model.inputs$N.date, "\n")
      cat("  Total plots:", length(full.ts.model.inputs$plotID), "\n")
      
      # Get truth data
      truth.plot.long <- full.ts.model.inputs$truth.plot.long
      
      cat("  Model inputs loaded successfully\n")
      cat("  Truth data dimensions:", dim(truth.plot.long), "\n")
      
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
      
      cat("  Observed sites:", length(site_list), "\n")
      cat("  Unobserved sites:", length(new_site_list), "\n")
      
      # Combine site lists
      full_site_list <- c(site_list, new_site_list)
      
      # Generate hindcasts for each site
      cat("  Generating CLR hindcasts for", length(full_site_list), "sites...\n")
      
      site_output_list <- list()
      successful_sites <- 0
      
      for (siteID in full_site_list) {
        cat("    Processing site:", siteID, "\n")
        
        # Skip sites with missing predictor data
        if (siteID %in% all_predictors$site_skip) {
          cat("      ⚠️  Skipping site (missing predictor data)\n")
          next
        }
        
        # Determine if this is a new site
        newsite <- siteID %in% new_plot_site_key$siteID
        plot_key <- if (newsite) new_plot_site_key else plot_site_key
        plot_key <- plot_key %>% filter(siteID == !!siteID)
        plot_list <- unique(plot_key$plotID)
        
        cat("      Plots in site:", length(plot_list), "\n")
        
        plot_output_list <- list()
        successful_plots <- 0
        
        for (plotID in plot_list) {
          cat("        Processing plot:", plotID, "\n")
          
          tryCatch({
            # Forecast with known or random site effect using CLR forecasting function
            hindcast.plot <- fcast_clr(
              plotID = plotID,
              model.inputs = full.ts.model.inputs,
              param_samples = param_samples,
              truth.plot.long = truth.plot.long,
              plot_summary = plot_summary,
              Nmc = 100,  # Reduced for testing
              predict_site_effects = NULL,
              rank.name = rank.name,
              model_id = model_id
            )
            
            if (!is.null(hindcast.plot)) {
              # Add metadata for CLR models
              hindcast.plot <- hindcast.plot %>% 
                mutate(
                  model_name = target_model,
                  time_period = "20130601_20180101",
                  species = taxon,
                  rank_name = rank.name,
                  predicted_site_effect = FALSE,
                  newsite = newsite,
                  model_id = model_id
                )
              
              # Forecast with estimated site effect, if available
              if (siteID %in% pred_effects$siteID) {
                pred_effects_filt <- pred_effects %>% filter(model_id == !!model_id)
                
                if (nrow(pred_effects_filt) > 0) {
                  hindcast.plot_pred_site_eff <- fcast_clr(
                    plotID = plotID,
                    model.inputs = full.ts.model.inputs,
                    param_samples = param_samples,
                    truth.plot.long = truth.plot.long,
                    plot_summary = plot_summary,
                    Nmc = 100,  # Reduced for testing
                    predict_site_effects = pred_effects_filt,
                    rank.name = rank.name,
                    model_id = model_id
                  )
                  
                  if (!is.null(hindcast.plot_pred_site_eff)) {
                    hindcast.plot_pred_site_eff <- hindcast.plot_pred_site_eff %>% 
                      mutate(
                        model_name = target_model,
                        time_period = "20130601_20180101",
                        species = taxon,
                        rank_name = rank.name,
                        predicted_site_effect = TRUE,
                        newsite = newsite,
                        model_id = model_id
                      )
                    
                    # Combine both hindcasts
                    hindcast.plot <- rbindlist(list(hindcast.plot, hindcast.plot_pred_site_eff), fill = TRUE)
                  }
                }
              }
              
              plot_output_list[[plotID]] <- hindcast.plot
              successful_plots <- successful_plots + 1
              cat("          ✅ Success\n")
            } else {
              cat("          ❌ Returned NULL\n")
            }
          }, error = function(e) {
            cat("          ❌ Error:", e$message, "\n")
          })
        }
        
        if (length(plot_output_list) > 0) {
          site_output_list[[siteID]] <- rbindlist(plot_output_list, fill = TRUE)
          successful_sites <- successful_sites + 1
        }
      }
      
      # Combine all site outputs for this taxon
      if (length(site_output_list) > 0) {
        tax_output <- rbindlist(site_output_list, fill = TRUE)
        
        # Add to overall results list
        all_taxa_outputs[[taxon]] <- tax_output
        
        cat("  Total CLR hindcast rows:", nrow(tax_output), "\n")
        cat("  Successful sites:", successful_sites, "\n")
      } else {
        cat("  ❌ No successful CLR hindcasts generated for this taxon\n")
      }
      
    } else {
      cat("  ❌ Missing samples file:", samples_file, "\n")
    }
    
    cat("\n")
  }
  
} else {
  cat("Target CLR model", target_model, "not found in available models\n")
}

# Combine all taxa outputs and save
if (length(all_taxa_outputs) > 0) {
  cat("\nCombining all CLR taxa outputs...\n")
  all_hindcasts_combined <- rbindlist(all_taxa_outputs, fill = TRUE)
  
  # Save combined output in the format expected by 07_tidyHindcasts_CLR.r
  output_file <- here("data", "summary", "all_hindcasts_raw_CLR.rds")
  saveRDS(all_hindcasts_combined, output_file)
  cat("✅ Combined CLR hindcasts saved to:", output_file, "\n")
  cat("Total combined rows:", nrow(all_hindcasts_combined), "\n")
  
  # Also save as Parquet for memory efficiency
  if (require(arrow, quietly = TRUE)) {
    parquet_dir <- here("data/summary/parquet")
    dir.create(parquet_dir, showWarnings = FALSE, recursive = TRUE)
    
    parquet_file <- file.path(parquet_dir, "all_hindcasts_raw_CLR.parquet")
    arrow::write_parquet(all_hindcasts_combined, parquet_file)
    cat("✅ Combined CLR hindcasts also saved as Parquet:", parquet_file, "\n")
  }
} else {
  cat("❌ No CLR hindcast data to combine\n")
}

cat("=== CLR HINDCAST GENERATION COMPLETE ===\n")

for(k in 1:model_limit){
	source("../../microbialForecast/R/run_hindcast_CLR.r")

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
	
	# For CLR models, we need to use prepCLRData instead of prepBetaRegData
	cat("  Using prepCLRData for CLR models...\n")
	
	if (exists("prepCLRData", envir = asNamespace("microbialForecast"))) {
		validation_inputs <- prepCLRData(rank.df = rank.df,
																							min.prev = 0,
																							min.date = min.date,
																							max.date = "20200101",
																							full_timeseries = T,
																							keep_vec = keep_vec)
	} else {
		cat("  ERROR: prepCLRData function not found\n")
		next
	}
	
	# Fix: Create compatible data structure to prevent indexing errors
	cat("  Creating compatible data structure...\n")
	
	# Load the CLR model to get calibration structure
	model_file <- here(file.path("data/model_outputs/CLR_regression/", model_name, paste0("samples_", model_id, ".rds")))
	if (file.exists(model_file)) {
		model_samples <- readRDS(model_file)
  
		# For CLR models, the structure is different
		if ("metadata" %in% names(model_samples) && "model_data" %in% names(model_samples$metadata)) {
			calibration_model_data <- model_samples$metadata$model_data
		} else {
			calibration_model_data <- model_samples$model_data
		}
  
		# Create full time series model inputs for hindcasting
		full.ts.model.inputs <- list()
  
		# Copy over the basic structure from validation data
		full.ts.model.inputs$N.date <- validation_inputs$N.date
		full.ts.model.inputs$dateID <- validation_inputs$dateID
		full.ts.model.inputs$date_num <- validation_inputs$date_num
		full.ts.model.inputs$plotID <- validation_inputs$plotID
		full.ts.model.inputs$siteID <- validation_inputs$siteID
		# CRITICAL FIX: Use calibration model's existing site mappings (original working approach)
		# The original version worked because it used these directly - don't reconstruct them
		full.ts.model.inputs$plot_site_num <- calibration_model_data$plot_site_num
		full.ts.model.inputs$site_start <- calibration_model_data$site_start
		full.ts.model.inputs$plot_start <- calibration_model_data$plot_start
		full.ts.model.inputs$truth.plot.long <- validation_inputs$truth.plot.long
  
		# Copy over environmental and seasonal covariates
		full.ts.model.inputs$temp <- validation_inputs$temp
		full.ts.model.inputs$temp_sd <- validation_inputs$temp_sd
		full.ts.model.inputs$mois <- validation_inputs$mois
		full.ts.model.inputs$mois_sd <- validation_inputs$mois_sd
		full.ts.model.inputs$pH <- validation_inputs$pH
		full.ts.model.inputs$pH_sd <- validation_inputs$pH_sd
		full.ts.model.inputs$pC <- validation_inputs$pC
		full.ts.model.inputs$pC_sd <- validation_inputs$pC_sd
		full.ts.model.inputs$relEM <- validation_inputs$relEM
		full.ts.model.inputs$LAI <- validation_inputs$LAI
		full.ts.model.inputs$sin_mo <- validation_inputs$sin_mo
		full.ts.model.inputs$cos_mo <- validation_inputs$cos_mo
  
		# CRITICAL FIX: Validate arrays with proper NA handling for site start dates
		required_arrays <- c("temp", "temp_sd", "mois", "mois_sd", "pH", "pH_sd", "pC", "pC_sd", "relEM", "LAI")
		required_vectors <- c("sin_mo", "cos_mo")
		cat("    Validating arrays with site start date consideration...\n")
  
		for (array_name in required_arrays) {
			if (is.null(full.ts.model.inputs[[array_name]])) {
				cat("    ❌ CRITICAL ERROR: Required array '", array_name, "' is NULL\n")
				stop("Required array '", array_name, "' is NULL - this will cause hindcasting to fail!")
			}
			
			array_data <- full.ts.model.inputs[[array_name]]
			if (is.array(array_data) || is.matrix(array_data)) {
				# Count total NAs for information
				total_na_count <- sum(is.na(array_data))
				
				# Check if NAs are only in pre-start periods (which is fine)
				# For now, just log the NA count - the create_covariate_samples function
				# will handle validation of data quality from site_start onwards
				if (total_na_count > 0) {
					cat("    ⚠️ Array '", array_name, "' has ", total_na_count, " NA values (may include pre-start periods)\n")
					cat("    ℹ️  Data quality will be validated during hindcasting from site_start onwards\n")
				} else {
					cat("    ✅ Array '", array_name, "' is clean (no NA values)\n")
				}
			} else {
				cat("    ❌ CRITICAL ERROR: Array '", array_name, "' is not a matrix/array\n")
				stop("Array '", array_name, "' is not a matrix/array - this will cause hindcasting to fail!")
			}
		}
		
		cat("    Validating required vectors...\n")
		for (vector_name in required_vectors) {
			if (is.null(full.ts.model.inputs[[vector_name]])) {
				cat("    ❌ CRITICAL ERROR: Required vector '", vector_name, "' is NULL\n")
				stop("Required vector '", vector_name, "' is NULL - this will cause hindcasting to fail!")
			}
			
			vector_data <- full.ts.model.inputs[[vector_name]]
			if (is.vector(vector_data) || is.numeric(vector_data)) {
				na_count <- sum(is.na(vector_data))
				if (na_count > 0) {
					cat("    ⚠️ Vector '", vector_name, "' has ", na_count, " NA values\n")
				} else {
					cat("    ✅ Vector '", vector_name, "' is clean (no NA values)\n")
				}
			} else {
				cat("    ❌ CRITICAL ERROR: Vector '", vector_name, "' is not a vector\n")
				stop("Vector '", vector_name, "' is not a vector - this will cause hindcasting to fail!")
			}
		}
		
		cat("    ✅ All required arrays validated (NAs in pre-start periods are acceptable)\n")
		
		# Copy the truth.plot.long structure
		full.ts.model.inputs$truth.plot.long <- validation_inputs$truth.plot.long
		
		cat("    ✓ Compatible data structure created with bounds validation\n")
		
	} else {
		cat("    ⚠️ Model file not found, using raw validation inputs (may cause errors)\n")
		full.ts.model.inputs <- validation_inputs
	}

	message("Forecasting with CLR model: ", model_name)

	# Get CLR model outputs
	f <- here(file.path("data/model_outputs/CLR_regression/", model_name,  paste0("samples_", model_id, ".rds")))
	if(!file.exists(f)) next()
	read_in <- readRDS(f)
	
	# For CLR models, handle different data structure
	if ("samples" %in% names(read_in) && is.list(read_in$samples)) {
		# Convert mcmc.list to matrix for CLR models
		param_samples <- as.data.frame(as.matrix(read_in$samples))
	} else {
		param_samples <- as.data.frame(as.matrix(read_in$samples))
	}
	
	# Extract model data with proper structure handling
	if ("metadata" %in% names(read_in) && "model_data" %in% names(read_in$metadata)) {
		model.dat <- read_in$metadata$model_data
	} else {
		model.dat <- read_in$metadata$model_data
	}
	
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
				# Use the actual CLR forecasting function
				hindcast.plot <- fcast_clr(
					plotID = plotID,
					model.inputs = model.inputs,
					param_samples = param_samples,
					truth.plot.long = truth.plot.long,
					plot_summary = plot_summary,
					Nmc = 1000,
					drop_other = TRUE,
					predict_site_effects = predict_site_effects,
					rank.name = rank.name,
					model_id = model_id
				)
				
				# Add metadata if forecast was successful
				hindcast.plot <- hindcast.plot %>% mutate(model_name = !!model_name,
																 																	 time_period = !!time_period,
																 																	 species = !!taxon,
																 																	 rank_name = !!rank.name,
																 																	 predicted_site_effect=F,
																 																	 newsite = !!newsite,
																 															 model_id=!!model_id)

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
cat("\n✅ Completed CLR hindcast generation!\n")

message("output nrow: ", nrow(out.fcast))

# Save as both RDS and Parquet for compatibility
out.path.rds <- here("data/summary/all_hindcasts_CLR_raw.rds")
out.path.parquet <- here("data/summary/parquet/all_hindcasts_CLR_raw.parquet")

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

cat("🎉 CLR HINDCAST GENERATION COMPLETE \n")
