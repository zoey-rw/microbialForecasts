# Create hindcasts for Dirichlet microbial forecasts
# Updated to follow the working version's approach but adapted for Dirichlet models

library(dplyr)
library(here)
library(foreach)
library(doParallel)
library(data.table)
source(here("source.R"))

cat("=== CREATING DIRICHLET HINDCASTS FOR ALL SITES ===\n")

# Load required data files (following working version)
cat("Loading required data files...\n")

# Read in microbial abundances
all.ranks <- c(readRDS("./data/clean/groupAbundances_16S_2023.rds"),
               readRDS("./data/clean/groupAbundances_ITS_2023.rds"))

# Read in Dirichlet model outputs to grab parameter estimates
calibration_model_summaries <- readRDS(here("data/summary/dirichlet_regression_summaries.rds"))

# Read in predicted site effects for Dirichlet models
dirichlet_site_effects_file <- here("data/summary/site_effects_unobserved_dirichlet.rds")
if (file.exists(dirichlet_site_effects_file)) {
  cat("Using Dirichlet-specific site effects file:", dirichlet_site_effects_file, "\n")
  read_in <- readRDS(dirichlet_site_effects_file)
  pred_effects <- read_in[[2]]
  pred_effects$fit = pred_effects$Median
  pred_effects$se.fit = pred_effects$se_fit
} else {
  cat("Dirichlet-specific site effects not found, using fallback\n")
  # Fall back to regular site effects if Dirichlet-specific doesn't exist
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
all_predictors = readRDS("./data/clean/all_predictor_data.rds")

cat("✅ Data files loaded successfully\n")

# Get available Dirichlet models
model_outputs_dir <- here("data", "model_outputs", "dirichlet_regression")
available_models <- list.dirs(model_outputs_dir, full.names = FALSE, recursive = FALSE)
cat("Available Dirichlet models:", paste(available_models, collapse = ", "), "\n\n")

# Focus on env_cycl model for now
target_model <- "env_cycl"
if (target_model %in% available_models) {
  
  # Get model files for this model type
  model_files <- list.files(file.path(model_outputs_dir, target_model), 
                           pattern = "\\.rds$", full.names = TRUE)
  
  cat("Found", length(model_files), "model files for", target_model, "\n")
  
  # Initialize output list
  tax_output_list <- list()
  
  # Progress tracking
  pb <- txtProgressBar(min = 0, max = length(model_files), style = 3)
  
  for (k in 1:length(model_files)) {
    setTxtProgressBar(pb, k)
    
    model_file <- model_files[k]
    model_id <- basename(model_file)
    model_id <- gsub("\\.rds$", "", model_id)
    
    cat("\nProcessing model:", model_id, "\n")
    
    # Parse model ID to get taxon and time period
    model_parts <- strsplit(model_id, "_")[[1]]
    if (length(model_parts) >= 4) {
      taxon <- model_parts[3]
      time_period <- paste(model_parts[4], model_parts[5], sep = "_")
    } else {
      cat("Warning: Could not parse model ID:", model_id, "\n")
      next
    }
    
    # Load model data
    tryCatch({
      model_data <- readRDS(model_file)
      
      # Extract components
      param_samples <- model_data$samples
      model.inputs <- model_data$model.inputs
      truth.plot.long <- model_data$model.inputs$truth.plot.long
      
      # Get plot summary from calibration summaries
      plot_summary <- calibration_model_summaries$plot_est[[model_id]]
      
      if (is.null(plot_summary)) {
        cat("Warning: No plot summary found for model:", model_id, "\n")
        next
      }
      
    }, error = function(e) {
      cat("Error loading model file:", model_file, "\n")
      cat("Error:", e$message, "\n")
      return(NULL)
    })
    
    # Get site lists
    plot_site_key <- model.inputs$plot_site_key
    new_plot_site_key <- model.inputs$new_plot_site_key
    
    # Combine all sites
    full_site_list <- unique(c(plot_site_key$siteID, new_plot_site_key$siteID))
    
    # Initialize site output list
    site_output_list = vector("list", length(full_site_list))
    site_counter = 1
    
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
        pred_rank = ifelse(fcast_type=="Taxonomic", rank.name, taxon)

        # Forecast with known or random site effect
        tryCatch({
          # Use the actual Dirichlet forecasting function
          hindcast.plot <- fcast_dirichlet(
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
          if (!is.null(hindcast.plot)) {
            hindcast.plot <- hindcast.plot %>% mutate(model_name = !!model_name,
                                                     time_period = !!time_period,
                                                     species = !!taxon,
                                                     rank_name = !!rank.name,
                                                     predicted_site_effect=F,
                                                     newsite = !!newsite,
                                                     model_id=!!model_id)
          }

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
  cat("\n✅ Completed Dirichlet hindcast generation!\n")

  message("output nrow: ", nrow(out.fcast))

  # Save as both RDS and Parquet for compatibility
  out.path.rds <- here("data/summary/all_hindcasts_dirichlet_raw.rds")
  out.path.parquet <- here("data/summary/parquet/all_hindcasts_dirichlet_raw.parquet")

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

  cat("🎉 DIRICHLET HINDCAST GENERATION COMPLETE \n")
  
} else {
  cat("Target model", target_model, "not found in available models\n")
  cat("Available models:", paste(available_models, collapse = ", "), "\n")
}
