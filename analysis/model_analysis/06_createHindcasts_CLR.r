#### Create forecasts for individual taxonomic groups using CLR models - CLR-SPECIFIC VERSION

testing=F
#### Reading in files ####

source("source_local.R")
# For CLR models, we'll use the CLR-specific forecasting function
source("microbialForecast/R/run_hindcast_CLR.r")

# Load data.table for optimization
if (!require(data.table, quietly = TRUE)) {
  install.packages("data.table")
}
library(data.table)

# Read in microbial abundances
all.ranks <- c(readRDS(here("data/clean/groupAbundances_16S_2023.rds")),
               readRDS(here("data/clean/groupAbundances_ITS_2023.rds")))

# Read in CLR model outputs to grab parameter estimates
# Try CLR-specific summaries first, then fall back to regular summaries
clr_summaries_file <- here("data/summary/clr_regression_summaries.rds")
if (file.exists(clr_summaries_file)) {
  calibration_model_summaries <- readRDS(clr_summaries_file)
  cat("Using CLR-specific model summaries\n")
} else {
  # Fall back to regular summaries if CLR-specific doesn't exist
  calibration_model_summaries <- readRDS(here("data/summary/logit_beta_regression_summaries.rds"))
  cat("Using regular model summaries (CLR-specific not found)\n")
}

# Read in predicted site effects for CLR models
# Try CLR-specific site effects first, then fall back to regular
clr_site_effects_file <- here("data/summary/site_effects_unobserved_CLR.rds")
if (file.exists(clr_site_effects_file)) {
  read_in = readRDS(clr_site_effects_file)
  pred_effects <- read_in[[2]]
  cat("Using CLR-specific site effects\n")
} else {
  # Fall back to regular site effects if CLR-specific doesn't exist
  read_in = readRDS(here("data/summary/site_effects_unobserved.rds"))
  pred_effects <- read_in[[2]]
  cat("Using regular site effects (CLR-specific not found)\n")
}

pred_effects$fit = pred_effects$Median
pred_effects$se.fit = pred_effects$se_fit

# Read in predictor data, just to get the list of sites missing pC data
all_predictors = readRDS(here("data/clean/all_predictor_data.rds"))

# Read in list of taxa that passed convergence
# Try CLR-specific convergence list first
clr_convergence_file <- here("data/summary/converged_taxa_list_CLR.rds")
if (file.exists(clr_convergence_file)) {
  keep_list <- readRDS(clr_convergence_file)
  cat("Using CLR-specific convergence list\n")
} else {
  # Fall back to regular convergence list if CLR-specific doesn't exist
  keep_list <- readRDS(here("data/summary/converged_taxa_list.rds"))
  cat("Using regular convergence list (CLR-specific not found)\n")
}

# Look for actual CLR models that exist in the output directory
clr_output_dir <- here("data/model_outputs/CLR_regression_expanded/")
if (dir.exists(clr_output_dir)) {
  # Find all CLR sample files (excluding individual chain files)
  clr_sample_files <- list.files(clr_output_dir, pattern = "samples_CLR_expanded_.*\\.rds$", 
                                 recursive = TRUE, full.names = FALSE)
  
  # Filter out individual chain files, keep only combined files
  clr_sample_files <- clr_sample_files[!grepl("_chain[0-9]+", clr_sample_files)]
  
  # Extract model IDs from file names
  clr_model_ids <- unique(gsub("samples_CLR_expanded_(.*)\\.rds", "\\1", clr_sample_files))
  clr_model_ids <- gsub(".*/", "", clr_model_ids)  # Remove directory path
  
  if (length(clr_model_ids) > 0) {
    model_id_list <- clr_model_ids
    cat("Found CLR models:", length(model_id_list), "models\n")
    cat("CLR models:", paste(clr_model_ids, collapse = ", "), "\n")
  } else {
    cat("No CLR sample files found in", clr_output_dir, "\n")
    cat("Processing all models from pred_effects:", length(unique(pred_effects$model_id)), "models\n")
    model_id_list = unique(pred_effects$model_id)
  }
} else {
  cat("CLR output directory not found:", clr_output_dir, "\n")
  cat("Processing all models from pred_effects:", length(unique(pred_effects$model_id)), "models\n")
  model_id_list = unique(pred_effects$model_id)
}

min.date = "20151101"

# Pre-allocate output list for better memory management
tax_output_list = vector("list", length(model_id_list))

# Add progress bar
cat("🚀 Starting CLR hindcast generation for", length(model_id_list), "models...\n")
pb <- txtProgressBar(min = 0, max = length(model_id_list), style = 3)

# for testing - run sequentially instead of in parallel
for(k in 1:length(model_id_list)){
  # For CLR models, we'll need to create a CLR-specific forecasting function
  # source("../../microbialForecast/R/run_hindcast.r")

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
  
  # For CLR models, we don't need to create "other" column since we're working with single-species data
  # rank.df_spec$other <- 1-rank.df_spec[[taxon]]
  
  # Prep validation data using full time series
  # For CLR models, we need to use prepCLRData instead of prepBetaRegData
  if (require(microbialForecast, quietly = TRUE)) {
    if (exists("prepCLRData", envir = asNamespace("microbialForecast"))) {
      full.ts.model.inputs <- prepCLRData(rank.df = rank.df,
                                          min.prev = 0,
                                          min.date = min.date,
                                          max.date = "20200101",
                                          s = taxon)
      cat("Using CLR data preparation\n")
    } else {
      # Fall back to beta regression data preparation if CLR function not available
      full.ts.model.inputs <- prepBetaRegData(rank.df = rank.df,
                                             min.prev = 0,
                                             min.date = min.date,
                                             max.date = "20200101",
                                             full_timeseries = T,
                                             keep_vec = keep_vec)
      cat("Using beta regression data preparation (CLR function not available)\n")
    }
  } else {
    # Fall back to beta regression data preparation if package not available
    full.ts.model.inputs <- prepBetaRegData(rank.df = rank.df,
                                           min.prev = 0,
                                           min.date = min.date,
                                           max.date = "20200101",
                                           full_timeseries = T,
                                           keep_vec = keep_vec)
    cat("Using beta regression data preparation (microbialForecast package not available)\n")
  }

  message("Forecasting with CLR model: ", model_name)

      # Get model outputs - look for CLR-specific output directory first
    clr_output_dir <- here("data/model_outputs/CLR_regression_expanded/")
    regular_output_dir <- here("data/model_outputs/logit_beta_regression/")
  
  if (dir.exists(clr_output_dir)) {
    # For CLR models, the file structure is: CLR_regression_expanded/model_type/samples_CLR_expanded_model_type_taxon_timeperiod.rds
    # Extract model type from model_id
    if (grepl("^cycl_only", model_id)) {
      model_type <- "cycl_only"
    } else if (grepl("^env_cov", model_id)) {
      model_type <- "env_cov"
    } else if (grepl("^env_cycl", model_id)) {
      model_type <- "env_cycl"
    } else {
      model_type <- "cycl_only"  # Default
    }
    
    f <- file.path(clr_output_dir, model_type, paste0("samples_CLR_expanded_", model_id, ".rds"))
  } else {
    f <- file.path(regular_output_dir, model_name, paste0("samples_", model_id, ".rds"))
  }
  
  if(!file.exists(f)) {
    cat("Model output file not found:", f, "\n")
    next()
  }
  
  read_in <- readRDS(f)
  
  # Handle different output structures for CLR vs regular models
  if ("samples" %in% names(read_in)) {
    # CLR model output structure
    if (is.list(read_in$samples) && length(read_in$samples) == 1) {
      # Single-chain model: extract samples from the first (and only) mcmc object
      param_samples <- as.data.frame(as.matrix(read_in$samples[[1]]))
    } else {
      # Multi-chain model: combine all chains
      param_samples <- as.data.frame(as.matrix(read_in$samples))
    }
    model.dat <- read_in$metadata$model_data
  } else {
    # Regular model output structure
    param_samples <- as.data.frame(as.matrix(read_in$samples))
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

  # For testing, limit the number of sites
  testing <- FALSE  # Set to TRUE to limit sites for testing
  if (testing == TRUE) {
    full_site_list <- c(head(site_list, 5), head(new_site_list, 5))
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

    for (plotID in plot_list) {
      # For CLR models, we now have a CLR-specific forecasting function
      # Call the CLR forecasting function
      tryCatch({
        cat("Calling CLR forecasting function for model:", model_id, "plot:", plotID, "\n")
        
        # Call the CLR forecasting function
        forecast_result <- fcast_clr(
          plotID = plotID,
          model.inputs = full.ts.model.inputs,
          param_samples = param_samples,
          truth.plot.long = truth.plot.long,
          plot_summary = plot_summary,
          Nmc = 1000,  # Number of MCMC samples to use
          drop_other = T,
          predict_site_effects = pred_effects,
          rank.name = rank.name,
          model_id = model_id
        )
        
        if (!is.null(forecast_result)) {
          # Add model metadata
          forecast_result$model_id <- model_id
          forecast_result$rank.name <- rank.name
          forecast_result$taxon <- taxon
          forecast_result$model_name <- model_name
          forecast_result$fcast_type <- fcast_type
          forecast_result$newsite <- newsite
          
          # Store the forecast result
          plot_output_list[[plot_counter]] <- forecast_result
          plot_counter = plot_counter + 1
          cat("✅ CLR forecast completed successfully for model:", model_id, "plot:", plotID, "\n")
        } else {
          cat("⚠️  CLR forecast returned NULL for model:", model_id, "plot:", plotID, "\n")
        }
        
      }, error = function(e) {
        cat("❌ Error in CLR forecasting for model", model_id, "plot", plotID, ":", e$message, "\n")
        cat("Error call:", paste(deparse(e$call), collapse="\n"), "\n")
      })
    }
    
    # Combine plot results for this site
    if (plot_counter > 1) {
      site_output_list[[site_counter]] <- rbindlist(plot_output_list[1:(plot_counter-1)], fill=TRUE)
      site_counter = site_counter + 1
    }
  }
  
  # Combine site results for this model
  if (site_counter > 1) {
    tax_output_list[[k]] <- rbindlist(site_output_list[1:(site_counter-1)], fill=TRUE)
  }
}

# Final output binding with data.table
out.fcast <- rbindlist(tax_output_list, fill = TRUE)

# Close progress bar
close(pb)
cat("\n✅ Completed CLR hindcast generation!\n")

message("output nrow: ", nrow(out.fcast))

# Save as both RDS and Parquet for compatibility
out.path.rds <- here("data/summary/CLR_hindcasts.rds")
out.path.parquet <- here("data/summary/parquet/CLR_hindcasts.parquet")

# Ensure parquet directory exists
dir.create(here("data/summary/parquet"), showWarnings = FALSE, recursive = TRUE)

# Save as RDS
saveRDS(out.fcast, out.path.rds)
cat("✓ Saved RDS file:", out.path.rds, "\n")

# Save as Parquet if arrow package is available
if (require(arrow, quietly = TRUE)) {
  arrow::write_parquet(out.fcast, out.path.parquet)
  cat("✓ Saved Parquet file:", out.path.parquet, "\n")
} else {
  cat("⚠️  Arrow package not available - skipping Parquet export\n")
}

cat("CLR hindcast generation completed!\n")

# Filter out NULL results and combine all forecasts
valid_forecasts <- tax_output_list[!sapply(tax_output_list, is.null)]
cat("Valid forecasts generated:", length(valid_forecasts), "out of", length(tax_output_list), "models\n")

if (length(valid_forecasts) > 0) {
  # Combine all forecasts into a single dataframe
  all_forecasts <- do.call(rbind, valid_forecasts)
  
  # Save the combined forecasts
  saveRDS(all_forecasts, here("data/model_outputs/CLR_hindcasts.rds"))
  cat("✅ CLR hindcasts saved to: CLR_hindcasts.rds\n")
  cat("Total forecast rows:", nrow(all_forecasts), "\n")
  
  # Also save individual model results
  saveRDS(valid_forecasts, here("data/model_outputs/CLR_hindcasts_by_model.rds"))
  cat("✅ Individual model forecasts saved to: CLR_hindcasts_by_model.rds\n")
  
  # Summary statistics
  cat("\nForecast Summary:\n")
  cat("Models with successful forecasts:", length(unique(all_forecasts$model_id)), "\n")
  cat("Total forecast timepoints:", nrow(all_forecasts), "\n")
  cat("Date range:", min(all_forecasts$dates, na.rm=TRUE), "to", max(all_forecasts$dates, na.rm=TRUE), "\n")
  
} else {
  cat("❌ No valid forecasts were generated\n")
  cat("Check the error messages above for issues\n")
}
