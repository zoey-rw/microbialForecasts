# Calculate forecast horizons

tryCatch({
  mem.maxVSize(Inf)
  cat("Memory limit increased to unlimited\n")
}, error = function(e) {
  cat("Note: Could not increase memory limit:", e$message, "\n")
})

MAKE_PLOTS <- FALSE

# Limit thread usage to prevent system overload
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1") 
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(VECLIB_MAXIMUM_THREADS = "1")
Sys.setenv(NUMEXPR_NUM_THREADS = "1")
Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1")
Sys.setenv(RCPP_PARALLEL_MIN_THREADS = "1")
Sys.setenv(RCPP_PARALLEL_MAX_THREADS = "1")

# Disable parallel processing in ggplot2 and other packages
options(ggplot2.continuous.colour = "viridis")
options(ggplot2.continuous.fill = "viridis")

# Visualize the forecast horizon
library(here)
source(here("source.R"))
# install.packages('egg', dependencies = TRUE)
library(egg)
library(ggpmisc) # for polynomial plotting function
library(data.table)
library(ggpubr)
if (requireNamespace("mgcv", quietly = TRUE)) {
  library(mgcv)  # For GAM models
  cat("mgcv package loaded for GAM support\n")
} else {
  cat("Warning: mgcv package not available. GAM models will not be used.\n")
}

# Re-apply thread limits after loading packages
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1") 
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(VECLIB_MAXIMUM_THREADS = "1")
Sys.setenv(NUMEXPR_NUM_THREADS = "1")
Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1")
Sys.setenv(RCPP_PARALLEL_MIN_THREADS = "1")
Sys.setenv(RCPP_PARALLEL_MAX_THREADS = "1")

# Disable parallel processing in R itself
options(mc.cores = 1)
options(cores = 1)

# Force single-threaded execution
if (requireNamespace("parallel", quietly = TRUE)) {
  options(mc.cores = 1)
}

# Use Parquet format if available for memory efficiency, otherwise RDS
# Try nanoparquet first (lightweight), then arrow, then fallback to RDS
parquet_file <- here("data/summary/parquet/all_hindcasts_plsr2.parquet")
rds_file <- here("data/summary/all_hindcasts_plsr2.rds")

if (file.exists(parquet_file)) {
  # Try nanoparquet first (lightweight, no heavy dependencies)
  parquet_success <- FALSE
  if (requireNamespace("nanoparquet", quietly = TRUE)) {
    cat("Attempting to use Parquet file with nanoparquet (memory-efficient)...\n")
    tryCatch({
    hindcast_in <- nanoparquet::read_parquet(parquet_file)
    cat("Loaded", nrow(hindcast_in), "rows from Parquet (nanoparquet)\n")
      parquet_success <- TRUE
    }, error = function(e) {
      cat("WARNING: Failed to read Parquet file with nanoparquet:", e$message, "\n")
      cat("  Falling back to RDS format...\n")
    })
  }
  
  if (!parquet_success && requireNamespace("arrow", quietly = TRUE)) {
    cat("Attempting to use Parquet file with arrow (memory-efficient)...\n")
    tryCatch({
    hindcast_in <- arrow::read_parquet(parquet_file)
    cat("Loaded", nrow(hindcast_in), "rows from Parquet (arrow)\n")
      parquet_success <- TRUE
    }, error = function(e) {
      cat("WARNING: Failed to read Parquet file with arrow:", e$message, "\n")
      cat("  Falling back to RDS format...\n")
    })
  }
  
  if (!parquet_success) {
    cat("WARNING: Parquet file exists but could not be read or no parquet reader available.\n")
    cat("  Falling back to RDS format (2.7GB)...\n")
    if (file.exists(rds_file)) {
      hindcast_in <- readRDS(rds_file)
      cat("Loaded", nrow(hindcast_in), "rows from RDS\n")
    } else {
      stop("No hindcast file found! Please run script 07_tidyHindcasts_v2.r first.")
    }
  }
} else if (file.exists(rds_file)) {
  cat("Using RDS file (2.7GB) - loading into memory...\n")
  cat("WARNING: Large file may hit memory limits. Consider converting to Parquet format.\n")
  hindcast_in <- readRDS(rds_file)
  cat("Loaded", nrow(hindcast_in), "rows from RDS\n")
} else {
  stop("No hindcast file found! Please run script 07_tidyHindcasts_v2.r first.")
}

# Convert to data.table immediately for memory efficiency
# If from Parquet, it may already be a data.frame - convert carefully
if (!inherits(hindcast_in, "data.table")) {
  if (inherits(hindcast_in, "data.frame")) {
    setDT(hindcast_in)
  } else {
    hindcast_in <- as.data.table(hindcast_in)
  }
}

# Fix list columns that should be numeric
if (is.list(hindcast_in$mean)) {
  cat("Converting mean column from list to numeric...\n")
  hindcast_in[, mean := sapply(mean, function(x) {
    if (length(x) == 0) NA_real_ else as.numeric(x)
  })]
}
if (is.list(hindcast_in$sd)) {
  cat("Converting sd column from list to numeric...\n")
  hindcast_in[, sd := sapply(sd, function(x) {
    if (length(x) == 0) NA_real_ else as.numeric(x)
  })]
}
if (is.list(hindcast_in$med)) {
  cat("Converting med column from list to numeric...\n")
  hindcast_in[, med := sapply(med, function(x) {
    if (length(x) == 0) NA_real_ else as.numeric(x)
  })]
}

# Convert truth column from character to numeric
cat("Converting truth column from character to numeric...\n")
hindcast_in[, truth := as.numeric(truth)]

# Assign pretty_group if missing or any NA (fill in missing values)
if (!"pretty_group" %in% names(hindcast_in) || any(is.na(hindcast_in$pretty_group))) {
  cat("Assigning pretty_group based on species taxonomy...\n")
  fg_names <- microbialForecast:::keep_fg_names
  rank_spec_names <- microbialForecast:::rank_spec_names
  all_bacteria <- unlist(rank_spec_names[grepl("_bac$", names(rank_spec_names))])
  all_fungi <- unlist(rank_spec_names[grepl("_fun$", names(rank_spec_names))])
  
  # Initialize pretty_group column if missing
  if (!"pretty_group" %in% names(hindcast_in)) {
    hindcast_in[, pretty_group := NA_character_]
  }
  
  # Only assign for rows where pretty_group is NA
  na_rows <- is.na(hindcast_in$pretty_group)
  if (any(na_rows)) {
    cat("Filling", sum(na_rows), "rows with missing pretty_group\n")
    
    # Assign pretty_group based on species - use data.table syntax for efficiency
    # First, assign for functional groups
    fg_species <- hindcast_in$species[na_rows] %in% fg_names
    if (any(fg_species)) {
      fg_indices <- which(na_rows)[fg_species]
      fg_kingdoms <- microbialForecast:::assign_fg_kingdoms(microbialForecast:::assign_fg_categories(hindcast_in$species[fg_indices]))
      hindcast_in[fg_indices, pretty_group := ifelse(fg_kingdoms == "16S", "Bacteria", "Fungi")]
      na_rows[fg_indices] <- FALSE  # Update na_rows to exclude these
    }
    
    # Then assign for taxonomic bacteria
    bac_species <- hindcast_in$species[na_rows] %in% all_bacteria
    if (any(bac_species)) {
      hindcast_in[which(na_rows)[bac_species], pretty_group := "Bacteria"]
      na_rows[which(na_rows)[bac_species]] <- FALSE
    }
    
    # Then assign for taxonomic fungi
    fun_species <- hindcast_in$species[na_rows] %in% all_fungi
    if (any(fun_species)) {
      hindcast_in[which(na_rows)[fun_species], pretty_group := "Fungi"]
      na_rows[which(na_rows)[fun_species]] <- FALSE
    }
    
    # Fallback: use rank_name pattern if species-based assignment failed
    still_na <- na_rows
    if (any(still_na) && "rank_name" %in% names(hindcast_in)) {
      hindcast_in[which(still_na), pretty_group := ifelse(grepl("_bac$|16S", rank_name), "Bacteria",
                                                           ifelse(grepl("_fun$|ITS", rank_name), "Fungi", NA_character_))]
    }
  }
  
  cat("Final pretty_group distribution:", table(hindcast_in$pretty_group, useNA = "ifany"), "\n")
}

# Add rank_only column if missing (first part of rank_name before underscore)
if (!"rank_only" %in% names(hindcast_in) && "rank_name" %in% names(hindcast_in)) {
  hindcast_in[, rank_only := sapply(strsplit(rank_name, "_", fixed = TRUE), function(x) {
    if (length(x) > 0 && x[1] != "") {
      # Handle special case: "functional" groups don't have _bac/_fun suffix
      if (x[1] %in% c("functional", "diversity")) {
        return(x[1])
      }
      # For taxonomic ranks, return first part (e.g., "phylum_bac" -> "phylum")
      return(x[1])
    } else {
      return(NA_character_)
    }
  })]
  # Set as ordered factor with standard levels
  rank_levels <- c("genus", "family", "order", "class", "phylum", "functional", "diversity")
  hindcast_in[, rank_only := factor(rank_only, levels = rank_levels, ordered = TRUE)]
  cat("Added rank_only column from rank_name\n")
}

converged <- readRDS(here("data/summary/weak_converged_taxa_list.rds"))
# converged_strict <- readRDS(here("data/summary/converged_taxa_list.rds"))
# 
# converged = converged_strict

# Normalize converged model IDs to match hindcast format
converged <- gsub("_beta_regression$", "", converged)

# FIX 1: Normalize model_id early (and only once) - use data.table syntax
hindcast_in[, model_id := gsub("_beta_regression$", "", model_id)]

# Use dateID as timepoint when timepoint is NA
hindcast_in[, corrected_timepoint := ifelse(is.na(timepoint) & !is.na(dateID), dateID, timepoint)]

# Rebuild last_obs from normalized hindcast_data using dateID (not timepoint)
# Use dateID when available (YYYYMM format), as it's needed for proper months calculation
# CRITICAL: Only use calibration data - models without calibration data are excluded
hindcast_cal <- hindcast_in[fcast_period != "hindcast" & !is.na(truth)]
if (nrow(hindcast_cal) > 0) {
  # Use calibration data if available
  last_obs <- hindcast_cal[, .(last_obs = ifelse(any(!is.na(dateID) & dateID >= 10000), 
                                                 max(dateID[!is.na(dateID) & dateID >= 10000], na.rm = TRUE),
                                                 max(corrected_timepoint, na.rm = TRUE))), 
                           by = .(model_id, plotID)]
  rm(hindcast_cal); gc(verbose = FALSE)
  
  # Identify models without calibration data
  all_models <- unique(hindcast_in$model_id)
  models_with_cal <- unique(last_obs$model_id)
  models_without_cal <- setdiff(all_models, models_with_cal)

  if (length(models_without_cal) > 0) {
    cat("Warning:", length(models_without_cal), "models have no calibration data and will be excluded from horizon calculation.\n")
    cat("  Sample excluded models:", paste(head(models_without_cal, 3), collapse=", "), if(length(models_without_cal) > 3) " ..." else "", "\n")
    cat("  Run analysis/model_analysis/diagnose_missing_calibration.r to investigate.\n")
  }
} else {
  # If no calibration data at all, create empty last_obs
  cat("Warning: No calibration data found for any model. No forecast horizons will be calculated.\n")
  last_obs <- data.table(model_id = character(), plotID = character(), last_obs = numeric())
  models_without_cal <- unique(hindcast_in$model_id)
}

hindcast_in[truth == 0, truth := 0.0001]

hindcast_data <- hindcast_in
rm(hindcast_in); gc(verbose = FALSE)


# Debug: Check data before filtering
cat("Debug: hindcast_data rows:", nrow(hindcast_data), "\n")
cat("Debug: fcast_period values:", unique(hindcast_data$fcast_period), "\n")
cat("Debug: new_site values:", unique(hindcast_data$new_site), "\n")
cat("Debug: truth NA count:", sum(is.na(hindcast_data$truth)), "\n")

# Calculate "time since final observation" for forecasts
cat("Starting data processing pipeline...\n")
setDT(hindcast_data)
setDT(last_obs)
setkey(hindcast_data, plotID, model_id)
setkey(last_obs, plotID, model_id)

# Filter to hindcast period with truth values
cat("Debug: Before filter - hindcast_data rows:", nrow(hindcast_data), "\n")
cat("Debug: Rows with !is.na(truth):", sum(!is.na(hindcast_data$truth)), "\n")
cat("Debug: Rows with fcast_period == 'hindcast':", sum(hindcast_data$fcast_period == "hindcast", na.rm=TRUE), "\n")
fcast_horizon_df <- hindcast_data[!is.na(truth) & fcast_period == "hindcast"]
cat("Debug: After filter - fcast_horizon_df rows:", nrow(fcast_horizon_df), "\n")

# CRITICAL: Only include models that have last_obs (i.e., have calibration data)
# Models without calibration data are excluded via nomatch = 0
fcast_horizon_df <- fcast_horizon_df[last_obs, on = c("plotID", "model_id"), nomatch = 0]
cat("Debug: After join with last_obs - fcast_horizon_df rows:", nrow(fcast_horizon_df), "\n")
cat("Debug: Models included after join:", length(unique(fcast_horizon_df$model_id)), "\n")
if ("pretty_group" %in% names(fcast_horizon_df)) {
  pg_table <- fcast_horizon_df[, .(n_models = length(unique(model_id))), by = pretty_group]
  cat("Debug: pretty_group in fcast_horizon_df after join (unique model_id per group):\n")
  print(pg_table)
  if (nrow(fcast_horizon_df[!is.na(pretty_group) & pretty_group == "Fungi"]) == 0) {
    cat("WARNING: No Fungi rows after join - Fungi model_ids may lack (plotID, model_id) in last_obs.\n")
    last_obs_models <- unique(last_obs$model_id)
    hindcast_model_pg <- unique(hindcast_data[!is.na(pretty_group), .(model_id, pretty_group)], by = "model_id")
    in_last_obs <- hindcast_model_pg[model_id %in% last_obs_models]
    cat("  model_ids in last_obs by pretty_group:\n")
    print(in_last_obs[, .N, by = pretty_group])
  }
}

# Convert dateIDs to months difference properly
# dateIDs are in YYYYMM format (e.g., 201306 = June 2013, 201801 = January 2018)
if ("dateID" %in% names(fcast_horizon_df)) {
  fcast_horizon_df[, timepoint_val := ifelse(!is.na(dateID) & dateID >= 10000, dateID, 
                                              ifelse(corrected_timepoint >= 10000, corrected_timepoint, NA_real_))]
} else {
  fcast_horizon_df[, timepoint_val := ifelse(corrected_timepoint >= 10000, corrected_timepoint, NA_real_)]
}

# Calculate months difference (vectorized)
valid_rows <- !is.na(fcast_horizon_df$timepoint_val) & fcast_horizon_df$timepoint_val >= 10000 &
              !is.na(fcast_horizon_df$last_obs) & fcast_horizon_df$last_obs >= 10000

fcast_horizon_df[, months_since_obs := NA_real_]
if (sum(valid_rows) > 0) {
  timepoint_year <- floor(fcast_horizon_df$timepoint_val[valid_rows] / 100)
  timepoint_month <- fcast_horizon_df$timepoint_val[valid_rows] %% 100
  last_obs_year <- floor(fcast_horizon_df$last_obs[valid_rows] / 100)
  last_obs_month <- fcast_horizon_df$last_obs[valid_rows] %% 100
  fcast_horizon_df[valid_rows, months_since_obs := (timepoint_year - last_obs_year) * 12 + (timepoint_month - last_obs_month)]
}

fcast_horizon_df[, timepoint_val := NULL]

cat("Completed initial data processing...\n")
cat("months_since_obs range:", range(fcast_horizon_df$months_since_obs, na.rm=TRUE), "\n")

# Calculate CRPS scores for each row before summarizing
cat("Starting CRPS calculation...\n")
if (!inherits(fcast_horizon_df, "data.table")) {
  setDT(fcast_horizon_df)
}
# Initialize both CRPS columns (truncated is primary for horizon calc, norm kept for reference)
fcast_horizon_df[, crps := NA_real_]
fcast_horizon_df[, crps_norm := NA_real_]
valid_rows <- !is.na(fcast_horizon_df$truth) & !is.na(fcast_horizon_df$mean) & !is.na(fcast_horizon_df$sd)
cat("Valid rows for CRPS:", sum(valid_rows), "out of", nrow(fcast_horizon_df), "\n")
if (sum(valid_rows) > 0) {
  cat("Computing CRPS for", sum(valid_rows), "rows...\n")
  chunk_size <- 10000
  n_chunks <- ceiling(sum(valid_rows) / chunk_size)
  valid_indices <- which(valid_rows)
  
  for (i in 1:n_chunks) {
    start_idx <- (i-1) * chunk_size + 1
    end_idx <- min(i * chunk_size, sum(valid_rows))
    chunk_indices <- valid_indices[start_idx:end_idx]
    
    cat("Processing CRPS chunk", i, "of", n_chunks, "(", length(chunk_indices), "rows)\n")
    # Filter out rows with non-positive SD before CRPS calculation
    valid_chunk <- chunk_indices[fcast_horizon_df[chunk_indices, sd > 0 & is.finite(sd)]]
    if (length(valid_chunk) > 0) {
      fcast_horizon_df[valid_chunk, crps := crps(truth, family = "tnorm",
                                                  location = mean,
                                                  scale = sd,
                                                  lower = 0, upper = 1)]
      fcast_horizon_df[valid_chunk, crps_norm := crps_norm(truth, mean, sd)]
    }
  }
  cat("CRPS calculation completed\n")
  crps_populated <- sum(!is.na(fcast_horizon_df$crps))
  cat("CRPS populated for", crps_populated, "rows (", round(100 * crps_populated / nrow(fcast_horizon_df), 1), "%)\n")
}

# Calculate average scores by model
cat("Calculating average scores by model...\n")
setDT(fcast_horizon_df)

# Helper function to safely calculate RSQ
safe_rsq <- function(truth, mean) {
  valid <- !is.na(truth) & !is.na(mean)
  if (sum(valid) > 1) {
    tryCatch({
      summary(lm(truth[valid] ~ mean[valid]))$r.squared
    }, error = function(e) NA_real_)
  } else {
    NA_real_
  }
}

# Determine grouping columns - include rank_only only if it exists
grouping_cols <- c("species", "model_name", "pretty_group", "model_id", "siteID", "months_since_obs")
if ("rank_only" %in% names(fcast_horizon_df)) {
  grouping_cols <- c("species", "model_name", "pretty_group", "rank_only", "model_id", "siteID", "months_since_obs")
}

fcast_horizon_model_mean <- fcast_horizon_df[, {
  valid <- !is.na(truth) & !is.na(mean)
  n_valid <- sum(valid)
  if (n_valid > 0) {
    rmse_val <- rmse(actual = truth[valid], predicted = mean[valid])
    rsq1_val <- if (!is.na(rmse_val)) {
      var_val <- var(truth[valid], na.rm = TRUE)
      if (!is.na(var_val) && var_val > 0) 1 - (rmse_val^2)/var_val else NA_real_
    } else NA_real_
    rsq_val <- if (n_valid > 1) safe_rsq(truth, mean) else NA_real_
    mape_val <- mape(actual = truth[valid], predicted = mean[valid])
    abund <- mean(truth, na.rm = TRUE)
    # FIX 3: Apply 0.005 abundance floor for RMSE.norm (consistent with script 08)
    abund_denom <- max(abund, 0.005)
    rmse_norm_val <- if (!is.na(rmse_val)) rmse_val/abund_denom else NA_real_

    list(
      # FIX 2: Use truncated CRPS as primary mean_crps (consistent with script 08)
      mean_crps = mean(crps, na.rm = TRUE),
      mean_crps_norm = mean(crps_norm, na.rm = TRUE),
      RMSE = rmse_val,
      # FIX 1: RSQ uses Nash-Sutcliffe (for horizon crossing); regression R² kept as RSQ.reg
      RSQ.1 = rsq1_val,
      RSQ = rsq1_val,
      RSQ.reg = rsq_val,
      MAPE = mape_val,
      abundance = abund,
      RMSE.norm = rmse_norm_val
    )
  } else {
    list(
      mean_crps = NA_real_,
      mean_crps_norm = NA_real_,
      RMSE = NA_real_,
      RSQ.1 = NA_real_,
      RSQ = NA_real_,
      RSQ.reg = NA_real_,
      MAPE = NA_real_,
      abundance = NA_real_,
      RMSE.norm = NA_real_
    )
  }
}, by = grouping_cols]

# NOTE: Do NOT clamp RSQ.1 to 0 here — negative NSE values are informative
# for the loess fit in pooled_loess_horizon. Clamping creates a false floor
# at the null line (NSE=0), causing the loess to cross the null sooner.
# A separate RSQ_plot column is used where clamping is needed for display.
fcast_horizon_model_mean[, RMSE.norm := ifelse(RMSE.norm > 5, 5, RMSE.norm)]
cat("Completed average scores calculation\n")

# Calculate number of observations per months out
fcast_n_obs <- fcast_horizon_df[, .(n_obs = .N), by = .(model_id, months_since_obs)]
setkey(fcast_horizon_model_mean, model_id, months_since_obs)
setkey(fcast_n_obs, model_id, months_since_obs)
fcast_horizon_model_mean <- fcast_n_obs[fcast_horizon_model_mean, on = c("model_id", "months_since_obs")] 


# Calculate null baselines from calibration data
# Three null types: site_mean, persistence, climatological
cat("Calculating null baselines (site mean, persistence, climatological)...\n")
hindcast_cal_data <- hindcast_data[fcast_period != "hindcast" & !is.na(truth)]
if (nrow(hindcast_cal_data) > 0) {
  # NULL 1: Site mean (existing) - mean truth per site across all calibration
  site_mean <- hindcast_cal_data[, .(site_mean = mean(truth, na.rm = TRUE),
                                      site_sd = sd(truth, na.rm = TRUE)),
                                  by = .(species, model_id, model_name, siteID)]
  overall_mean <- hindcast_cal_data[, .(overall_mean = mean(truth, na.rm = TRUE),
                                        overall_sd = sd(truth, na.rm = TRUE)),
                                    by = .(species, model_id)]

  # NULL 2: Persistence - last observed calibration truth per plot
  # For each plot, take the truth value at the last calibration timepoint
  if ("dateID" %in% names(hindcast_cal_data)) {
    hindcast_cal_data[, dateID_int := as.integer(dateID)]
    persistence_null <- hindcast_cal_data[!is.na(truth), .SD[which.max(dateID_int)],
                                          by = .(model_id, plotID)][, .(model_id, plotID, persistence_val = truth)]
    cat("  Persistence null: computed for", nrow(persistence_null), "model-plot combinations\n")
  } else {
    persistence_null <- data.table(model_id = character(), plotID = character(), persistence_val = numeric())
    cat("  Warning: No dateID column - persistence null unavailable\n")
  }

  # NULL 3: Climatological mean - mean truth per site per month-of-year
  if ("dateID" %in% names(hindcast_cal_data)) {
    hindcast_cal_data[, month_of_year := as.integer(dateID) %% 100]
    climatology_null <- hindcast_cal_data[!is.na(truth), .(clim_mean = mean(truth, na.rm = TRUE),
                                                            clim_sd = sd(truth, na.rm = TRUE),
                                                            clim_n = .N),
                                           by = .(model_id, siteID, month_of_year)]
    # Fallback: use site_sd when per-month sample is too small
    climatology_null[clim_n < 3, clim_sd := NA_real_]
    cat("  Climatological null: computed for", nrow(climatology_null), "model-site-month combinations\n")
  } else {
    climatology_null <- data.table(model_id = character(), siteID = character(),
                                   month_of_year = integer(), clim_mean = numeric(),
                                   clim_sd = numeric(), clim_n = integer())
    cat("  Warning: No dateID column - climatological null unavailable\n")
  }

  rm(hindcast_cal_data); gc(verbose = FALSE)
} else {
  cat("Warning: No calibration data found. Using hindcast data for site means.\n")
  hindcast_for_mean <- hindcast_data[fcast_period == "hindcast" & !is.na(truth)]
  if (nrow(hindcast_for_mean) > 0) {
    site_mean <- hindcast_for_mean[, .(site_mean = mean(truth, na.rm = TRUE),
                                        site_sd = sd(truth, na.rm = TRUE)),
                                    by = .(species, model_id, model_name, siteID)]
    overall_mean <- hindcast_for_mean[, .(overall_mean = mean(truth, na.rm = TRUE),
                                          overall_sd = sd(truth, na.rm = TRUE)),
                                      by = .(species, model_id)]
    rm(hindcast_for_mean); gc(verbose = FALSE)
  } else {
    site_mean <- data.table(species = character(), model_id = character(), model_name = character(),
                           siteID = character(), site_mean = numeric(), site_sd = numeric())
    overall_mean <- data.table(species = character(), model_id = character(),
                              overall_mean = numeric(), overall_sd = numeric())
  }
  persistence_null <- data.table(model_id = character(), plotID = character(), persistence_val = numeric())
  climatology_null <- data.table(model_id = character(), siteID = character(),
                                 month_of_year = integer(), clim_mean = numeric(),
                                 clim_sd = numeric(), clim_n = integer())
  rm(hindcast_cal_data); gc(verbose = FALSE)
}
cat("  Null baselines computed\n")

setkey(site_mean, species, model_id)
setkey(overall_mean, species, model_id)
historical_mean <- merge(site_mean, overall_mean, by = c("species", "model_id"), all.x = TRUE)
rm(site_mean, overall_mean); gc(verbose = FALSE)
cat("Completed historical means calculation (", nrow(historical_mean), "combinations)\n")

# Debug: Check data before merge
cat("Debug: fcast_horizon_df rows:", nrow(fcast_horizon_df), "\n")
cat("Debug: historical_mean rows:", nrow(historical_mean), "\n")
cat("Debug: fcast_horizon_df columns:", colnames(fcast_horizon_df), "\n")
cat("Debug: historical_mean columns:", colnames(historical_mean), "\n")

# Merge historical means into main df, and use to calculate null scores
setDT(historical_mean)
setkey(historical_mean, species, model_id, siteID, model_name)
setkey(fcast_horizon_df, species, model_id, siteID, model_name)
fcast_horizon_df2 <- historical_mean[fcast_horizon_df, on = c("species", "model_id", "siteID", "model_name")]


# Debug: Check site_sd values before filtering
cat("Debug: fcast_horizon_df2 rows:", nrow(fcast_horizon_df2), "\n")
cat("Debug: site_sd unique count:", length(unique(fcast_horizon_df2$site_sd)), "\n")
cat("Debug: site_sd summary: ",
    "min=", suppressWarnings(min(fcast_horizon_df2$site_sd, na.rm=TRUE)),
    " median=", suppressWarnings(median(fcast_horizon_df2$site_sd, na.rm=TRUE)),
    " max=", suppressWarnings(max(fcast_horizon_df2$site_sd, na.rm=TRUE)), "\n")

# Calculate null site scores
# Null = one number per (model_id, siteID): mean CRPS/RSQ/RMSE when "predicting" with
# site_mean (historical mean abundance at that site) across *all* forecast timepoints in the
# hindcast window. So we compare: forecast curve (performance at 0, 3, 6, 9, ... months) vs
# this single pooled null. Short horizons occur when: (1) the null is strong (site persistence
# is a good baseline), so the model only stays better for a short time; (2) evaluation is
# sparse (e.g. only 0, 3, 6, 9, 12 months), so the fitted curve crosses soon after 0; (3) the
# model degrades quickly with lead time. The null is NOT lead-time-specific (we do not compare
# "model at 3 months" to "null at 3 months" only).
cat("Calculating null site scores...\n")
setDT(fcast_horizon_df2)

# Determine grouping columns for null calculation
null_grouping_cols <- c("species", "model_name", "pretty_group", "model_id", "siteID")
if ("rank_only" %in% names(fcast_horizon_df2)) {
  null_grouping_cols <- c("species", "model_name", "pretty_group", "rank_only", "model_id", "siteID")
}

fcast_horizon_null_site <- fcast_horizon_df2[site_sd != 0, {
  valid_crps_trunc <- !is.na(truth) & !is.na(site_mean) & !is.na(site_sd)
  valid_rmse <- !is.na(truth) & !is.na(site_mean)

  # FIX 2: Use truncated CRPS as primary null (consistent with script 08)
  null_crps_trunc_val <- if (sum(valid_crps_trunc) > 0) {
    tryCatch({
      mean(crps(truth[valid_crps_trunc], family = "tnorm",
                location = site_mean[valid_crps_trunc],
                scale = site_sd[valid_crps_trunc],
                lower = 0, upper = 1), na.rm = TRUE)
    }, error = function(e) NA_real_)
  } else NA_real_

  # Keep unbounded CRPS for reference
  null_crps_norm_val <- if (sum(valid_crps_trunc) > 0) {
    mean(crps_norm(truth[valid_crps_trunc], site_mean[valid_crps_trunc], site_sd[valid_crps_trunc]), na.rm = TRUE)
  } else NA_real_

  null_rmse_val <- if (sum(valid_rmse) > 0) {
    rmse(actual = truth[valid_rmse], predicted = site_mean[valid_rmse])
  } else NA_real_

  # FIX 1: RSQ.1 (Nash-Sutcliffe) is the meaningful null; lm(truth ~ constant) always gives R²=0
  null_rsq1_val <- if (sum(valid_rmse) > 0 && !is.na(null_rmse_val)) {
    var_val <- var(truth[valid_rmse], na.rm = TRUE)
    if (!is.na(var_val) && var_val > 0) 1 - (null_rmse_val^2)/var_val else NA_real_
  } else NA_real_

  null_mape_val <- if (sum(valid_rmse) > 0) {
    mape(actual = truth[valid_rmse], predicted = site_mean[valid_rmse])
  } else NA_real_

  abund <- mean(truth, na.rm = TRUE)
  # FIX 3: Apply 0.005 abundance floor for RMSE.norm (consistent with script 08)
  abund_denom <- max(abund, 0.005)
  rmse_norm_val <- if (!is.na(null_rmse_val)) null_rmse_val/abund_denom else NA_real_

  list(
    null_CRPS = null_crps_norm_val,
    null_CRPS_truncated = null_crps_trunc_val,
    # FIX 2: mean_crps now uses truncated CRPS (was unbounded)
    null_mean_crps = null_crps_trunc_val,
    null_RMSE = null_rmse_val,
    null_RSQ.1 = null_rsq1_val,
    null_MAPE = null_mape_val,
    # FIX 1: Set null_RSQ = null_RSQ.1 (Nash-Sutcliffe) instead of broken lm R²
    null_RSQ = null_rsq1_val,
    abundance = abund,
    RMSE.norm = rmse_norm_val
  )
}, by = null_grouping_cols]

fcast_horizon_null_site[, null_RSQ.1 := ifelse(null_RSQ.1 < 0, 0, null_RSQ.1)]
fcast_horizon_null_site[, null_RSQ := ifelse(null_RSQ < 0, 0, null_RSQ)]
fcast_horizon_null_site[, months_since_obs := Inf]
fcast_horizon_null_site[, null_type := "site_mean"]
cat("Completed site mean null scores\n")

# ── Persistence null scores ──────────────────────────────────────────────────
# Null prediction = last calibration truth per plot, applied to all future timepoints
cat("Calculating persistence null scores...\n")
if (nrow(persistence_null) > 0 && nrow(fcast_horizon_df2) > 0) {
  # Join persistence predictions into the hindcast data
  fcast_horizon_df2_persist <- merge(fcast_horizon_df2, persistence_null, by = c("model_id", "plotID"), all.x = TRUE)

  # Compute persistence null per (model_id, siteID) — aggregate across plots and timepoints
  persist_null_site <- fcast_horizon_df2_persist[!is.na(persistence_val) & !is.na(truth), {
    valid <- !is.na(truth) & !is.na(persistence_val)
    if (sum(valid) < 2) return(NULL)

    p_rmse <- rmse(actual = truth[valid], predicted = persistence_val[valid])
    p_rsq1 <- if (!is.na(p_rmse)) {
      var_val <- var(truth[valid], na.rm = TRUE)
      if (!is.na(var_val) && var_val > 0) max(0, 1 - (p_rmse^2) / var_val) else NA_real_
    } else NA_real_

    # Use site_sd for CRPS uncertainty (persistence has no natural uncertainty estimate)
    p_crps <- if (!is.na(site_sd[1]) && site_sd[1] > 0) {
      tryCatch(mean(crps(truth[valid], family = "tnorm",
                         location = persistence_val[valid],
                         scale = site_sd[valid],
                         lower = 0, upper = 1), na.rm = TRUE),
               error = function(e) NA_real_)
    } else NA_real_

    abund <- mean(truth, na.rm = TRUE)
    abund_denom <- max(abund, 0.005)

    list(null_CRPS = NA_real_, null_CRPS_truncated = p_crps, null_mean_crps = p_crps,
         null_RMSE = p_rmse, null_RSQ.1 = p_rsq1, null_MAPE = NA_real_,
         null_RSQ = p_rsq1, abundance = abund, RMSE.norm = p_rmse / abund_denom)
  }, by = null_grouping_cols]

  if (nrow(persist_null_site) > 0) {
    persist_null_site[, null_RSQ.1 := ifelse(null_RSQ.1 < 0, 0, null_RSQ.1)]
    persist_null_site[, null_RSQ := ifelse(null_RSQ < 0, 0, null_RSQ)]
    persist_null_site[, months_since_obs := Inf]
    persist_null_site[, null_type := "persistence"]
    cat("  Persistence null scores:", nrow(persist_null_site), "site-model combinations\n")
  }
  rm(fcast_horizon_df2_persist)
} else {
  persist_null_site <- data.table()
  cat("  Persistence null: skipped (no data)\n")
}

# ── Climatological null scores ───────────────────────────────────────────────
# Null prediction = site-month mean from calibration, matched to hindcast month
cat("Calculating climatological null scores...\n")
if (nrow(climatology_null) > 0 && nrow(fcast_horizon_df2) > 0 && "dateID" %in% names(fcast_horizon_df2)) {
  fcast_horizon_df2[, month_of_year := as.integer(dateID) %% 100]
  fcast_horizon_df2_clim <- merge(fcast_horizon_df2, climatology_null,
                                   by = c("model_id", "siteID", "month_of_year"), all.x = TRUE)
  # Fallback clim_sd to site_sd when per-month estimate unavailable
  fcast_horizon_df2_clim[is.na(clim_sd), clim_sd := site_sd]

  clim_null_site <- fcast_horizon_df2_clim[!is.na(clim_mean) & !is.na(truth), {
    valid <- !is.na(truth) & !is.na(clim_mean)
    if (sum(valid) < 2) return(NULL)

    c_rmse <- rmse(actual = truth[valid], predicted = clim_mean[valid])
    c_rsq1 <- if (!is.na(c_rmse)) {
      var_val <- var(truth[valid], na.rm = TRUE)
      if (!is.na(var_val) && var_val > 0) max(0, 1 - (c_rmse^2) / var_val) else NA_real_
    } else NA_real_

    c_crps <- if (any(!is.na(clim_sd) & clim_sd > 0)) {
      valid_sd <- valid & !is.na(clim_sd) & clim_sd > 0
      if (sum(valid_sd) > 0) {
        tryCatch(mean(crps(truth[valid_sd], family = "tnorm",
                           location = clim_mean[valid_sd],
                           scale = clim_sd[valid_sd],
                           lower = 0, upper = 1), na.rm = TRUE),
                 error = function(e) NA_real_)
      } else NA_real_
    } else NA_real_

    abund <- mean(truth, na.rm = TRUE)
    abund_denom <- max(abund, 0.005)

    list(null_CRPS = NA_real_, null_CRPS_truncated = c_crps, null_mean_crps = c_crps,
         null_RMSE = c_rmse, null_RSQ.1 = c_rsq1, null_MAPE = NA_real_,
         null_RSQ = c_rsq1, abundance = abund, RMSE.norm = c_rmse / abund_denom)
  }, by = null_grouping_cols]

  if (nrow(clim_null_site) > 0) {
    clim_null_site[, null_RSQ.1 := ifelse(null_RSQ.1 < 0, 0, null_RSQ.1)]
    clim_null_site[, null_RSQ := ifelse(null_RSQ < 0, 0, null_RSQ)]
    clim_null_site[, months_since_obs := Inf]
    clim_null_site[, null_type := "climatological"]
    cat("  Climatological null scores:", nrow(clim_null_site), "site-model combinations\n")
  }
  rm(fcast_horizon_df2_clim)
} else {
  clim_null_site <- data.table()
  cat("  Climatological null: skipped (no data)\n")
}

# Combine all null types into a single table
fcast_horizon_null_all <- rbindlist(list(fcast_horizon_null_site, persist_null_site, clim_null_site), fill = TRUE)
cat("Total null scores:", nrow(fcast_horizon_null_all), "across",
    length(unique(fcast_horizon_null_all$null_type)), "null types\n")


# Calculate last observation scores
cat("Calculating last observation scores...\n")
setDT(last_obs)
cat("Debug: last_obs structure - rows:", nrow(last_obs), ", sample:\n")
print(head(last_obs, 3))
cat("Debug: last_obs values range:", range(last_obs$last_obs, na.rm=TRUE), "\n")

# Use dateID for matching if available, otherwise use corrected_timepoint
if ("dateID" %in% names(hindcast_data)) {
  setkey(hindcast_data, plotID, model_id, dateID)
  last_obs_for_join <- copy(last_obs)
  setnames(last_obs_for_join, "last_obs", "dateID")
  setkey(last_obs_for_join, plotID, model_id, dateID)
  last_obs_values <- hindcast_data[last_obs_for_join, on = c("plotID", "model_id", "dateID"), nomatch = 0]
  rm(last_obs_for_join)
} else {
  last_obs[, corrected_timepoint := last_obs]
  setkey(last_obs, plotID, model_id, corrected_timepoint)
  setkey(hindcast_data, plotID, model_id, corrected_timepoint)
  last_obs_values <- hindcast_data[last_obs, on = c("plotID", "model_id", "corrected_timepoint"), nomatch = 0]
}

cat("Debug: last_obs_values rows after join:", nrow(last_obs_values), "\n")
if (nrow(last_obs_values) == 0) {
  cat("WARNING: No matches found in last_obs join! This will cause all metrics to be NA.\n")
  cat("  Checking hindcast_data columns:", paste(head(names(hindcast_data), 10), collapse=", "), "...\n")
  cat("  Sample hindcast_data plotID/model_id:", 
      paste(head(unique(paste(hindcast_data$plotID, hindcast_data$model_id, sep=":")), 3), collapse=", "), "\n")
  cat("  Sample last_obs plotID/model_id:", 
      paste(head(unique(paste(last_obs$plotID, last_obs$model_id, sep=":")), 3), collapse=", "), "\n")
}

rm(hindcast_data); gc(verbose = FALSE)

last_obs_values[, months_since_obs := 0]

# Determine grouping columns for last_obs_score
# CRITICAL: Include siteID if available so month 0 data is site-specific (matches fcast_horizon_model_mean structure)
last_obs_grouping_cols <- c("species", "model_name", "pretty_group", "model_id", "months_since_obs")
if ("rank_only" %in% names(last_obs_values)) {
  last_obs_grouping_cols <- c("species", "model_name", "pretty_group", "rank_only", "model_id", "months_since_obs")
}
# Add siteID if available in last_obs_values (needed for site-specific horizon calculations)
if ("siteID" %in% names(last_obs_values)) {
  last_obs_grouping_cols <- c(last_obs_grouping_cols, "siteID")
} else if ("plotID" %in% names(last_obs_values)) {
  # Extract siteID from plotID if available
  last_obs_values[, siteID := gsub("_.*$", "", plotID)]
  last_obs_grouping_cols <- c(last_obs_grouping_cols, "siteID")
}

last_obs_score <- last_obs_values[, {
  valid_crps <- !is.na(truth) & !is.na(mean) & !is.na(sd)
  valid_rmse <- !is.na(truth) & !is.na(mean)

  # FIX 2: Truncated CRPS as primary
  crps_trunc_val <- if (sum(valid_crps) > 0) {
    tryCatch({
      mean(crps(truth[valid_crps], family = "tnorm",
                location = mean[valid_crps],
                scale = sd[valid_crps],
                lower = 0, upper = 1), na.rm = TRUE)
    }, error = function(e) NA_real_)
  } else NA_real_

  crps_norm_val <- if (sum(valid_crps) > 0) {
    mean(crps_norm(truth[valid_crps], mean[valid_crps], sd[valid_crps]), na.rm = TRUE)
  } else NA_real_

  rmse_val <- if (sum(valid_rmse) > 0) {
    rmse(actual = truth[valid_rmse], predicted = mean[valid_rmse])
  } else NA_real_

  rsq1_val <- if (sum(valid_rmse) > 0 && !is.na(rmse_val)) {
    var_val <- var(truth[valid_rmse], na.rm = TRUE)
    if (!is.na(var_val) && var_val > 0) 1 - (rmse_val^2)/var_val else NA_real_
  } else NA_real_

  mape_val <- if (sum(valid_rmse) > 0) {
    mape(actual = truth[valid_rmse], predicted = mean[valid_rmse])
  } else NA_real_

  rsq_val <- if (sum(valid_rmse) > 1) {
    tryCatch({
      summary(lm(truth[valid_rmse] ~ mean[valid_rmse]))$r.squared
    }, error = function(e) NA_real_)
  } else NA_real_

  abund <- mean(truth, na.rm = TRUE)
  # FIX 3: Apply 0.005 abundance floor for RMSE.norm (consistent with script 08)
  abund_denom <- max(abund, 0.005)
  rmse_norm_val <- if (!is.na(rmse_val)) rmse_val/abund_denom else NA_real_

  list(
    CRPS = crps_norm_val,
    CRPS_truncated = crps_trunc_val,
    # FIX 2: mean_crps uses truncated CRPS (used for horizon crossing)
    mean_crps = crps_trunc_val,
    mean_crps_norm = crps_norm_val,
    RMSE = rmse_val,
    RSQ.1 = rsq1_val,
    MAPE = mape_val,
    # FIX 1: RSQ uses Nash-Sutcliffe for horizon crossing; regression R² kept as RSQ.reg
    RSQ = rsq1_val,
    RSQ.reg = rsq_val,
    abundance = abund,
    RMSE.norm = rmse_norm_val
  )
}, by = last_obs_grouping_cols]

last_obs_score[, RSQ.1 := ifelse(RSQ.1 < 0, 0, RSQ.1)]
last_obs_score[, RSQ := ifelse(RSQ < 0, 0, RSQ)]
last_obs_score[, RMSE.norm := ifelse(RMSE.norm > 5, 5, RMSE.norm)]

# Debug: Check if last_obs_score has valid data
cat("Debug: last_obs_score rows:", nrow(last_obs_score), "\n")
if (nrow(last_obs_score) > 0) {
  cat("Debug: last_obs_score sample (first 5 rows):\n")
  print(head(last_obs_score[, .(model_id, months_since_obs, RSQ, CRPS, CRPS_truncated, mean_crps, RMSE.norm)], 5))
  cat("Debug: last_obs_score - RSQ non-NA:", sum(!is.na(last_obs_score$RSQ)), 
      ", CRPS_truncated non-NA:", sum(!is.na(last_obs_score$CRPS_truncated)),
      ", mean_crps non-NA:", sum(!is.na(last_obs_score$mean_crps)), 
      ", RMSE.norm non-NA:", sum(!is.na(last_obs_score$RMSE.norm)), "\n")
  
  # Check why metrics might be NA - look at source data
  if (nrow(last_obs_values) > 0) {
    cat("Debug: last_obs_values sample (checking if truth/mean/sd exist):\n")
    sample_vals <- head(last_obs_values[, .(model_id, truth, mean, sd, months_since_obs)], 5)
    print(sample_vals)
    cat("Debug: last_obs_values - truth non-NA:", sum(!is.na(last_obs_values$truth)),
        ", mean non-NA:", sum(!is.na(last_obs_values$mean)),
        ", sd non-NA:", sum(!is.na(last_obs_values$sd)), "\n")
  }
  
  if (!"mean_crps" %in% names(last_obs_score)) {
    stop("last_obs_score must have mean_crps; check last_obs_score aggregation block.")
  }
} else {
  cat("WARNING: last_obs_score is empty! This means no last observation data was found.\n")
}

# Keep a site-level version of last_obs_values for per-site horizon calculation
# This contains the raw data at months_since_obs = 0 for each site
if (nrow(last_obs_values) > 0 && "siteID" %in% names(last_obs_values)) {
  last_obs_values_site <- copy(last_obs_values[, .(model_id, siteID, plotID, truth, mean, sd, dateID, corrected_timepoint)])
  # Calculate crps if needed
  if (!"crps" %in% names(last_obs_values_site) && "sd" %in% names(last_obs_values_site)) {
    valid_crps <- !is.na(last_obs_values_site$truth) & 
                  !is.na(last_obs_values_site$mean) & 
                  !is.na(last_obs_values_site$sd)
    if (sum(valid_crps) > 0) {
      last_obs_values_site[valid_crps, crps := crps_norm(truth, mean, sd)]
    }
  }
} else {
  last_obs_values_site <- data.table()
}

rm(last_obs_values); gc(verbose = FALSE)
cat("Completed last observation scores calculation\n")


# Combine for plotting! - use data.table
cat("Combining data for plotting...\n")
setDT(fcast_horizon_null_all)
setDT(last_obs_score)
setDT(fcast_horizon_model_mean)

# Rename columns by removing "null_" prefix (for backward compat, use site_mean null in to_plot)
fcast_horizon_null_site_compat <- fcast_horizon_null_all[null_type == "site_mean"]
null_cols <- grep("^null_", names(fcast_horizon_null_site_compat), value = TRUE)
new_names <- gsub("^null_", "", null_cols)
setnames(fcast_horizon_null_site_compat, null_cols, new_names)

comparison_scores <- rbindlist(list(fcast_horizon_null_site_compat, last_obs_score), fill = TRUE)
rm(fcast_horizon_null_site_compat)

# Merge n_obs from fcast_n_obs into fcast_horizon_model_mean if not already present
if (!"n_obs" %in% names(fcast_horizon_model_mean)) {
  n_obs_by_model <- fcast_n_obs[, .(n_obs = sum(n_obs, na.rm = TRUE)), by = model_id]
  setkey(n_obs_by_model, model_id)
  setkey(fcast_horizon_model_mean, model_id)
  fcast_horizon_model_mean <- n_obs_by_model[fcast_horizon_model_mean, on = "model_id"]
  rm(n_obs_by_model); gc(verbose = FALSE)
}

# Filter for models with sufficient observations (lower threshold to include more models)
# Note: If hindcast data have few discrete months_since_obs (e.g. only 9 and 12), horizon values
# will cluster at those months; continuous variation per taxon requires more evaluation timepoints.
fcast_horizon_model_mean_filtered <- fcast_horizon_model_mean[is.na(n_obs) | n_obs > 5]
if ("pretty_group" %in% names(fcast_horizon_model_mean_filtered)) {
  pg_after <- fcast_horizon_model_mean_filtered[, .(n_models = length(unique(model_id))), by = pretty_group]
  cat("Debug: pretty_group in fcast_horizon_model_mean_filtered (after n_obs > 5):\n")
  print(pg_after)
  n_dropped <- nrow(fcast_horizon_model_mean[!is.na(n_obs) & n_obs <= 5])
  if (n_dropped > 0) {
    dropped_pg <- fcast_horizon_model_mean[!is.na(n_obs) & n_obs <= 5]
    if ("pretty_group" %in% names(dropped_pg)) {
      cat("  Dropped (n_obs <= 5) by pretty_group:\n")
      print(dropped_pg[, .(n_models = length(unique(model_id))), by = pretty_group])
    }
  }
}

# Check column availability before combining
comp_cols <- names(comparison_scores)
model_cols <- names(fcast_horizon_model_mean_filtered)

# mean_crps must exist in both: null block and last_obs_score produce it; fcast_horizon_model_mean calculates it
if (!"mean_crps" %in% comp_cols) {
  stop("comparison_scores must have mean_crps. Ensure fcast_horizon_null_site includes null_mean_crps and last_obs_score includes mean_crps.")
}
if (!"mean_crps" %in% model_cols) {
  stop("fcast_horizon_model_mean_filtered must have mean_crps. Check fcast_horizon_model_mean aggregation.")
}

to_plot <- rbindlist(list(comparison_scores, fcast_horizon_model_mean_filtered), fill = TRUE)
rm(fcast_horizon_model_mean_filtered); gc(verbose = FALSE)

# CRITICAL FIX: Add month 0 data (last_obs_score) to fcast_horizon_model_mean
# This ensures month 0 data is available for the improved horizon calculation methods
# last_obs_score has months_since_obs = 0 and contains the calibration period performance
if (nrow(last_obs_score) > 0) {
  # Check if last_obs_score has siteID (it should if we added it to grouping columns)
  if (!"siteID" %in% names(last_obs_score)) {
    # If no siteID, we need to expand to all sites for each model
    # Get all siteIDs from fcast_horizon_model_mean for models in last_obs_score
    model_sites <- fcast_horizon_model_mean[model_id %in% unique(last_obs_score$model_id), 
                                           .(siteID = unique(siteID)), by = model_id]
    if (nrow(model_sites) > 0) {
      # Expand last_obs_score to include all sites
      last_obs_score <- last_obs_score[model_sites, on = "model_id", allow.cartesian = TRUE]
      last_obs_score[, siteID := i.siteID]
      last_obs_score[, i.siteID := NULL]
    } else {
      cat("WARNING: Cannot expand last_obs_score to sites - no matching models in fcast_horizon_model_mean\n")
    }
  }
  
  if ("siteID" %in% names(last_obs_score)) {
    # Filter last_obs_score to only include models that are in fcast_horizon_model_mean
    last_obs_score_filtered <- last_obs_score[model_id %in% unique(fcast_horizon_model_mean$model_id)]
    
    # Ensure last_obs_score has the same structure as fcast_horizon_model_mean
    # Add missing columns if needed
    model_mean_cols <- names(fcast_horizon_model_mean)
    last_obs_cols <- names(last_obs_score_filtered)
    missing_cols <- setdiff(model_mean_cols, last_obs_cols)
    
    # Add missing columns with NA values
    for (col in missing_cols) {
      if (!col %in% c("model_id", "siteID", "months_since_obs", "species", "model_name", "pretty_group", "rank_only")) {
        last_obs_score_filtered[, (col) := NA_real_]
      }
    }
    
    # Ensure RSQ column exists (may be RSQ.1 in last_obs_score)
    if (!"RSQ" %in% names(last_obs_score_filtered) && "RSQ.1" %in% names(last_obs_score_filtered)) {
      last_obs_score_filtered[, RSQ := RSQ.1]
    }
    
    # Combine month 0 data with hindcast data
    fcast_horizon_model_mean <- rbindlist(list(fcast_horizon_model_mean, last_obs_score_filtered), fill = TRUE)
    cat("Added", nrow(last_obs_score_filtered), "month 0 data points to fcast_horizon_model_mean\n")
  }
}

to_plot_null <- copy(fcast_horizon_null_all[null_type == "site_mean" & model_id %in% unique(to_plot$model_id)])
# Strip null_ prefix for backward compat with downstream code
tp_null_save_cols <- grep("^null_", names(to_plot_null), value = TRUE)
if (length(tp_null_save_cols) > 0) setnames(to_plot_null, tp_null_save_cols, gsub("^null_", "", tp_null_save_cols))
cat("Completed data combination\n")


saveRDS(list(to_plot, to_plot_null, fcast_horizon_null_all, fcast_horizon_model_mean),
				here("data/summary/fcast_horizon_input.rds"))

# Clean up large objects before loading back
# Keep fcast_horizon_df for site-level aggregation in horizon calculation
rm(fcast_horizon_df2, comparison_scores, historical_mean)
suppressWarnings(rm(persist_null_site, clim_null_site, persistence_null, climatology_null))
gc(verbose = FALSE, full = TRUE)

in_list <- readRDS(here("data/summary/fcast_horizon_input.rds"))
to_plot <- as.data.table(in_list[[1]])
to_plot_null <- as.data.table(in_list[[2]])
# Strip null_ prefix from to_plot_null for backward compatibility with plot code
tp_null_cols <- grep("^null_", names(to_plot_null), value = TRUE)
if (length(tp_null_cols) > 0) setnames(to_plot_null, tp_null_cols, gsub("^null_", "", tp_null_cols))
fcast_horizon_null_all <- as.data.table(in_list[[3]])
fcast_horizon_model_mean <- as.data.table(in_list[[4]])
rm(in_list); gc(verbose = FALSE)

# For backward compat, fcast_horizon_null_site has null_ prefix stripped (used in per-model loop)
fcast_horizon_null_site <- copy(fcast_horizon_null_all[null_type == "site_mean"])
null_cols_to_rename <- grep("^null_", names(fcast_horizon_null_site), value = TRUE)
setnames(fcast_horizon_null_site, null_cols_to_rename, gsub("^null_", "", null_cols_to_rename))

model_id_list <- unique(to_plot$model_id)

# Filter to only weakly converged models (convergence IDs already had
# _beta_regression suffix stripped at line 224 to match hindcast model_id format)
if (!is.null(converged) && length(converged) > 0) {
  model_id_list <- model_id_list[model_id_list %in% converged]
  cat("Filtered to", length(model_id_list), "converged models\n")
}

# Configuration: Model selection for forecast horizon calculation
# Options: "loess" (default), "exponential" (GLM with log link), "gam" (GAM with mgcv), or "all" (test all methods)
FCAST_HORIZON_METHOD <- "all"  # Options: "loess", "exponential", "gam", "all"
USE_EXPONENTIAL_DECAY <- TRUE  # Enable exponential decay as an option (used when METHOD includes it)

# Horizon pooling: "pooled" (pool across sites, default), "per_site" (per-site then median), "both"
HORIZON_POOLING <- "pooled"

# Test mode: Only process a few models for testing (set FALSE for full runs with both kingdoms)
TEST_MODE <- FALSE
TEST_MODEL_COUNT <- 5  # Number of models to test per group when TEST_MODE is TRUE

# Test mode: Limit to a few models for testing; sample from both kingdoms so Fungi is included
if (TEST_MODE) {
  cat("\n=== TEST MODE ENABLED ===\n")
  cat("Processing only", TEST_MODEL_COUNT, "models per group for testing\n")
  model_pg <- unique(to_plot[, .(model_id, pretty_group)], by = "model_id")
  model_types <- sapply(model_id_list, function(x) strsplit(x, "_")[[1]][1:2])
  model_type_str <- sapply(model_types, function(x) paste(x, collapse = "_"))
  unique_types <- unique(model_type_str)
  test_models <- c()
  n_per_kingdom <- max(1L, floor(TEST_MODEL_COUNT / 2))
  for (type in unique_types) {
    type_models <- model_id_list[model_type_str == type]
    type_bac <- intersect(type_models, model_pg[pretty_group == "Bacteria"]$model_id)
    type_fun <- intersect(type_models, model_pg[pretty_group == "Fungi"]$model_id)
    # Take from both kingdoms, then fill to TEST_MODEL_COUNT
    picked <- c(head(type_bac, n_per_kingdom), head(type_fun, n_per_kingdom))
    remaining <- setdiff(type_models, picked)
    picked <- c(picked, head(remaining, max(0L, TEST_MODEL_COUNT - length(picked))))
    test_models <- c(test_models, picked)
  }
  model_id_list <- unique(test_models)
  cat("Selected", length(model_id_list), "test models\n")
  if (nrow(model_pg) > 0) {
    sel_pg <- model_pg[model_id %in% model_id_list, .N, by = pretty_group]
    cat("  By pretty_group:", paste(sel_pg[, paste(pretty_group, N, sep = "=")], collapse = ", "), "\n")
  }
  cat("Test models:", paste(head(model_id_list, 10), collapse = ", "), if (length(model_id_list) > 10) " ...\n" else "\n")
}

# Use all available models (or test subset)
cat("\n=== MODEL PROCESSING SUMMARY ===\n")
cat("Total models to process:", length(model_id_list), "\n")
if (length(model_id_list) > 0) {
cat("First 10 model_ids:", paste(head(model_id_list, 10), collapse=", "), if (length(model_id_list)>10) " ...\n" else "\n")
  
  # Check species diversity
  if ("species" %in% names(to_plot)) {
    unique_species <- unique(to_plot$species[!is.na(to_plot$species)])
    cat("Unique species in to_plot:", length(unique_species), "\n")
    if (length(unique_species) <= 5) {
      cat("  Species:", paste(unique_species, collapse=", "), "\n")
    } else {
      cat("  First 10 species:", paste(head(unique_species, 10), collapse=", "), " ...\n")
    }
  }
  
  # Check pretty_group diversity
  if ("pretty_group" %in% names(to_plot)) {
    unique_groups <- unique(to_plot$pretty_group[!is.na(to_plot$pretty_group)])
    cat("Pretty groups in to_plot:", paste(unique_groups, collapse=", "), "\n")
    cat("  Count by group:\n")
    print(table(to_plot$pretty_group, useNA="ifany"))
  }
} else {
  cat("⚠️  WARNING: No models found in to_plot! This will cause the script to fail.\n")
  cat("  Checking data sources...\n")
  if (exists("comparison_scores")) {
    cat("  comparison_scores rows:", nrow(comparison_scores), "\n")
    if (nrow(comparison_scores) > 0 && "model_id" %in% names(comparison_scores)) {
      cat("  comparison_scores unique models:", length(unique(comparison_scores$model_id)), "\n")
    }
  }
  if (exists("fcast_horizon_model_mean")) {
    cat("  fcast_horizon_model_mean rows:", nrow(fcast_horizon_model_mean), "\n")
    if (nrow(fcast_horizon_model_mean) > 0 && "model_id" %in% names(fcast_horizon_model_mean)) {
      cat("  fcast_horizon_model_mean unique models:", length(unique(fcast_horizon_model_mean$model_id)), "\n")
    }
  }
}
cat("================================\n\n")

# Group models by type (e.g., env_cycl, cycl_only, env_cov) for memory-efficient batch processing
extract_model_type <- function(model_id) {
  # Extract model type from model_id (first part before first underscore or number)
  parts <- strsplit(model_id, "_")[[1]]
  if (length(parts) >= 2 && parts[1] == "env" && parts[2] == "cycl") {
    return("env_cycl")
  } else if (length(parts) >= 2 && parts[1] == "env" && parts[2] == "cov") {
    return("env_cov")
  } else if (length(parts) >= 2 && parts[1] == "cycl" && parts[2] == "only") {
    return("cycl_only")
  } else if (length(parts) >= 1) {
    return(parts[1])
  } else {
    return("other")
  }
}

# Group model_ids by type
model_groups <- split(model_id_list, sapply(model_id_list, extract_model_type))
cat("\nGrouping models by type for batch processing:\n")
for (group_name in names(model_groups)) {
  cat("  ", group_name, ":", length(model_groups[[group_name]]), "models\n")
}
cat("Total groups:", length(model_groups), "\n\n")

out_df_list = list()
out_figure_list = list()
out_parallel_list = list()
setDT(to_plot)
setDT(fcast_horizon_null_site)

# Process each model group separately to reduce memory usage
for (group_idx in seq_along(model_groups)) {
  group_name <- names(model_groups)[group_idx]
  group_model_ids <- model_groups[[group_name]]
  
  cat("\n=== Processing group", group_idx, "of", length(model_groups), ":", group_name, "===\n")
  cat("  Models in this group:", length(group_model_ids), "\n")
  
  # Load hindcast data for this group only
  cat("  Loading hindcast data for this group...\n")
  parquet_path <- here('data/summary/parquet/all_hindcasts_plsr2.parquet')
  rds_path <- here('data/summary/all_hindcasts_plsr2.rds')
  
  if (file.exists(parquet_path)) {
    if (requireNamespace("nanoparquet", quietly = TRUE)) {
      cat("    Using Parquet file with nanoparquet...\n")
      hindcast_data_all <- as.data.table(nanoparquet::read_parquet(parquet_path))
    } else if (requireNamespace("arrow", quietly = TRUE)) {
      cat("    Using Parquet file with arrow...\n")
      hindcast_data_all <- as.data.table(arrow::read_parquet(parquet_path))
    } else {
      cat("    WARNING: Parquet file exists but no reader available. Using RDS...\n")
      if (file.exists(rds_path)) {
        hindcast_data_all <- readRDS(rds_path)
        setDT(hindcast_data_all)
      } else {
        stop("No hindcast data file found!")
      }
    }
  } else if (file.exists(rds_path)) {
    cat("    WARNING: Using full RDS file (2.7GB) - may hit memory limits...\n")
    hindcast_data_all <- readRDS(rds_path)
    setDT(hindcast_data_all)
  } else {
    stop("No hindcast data file found!")
  }
  
  setDT(hindcast_data_all)
  hindcast_data_all[, model_id := gsub('_beta_regression$', '', model_id)]
  hindcast_data_all[, corrected_timepoint := ifelse(is.na(timepoint) & !is.na(dateID), dateID, timepoint)]
  
  # Filter to only models in this group
  hindcast_data_all <- hindcast_data_all[model_id %in% group_model_ids]
  cat("    Filtered to", nrow(hindcast_data_all), "rows for", length(group_model_ids), "models\n")
  
  # Calculate last observation dates for models in this group
  # CRITICAL: Only use calibration data - models without calibration are excluded
  cat("  Calculating last observation dates for this group...\n")
  hindcast_cal_all <- hindcast_data_all[fcast_period != 'hindcast' & !is.na(truth)]
if (nrow(hindcast_cal_all) > 0) {
  # Use calibration data if available
  last_obs_all <- hindcast_cal_all[, .(last_obs = ifelse(any(!is.na(dateID) & dateID >= 10000), 
                                                         max(dateID[!is.na(dateID) & dateID >= 10000], na.rm = TRUE),
                                                         max(corrected_timepoint, na.rm = TRUE))), 
                                   by = .(model_id, plotID, siteID)]
} else {
  # If no calibration data, use the earliest hindcast point as last_obs
  cat("  Warning: No calibration data found. Using earliest hindcast point as last_obs.\n")
  hindcast_only_all <- hindcast_data_all[fcast_period == 'hindcast' & !is.na(truth)]
  if (nrow(hindcast_only_all) > 0) {
    last_obs_all <- hindcast_only_all[, .(last_obs = ifelse(any(!is.na(dateID) & dateID >= 10000), 
                                                            min(dateID[!is.na(dateID) & dateID >= 10000], na.rm = TRUE),
                                                            min(corrected_timepoint, na.rm = TRUE))), 
                                         by = .(model_id, plotID, siteID)]
  } else {
    last_obs_all <- data.table(model_id = character(), plotID = character(), siteID = character(), last_obs = numeric())
  }
}

# Get last observation data points
setkey(hindcast_data_all, plotID, model_id, dateID)
last_obs_all[, dateID := last_obs]
setkey(last_obs_all, plotID, model_id, dateID)
last_obs_values_all <- hindcast_data_all[last_obs_all, on = c('plotID', 'model_id', 'dateID'), nomatch = 0]

# Calculate CRPS if needed — use truncated CRPS (FIX 2)
if (nrow(last_obs_values_all) > 0 && "sd" %in% names(last_obs_values_all)) {
  cat("  Calculating truncated CRPS for last observation points...\n")
  last_obs_values_all[, crps := NA_real_]

  valid_crps <- !is.na(last_obs_values_all$truth) &
                !is.na(last_obs_values_all$mean) &
                !is.na(last_obs_values_all$sd)
  if (sum(valid_crps) > 0) {
    chunk_size <- 10000
    n_chunks <- ceiling(sum(valid_crps) / chunk_size)
    valid_indices <- which(valid_crps)

    for (i in 1:n_chunks) {
      start_idx <- (i-1) * chunk_size + 1
      end_idx <- min(i * chunk_size, sum(valid_crps))
      chunk_indices <- valid_indices[start_idx:end_idx]

      if (i %% 10 == 0) cat("    Processing CRPS chunk", i, "of", n_chunks, "\n")
      # FIX 2: Use truncated CRPS [0,1] instead of unbounded crps_norm
      valid_chunk <- chunk_indices[last_obs_values_all[chunk_indices, sd > 0 & is.finite(sd)]]
      if (length(valid_chunk) > 0) {
        last_obs_values_all[valid_chunk, crps := crps(truth, family = "tnorm",
                                                       location = mean, scale = sd,
                                                       lower = 0, upper = 1)]
      }
    }
  }
}

# Pre-calculate site-level metrics for last observation points
if (nrow(last_obs_values_all) > 0 && "siteID" %in% names(last_obs_values_all)) {
  cat("  Aggregating last observation metrics by site...\n")
  last_obs_values_site_precalc <- last_obs_values_all[, {
    valid <- !is.na(truth) & !is.na(mean)
    if (sum(valid) > 0) {
      rmse_val <- rmse(actual = truth[valid], predicted = mean[valid])
      # FIX 1: Use RSQ.1 (Nash-Sutcliffe) instead of lm R²
      rsq1_val <- if (!is.na(rmse_val)) {
        var_val <- var(truth[valid], na.rm = TRUE)
        if (!is.na(var_val) && var_val > 0) max(0, 1 - (rmse_val^2)/var_val) else NA_real_
      } else NA_real_
      crps_val <- if ("crps" %in% names(last_obs_values_all)) {
        mean(crps[valid], na.rm = TRUE)
      } else NA_real_
      abund <- mean(truth[valid], na.rm=TRUE)
      # FIX 3: Apply 0.005 abundance floor
      abund_denom <- max(abund, 0.005)
      rmse_norm_val <- if (!is.na(rmse_val)) rmse_val/abund_denom else NA_real_

      list(
        RSQ = rsq1_val,
        mean_crps = crps_val,
        RMSE.norm = rmse_norm_val,
        months_since_obs = 0
      )
    } else {
      list(RSQ = NA_real_, mean_crps = NA_real_, RMSE.norm = NA_real_, months_since_obs = 0)
    }
  }, by = .(model_id, siteID)]
  
  cat("  Pre-calculated last observation points for", length(unique(last_obs_values_site_precalc$model_id)), 
      "models and", length(unique(last_obs_values_site_precalc$siteID)), "sites\n")
  }
  
  # Clean up large objects before processing individual models
  rm(hindcast_data_all, hindcast_cal_all, last_obs_all, last_obs_values_all)
  gc(verbose = FALSE)
  
  # Process each model in this group
  # CRITICAL: Filter to only models that have last_obs (calibration data)
  # Models without calibration data are excluded from horizon calculation
  models_with_last_obs <- if (exists("last_obs") && nrow(last_obs) > 0) {
    unique(last_obs$model_id)
  } else {
    character(0)
  }
  
  # Filter group_model_ids to only include models with last_obs
  group_model_ids_valid <- intersect(group_model_ids, models_with_last_obs)
  group_model_ids_excluded <- setdiff(group_model_ids, models_with_last_obs)
  
  if (length(group_model_ids_excluded) > 0) {
    cat("  Excluding", length(group_model_ids_excluded), "models without calibration data from horizon calculation\n")
  }
  
  for(current_model_id in group_model_ids_valid) {
  cat("Processing:", current_model_id, "\n")
  cat("Starting processing for:", current_model_id, "\n")
  gam_skip_warned_this_model <- FALSE
	
single_tax = to_plot[model_id == current_model_id]

# Skip if no data for this model_id
if (nrow(single_tax) == 0) {
  cat("Warning: No data for model_id", current_model_id, "- skipping\n")
  next
}

# Skip if insufficient data for reliable curve fitting
if (nrow(single_tax) < 4) {
  cat("Warning: Insufficient data for model_id", current_model_id, "(", nrow(single_tax), "rows) - skipping\n")
  next
}

cat("Found", nrow(single_tax), "rows for:", current_model_id, "\n")

setDT(single_tax)
if (!"mean_crps" %in% names(single_tax)) {
  stop("single_tax must have mean_crps. Ensure fcast_horizon_model_mean and last_obs_score produce mean_crps.")
}

single_tax_null <- fcast_horizon_null_site[fcast_horizon_null_site$model_id == current_model_id, ]

# Skip if no null data for this model_id
if (nrow(single_tax_null) == 0) {
  cat("Warning: No null data for model_id", current_model_id, "- skipping\n")
  next
}

setDT(single_tax_null)

# Collapse null baselines to single values (median across sites)
# After FIX 1, RSQ and RSQ.1 are both Nash-Sutcliffe in the null; RSQ is used for horizon crossing
single_tax_null_collapse <- single_tax_null[, .(RSQ = median(RSQ, na.rm = TRUE),
                                                RSQ.1 = median(RSQ.1, na.rm = TRUE),
                                                mean_crps = median(mean_crps, na.rm = TRUE),
                                                RMSE.norm = median(RMSE.norm, na.rm = TRUE))]

# FIX 8: Only adjust horizon values for modeling/plotting, not storage
# Exclude null scores (Inf months_since_obs) from curve fitting - they're only for baseline
single_tax_for_fit <- copy(single_tax[is.finite(months_since_obs)])

# Debug: Check data at months_since_obs = 0 before conversion (for first model)
debug_this_model_now <- (current_model_id == model_id_list[1])
if (debug_this_model_now) {
  cat("\n=== CHECKING DATA AT months_since_obs = 0 ===\n")
  zero_data_before <- single_tax_for_fit[months_since_obs == 0]
  cat("Rows at months_since_obs = 0 (before conversion):", nrow(zero_data_before), "\n")
  if (nrow(zero_data_before) > 0) {
    cat("Sample data at 0:\n")
    print(head(zero_data_before[, .(months_since_obs, RSQ, mean_crps, RMSE.norm)], 3))
  }
}

single_tax_for_fit[, months_since_obs := ifelse(months_since_obs == 0, 1e-3, months_since_obs)]

# Debug: Check data after conversion
if (debug_this_model_now) {
  zero_data_after <- single_tax_for_fit[abs(months_since_obs - 0.001) < 0.0001]
  cat("Rows at months_since_obs ~ 0.001 (after conversion):", nrow(zero_data_after), "\n")
  if (nrow(zero_data_after) > 0) {
    cat("Sample data at 0.001:\n")
    print(head(zero_data_after[, .(months_since_obs, RSQ, mean_crps, RMSE.norm)], 3))
  }
  cat("==================================================\n\n")
}

# ---- HEADLESS SMOOTH & HORIZON ----
# Check available columns and investigate if missing
available_cols <- names(single_tax_for_fit)
missing_cols <- c()
if (!"RSQ" %in% available_cols) missing_cols <- c(missing_cols, "RSQ")
if (!"mean_crps" %in% available_cols) missing_cols <- c(missing_cols, "mean_crps")
if (!"RMSE.norm" %in% available_cols) missing_cols <- c(missing_cols, "RMSE.norm")

if (length(missing_cols) > 0) {
  cat("INVESTIGATING: Missing required columns for", current_model_id, ":", paste(missing_cols, collapse=", "), "\n")
  cat("  Available columns:", paste(available_cols, collapse=", "), "\n")
  cat("  Checking source data (to_plot) columns:", paste(names(to_plot), collapse=", "), "\n")
  cat("  Checking fcast_horizon_model_mean columns:", paste(names(fcast_horizon_model_mean), collapse=", "), "\n")
  cat("  Checking comparison_scores columns:", 
      if(exists("comparison_scores")) paste(names(comparison_scores), collapse=", ") else "not available", "\n")
  
  # Check for alternative column names
  if ("RSQ.1" %in% available_cols && "RSQ" %in% missing_cols) {
    cat("  Found RSQ.1 as alternative to RSQ, using it\n")
    single_tax_for_fit[, RSQ := RSQ.1]
    missing_cols <- setdiff(missing_cols, "RSQ")
  }
  
  if (length(missing_cols) > 0) {
    cat("  ERROR: Cannot proceed without:", paste(missing_cols, collapse=", "), "\n")
    cat("  This suggests a problem in data combination or column naming\n")
    next
  }
}

# Filter to reasonable forecast horizons (0-20 months) for fitting
# Keep plot-level data (each plot contributes one point per timepoint)
# We need plot-level data with siteID to fit separate loess curves per site
if (!exists("fcast_horizon_df") || is.null(fcast_horizon_df)) {
  # fcast_horizon_df not available - need to get plot-level data from to_plot
  # But to_plot is already aggregated, so we can't get plot-level data
  cat("Warning: fcast_horizon_df not available and to_plot is aggregated - cannot fit per-site curves for", current_model_id, "\n")
  cat("  Falling back to single curve across all sites\n")
  single_tax_plot_level <- single_tax_for_fit[months_since_obs <= 20 & 
                                                is.finite(months_since_obs) & 
                                                months_since_obs >= 0]
} else {
  # Use fcast_horizon_df to get plot-level data with siteID
  fcast_horizon_model_subset <- fcast_horizon_df[model_id == current_model_id & 
                                                  months_since_obs <= 20 & 
                                                  is.finite(months_since_obs) & 
                                                  months_since_obs >= 0]
  
  if (nrow(fcast_horizon_model_subset) == 0) {
    cat("Warning: No data in fcast_horizon_df for", current_model_id, "- skipping\n")
    next
  }
  
  # Calculate metrics at plot level (each plot contributes one point per timepoint)
  if (!"siteID" %in% names(fcast_horizon_model_subset)) {
    # Extract siteID from plotID if available
    if ("plotID" %in% names(fcast_horizon_model_subset)) {
      fcast_horizon_model_subset[, siteID := gsub("_.*$", "", plotID)]
    } else {
      cat("Warning: No siteID or plotID in fcast_horizon_df for", current_model_id, "- skipping\n")
      next
    }
  }
  
  # Calculate truncated CRPS if not available (FIX 2)
  if (!"crps" %in% names(fcast_horizon_model_subset) && "sd" %in% names(fcast_horizon_model_subset)) {
    valid_crps <- !is.na(fcast_horizon_model_subset$truth) &
                  !is.na(fcast_horizon_model_subset$mean) &
                  !is.na(fcast_horizon_model_subset$sd) &
                  fcast_horizon_model_subset$sd > 0
    if (sum(valid_crps) > 0) {
      fcast_horizon_model_subset[valid_crps, crps := crps(truth, family = "tnorm",
                                                           location = mean, scale = sd,
                                                           lower = 0, upper = 1)]
    }
  }

  # Calculate metrics per site per timepoint (aggregate across all plots at that site/timepoint)
  single_tax_plot_level <- fcast_horizon_model_subset[, {
    valid <- !is.na(truth) & !is.na(mean)
    if (sum(valid) > 0) {
      rmse_val <- rmse(actual = truth[valid], predicted = mean[valid])
      # FIX 1: Use RSQ.1 (Nash-Sutcliffe) as RSQ
      rsq1_val <- if (!is.na(rmse_val)) {
        var_val <- var(truth[valid], na.rm = TRUE)
        if (!is.na(var_val) && var_val > 0) max(0, 1 - (rmse_val^2)/var_val) else NA_real_
      } else NA_real_
      abund <- mean(truth, na.rm = TRUE)
      # FIX 3: Apply 0.005 abundance floor
      abund_denom <- max(abund, 0.005)
      rmse_norm_val <- if (!is.na(rmse_val)) rmse_val/abund_denom else NA_real_

      crps_val <- if ("crps" %in% names(fcast_horizon_model_subset)) {
        mean(crps[valid], na.rm = TRUE)
      } else NA_real_

      list(
        RSQ = rsq1_val,
        mean_crps = crps_val,
        RMSE.norm = rmse_norm_val
      )
    } else {
      list(RSQ = NA_real_, mean_crps = NA_real_, RMSE.norm = NA_real_)
    }
  }, by = .(siteID, months_since_obs)]
  
  # Add last observation scores per site from pre-calculated data
  # These represent the best performance at the last observation point (months_since_obs = 0)
  if (nrow(last_obs_values_site_precalc) > 0) {
    # Get pre-calculated last observation points for this model
    last_obs_site_level <- last_obs_values_site_precalc[model_id == current_model_id & 
                                                        siteID %in% unique(single_tax_plot_level$siteID)]
    
    if (nrow(last_obs_site_level) > 0) {
      # Only add sites that don't already have months_since_obs = 0
      existing_sites_at_zero <- single_tax_plot_level[months_since_obs == 0, unique(siteID)]
      last_obs_site_level <- last_obs_site_level[!siteID %in% existing_sites_at_zero & 
                                                 (is.finite(RSQ) | is.finite(mean_crps) | is.finite(RMSE.norm))]
      if (nrow(last_obs_site_level) > 0) {
        single_tax_plot_level <- rbind(single_tax_plot_level, last_obs_site_level[, .(siteID, RSQ, mean_crps, RMSE.norm, months_since_obs)], fill = TRUE)
        cat("Added", nrow(last_obs_site_level), "last observation points (months_since_obs = 0) from pre-calculated data\n")
      }
    }
  }
}

# Check if we have enough data points
if (nrow(single_tax_plot_level) < 2) {
  cat("Warning: Insufficient data points for horizon calculation (", nrow(single_tax_plot_level), "rows) - skipping\n")
  next
}

# Check if we have siteID
if (!"siteID" %in% names(single_tax_plot_level)) {
  cat("Warning: siteID not available for per-site fitting for", current_model_id, "- skipping\n")
  next
}

# Debug: Check data before fitting
if (debug_this_model_now) {
  cat("\n=== CHECKING single_tax_plot_level (AGGREGATED BY SITE) ===\n")
  cat("Total rows (siteID x months_since_obs):", nrow(single_tax_plot_level), "\n")
  cat("Unique siteIDs:", length(unique(single_tax_plot_level$siteID)), "\n")
  cat("Unique months_since_obs values:", length(unique(single_tax_plot_level$months_since_obs)), "\n")
  cat("months_since_obs range:", range(single_tax_plot_level$months_since_obs, na.rm=TRUE), "\n")
  cat("Valid data points:\n")
  cat("  RSQ finite:", sum(is.finite(single_tax_plot_level$RSQ)), "of", nrow(single_tax_plot_level), "\n")
  cat("  mean_crps finite:", sum(is.finite(single_tax_plot_level$mean_crps)), "of", nrow(single_tax_plot_level), "\n")
  cat("  RMSE.norm finite:", sum(is.finite(single_tax_plot_level$RMSE.norm)), "of", nrow(single_tax_plot_level), "\n")
  cat("Sample data:\n")
  print(head(single_tax_plot_level[order(months_since_obs), .(siteID, months_since_obs, RSQ, mean_crps, RMSE.norm)], 15))
  cat("==================================================\n\n")
}

# baselines (scalars from your collapsed null) - calculate before per-site fitting
if (nrow(single_tax_null_collapse) == 0) {
  cat("Warning: single_tax_null_collapse is empty for", current_model_id, "\n")
  next
}

rsq_null_line  <- if (!is.na(single_tax_null_collapse$RSQ[1]) && is.finite(single_tax_null_collapse$RSQ[1])) {
  single_tax_null_collapse$RSQ[1]
} else {
  if (nrow(single_tax_null) > 0 && "RSQ" %in% names(single_tax_null)) {
    median(single_tax_null$RSQ, na.rm = TRUE)
  } else {
    NA_real_
  }
}

crps_null_line <- if (!is.na(single_tax_null_collapse$mean_crps[1]) && is.finite(single_tax_null_collapse$mean_crps[1])) {
  single_tax_null_collapse$mean_crps[1]
} else {
  if (nrow(single_tax_null) > 0 && "mean_crps" %in% names(single_tax_null)) {
    median(single_tax_null$mean_crps, na.rm = TRUE)
  } else {
    NA_real_
  }
}

rmse_null_line <- if (!is.na(single_tax_null_collapse$RMSE.norm[1]) && is.finite(single_tax_null_collapse$RMSE.norm[1])) {
  single_tax_null_collapse$RMSE.norm[1]
} else {
  if (nrow(single_tax_null) > 0 && "RMSE.norm" %in% names(single_tax_null)) {
    median(single_tax_null$RMSE.norm, na.rm = TRUE)
  } else {
    NA_real_
  }
}

# ── POOLED HORIZON COMPUTATION ───────────────────────────────────────────────
# Pool all sites for this model, fit single loess/GAM on pooled data
if (HORIZON_POOLING %in% c("pooled", "both")) {

  # Get null lines for each null type (median across sites)
  get_null_lines <- function(null_type_str) {
    nt <- fcast_horizon_null_all[null_type == null_type_str & model_id == current_model_id]
    if (nrow(nt) == 0) return(list(rsq = NA_real_, crps = NA_real_, rmse = NA_real_))
    # Handle both null_-prefixed and plain column names
    get_col <- function(dt, candidates) {
      for (c in candidates) {
        if (c %in% names(dt)) return(median(dt[[c]], na.rm = TRUE))
      }
      NA_real_
    }
    list(
      rsq = get_col(nt, c("null_RSQ", "RSQ")),
      crps = get_col(nt, c("null_mean_crps", "null_CRPS_truncated", "mean_crps")),
      rmse = get_col(nt, c("RMSE.norm"))
    )
  }

  null_sm <- get_null_lines("site_mean")
  null_persist <- get_null_lines("persistence")
  null_clim <- get_null_lines("climatological")

  # Use fcast_horizon_model_mean which has per-site-per-timepoint aggregated metrics
  pooled_data <- fcast_horizon_model_mean[model_id == current_model_id &
                                           is.finite(months_since_obs) &
                                           months_since_obs >= 0 & months_since_obs <= 20]

  # Filter out months with fewer than 3 sites to avoid single-site artifacts.
  # E.g., month 3 may only have data from 1 tropical site, which can dominate
  # the loess fit and produce artificially short horizons.
  if (nrow(pooled_data) > 0 && "siteID" %in% names(pooled_data)) {
    sites_per_month <- pooled_data[months_since_obs > 0, .(n_sites = uniqueN(siteID)), by = months_since_obs]
    sparse_months <- sites_per_month[n_sites < 3]$months_since_obs
    if (length(sparse_months) > 0) {
      pooled_data <- pooled_data[!(months_since_obs %in% sparse_months)]
    }
  }

  # Add month 0 from last_obs if not present
  if (exists("last_obs_values_site_precalc") && nrow(last_obs_values_site_precalc) > 0) {
    m0 <- last_obs_values_site_precalc[model_id == current_model_id]
    if (nrow(m0) > 0 && nrow(pooled_data[months_since_obs == 0]) == 0) {
      pooled_data <- rbind(pooled_data, m0[, .(siteID, RSQ, mean_crps, RMSE.norm, months_since_obs)], fill = TRUE)
    }
  }

  pooled_loess_horizon <- function(metric_vals, time_vals, null_val, direction) {
    ok <- is.finite(metric_vals) & is.finite(time_vals)
    if (sum(ok) < 4 || is.na(null_val) || !is.finite(null_val)) return(NA_real_)

    # Only search for crossings at or after the first real hindcast month.
    # Month 0 is the last calibration observation (anchor for the curve fit)
    # but crossings in the interpolation gap (months 0-3) are artifacts.
    first_hindcast <- min(time_vals[ok & time_vals > 0], na.rm = TRUE)
    if (!is.finite(first_hindcast)) return(NA_real_)

    fit <- tryCatch(loess(metric_vals[ok] ~ time_vals[ok], span = 0.75,
                          control = loess.control(surface = "direct")),
                    error = function(e) NULL)
    if (is.null(fit)) {
      # Discrete crossing fallback
      if (direction == "below") {
        cr <- which(metric_vals[ok] < null_val & time_vals[ok] >= first_hindcast)
      } else {
        cr <- which(metric_vals[ok] > null_val & time_vals[ok] >= first_hindcast)
      }
      return(if (length(cr) > 0) time_vals[ok][min(cr)] else NA_real_)
    }
    xg <- seq(min(time_vals[ok]), max(time_vals[ok]), by = 0.5)
    pred <- predict(fit, xg)
    if (direction == "below") {
      cr <- which(pred < null_val & xg >= first_hindcast)
    } else {
      cr <- which(pred > null_val & xg >= first_hindcast)
    }
    if (length(cr) > 0) xg[min(cr)] else NA_real_
  }

  if (nrow(pooled_data) >= 4) {
    # Compute pooled horizons for each null type
    compute_pooled_for_null <- function(nl) {
      list(
        rsq  = pooled_loess_horizon(pooled_data$RSQ, pooled_data$months_since_obs, nl$rsq, "below"),
        crps = pooled_loess_horizon(pooled_data$mean_crps, pooled_data$months_since_obs, nl$crps, "above"),
        rmse = pooled_loess_horizon(pooled_data$RMSE.norm, pooled_data$months_since_obs, nl$rmse, "above")
      )
    }

    ph_sm <- compute_pooled_for_null(null_sm)
    ph_persist <- compute_pooled_for_null(null_persist)
    ph_clim <- compute_pooled_for_null(null_clim)
  } else {
    # Not enough data points for loess fit; return NA instead of a misleading default
    ph_sm <- ph_persist <- ph_clim <- list(rsq = NA_real_, crps = NA_real_, rmse = NA_real_)
  }
}

# If pooled-only, set the main horizon variables and skip per-site loop
if (HORIZON_POOLING == "pooled") {
  # Primary horizons use site_mean null (backward compatible)
  rsq_fcast_horizon <- if (is.na(ph_sm$rsq)) NA_real_ else min(ph_sm$rsq, 20)
  crps_fcast_horizon <- if (is.na(ph_sm$crps)) NA_real_ else min(ph_sm$crps, 20)
  rmse_fcast_horizon <- if (is.na(ph_sm$rmse)) NA_real_ else min(ph_sm$rmse, 20)

  # GAM/LM columns: use persistence and climatological null horizons
  rsq_fcast_horizon_gam <- if (is.na(ph_persist$rsq)) NA_real_ else min(ph_persist$rsq, 20)
  crps_fcast_horizon_gam <- if (is.na(ph_persist$crps)) NA_real_ else min(ph_persist$crps, 20)
  rmse_fcast_horizon_gam <- if (is.na(ph_persist$rmse)) NA_real_ else min(ph_persist$rmse, 20)
  rsq_fcast_horizon_lm <- if (is.na(ph_clim$rsq)) NA_real_ else min(ph_clim$rsq, 20)
  crps_fcast_horizon_lm <- if (is.na(ph_clim$crps)) NA_real_ else min(ph_clim$crps, 20)
  rmse_fcast_horizon_lm <- if (is.na(ph_clim$rmse)) NA_real_ else min(ph_clim$rmse, 20)

  # Store null lines (site_mean for primary)
  rsq_null_line <- null_sm$rsq
  crps_null_line <- null_sm$crps
  rmse_null_line <- null_sm$rmse

  # Skip straight to output assembly
  goto_output <- TRUE
} else {
  goto_output <- FALSE
}

if (!goto_output) {
# ── PER-SITE HORIZON COMPUTATION (original approach) ────────────────────────

# crossing helper: first x where y crosses thresh (dir: below = y < thresh, above = y > thresh).
# min_x: only consider crossings at x > min_x (avoids horizon 0; use 0 to require first month after last obs).
first_crossing <- function(x, y, thresh, dir = c("below","above"), min_x = 0) {
  dir <- match.arg(dir)
  if (is.na(thresh) || !is.finite(thresh)) return(Inf)
  ok  <- is.finite(x) & is.finite(y)
  if (!any(ok)) return(Inf)
  x <- x[ok]; y <- y[ok]
  if (length(x) == 0 || length(y) == 0) return(Inf)
  idx <- if (dir == "below") {
    which(y < thresh)
  } else {
    which(y > thresh)
  }
  if (length(idx) == 0) return(Inf)
  candidates <- x[idx]
  if (is.finite(min_x) && min_x >= 0) candidates <- candidates[candidates > min_x]
  if (length(candidates) == 0) return(Inf)
  result <- min(candidates, na.rm = TRUE)
  if (!is.finite(result)) return(Inf)
  result
}

# define a common prediction grid
xgrid <- seq(0, 20, by = 0.1)

safe_loess <- function(y, x) {
  # need at least 4 unique x for loess to behave; fallback to linear if not
  valid <- is.finite(x) & is.finite(y)
  if (sum(valid) < 2) return(NULL)
  ux <- unique(x[valid])
  if (length(ux) >= 4 && sum(valid) >= 4 && sd(x[valid], na.rm=TRUE) > 0) {
    tryCatch({
      suppressWarnings(stats::loess(y[valid] ~ x[valid], span = 0.8, control = loess.control(surface = "direct")))
    }, error = function(e) {
      if (length(ux) >= 2 && sd(x[valid], na.rm=TRUE) > 0) {
        stats::lm(y[valid] ~ x[valid])  # fallback
      } else {
        NULL
      }
    })
  } else if (length(ux) >= 2 && sd(x[valid], na.rm=TRUE) > 0) {
    stats::lm(y[valid] ~ x[valid])  # fallback
  } else {
    NULL
  }
}

# GAM fit using mgcv
safe_gam <- function(y, x) {
  if (!requireNamespace("mgcv", quietly = TRUE)) return(NULL)
  valid <- is.finite(x) & is.finite(y)
  if (sum(valid) < 3) return(NULL)  # GAM needs at least 3 points
  ux <- unique(x[valid])
  if (length(ux) < 2 || sd(x[valid], na.rm=TRUE) == 0) return(NULL)
  
  y_vals <- y[valid]
  x_vals <- x[valid]
  
  # Check for sufficient data variation
  if (sd(y_vals, na.rm=TRUE) == 0) {
    # All y values are the same, GAM won't work well
    return(NULL)
  }
  
  tryCatch({
    # Use a smooth term with appropriate basis dimension
    # k controls the basis dimension (max degrees of freedom)
    # For small datasets, use fewer basis functions
    k_val <- min(5, length(ux) - 1)  # At least one less than unique x values
    if (k_val < 2) k_val <- 2  # Minimum of 2
    
    # Create data frame with explicit column name for GAM
    df_gam <- data.frame(y_vals = y_vals, x_vals = x_vals)
    
    # Fit GAM with explicit variable names
    # Use "tp" (thin plate) basis which is more robust than "cs" for small datasets
    fit <- suppressWarnings(
      mgcv::gam(y_vals ~ s(x_vals, bs = "tp", k = k_val), data = df_gam, method = "REML")
    )
    
    # Verify the fit is valid and converged
    if (is.null(fit) || !inherits(fit, "gam")) {
      return(NULL)
    }
    
    # Check if model converged
    if (!is.null(fit$converged) && !fit$converged) {
      return(NULL)
    }
    
    return(fit)
  }, error = function(e) {
    # If GAM fails, try simpler model with fewer basis functions
    tryCatch({
      df_gam <- data.frame(y_vals = y_vals, x_vals = x_vals)
      # Try with even fewer basis functions
      k_simple <- min(3, length(ux) - 1)
      if (k_simple < 2) k_simple <- 2
      fit <- suppressWarnings(
        mgcv::gam(y_vals ~ s(x_vals, bs = "tp", k = k_simple), data = df_gam, method = "REML")
      )
      if (is.null(fit) || !inherits(fit, "gam")) {
        return(NULL)
      }
      if (!is.null(fit$converged) && !fit$converged) {
        return(NULL)
      }
      return(fit)
    }, error = function(e2) {
      NULL
    })
  })
}

# Exponential decay fit using glm with log link
# This forces an exponential decay shape: performance starts high and decays smoothly
# For RSQ: decreases over time (high to low/null)
# For CRPS/RMSE: increases over time (low to high/null)
safe_glm_exp <- function(y, x, null_val) {
  valid <- is.finite(x) & is.finite(y)
  if (sum(valid) < 2) return(NULL)
  ux <- unique(x[valid])
  if (length(ux) < 2 || sd(x[valid], na.rm=TRUE) == 0) return(NULL)
  
  y_vals <- y[valid]
  x_vals <- x[valid]
  
  # For RSQ: values can be 0 or negative, need to shift to positive for log link
  # For CRPS/RMSE: should be positive, but add small epsilon to be safe
  # Strategy: shift by (min - 1) to ensure all values are > 1 (safer for log link)
  y_min <- min(y_vals, na.rm = TRUE)
  if (y_min <= 0) {
    # Shift so minimum becomes 1.0 (safer for log link than just above 0)
    shift_amount <- 1.0 - y_min
  } else {
    # Already positive, but ensure minimum is at least 1.0 for numerical stability
    shift_amount <- max(0, 1.0 - y_min)
  }
  y_shifted <- y_vals + shift_amount
  
  # Check if all values are positive after shift
  if (!all(y_shifted > 0)) {
    return(NULL)
  }
  
  # Store shift amount as an attribute for later use in prediction
  tryCatch({
    # Fit exponential decay: log(y) = a + b*x
    # For RSQ: b should be negative (decreasing)
    # For CRPS/RMSE: b should be positive (increasing)
    df <- data.frame(y = y_shifted, x = x_vals)
    fit <- suppressWarnings(
      stats::glm(y ~ x, data = df, family = stats::gaussian(link = "log"))
    )
    # Store shift amount as attribute
    attr(fit, "shift_amount") <- shift_amount
    return(fit)
  }, error = function(e) {
    NULL
  })
}

predict_on_grid <- function(fit, newx) {
  if (is.null(fit)) return(rep(NA_real_, length(newx)))
  tryCatch({
    # Handle different model types
    if (inherits(fit, "glm")) {
      # For glm with log link, predict returns on log scale, need to exponentiate
      # The variable name in glm model is "x" (from data.frame(y = ..., x = ...))
      pred_log <- stats::predict(fit, data.frame(x = newx), type = "link")
      pred_vals <- as.numeric(exp(pred_log))
      # Subtract shift amount that was added during fitting
      shift_amount <- attr(fit, "shift_amount")
      if (is.null(shift_amount)) {
        # Fallback: assume old epsilon method
        shift_amount <- 1e-6
      }
      pred_vals - shift_amount
    } else if (inherits(fit, "gam")) {
      # For GAM, predict with newdata
      # The GAM model was fit with data.frame(y_vals, x_vals), so use x_vals
      pred_df <- data.frame(x_vals = newx)
      # Use predict.gam with newdata
      tryCatch({
        pred_vals <- as.numeric(mgcv::predict.gam(fit, newdata = pred_df, type = "response"))
        # Check for any invalid predictions
        if (any(!is.finite(pred_vals))) {
          # If any predictions are invalid, try to fix them
          pred_vals[!is.finite(pred_vals)] <- NA_real_
        }
        pred_vals
      }, error = function(e) {
        # If prediction fails, return NA
        rep(NA_real_, length(newx))
      })
    } else if (inherits(fit, "loess")) {
      # For loess, use x as variable name
      as.numeric(stats::predict(fit, data.frame(x = newx)))
    } else {
      # For lm, try to infer variable name from model
      x_var_name <- names(fit$model)[2]  # Usually the second column is x
      pred_df <- data.frame(x = newx)
      names(pred_df) <- x_var_name
      as.numeric(stats::predict(fit, pred_df))
    }
  }, error = function(e) {
    # Try alternative prediction with generic x
    tryCatch({
      if (inherits(fit, "gam")) {
        # For GAM, try with x_vals (the variable name used in safe_gam)
        pred_df <- data.frame(x_vals = newx)
        pred_vals <- as.numeric(mgcv::predict.gam(fit, newdata = pred_df, type = "response"))
        pred_vals[!is.finite(pred_vals)] <- NA_real_
        pred_vals
      } else {
      as.numeric(stats::predict(fit, newdata = data.frame(x = newx)))
      }
    }, error = function(e2) {
      rep(NA_real_, length(newx))
    })
  })
}

# Fit loess curves PER SITE, then aggregate horizons across sites
# Each plot contributes one point per timepoint, but curves are calculated per-site
# Use SITE-SPECIFIC null baselines for comparison
unique_sites <- unique(single_tax_plot_level$siteID[!is.na(single_tax_plot_level$siteID)])
site_horizons <- data.table(siteID = character(),
                            rsq_horizon = numeric(),
                            crps_horizon = numeric(),
                            rmse_horizon = numeric(),
                            method = character(),
                            rsq_horizon_gam = numeric(),
                            rsq_horizon_lm = numeric(),
                            crps_horizon_gam = numeric(),
                            crps_horizon_lm = numeric(),
                            rmse_horizon_gam = numeric(),
                            rmse_horizon_lm = numeric())

# Get site-specific null baselines from fcast_horizon_null_site
# Filter for this model and extract site-level nulls
if ("siteID" %in% names(fcast_horizon_null_site)) {
  model_null_site <- fcast_horizon_null_site[fcast_horizon_null_site$model_id == current_model_id & !is.na(siteID), ]
} else {
  # Fallback: if siteID not in null data, calculate from raw data
  model_null_site <- data.table()
}

for (site in unique_sites) {
  site_data <- single_tax_plot_level[siteID == site & 
                                      is.finite(months_since_obs) & 
                                      months_since_obs >= 0 & 
                                      months_since_obs <= 20]
  
  if (nrow(site_data) < 2) {
    next  # Skip sites with insufficient data
  }
  
  # GET SITE-SPECIFIC NULL BASELINES
  # NO FALLBACKS - site-specific nulls are required
  if (!("siteID" %in% names(fcast_horizon_null_site))) {
    cat("ERROR: siteID not found in fcast_horizon_null_site - cannot calculate site-specific horizons\n")
    next
  }
  
  if (nrow(model_null_site) == 0) {
    cat("ERROR: No site-specific null data found for model", current_model_id, "- skipping\n")
    next
  }
  
  if (!(site %in% model_null_site$siteID)) {
    cat("WARNING: Site", site, "not found in site-specific null data for model", current_model_id, "- skipping this site\n")
    next
  }
  
  site_null_row <- model_null_site[siteID == site]
  if (nrow(site_null_row) == 0) {
    cat("WARNING: No null data row found for site", site, "in model", current_model_id, "- skipping this site\n")
    next
  }
  
  # Extract null values for this site
  # Handle both renamed (without "null_" prefix) and original (with "null_" prefix) column names
  rsq_col <- if ("RSQ" %in% names(site_null_row)) "RSQ" 
             else if ("RSQ.1" %in% names(site_null_row)) "RSQ.1"
             else if ("null_RSQ" %in% names(site_null_row)) "null_RSQ"
             else if ("null_RSQ.1" %in% names(site_null_row)) "null_RSQ.1"
             else NULL
  rmse_col <- if ("RMSE.norm" %in% names(site_null_row)) "RMSE.norm"
              else if ("null_RMSE.norm" %in% names(site_null_row)) "null_RMSE.norm"
              else NULL
  
  this_site_null_rsq <- if (!is.null(rsq_col) && !is.na(site_null_row[[rsq_col]][1]) && is.finite(site_null_row[[rsq_col]][1])) {
    site_null_row[[rsq_col]][1]
  } else {
    cat("ERROR: Cannot extract RSQ null for site", site, "in model", current_model_id, "- skipping this site\n")
    next
  }
  
  this_site_null_crps <- if ("mean_crps" %in% names(site_null_row) && !is.na(site_null_row$mean_crps[1]) && is.finite(site_null_row$mean_crps[1])) {
    site_null_row$mean_crps[1]
  } else {
    cat("ERROR: Cannot extract mean_crps null for site", site, "in model", current_model_id, "- skipping this site\n")
    next
  }
  
  this_site_null_rmse <- if (!is.null(rmse_col) && !is.na(site_null_row[[rmse_col]][1]) && is.finite(site_null_row[[rmse_col]][1])) {
    site_null_row[[rmse_col]][1]
  } else {
    cat("ERROR: Cannot extract RMSE null for site", site, "in model", current_model_id, "- skipping this site\n")
    next
  }
  
  # IMPROVED HORIZON CALCULATION METHODS
  # Method 1: Null-Intersection (Binary Crossing) - improved loess
  # Method 2: GAM-Threshold (Smoothing) 
  # Method 3: CRPS-Crossing (Relative Improvement) 
  
  methods_to_try <- if (FCAST_HORIZON_METHOD == "all") {
    c("loess", "exponential", "gam")
  } else if (is.character(FCAST_HORIZON_METHOD) && length(FCAST_HORIZON_METHOD) == 1) {
    FCAST_HORIZON_METHOD
  } else {
    FCAST_HORIZON_METHOD
  }
  
  # Store results for each method
  site_method_results <- list()
  
  # Get model-level data for this site (from fcast_horizon_model_mean)
  # NO FALLBACKS - must have proper data structure
  if (!"siteID" %in% names(fcast_horizon_model_mean)) {
    cat("ERROR: fcast_horizon_model_mean missing siteID column - cannot calculate site-specific horizons for", site, "\n")
    next
  }
  
  site_model_data <- fcast_horizon_model_mean[model_id == current_model_id & siteID == site & 
                                              is.finite(months_since_obs) & months_since_obs >= 0 & months_since_obs <= 20]
  
  # Add month 0 data from last_obs_values_site_precalc if available
  if (exists("last_obs_values_site_precalc") && nrow(last_obs_values_site_precalc) > 0) {
    m0_data <- last_obs_values_site_precalc[model_id == current_model_id & siteID == site]
    if (nrow(m0_data) > 0) {
      # Check if month 0 already exists
      if (nrow(site_model_data[months_since_obs == 0 | abs(months_since_obs) < 0.01]) == 0) {
        m0_row <- m0_data[, .(months_since_obs = 0, RSQ, mean_crps, RMSE.norm)]
        # Ensure RSQ exists
        if (!"RSQ" %in% names(m0_row) && "RSQ.1" %in% names(m0_data)) {
          m0_row[, RSQ := m0_data$RSQ.1[1]]
        }
        site_model_data <- rbind(site_model_data, m0_row, fill = TRUE)
      }
    }
  }
  
  if (nrow(site_model_data) < 2) {
    cat("WARNING: Insufficient model-level data for", current_model_id, "site", site, "- skipping this site\n")
    next
  }
  
  # Ensure required columns exist - fail if missing
  if (!"RSQ" %in% names(site_model_data)) {
    if ("RSQ.1" %in% names(site_model_data)) {
      site_model_data[, RSQ := RSQ.1]
    } else {
      cat("ERROR: Missing RSQ column in site_model_data for", current_model_id, "site", site, "- skipping\n")
      next
    }
  }
  
  if (!"mean_crps" %in% names(site_model_data)) {
    cat("ERROR: Missing mean_crps column in site_model_data for", current_model_id, "site", site, "- skipping\n")
    next
  }
  
  if (!"RMSE.norm" %in% names(site_model_data)) {
    cat("ERROR: Missing RMSE.norm column in site_model_data for", current_model_id, "site", site, "- skipping\n")
    next
  }

  # Skip sites where ALL metric values are NA — no valid scoring data exists.
  # Without this check, loess/discrete fallbacks silently assign max_months,
  # producing bogus horizons (e.g. 9 or 12 months) for taxa that were never scored.
  n_valid_rsq  <- sum(is.finite(site_model_data$RSQ) & site_model_data$months_since_obs > 0, na.rm = TRUE)
  n_valid_crps <- sum(is.finite(site_model_data$mean_crps) & site_model_data$months_since_obs > 0, na.rm = TRUE)
  n_valid_rmse <- sum(is.finite(site_model_data$RMSE.norm) & site_model_data$months_since_obs > 0, na.rm = TRUE)
  if (n_valid_rsq == 0 && n_valid_crps == 0 && n_valid_rmse == 0) {
    next
  }

  for (method in methods_to_try) {
    if (method == "loess") {
      # Method 1: Fit loess/lm and find crossing on a fine grid (continuous horizon)
      # Fallback to discrete crossing if fit fails
      if (!is.finite(this_site_null_rsq) || !is.finite(this_site_null_crps) || !is.finite(this_site_null_rmse)) {
        cat("ERROR: Invalid null baselines for", current_model_id, "site", site, "- skipping method loess\n")
        next
      }
      max_months <- if (nrow(site_model_data) > 0) max(site_model_data$months_since_obs, na.rm = TRUE) else Inf
      min_months <- if (nrow(site_model_data) > 0) min(site_model_data$months_since_obs, na.rm = TRUE) else 0
      initial_data <- site_model_data[months_since_obs == min_months | abs(months_since_obs) < 0.01]
      if (nrow(initial_data) == 0 && nrow(site_model_data) > 0) {
        initial_data <- site_model_data[which.min(site_model_data$months_since_obs)]
      }
      # Grid within data range for interpolation (no extrapolation)
      grid_max <- min(20, max_months)
      xgrid_site <- seq(max(0, min_months), grid_max, by = 0.1)
      # RSQ: fit curve and find first x where predicted RSQ drops below null
      fit_rsq <- safe_loess(site_model_data$RSQ, site_model_data$months_since_obs)
      if (!is.null(fit_rsq)) {
        yhat_rsq <- predict_on_grid(fit_rsq, xgrid_site)
        rsq_horizon_site <- first_crossing(xgrid_site, yhat_rsq, this_site_null_rsq, "below", min_x = 0)
        if (!is.finite(rsq_horizon_site)) rsq_horizon_site <- max_months
      } else {
        rsq_crossings <- site_model_data[!is.na(RSQ) & is.finite(RSQ) & RSQ < this_site_null_rsq & months_since_obs > 0, months_since_obs]
        rsq_horizon_site <- if (length(rsq_crossings) > 0) min(rsq_crossings, na.rm = TRUE) else max_months
      }
      # CRPS: fit curve and find first x where predicted CRPS rises above null (only x > 0)
      fit_crps <- safe_loess(site_model_data$mean_crps, site_model_data$months_since_obs)
      if (!is.null(fit_crps)) {
        yhat_crps <- predict_on_grid(fit_crps, xgrid_site)
        crps_horizon_site <- first_crossing(xgrid_site, yhat_crps, this_site_null_crps, "above", min_x = 0)
        if (!is.finite(crps_horizon_site)) crps_horizon_site <- max_months
      } else {
        crps_crossings <- site_model_data[!is.na(mean_crps) & is.finite(mean_crps) & mean_crps > this_site_null_crps & months_since_obs > 0, months_since_obs]
        crps_horizon_site <- if (length(crps_crossings) > 0) min(crps_crossings, na.rm = TRUE) else max_months
      }
      # RMSE: fit curve and find first x where predicted RMSE rises above null (only x > 0)
      fit_rmse <- safe_loess(site_model_data$RMSE.norm, site_model_data$months_since_obs)
      if (!is.null(fit_rmse)) {
        yhat_rmse <- predict_on_grid(fit_rmse, xgrid_site)
        rmse_horizon_site <- first_crossing(xgrid_site, yhat_rmse, this_site_null_rmse, "above", min_x = 0)
        if (!is.finite(rmse_horizon_site)) rmse_horizon_site <- max_months
      } else {
        rmse_crossings <- site_model_data[!is.na(RMSE.norm) & is.finite(RMSE.norm) & RMSE.norm > this_site_null_rmse & months_since_obs > 0, months_since_obs]
        rmse_horizon_site <- if (length(rmse_crossings) > 0) min(rmse_crossings, na.rm = TRUE) else max_months
      }
      # Cap at max_months when curve never crosses
      if (rsq_horizon_site >= grid_max) rsq_horizon_site <- max_months
      if (crps_horizon_site >= grid_max) crps_horizon_site <- max_months
      if (rmse_horizon_site >= grid_max) rmse_horizon_site <- max_months
      site_method_results[["loess"]] <- list(rsq_horizon = rsq_horizon_site, crps_horizon = crps_horizon_site, rmse_horizon = rmse_horizon_site)
      
    } else if (method == "gam") {
      # Method 2: GAM-Threshold (Smoothing to find stable horizon)
      # NO FALLBACKS - GAM must succeed or return NA
      if (!is.finite(this_site_null_rsq) || !is.finite(this_site_null_crps) || !is.finite(this_site_null_rmse)) {
        cat("ERROR: Invalid null baselines for", current_model_id, "site", site, "- skipping method gam\n")
        next
      }
      
      calculate_gam_horizon <- function(dt, metric_col, null_val, dir = "above") {
        if (nrow(dt) < 3) return(NA_real_)
        valid_data <- dt[!is.na(get(metric_col)) & is.finite(get(metric_col))]
        if (nrow(valid_data) < 3) return(NA_real_)
        if (!is.finite(null_val)) return(NA_real_)
        min_months <- min(valid_data$months_since_obs, na.rm = TRUE)
        max_months <- max(valid_data$months_since_obs, na.rm = TRUE)
        min_val <- valid_data[months_since_obs == min_months, get(metric_col)][1]
        if (dir == "above" && !is.na(min_val) && is.finite(min_val) && min_val > null_val) {
          improves <- valid_data[get(metric_col) < null_val, months_since_obs]
          if (length(improves) > 0) {
            improve_point <- min(improves, na.rm = TRUE)
            worse_after_improve <- valid_data[months_since_obs > improve_point & get(metric_col) > null_val, months_since_obs]
            if (length(worse_after_improve) > 0) return(min(worse_after_improve, na.rm = TRUE))
          }
          if (min_months <= 0) return(NA_real_)
          return(min_months)
        } else if (dir == "below" && !is.na(min_val) && is.finite(min_val) && min_val < null_val) {
          improves <- valid_data[get(metric_col) > null_val, months_since_obs]
          if (length(improves) > 0) {
            improve_point <- min(improves, na.rm = TRUE)
            worse_after_improve <- valid_data[months_since_obs > improve_point & get(metric_col) < null_val, months_since_obs]
            if (length(worse_after_improve) > 0) return(min(worse_after_improve, na.rm = TRUE))
          }
          if (min_months <= 0) return(NA_real_)
          return(min_months)
        }
        tryCatch({
          k_edf <- min(5, nrow(valid_data) - 1)
          if (k_edf < 2) k_edf <- 2
          m <- mgcv::gam(get(metric_col) ~ s(months_since_obs, k = k_edf, bs = "tp"),
                        data = valid_data, method = "REML")
          if (!is.null(m$converged) && !m$converged) return(NA_real_)
          grid <- data.frame(months_since_obs = seq(min_months, max_months, by = 0.1))
          grid$pred <- as.numeric(predict(m, grid))
          if (dir == "above") {
            cross <- grid$months_since_obs[grid$pred >= null_val & is.finite(grid$pred) & grid$months_since_obs > 0]
          } else {
            cross <- grid$months_since_obs[grid$pred <= null_val & is.finite(grid$pred) & grid$months_since_obs > 0]
          }
          if (length(cross) == 0) max_months else min(cross, na.rm = TRUE)
        }, error = function(e) NA_real_)
      }
      
      calculate_lm_horizon <- function(dt, metric_col, null_val, dir = "above") {
        if (nrow(dt) < 3) return(NA_real_)
        valid_data <- dt[!is.na(get(metric_col)) & is.finite(get(metric_col))]
        if (nrow(valid_data) < 3) return(NA_real_)
        if (!is.finite(null_val)) return(NA_real_)
        min_months <- min(valid_data$months_since_obs, na.rm = TRUE)
        max_months <- max(valid_data$months_since_obs, na.rm = TRUE)
        tryCatch({
          form <- as.formula(paste(metric_col, "~ months_since_obs"))
          m <- lm(form, data = valid_data)
          grid <- data.frame(months_since_obs = seq(min_months, max_months, by = 0.1))
          grid$pred <- as.numeric(predict(m, grid))
          if (dir == "above") {
            cross <- grid$months_since_obs[grid$pred >= null_val & is.finite(grid$pred) & grid$months_since_obs > 0]
          } else {
            cross <- grid$months_since_obs[grid$pred <= null_val & is.finite(grid$pred) & grid$months_since_obs > 0]
          }
          if (length(cross) == 0) max_months else min(cross, na.rm = TRUE)
        }, error = function(e) NA_real_)
      }
      
      rsq_horizon_site_gam <- calculate_gam_horizon(site_model_data, "RSQ", this_site_null_rsq, "below")
      rsq_horizon_site_lm  <- calculate_lm_horizon(site_model_data, "RSQ", this_site_null_rsq, "below")
      crps_horizon_site_gam <- calculate_gam_horizon(site_model_data, "mean_crps", this_site_null_crps, "above")
      crps_horizon_site_lm  <- calculate_lm_horizon(site_model_data, "mean_crps", this_site_null_crps, "above")
      rmse_horizon_site_gam <- calculate_gam_horizon(site_model_data, "RMSE.norm", this_site_null_rmse, "above")
      rmse_horizon_site_lm  <- calculate_lm_horizon(site_model_data, "RMSE.norm", this_site_null_rmse, "above")
      rsq_horizon_site <- rsq_horizon_site_gam
      crps_horizon_site <- crps_horizon_site_gam
      rmse_horizon_site <- rmse_horizon_site_gam
      
      # Check which metrics failed and why
      failed_metrics <- c()
      if (is.na(rsq_horizon_site)) failed_metrics <- c(failed_metrics, "RSQ")
      if (is.na(crps_horizon_site)) failed_metrics <- c(failed_metrics, "CRPS")
      if (is.na(rmse_horizon_site)) failed_metrics <- c(failed_metrics, "RMSE")
      
      # If all three horizon metrics returned NA, skip this site (no horizon to aggregate)
      if (length(failed_metrics) == 3) {
        valid_rsq <- sum(!is.na(site_model_data$RSQ) & is.finite(site_model_data$RSQ))
        valid_crps <- sum(!is.na(site_model_data$mean_crps) & is.finite(site_model_data$mean_crps))
        valid_rmse <- sum(!is.na(site_model_data$RMSE.norm) & is.finite(site_model_data$RMSE.norm))
        if (!gam_skip_warned_this_model) {
          cat("WARNING: GAM horizon skipped for", current_model_id, "site", site,
              "- valid points: RSQ=", valid_rsq, ", CRPS=", valid_crps, ", RMSE=", valid_rmse, " (need ≥3)\n")
          gam_skip_warned_this_model <<- TRUE
        }
        next
      } else if (length(failed_metrics) > 0) {
        # Some metrics failed but not all - use available ones, set failed ones to NA
        if (is.na(rsq_horizon_site)) rsq_horizon_site <- NA_real_
        if (is.na(crps_horizon_site)) crps_horizon_site <- NA_real_
        if (is.na(rmse_horizon_site)) rmse_horizon_site <- NA_real_
      }
      
    } else if (method == "exponential") {
      # Method 3: CRPS-Crossing (Relative Improvement - compare to month 0 performance)
      # Work with available metrics - don't require all three
      m0_data <- site_model_data[months_since_obs == 0 | abs(months_since_obs) < 0.01]
      if (nrow(m0_data) == 0) {
        cat("ERROR: No month 0 data for", current_model_id, "site", site, "- skipping method exponential\n")
        next
      }
      
      start_rsq <- m0_data$RSQ[1]
      start_crps <- m0_data$mean_crps[1]
      start_rmse <- m0_data$RMSE.norm[1]
      
      # Check which metrics are available and valid
      has_valid_rsq <- is.finite(start_rsq)
      has_valid_crps <- is.finite(start_crps)
      has_valid_rmse <- is.finite(start_rmse)
      
      # Need at least one valid metric
      if (!has_valid_rsq && !has_valid_crps && !has_valid_rmse) {
        cat("ERROR: No valid month 0 performance values for", current_model_id, "site", site, "- skipping method exponential\n")
        next
      }
      
      max_months <- if (nrow(site_model_data) > 0) max(site_model_data$months_since_obs, na.rm = TRUE) else Inf
      
      # Calculate RSQ horizon if available
      if (has_valid_rsq) {
        rsq_threshold <- start_rsq * 0.67
        rsq_crossings <- site_model_data[!is.na(RSQ) & is.finite(RSQ) & RSQ < rsq_threshold & months_since_obs > 0, months_since_obs]
        rsq_horizon_site <- if (length(rsq_crossings) > 0) min(rsq_crossings, na.rm = TRUE) else max_months
      } else {
        rsq_horizon_site <- NA_real_
      }
      
      # Calculate CRPS horizon if available
      if (has_valid_crps) {
        crps_threshold <- start_crps * 1.5
        crps_crossings <- site_model_data[!is.na(mean_crps) & is.finite(mean_crps) & mean_crps > crps_threshold & months_since_obs > 0, months_since_obs]
        crps_horizon_site <- if (length(crps_crossings) > 0) min(crps_crossings, na.rm = TRUE) else max_months
      } else {
        crps_horizon_site <- NA_real_
      }
      
      # Calculate RMSE horizon if available
      if (has_valid_rmse) {
        rmse_threshold <- start_rmse * 1.5
        rmse_crossings <- site_model_data[!is.na(RMSE.norm) & is.finite(RMSE.norm) & RMSE.norm > rmse_threshold & months_since_obs > 0, months_since_obs]
        rmse_horizon_site <- if (length(rmse_crossings) > 0) min(rmse_crossings, na.rm = TRUE) else max_months
      } else {
        rmse_horizon_site <- NA_real_
      }
    }
    
    # Store results for this method (gam includes separate _gam and _lm columns)
    if (method == "gam") {
      site_method_results[[method]] <- list(
        rsq_horizon = rsq_horizon_site,
        crps_horizon = crps_horizon_site,
        rmse_horizon = rmse_horizon_site,
        rsq_horizon_gam = rsq_horizon_site_gam,
        rsq_horizon_lm = rsq_horizon_site_lm,
        crps_horizon_gam = crps_horizon_site_gam,
        crps_horizon_lm = crps_horizon_site_lm,
        rmse_horizon_gam = rmse_horizon_site_gam,
        rmse_horizon_lm = rmse_horizon_site_lm
      )
    } else {
      site_method_results[[method]] <- list(
        rsq_horizon = rsq_horizon_site,
        crps_horizon = crps_horizon_site,
        rmse_horizon = rmse_horizon_site
      )
    }
  }
  
  # Select method to use (if testing all, use loess as default for now, but store all)
  selected_method <- if (FCAST_HORIZON_METHOD == "all") "loess" else FCAST_HORIZON_METHOD
  if (!selected_method %in% names(site_method_results)) {
    selected_method <- names(site_method_results)[1]  # Fallback to first available
  }
  
  result <- site_method_results[[selected_method]]
  has_gam_lm <- selected_method == "gam" && "rsq_horizon_gam" %in% names(result)
  # Floor horizons in (0, 1) to 1 month so we don't report sub-month "0.1" as meaningful skill
  floor_min_horizon <- function(h) {
    if (!is.finite(h)) return(h)
    if (h > 0 && h < 1) return(1)
    h
  }
  site_row <- data.table(siteID = site,
                         rsq_horizon = floor_min_horizon(result$rsq_horizon),
                         crps_horizon = floor_min_horizon(result$crps_horizon),
                         rmse_horizon = floor_min_horizon(result$rmse_horizon),
                         method = selected_method,
                         rsq_horizon_gam = if (has_gam_lm) floor_min_horizon(result$rsq_horizon_gam) else NA_real_,
                         rsq_horizon_lm = if (has_gam_lm) floor_min_horizon(result$rsq_horizon_lm) else NA_real_,
                         crps_horizon_gam = if (has_gam_lm) floor_min_horizon(result$crps_horizon_gam) else NA_real_,
                         crps_horizon_lm = if (has_gam_lm) floor_min_horizon(result$crps_horizon_lm) else NA_real_,
                         rmse_horizon_gam = if (has_gam_lm) floor_min_horizon(result$rmse_horizon_gam) else NA_real_,
                         rmse_horizon_lm = if (has_gam_lm) floor_min_horizon(result$rmse_horizon_lm) else NA_real_)
  site_horizons <- rbindlist(list(site_horizons, site_row), fill = TRUE)
  
  # If testing all methods, also store comparison data for first model only
  if (FCAST_HORIZON_METHOD == "all" && current_model_id == model_id_list[1] && site == unique_sites[1]) {
    cat("\n=== METHOD COMPARISON FOR FIRST MODEL/SITE ===\n")
    cat("Site:", site, "Model:", current_model_id, "\n")
    for (m in names(site_method_results)) {
      r <- site_method_results[[m]]
      cat("  Method:", m, "\n")
      cat("    RSQ horizon:", r$rsq_horizon, "months\n")
      cat("    CRPS horizon:", r$crps_horizon, "months\n")
      cat("    RMSE horizon:", r$rmse_horizon, "months\n")
    }
  }
}

# Aggregate horizons across sites (use median to be robust to outliers)
if (nrow(site_horizons) > 0) {
  # Filter finite values for each metric
  rsq_finite <- site_horizons$rsq_horizon[is.finite(site_horizons$rsq_horizon)]
  crps_finite <- site_horizons$crps_horizon[is.finite(site_horizons$crps_horizon)]
  rmse_finite <- site_horizons$rmse_horizon[is.finite(site_horizons$rmse_horizon)]
  
  # Calculate median across sites. NA when no sites produced a valid horizon
  # (no data or no null crossing). Do NOT substitute a fallback value — that
  # would produce misleading horizons for taxa with missing scoring data.
  rsq_fcast_horizon  <- if (length(rsq_finite)  > 0) median(rsq_finite,  na.rm = TRUE) else NA_real_
  crps_fcast_horizon <- if (length(crps_finite) > 0) median(crps_finite, na.rm = TRUE) else NA_real_
  rmse_fcast_horizon <- if (length(rmse_finite) > 0) median(rmse_finite, na.rm = TRUE) else NA_real_
  # Aggregate GAM and LM horizons separately when present
  rsq_fcast_horizon_gam <- if ("rsq_horizon_gam" %in% names(site_horizons)) {
    x <- site_horizons$rsq_horizon_gam[is.finite(site_horizons$rsq_horizon_gam)]
    if (length(x) > 0) median(x, na.rm = TRUE) else NA_real_
  } else NA_real_
  rsq_fcast_horizon_lm <- if ("rsq_horizon_lm" %in% names(site_horizons)) {
    x <- site_horizons$rsq_horizon_lm[is.finite(site_horizons$rsq_horizon_lm)]
    if (length(x) > 0) median(x, na.rm = TRUE) else NA_real_
  } else NA_real_
  crps_fcast_horizon_gam <- if ("crps_horizon_gam" %in% names(site_horizons)) {
    x <- site_horizons$crps_horizon_gam[is.finite(site_horizons$crps_horizon_gam)]
    if (length(x) > 0) median(x, na.rm = TRUE) else NA_real_
  } else NA_real_
  crps_fcast_horizon_lm <- if ("crps_horizon_lm" %in% names(site_horizons)) {
    x <- site_horizons$crps_horizon_lm[is.finite(site_horizons$crps_horizon_lm)]
    if (length(x) > 0) median(x, na.rm = TRUE) else NA_real_
  } else NA_real_
  rmse_fcast_horizon_gam <- if ("rmse_horizon_gam" %in% names(site_horizons)) {
    x <- site_horizons$rmse_horizon_gam[is.finite(site_horizons$rmse_horizon_gam)]
    if (length(x) > 0) median(x, na.rm = TRUE) else NA_real_
  } else NA_real_
  rmse_fcast_horizon_lm <- if ("rmse_horizon_lm" %in% names(site_horizons)) {
    x <- site_horizons$rmse_horizon_lm[is.finite(site_horizons$rmse_horizon_lm)]
    if (length(x) > 0) median(x, na.rm = TRUE) else NA_real_
  } else NA_real_
} else {
  cat("Warning: No valid site horizons calculated for", current_model_id, "- skipping\n")
  next
}

# Horizons are already calculated from site-level data above (lines 1404-1417)
# The site-level horizons are the correct values to use
# DO NOT recalculate using yhat vectors that were never computed (they're all NA)
# The previous code at lines 1543-1545 was overwriting correct calculations with Inf

# Debug: Print values for first model to understand what's happening
debug_this_model <- (current_model_id == model_id_list[1])
if (debug_this_model) {
  cat("\n=== DETAILED HORIZON CALCULATION DEBUG FOR FIRST MODEL ===\n")
  cat("Model ID:", current_model_id, "\n\n")
  
  cat("1. RAW DATA (single_tax_for_fit):\n")
  cat("   Total rows:", nrow(single_tax_for_fit), "\n")
  cat("   months_since_obs range:", range(single_tax_for_fit$months_since_obs, na.rm=TRUE), "\n")
  cat("   Sample of raw data:\n")
  print(head(single_tax_for_fit[order(months_since_obs), .(months_since_obs, RSQ, mean_crps, RMSE.norm)], 10))
  
  cat("\n2. DATA FOR LOESS FITTING (single_tax_plot_level - AGGREGATED BY SITE):\n")
  cat("   Total rows (siteID x months_since_obs):", nrow(single_tax_plot_level), "\n")
  cat("   Unique months_since_obs values:", length(unique(single_tax_plot_level$months_since_obs)), "\n")
  cat("   Sample of data used for fitting:\n")
  print(head(single_tax_plot_level[order(months_since_obs), .(months_since_obs, RSQ, mean_crps, RMSE.norm)], 15))
  
  cat("\n3. NULL BASELINES:\n")
  cat("   RSQ null:", rsq_null_line, "\n")
  cat("   CRPS null:", crps_null_line, "\n")
  cat("   RMSE null:", rmse_null_line, "\n")
  cat("   Source (single_tax_null_collapse):\n")
  print(single_tax_null_collapse)
  
  cat("\n4. SITE-LEVEL HORIZONS:\n")
  cat("   Number of sites:", nrow(site_horizons), "\n")
  if (nrow(site_horizons) > 0) {
    cat("   RSQ horizons by site:", paste(round(site_horizons$rsq_horizon, 2), collapse=", "), "\n")
    cat("   CRPS horizons by site:", paste(round(site_horizons$crps_horizon, 2), collapse=", "), "\n")
    cat("   RMSE horizons by site:", paste(round(site_horizons$rmse_horizon, 2), collapse=", "), "\n")
  }
  
  cat("\n5. AGGREGATED HORIZONS (median across sites - FINAL VALUES):\n")
  cat("   RSQ horizon:", rsq_fcast_horizon, "months\n")
  cat("   CRPS horizon:", crps_fcast_horizon, "months\n")
  cat("   RMSE horizon:", rmse_fcast_horizon, "months\n")
  
  cat("\n==================================================\n\n")
}

# Horizons are already calculated from site-level data (lines 1404-1417) - DO NOT overwrite them!
# The old code below would overwrite correct calculations with Inf:
# rsq_fcast_horizon  <- first_crossing(xgrid, yhat_rsq,  rsq_null_line,  dir = "below")
# crps_fcast_horizon <- first_crossing(xgrid, yhat_crps, crps_null_line, dir = "above")
# rmse_fcast_horizon <- first_crossing(xgrid, yhat_rmse, rmse_null_line, dir = "above")
# This is removed because yhat_* vectors are all NA and would produce Inf

# Cap finite horizons at 20; preserve NA for taxa with no valid scoring data
max_horizon <- 20
cap <- function(x) if (is.finite(x)) min(x, max_horizon) else x
rsq_fcast_horizon  <- cap(rsq_fcast_horizon)
crps_fcast_horizon <- cap(crps_fcast_horizon)
rmse_fcast_horizon <- cap(rmse_fcast_horizon)
rsq_fcast_horizon_gam  <- cap(rsq_fcast_horizon_gam)
rsq_fcast_horizon_lm   <- cap(rsq_fcast_horizon_lm)
crps_fcast_horizon_gam <- cap(crps_fcast_horizon_gam)
crps_fcast_horizon_lm  <- cap(crps_fcast_horizon_lm)
rmse_fcast_horizon_gam <- cap(rmse_fcast_horizon_gam)
rmse_fcast_horizon_lm  <- cap(rmse_fcast_horizon_lm)
# ------------------------------------

if (MAKE_PLOTS) {
fig_rsq <- ggplot(single_tax_for_fit,
									aes(x = months_since_obs, y = RSQ, color = pretty_group)) +
	xlim(c(0,20)) +
	geom_point() +
	theme_minimal(base_size = 16) +
	ylab("Mean forecast accuracy (RSQ 1:1)") +
	xlab("Forecast horizon (months since last observation)") +
	labs(color="Kingdom") +
	geom_hline(yintercept = single_tax_null_collapse$RSQ[1], na.rm = TRUE, alpha=.8, linetype=2) + 
	geom_smooth(method = "loess", span = 0.8, se = FALSE, fullrange=TRUE)
		#  stat_poly_line(formula = y ~x, se = F, method = 'glm',	
		#  							 method.args = list(family = gaussian(link = 'log'))) +
		# stat_poly_eq(formula = y ~x, use_label("R2")) + guides(color="none")


fig_crps <- ggplot(single_tax_for_fit,
									aes(x = months_since_obs, y = mean_crps, color = pretty_group)) +
	xlim(c(0,20)) +
	geom_point() +
	theme_minimal(base_size = 16) +
	ylab("Mean forecast error (CRPS)") +
	xlab("Forecast horizon (months since last observation)") +
	labs(color="Kingdom") +
	geom_hline(yintercept = single_tax_null_collapse$mean_crps[1], na.rm = TRUE, alpha=.8, linetype=2) + 
	geom_smooth(method = "loess", span = 0.8, se = FALSE, fullrange=TRUE)
# +
# 	 stat_poly_line(formula = y ~ poly(x, 2, raw = TRUE), se = F) +
# 	stat_poly_eq(formula = y ~ poly(x, 2, raw = TRUE), use_label("R2")) + guides(color="none")


fig_rmse <- ggplot(single_tax_for_fit,
									 aes(x = months_since_obs, y = RMSE.norm, color = pretty_group)) +
	xlim(c(0,20)) +
	geom_point() +
	theme_minimal(base_size = 16) +
	ylab("Mean forecast error (RMSE)") +
	xlab("Forecast horizon (months since last observation)") +
	labs(color="Kingdom") +
	geom_hline(yintercept = single_tax_null_collapse$RMSE.norm[1], na.rm = TRUE, alpha=.8, linetype=2) + 
	geom_smooth(method = "loess", span = 0.8, se = FALSE, fullrange=TRUE)
# +
# 	stat_poly_line(formula = y ~ poly(x, 2, raw = TRUE), se = F) +
# 	stat_poly_eq(formula = y ~ poly(x, 2, raw = TRUE), use_label("R2")) + guides(color="none")

  save_fig = ggarrange(fig_rsq + rremove("xlab"), 
  										 fig_crps + rremove("xlab"), 
  										 fig_rmse, labels = c("A","B","C"), top = current_model_id)
  out_figure_list[[current_model_id]] <- save_fig
} else {
  out_figure_list[[current_model_id]] <- NULL
}

# Horizon calculation now done above with headless computation
# Just to see how close it got, if infinite
# lowest_y_val = loess_line$y[which.min(loess_line$y)]
# rsq_fit_info = horizon_rsq_data[[4]][,c("adj.r.squared","r.squared","p.value")] %>% 
# 	setNames(., paste0("rsq_",names(.)))


# RMSE horizon calculation now done above with headless computation
# Just to see how close it got, if infinite
# lowest_y_val = loess_line$y[which.min(loess_line$y)]
# rmse_fit_info = horizon_rmse_data[[4]][,c("adj.r.squared","r.squared","p.value")]  %>% 
# 	setNames(., paste0("rmse_",names(.)))



# CRPS horizon calculation now done above with headless computation
# Just to see how close it got, if infinite
# lowest_y_val = loess_line$y[which.min(loess_line$y)]
# crps_fit_info = horizon_crps_data[[4]][,c("adj.r.squared","r.squared","p.value")] %>% 
# 	setNames(., paste0("crps_",names(.)))


# Horizon capping now done above with headless computation

} # end if (!goto_output) — per-site computation

# For pooled mode, skip plots
if (HORIZON_POOLING == "pooled") {
  out_figure_list[[current_model_id]] <- NULL
}

# Output columns use three null baselines for horizon crossing:
#   *_fcast_horizon      = site-mean null (primary; backward compatible)
#   *_fcast_horizon_gam  = persistence null (last observed value carried forward)
#   *_fcast_horizon_lm   = climatological null (monthly mean from calibration period)
# The *_horizon_persist and *_horizon_clim columns duplicate _gam/_lm for clarity.
# RSQ throughout uses Nash-Sutcliffe efficiency (1 - SS_res/SS_tot), not regression R².
# Regression R² is computed as RSQ.reg in scoring tables but is NOT used for horizon
# crossing because Nash-Sutcliffe is the standard forecast skill benchmark.
out_df = cbind.data.frame(single_tax[1,1:5],
													rsq_fcast_horizon = rsq_fcast_horizon,
													rsq_fcast_horizon_gam = rsq_fcast_horizon_gam,
													rsq_fcast_horizon_lm = rsq_fcast_horizon_lm,
													rsq_null_line = rsq_null_line,
													rmse_fcast_horizon = rmse_fcast_horizon,
													rmse_fcast_horizon_gam = rmse_fcast_horizon_gam,
													rmse_fcast_horizon_lm = rmse_fcast_horizon_lm,
													rmse_null_line = rmse_null_line,
													crps_fcast_horizon = crps_fcast_horizon,
													crps_fcast_horizon_gam = crps_fcast_horizon_gam,
													crps_fcast_horizon_lm = crps_fcast_horizon_lm,
													crps_null_line = crps_null_line,
													# Pooled horizons with alternative null baselines
													rsq_horizon_persist = if (exists("ph_persist") && !is.na(ph_persist$rsq)) min(ph_persist$rsq, 20) else NA_real_,
													crps_horizon_persist = if (exists("ph_persist") && !is.na(ph_persist$crps)) min(ph_persist$crps, 20) else NA_real_,
													rmse_horizon_persist = if (exists("ph_persist") && !is.na(ph_persist$rmse)) min(ph_persist$rmse, 20) else NA_real_,
													rsq_horizon_clim = if (exists("ph_clim") && !is.na(ph_clim$rsq)) min(ph_clim$rsq, 20) else NA_real_,
													crps_horizon_clim = if (exists("ph_clim") && !is.na(ph_clim$crps)) min(ph_clim$crps, 20) else NA_real_,
													rmse_horizon_clim = if (exists("ph_clim") && !is.na(ph_clim$rmse)) min(ph_clim$rmse, 20) else NA_real_)
out_df_list[[current_model_id]] = out_df

cat("Completed processing for:", current_model_id, "\n")
cat("  - RSQ horizon:", rsq_fcast_horizon, "months\n")
cat("  - RMSE horizon:", rmse_fcast_horizon, "months\n") 
cat("  - CRPS horizon:", crps_fcast_horizon, "months\n")
    out_parallel_list[[current_model_id]] <- list(out_df, out_figure_list[[current_model_id]])

    # Clean up memory between models
    suppressWarnings(rm(single_tax, single_tax_null, single_tax_null_collapse,
                        single_tax_for_fit, single_tax_plot_level, pooled_data,
                        ph_sm, ph_persist, ph_clim, null_sm, null_persist, null_clim,
                        goto_output))
    if (MAKE_PLOTS) {
      suppressWarnings(rm(fig_rsq, fig_crps, fig_rmse))
    }
    suppressWarnings(rm(rsq_fcast_horizon, crps_fcast_horizon, rmse_fcast_horizon))
    suppressWarnings(rm(fit_rsq, fit_crps, fit_rmse, yhat_rsq, yhat_crps, yhat_rmse,
                        rsq_fcast_horizon_gam, rsq_fcast_horizon_lm,
                        crps_fcast_horizon_gam, crps_fcast_horizon_lm,
                        rmse_fcast_horizon_gam, rmse_fcast_horizon_lm,
                        site_horizons))
    gc(verbose = FALSE, full = TRUE)  # Force garbage collection
  }
  
  # Clean up group-level data after processing all models in this group
  cat("  Completed group", group_name, "- cleaning up memory...\n")
  if (exists("last_obs_values_site_precalc")) rm(last_obs_values_site_precalc)
  gc(verbose = FALSE, full = TRUE)
  cat("  Group", group_name, "complete. Processed", length(group_model_ids), "models.\n")
}

# Clean up fcast_horizon_df now that we're done with the loop
rm(fcast_horizon_df)
gc(verbose = FALSE)

out_df_list = lapply(out_parallel_list, "[[", 1)
out_figure_list = lapply(out_parallel_list, "[[", 2)

cat("Debug: Number of parallel results:", length(out_parallel_list), "\n")
cat("Debug: Number of df results:", length(out_df_list), "\n")

# Filter out empty data frames
out_df_list = out_df_list[vapply(out_df_list, function(x) nrow(x) > 0, FUN.VALUE = logical(1))]
cat("Debug: Number of non-empty df results:", length(out_df_list), "\n")

if (length(out_df_list) == 0) {
  stop("No results from parallel processing - all iterations may have failed")
}

fcast_horizon_results = data.table::rbindlist(out_df_list) %>% as.data.frame()
cat("Debug: fcast_horizon_results columns:", colnames(fcast_horizon_results), "\n")
cat("Debug: fcast_horizon_results rows:", nrow(fcast_horizon_results), "\n")

# Summary of horizon results
cat("\n=== FORECAST HORIZON SUMMARY ===\n")
cat("Site mean null (primary):\n")
cat("  RSQ  - Mean:", round(mean(fcast_horizon_results$rsq_fcast_horizon, na.rm=T), 2),
    " Median:", median(fcast_horizon_results$rsq_fcast_horizon, na.rm=T), "\n")
cat("  RMSE - Mean:", round(mean(fcast_horizon_results$rmse_fcast_horizon, na.rm=T), 2),
    " Median:", median(fcast_horizon_results$rmse_fcast_horizon, na.rm=T), "\n")
cat("  CRPS - Mean:", round(mean(fcast_horizon_results$crps_fcast_horizon, na.rm=T), 2),
    " Median:", median(fcast_horizon_results$crps_fcast_horizon, na.rm=T), "\n")
if ("crps_horizon_persist" %in% names(fcast_horizon_results)) {
  cat("Persistence null:\n")
  cat("  RSQ  - Mean:", round(mean(fcast_horizon_results$rsq_horizon_persist, na.rm=T), 2),
      " Median:", median(fcast_horizon_results$rsq_horizon_persist, na.rm=T), "\n")
  cat("  RMSE - Mean:", round(mean(fcast_horizon_results$rmse_horizon_persist, na.rm=T), 2),
      " Median:", median(fcast_horizon_results$rmse_horizon_persist, na.rm=T), "\n")
  cat("  CRPS - Mean:", round(mean(fcast_horizon_results$crps_horizon_persist, na.rm=T), 2),
      " Median:", median(fcast_horizon_results$crps_horizon_persist, na.rm=T), "\n")
}
if ("crps_horizon_clim" %in% names(fcast_horizon_results)) {
  cat("Climatological null:\n")
  cat("  RSQ  - Mean:", round(mean(fcast_horizon_results$rsq_horizon_clim, na.rm=T), 2),
      " Median:", median(fcast_horizon_results$rsq_horizon_clim, na.rm=T), "\n")
  cat("  RMSE - Mean:", round(mean(fcast_horizon_results$rmse_horizon_clim, na.rm=T), 2),
      " Median:", median(fcast_horizon_results$rmse_horizon_clim, na.rm=T), "\n")
  cat("  CRPS - Mean:", round(mean(fcast_horizon_results$crps_horizon_clim, na.rm=T), 2),
      " Median:", median(fcast_horizon_results$crps_horizon_clim, na.rm=T), "\n")
}
cat("================================\n\n")



# Use data.table melt instead of pivot_longer for memory efficiency
setDT(fcast_horizon_results)
measure_vars <- c("rsq_fcast_horizon", "rmse_fcast_horizon", "crps_fcast_horizon",
                  "rsq_null_line", "rmse_null_line", "crps_null_line",
                  "rsq_horizon_persist", "crps_horizon_persist", "rmse_horizon_persist",
                  "rsq_horizon_clim", "crps_horizon_clim", "rmse_horizon_clim")
measure_vars <- intersect(measure_vars, names(fcast_horizon_results))
if (length(measure_vars) == 0) {
  stop("fcast_horizon_results has none of the expected horizon columns: ",
       paste(c("rsq_fcast_horizon", "rmse_fcast_horizon", "crps_fcast_horizon",
               "rsq_null_line", "rmse_null_line", "crps_null_line"), collapse = ", "))
}
id_cols <- setdiff(names(fcast_horizon_results), measure_vars)
fcast_horizon_long <- melt(fcast_horizon_results,
                          id.vars = id_cols,
                          measure.vars = measure_vars,
                          variable.name = "horizon_parameter",
                          value.name = "value",
                          variable.factor = FALSE)

# Extract metric and parameter_type using data.table syntax
fcast_horizon_long[, metric := sapply(strsplit(horizon_parameter, "_"), function(x) x[1])]
fcast_horizon_long[, parameter_type := sapply(strsplit(horizon_parameter, "_"), function(x) x[2])]
fcast_horizon_long[, parameter_type := ifelse(parameter_type == "fcast", "horizon", parameter_type)]


# Save as RDS
saveRDS(list(fcast_horizon_results, out_figure_list, fcast_horizon_long),
				here("data/summary/fcast_horizon_df.rds"))
cat("✓ Saved RDS file: data/summary/fcast_horizon_df.rds\n")

# Also save as Parquet for memory efficiency
parquet_dir <- here("data/summary/parquet")
if (!dir.exists(parquet_dir)) {
  dir.create(parquet_dir, showWarnings = FALSE, recursive = TRUE)
}

parquet_file <- here("data/summary/parquet/fcast_horizon_df.parquet")
parquet_saved <- FALSE

# Try nanoparquet first (lightweight)
if (requireNamespace("nanoparquet", quietly = TRUE)) {
  tryCatch({
    # Save each element separately as parquet files
    nanoparquet::write_parquet(fcast_horizon_results, here("data/summary/parquet/fcast_horizon_results.parquet"))
    nanoparquet::write_parquet(fcast_horizon_long, here("data/summary/parquet/fcast_horizon_long.parquet"))
    cat("✓ Saved Parquet files using nanoparquet:\n")
    cat("  - data/summary/parquet/fcast_horizon_results.parquet\n")
    cat("  - data/summary/parquet/fcast_horizon_long.parquet\n")
    parquet_saved <- TRUE
  }, error = function(e) {
    cat("⚠️  Failed to save Parquet with nanoparquet:", e$message, "\n")
  })
}

# Fallback to arrow if nanoparquet failed
if (!parquet_saved && requireNamespace("arrow", quietly = TRUE)) {
  tryCatch({
    arrow::write_parquet(fcast_horizon_results, here("data/summary/parquet/fcast_horizon_results.parquet"))
    arrow::write_parquet(fcast_horizon_long, here("data/summary/parquet/fcast_horizon_long.parquet"))
    cat("✓ Saved Parquet files using arrow:\n")
    cat("  - data/summary/parquet/fcast_horizon_results.parquet\n")
    cat("  - data/summary/parquet/fcast_horizon_long.parquet\n")
    parquet_saved <- TRUE
  }, error = function(e) {
    cat("⚠️  Failed to save Parquet with arrow:", e$message, "\n")
  })
}

if (!parquet_saved) {
  cat("⚠️  Neither nanoparquet nor arrow available - Parquet files not saved\n")
  cat("  Install 'nanoparquet' or 'arrow' package to enable Parquet format.\n")
}

# ============================================================================
# KINGDOM ANALYSIS: Compare forecast horizons between Bacteria and Fungi
# ============================================================================

cat("\n=== ANALYZING FORECAST HORIZONS BY KINGDOM ===\n")

# Function to classify species as bacteria or fungi based on taxonomic patterns
classify_kingdom <- function(species_name) {
  # Common fungal phyla and classes
  fungal_patterns <- c("ascomycota", "basidiomycota", "chytridiomycota", "glomeromycota", 
                       "mortierellomycota", "mucoromycota", "rozellomycota",
                       "ascomycetes", "basidiomycetes", "chytridiomycetes", 
                       "glomeromycetes", "mortierellomycetes", "mucoromycetes",
                       "eurotiomycetes", "sordariomycetes", "dothideomycetes",
                       "leotiomycetes", "pezizomycetes", "tremellomycetes",
                       "agaricomycetes", "microbotryomycetes", "spizellomycetes",
                       "geminibasidiomycetes", "umbelopsidomycetes",
                       "eurotiales", "sordariales", "pleosporales", "helotiales",
                       "pezizales", "tremellales", "agaricales", "capnodiales",
                       "chaetothyriales", "chaetosphaeriales", "hypocreales",
                       "thelebolales", "xylariales", "cantharellales", "russulales",
                       "thelephorales", "mortierellales", "umbelopsidales",
                       "aspergillus", "penicillium", "fusarium", "mortierella",
                       "trichoderma", "cladosporium", "cortinarius", "russula",
                       "entoloma", "talaromyces", "exophiala", "coniochaeta",
                       "cladophialophora", "metarhizium", "oidiodendron",
                       "cenococcum", "glomus", "ectomycorrhizal", "endophyte",
                       "saprotroph", "plant_pathogen", "animal_pathogen",
                       "lichenized", "lignolytic", "chitinolytic", "cellulolytic",
                       "n_fixation", "nitrification", "denitrification",
                       "acetogen", "anaerobic", "copiotroph", "oligotroph")
  
  # Check if species name contains any fungal patterns
  species_lower <- tolower(species_name)
  
  for (pattern in fungal_patterns) {
    if (grepl(pattern, species_lower, fixed = TRUE)) {
      return("Fungi")
    }
  }
  
  # If no fungal patterns found, classify as bacteria
  return("Bacteria")
}

# Apply kingdom classification to horizon results
fcast_horizon_results$kingdom <- sapply(fcast_horizon_results$species, classify_kingdom)
fcast_horizon_long$kingdom <- sapply(fcast_horizon_long$species, classify_kingdom)

# Check classification results
cat("Kingdom distribution in horizon results:\n")
print(table(fcast_horizon_results$kingdom, useNA = "ifany"))

# Filter for RSQ and RMSE horizons
rsq_horizons <- fcast_horizon_long[fcast_horizon_long$metric == "rsq" & fcast_horizon_long$horizon_parameter == "rsq_fcast_horizon", ]
rmse_horizons <- fcast_horizon_long[fcast_horizon_long$metric == "rmse" & fcast_horizon_long$horizon_parameter == "rmse_fcast_horizon", ]

# RSQ Horizon Analysis
cat("\nRSQ Forecast Horizons by Kingdom:\n")
# Filter out invalid values (NaN, Inf, -Inf)
rsq_horizons_valid <- rsq_horizons[is.finite(rsq_horizons$value) & !is.na(rsq_horizons$value), ]
if (nrow(rsq_horizons_valid) > 0) {
  rsq_summary <- aggregate(value ~ kingdom, data = rsq_horizons_valid, 
                          FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                            median = median(x, na.rm = TRUE),
                                            sd = sd(x, na.rm = TRUE),
                                            n = length(x)))
  print(rsq_summary)
  
  # Statistical test
  cat("\nRSQ horizon comparison (Bacteria vs Fungi):\n")
  if (nrow(rsq_horizons_valid) > 0 && length(unique(rsq_horizons_valid$kingdom)) == 2) {
    rsq_test <- t.test(value ~ kingdom, data = rsq_horizons_valid)
    print(rsq_test)
  } else {
    cat("Only", length(unique(rsq_horizons_valid$kingdom)), "kingdom(s) available - skipping comparison\n")
  }
  
  # Summary statistics
  bacteria_rsq <- rsq_horizons_valid$value[rsq_horizons_valid$kingdom == "Bacteria"]
  fungi_rsq <- rsq_horizons_valid$value[rsq_horizons_valid$kingdom == "Fungi"]
  
  cat("  Bacteria: Mean =", round(mean(bacteria_rsq, na.rm = TRUE), 2), 
      "months, Median =", round(median(bacteria_rsq, na.rm = TRUE), 2), 
      "months, N =", length(bacteria_rsq), "\n")
  cat("  Fungi: Mean =", round(mean(fungi_rsq, na.rm = TRUE), 2), 
      "months, Median =", round(median(fungi_rsq, na.rm = TRUE), 2), 
      "months, N =", length(fungi_rsq), "\n")
} else {
  cat("No RSQ horizon data found\n")
}

# RMSE Horizon Analysis
cat("\nRMSE Forecast Horizons by Kingdom:\n")
# Filter out invalid values (NaN, Inf, -Inf)
rmse_horizons_valid <- rmse_horizons[is.finite(rmse_horizons$value) & !is.na(rmse_horizons$value), ]
if (nrow(rmse_horizons_valid) > 0) {
  rmse_summary <- aggregate(value ~ kingdom, data = rmse_horizons_valid, 
                           FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                             median = median(x, na.rm = TRUE),
                                             sd = sd(x, na.rm = TRUE),
                                             n = length(x)))
  print(rmse_summary)
  
  # Statistical test
  cat("\nRMSE horizon comparison (Bacteria vs Fungi):\n")
  if (nrow(rmse_horizons_valid) > 0 && length(unique(rmse_horizons_valid$kingdom)) == 2) {
    rmse_test <- t.test(value ~ kingdom, data = rmse_horizons_valid)
    print(rmse_test)
  } else {
    cat("Only", length(unique(rmse_horizons_valid$kingdom)), "kingdom(s) available - skipping comparison\n")
  }
  
  # Summary statistics
  bacteria_rmse <- rmse_horizons_valid$value[rmse_horizons_valid$kingdom == "Bacteria"]
  fungi_rmse <- rmse_horizons_valid$value[rmse_horizons_valid$kingdom == "Fungi"]
  
  cat("  Bacteria: Mean =", round(mean(bacteria_rmse, na.rm = TRUE), 2), 
      "months, Median =", round(median(bacteria_rmse, na.rm = TRUE), 2), 
      "months, N =", length(bacteria_rmse), "\n")
  cat("  Fungi: Mean =", round(mean(fungi_rmse, na.rm = TRUE), 2), 
      "months, Median =", round(median(fungi_rmse, na.rm = TRUE), 2), 
      "months, N =", length(fungi_rmse), "\n")
} else {
  cat("No RMSE horizon data found\n")
}

# Overall horizon distribution
cat("\n=== OVERALL HORIZON DISTRIBUTION ===\n")
if (nrow(rsq_horizons) > 0) {
  cat("RSQ Horizon Distribution:\n")
  cat("  Overall: Min =", min(rsq_horizons$value, na.rm = TRUE), 
      ", Max =", max(rsq_horizons$value, na.rm = TRUE), 
      ", Mean =", round(mean(rsq_horizons$value, na.rm = TRUE), 2), "\n")
  
  # Count how many are at the maximum (20 months)
  max_horizon_count <- sum(rsq_horizons$value == 20, na.rm = TRUE)
  cat("  Models at maximum horizon (20 months):", max_horizon_count, 
      "out of", nrow(rsq_horizons), "(", round(100 * max_horizon_count / nrow(rsq_horizons), 1), "%)\n")
}

if (nrow(rmse_horizons) > 0) {
  cat("\nRMSE Horizon Distribution:\n")
  cat("  Overall: Min =", min(rmse_horizons$value, na.rm = TRUE), 
      ", Max =", max(rmse_horizons$value, na.rm = TRUE), 
      ", Mean =", round(mean(rmse_horizons$value, na.rm = TRUE), 2), "\n")
  
  # Count how many are at the maximum (20 months)
  max_horizon_count <- sum(rmse_horizons$value == 20, na.rm = TRUE)
  cat("  Models at maximum horizon (20 months):", max_horizon_count, 
      "out of", nrow(rmse_horizons), "(", round(100 * max_horizon_count / nrow(rmse_horizons), 1), "%)\n")
}

cat("\n=== KINGDOM ANALYSIS COMPLETE ===\n")

# Confirm that all Inf results are from a forecast horizon beyond 12 months:
inf_horizon = fcast_horizon_results[fcast_horizon_results$rsq_fcast_horizon==Inf,]
if (length(unique(to_plot$model_id)) > 1 && nrow(inf_horizon) > 0 && exists("to_plot")) {
  ggplot(to_plot %>% filter(model_id %in% inf_horizon$model_id),
  			 aes(x = months_since_obs, y = RSQ, color = pretty_group, group=model_id)) +
  	
  	facet_wrap( ~model_id, 
  						 labeller = labeller(model_name = model.labs), scales="free") +
  	geom_point(show.legend = F) +
	theme_minimal(base_size = 16) +
	ylab("Mean forecast accuracy (RSQ 1:1)") +
	xlab("Forecast horizon (months since last observation)") +
	labs(color="Kingdom") +
	geom_hline(data =to_plot_null %>% filter(model_id %in% inf_horizon$model_id & is.infinite(months_since_obs)), 
						 aes(yintercept = RSQ,
																				color = pretty_group, 
																				group=model_name), na.rm = T, alpha=.8, show.legend = F) + 
  	geom_smooth(method = "glm", formula = y~x,	method.args = list(family = gaussian(link = 'log')))
  # +
  # 	stat_poly_line(formula = y ~ poly(x, 2, raw = TRUE),se=F, show.legend = F) +
  # 	stat_poly_eq(formula = y ~ poly(x, 2, raw = TRUE), use_label("R2"), show.legend = F)
} else {
  cat("Skipping infinite horizon faceted plot - only 1 model available\n")
}

# Change infinite horizon to max (12 months)
# fcast_horizon_results <- fcast_horizon_results %>% mutate(horizon = ifelse(horizon==Inf, 12, horizon))

# View overall horizon results!
if (length(unique(fcast_horizon_results$model_name)) > 1 && 
    "rank_only" %in% names(fcast_horizon_results) &&
    "rsq_fcast_horizon" %in% names(fcast_horizon_results) &&
    sum(is.finite(fcast_horizon_results$rsq_fcast_horizon)) > 0) {
  tryCatch({
    ggplot(fcast_horizon_results,
    			 aes(x = rank_only, y = rsq_fcast_horizon, color = pretty_group)) +
    	facet_grid(pretty_group ~model_name, 
    						 labeller = labeller(model_name = model.labs), scales="free") + 
    	coord_flip() +
    	geom_boxplot() +
    	geom_point(position=position_jitter(), alpha=.3)
  }, error = function(e) {
    cat("Skipping faceted plot - error:", e$message, "\n")
  })
} else {
  cat("Skipping faceted plot - insufficient data or only 1 model\n")
}

if (length(unique(fcast_horizon_long$model_name)) > 1) {
  if ("rank_only" %in% names(fcast_horizon_long) && 
      nrow(fcast_horizon_long %>% filter(parameter_type=="null")) > 0 &&
      length(unique(fcast_horizon_long$metric)) > 0 &&
      length(unique(fcast_horizon_long$model_name)) > 0) {
    tryCatch({
      ggplot(fcast_horizon_long %>% filter(parameter_type=="null"),
      			 aes(x = rank_only, y = value, color = pretty_group)) +
      	facet_grid(metric ~model_name, 
      						 labeller = labeller(model_name = model.labs), scales="free") + 
      	#coord_flip() +
      	geom_boxplot() +
      	geom_point(position=position_jitter(), alpha=.3)
    }, error = function(e) {
      cat("Skipping null parameter faceted plot - error:", e$message, "\n")
    })
  } else {
    cat("Skipping null parameter faceted plot - insufficient data\n")
  }
} else {
  cat("Skipping null parameter faceted plot - only 1 model available\n")
}



to_plot



if (length(unique(to_plot$model_name)) > 1 && exists("to_plot") && exists("converged") && exists("to_plot_null")) {
  tryCatch({
    plot_data_check <- to_plot %>% filter(rank_only=="functional") %>% filter(model_id %in% converged)
    if (nrow(plot_data_check) > 0 && length(unique(plot_data_check$pretty_group)) > 0 && length(unique(plot_data_check$model_name)) > 0) {
      ggplot(plot_data_check,
      			 aes(x = months_since_obs, y = RSQ, color = pretty_group)) +
      	facet_grid(pretty_group ~model_name, 
      						 labeller = labeller(model_name = model.labs), scales="free") + 
      	#coord_flip() +
      #	geom_boxplot() +
      	geom_point(position=position_jitter(), alpha=.3) +
      	stat_poly_line(formula = y ~ poly(x, 2, raw = TRUE), se = F) +
      	stat_poly_eq(formula = y ~ poly(x, 2, raw = TRUE), use_label("R2")) + guides(color="none") +
      	geom_hline(data = to_plot_null, aes(yintercept = RSQ.1), alpha=.1)
    }
  }, error = function(e) {
    cat("Skipping functional faceted plot 1 - error:", e$message, "\n")
  })
} else {
  cat("Skipping functional faceted plot 1 - insufficient data or only 1 model\n")
}



if (length(unique(to_plot$model_name)) > 1 && exists("to_plot") && exists("converged") && exists("to_plot_null")) {
  tryCatch({
    plot_data_check <- to_plot %>% filter(rank_only=="functional") %>% filter(model_id %in% converged)
    if (nrow(plot_data_check) > 0 && length(unique(plot_data_check$pretty_group)) > 0 && length(unique(plot_data_check$model_name)) > 0) {
      ggplot(plot_data_check,
      			 aes(x = months_since_obs, y = RSQ, color = pretty_group)) +
      	facet_grid(pretty_group ~model_name, 
      						 labeller = labeller(model_name = model.labs), scales="free") +
      	geom_hline(data = to_plot_null, aes(yintercept = RSQ.1), alpha=.1)
    }
  }, error = function(e) {
    cat("Skipping functional faceted plot 2 - error:", e$message, "\n")
  })
} else {
  cat("Skipping functional faceted plot 2 - insufficient data or only 1 model\n")
}

	to_plot_null

to_plot_null


# Only filter if rank_only column exists
if (exists("to_plot_null") && "rank_only" %in% names(to_plot_null)) {
  to_plot_null = to_plot_null %>% group_by(pretty_group,model_name) %>% mutate(median_null = median(RSQ)) %>% 
    filter(rank_only %in% c("phylum","functional")) %>% filter(model_name %in% c("env_cycl") | grepl("driver_uncertainty", model_name))
} else if (exists("to_plot_null")) {
  # If rank_only doesn't exist, just filter by model_name
  to_plot_null = to_plot_null %>% group_by(pretty_group,model_name) %>% mutate(median_null = median(RSQ)) %>% 
    filter(model_name %in% c("env_cycl") | grepl("driver_uncertainty", model_name))
}
if (length(unique(to_plot$model_name)) > 1 && exists("to_plot")) {
  # Filter and validate data before plotting
  plot_data <- to_plot %>% 
    {if ("rank_only" %in% names(.)) filter(., rank_only %in% c("phylum","functional")) else .} %>%
    filter(model_name %in% c("env_cycl") | grepl("driver_uncertainty", model_name)) %>%
    filter(!is.na(months_since_obs) & !is.na(RSQ) & is.finite(months_since_obs) & is.finite(RSQ))
  
  if (nrow(plot_data) > 0) {
    p <- ggplot(plot_data,
      			 aes(x = months_since_obs, y = RSQ, color = pretty_group)) +
      	
      	facet_grid(rows=vars(pretty_group), 
      						 labeller = labeller(model_name = model.labs), scales="free") +
      	geom_point(position=position_jitter(), alpha=.1, show.legend = F) +
      	# Only add smooth if we have enough data points
      	{ if (nrow(plot_data) >= 3) geom_smooth(span=2, method="loess", show.legend = F, se=F) else NULL } +
    	stat_cor(
    		aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")),
    		p.digits = 1, label.x.npc = .15, label.y.npc = .9, show.legend = F
    	) +
	theme_bw(base_size = 20) +
	ylab("Forecast accuracy (RSQ)") +
	xlab("Months since last observation") +
	geom_hline(data = to_plot_null, aes(yintercept = RSQ , group=model_name), alpha=.05) +
	geom_hline(data = to_plot_null, aes(yintercept = median_null, group=model_name),linetype=2, alpha=.5) 

arrows <- 
	tibble(
		x1 = c(8, 12),
		x2 = c(7.2, 10),
		y1 = c(.5, .3),
		y2 = c(.4, .12), 
		pretty_group = c("Bacteria", "Fungi")
	)

dat_text <- data.frame(
	pretty_group = c("Bacteria", "Fungi","Bacteria", "Fungi"),
	months_since_obs     = c(9 , 11.5, 2.2, 2.2),
	RSQ     = c(.55,.42, .44, .18),
	label = c("Forecast horizon: \n7 months","Forecast horizon: \n10 months",
						"Null predictability","Null predictability")
)
    p + geom_curve(
    	data = arrows, aes(x = x1, y = y1, xend = x2, yend = y2),
    	arrow = arrow(length = unit(0.08, "inch")), size = 0.8, curvature = -.3, color = "gray20") + 
    	geom_text(data = dat_text, aes(label=label), color=1, size=5, show.legend = F) + guides(pretty_group=NULL)
  } else {
    cat("Skipping final faceted plot - insufficient data for plotting\n")
  }
} else {
  cat("Skipping final faceted plot - only 1 model available\n")
}

