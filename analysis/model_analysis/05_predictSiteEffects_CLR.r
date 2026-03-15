# OPTIMIZED VERSION: Determine important predictors of site effects (random effects) for CLR models
# Creates the "site_effects_dredged_CLR.rds" files
# Has to be run before creating any CLR forecasts (06_createHindcasts_CLR)
# 
# PERFORMANCE OPTIMIZATIONS:
# 1. Sequential processing (can be switched to parallel)
# 2. Pre-load all unobserved data files once
# 3. Optimized PLSR operations
# 4. Limited dredge search space
# 5. Streamlined data operations

source("../../source.R")
pacman::p_load(stringr, forestplot, gridExtra, ggpubr, MuMIn, pls, ggforce, reshape2, plotrix, spectratrait, 
               parallel, doParallel, foreach)

# Load the PLSR function - try fixed version first, then fallback
fixed_function_path <- here("microbialForecast/R/spectra_site_eff_permutation_fixed.r")
if (file.exists(fixed_function_path)) {
  source(fixed_function_path)
  cat("Loaded fixed spectra function\n")
} else {
  cat("WARNING: Fixed spectra function not found, using existing one\n")
  if (!exists("site_eff_uncertainties")) {
    cat("ERROR: site_eff_uncertainties function not available\n")
    stop("Required PLSR function not found")
  }
}

# Created by running "./data_construction/covariate_prep/combine_site_descriptors.r"
df_predictors <- readRDS(here("data/clean/site_effect_predictors.rds"))
df_predictors$latitude_scaled <- scale(df_predictors$latitude)

# FIXED: Read in CLR site effect estimates from model; split into calibration timeperiod and full dataset
# For CLR models, we need to check if site_effects_CLR.rds exists, otherwise use regular site_effects.rds
clr_site_effects_file <- here("data/summary/site_effects_CLR.rds")
if (file.exists(clr_site_effects_file)) {
  cat("Using CLR-specific site effects file:", clr_site_effects_file, "\n")
  site_effects_all <- readRDS(clr_site_effects_file) %>% filter(nchar(siteID) > 3) #remove numeric siteIDs
} else {
  # Fall back to regular site effects if CLR-specific doesn't exist
  fallback_file <- here("data/summary/site_effects.rds")
  if (file.exists(fallback_file)) {
    cat("CLR-specific site effects not found, using fallback:", fallback_file, "\n")
    site_effects_all <- readRDS(fallback_file) %>% filter(nchar(siteID) > 3)
  } else {
    cat("ERROR: No site effects files found!\n")
    cat("This script requires site effects data from upstream analysis\n")
    stop("Site effects data not available")
  }
}

cat("Loaded site effects data with", nrow(site_effects_all), "rows\n")

if (nrow(site_effects_all) == 0) {
  stop("No site effects data loaded")
}

site_effects_calibration <- site_effects_all %>%
  filter(time_period == "2015-11_2018-01") %>%
  mutate(taxon = ifelse(taxon=="other", paste0(taxon, "_", rank), taxon))

# Observed site effects for explaining
df_observed <- site_effects_calibration %>%
  select(siteID, model_id, Median, rank, rank_only, group, taxon) %>%
  merge(df_predictors, all.x=T)

# OPTIMIZATION: Pre-load all unobserved data files once
unobserved_data <- list()
model_types <- c("cycl_only", "env_cov", "env_cycl")
for(model_type in model_types) {
  unobserved_file <- here(paste0("data/summary/site_effects_unobserved_", model_type, "_CLR.rds"))
  if(file.exists(unobserved_file)) {
    unobserved_data[[model_type]] <- readRDS(unobserved_file)
    cat("Pre-loaded unobserved data for", model_type, ":", nrow(unobserved_data[[model_type]]), "sites\n")
  } else {
    # Try fallback to regular unobserved data
    fallback_file <- here(paste0("data/summary/site_effects_unobserved_", model_type, ".rds"))
    if(file.exists(fallback_file)) {
      unobserved_data[[model_type]] <- readRDS(fallback_file)
      cat("Pre-loaded fallback unobserved data for", model_type, ":", nrow(unobserved_data[[model_type]]), "sites\n")
    } else {
      cat("Warning: Unobserved data file not found:", unobserved_file, "\n")
    }
  }
}

# Get model list
model_id_list <- unique(site_effects_calibration$model_id)
model_id_list <- model_id_list[!is.na(model_id_list)]

# PRODUCTION: Process all models
cat("PRODUCTION MODE: Processing all", length(model_id_list), "CLR models\n")

cat("Processing", length(model_id_list), "CLR models sequentially\n")
cat("Models:", paste(model_id_list, collapse=", "), "\n")

# OPTIMIZATION: Pre-define predictor list to avoid repeated creation
# Mehlich III extractable metals replace former oxalate columns (megapit update)
pred_list <- c("MAT", "MAP", "latitude_scaled",
               "caNh4d", "kNh4d", "mgNh4d",
               "naNh4d", "cecdNh4", "feMjelm", "mnMjelm",
               "pMjelm", "siMjelm", "totalP")

# TESTING: Use sequential processing for debugging
cat("Using sequential processing for debugging\n")
dredge.out <- list()

for(s in model_id_list) {
  
  cat("Processing model:", s, "\n")
  
  tryCatch({
    # Determine model type
    if(grepl("^cycl_only", s)) {
      model_type <- "cycl_only"
    } else if(grepl("^env_cov", s)) {
      model_type <- "env_cov"
    } else if(grepl("^env_cycl", s)) {
      model_type <- "env_cycl"
    } else {
      cat("  Warning: Unknown model type for:", s, "\n")
      next
    }
    
    # Get pre-loaded unobserved data
    if(is.null(unobserved_data[[model_type]])) {
      cat("  Warning: No unobserved data for model type:", model_type, "\n")
      next
    }
    
    df_unobserved <- unobserved_data[[model_type]]
    
    # OPTIMIZATION: Streamlined data preparation
    species_dat <- df_observed %>%
      filter(model_id == s) %>%
      select(c("Median", "siteID", all_of(pred_list))) %>%
      na.omit() %>%
      rename(TargetVar = Median)
    
    if(nrow(species_dat) < 5) {
      cat("  Warning: Not enough data for model:", s, "(", nrow(species_dat), "sites)\n")
      next
    }
    
    # OPTIMIZATION: Pre-filter unobserved data
    new_data <- df_unobserved %>%
      select(all_of(pred_list)) %>%
      na.omit()
    
    if(nrow(new_data) == 0) {
      cat("  Warning: No unobserved data available for model:", s, "\n")
      next
    }
    
    # OPTIMIZATION: Use site_eff_uncertainties function with error handling
    cat("  Running PLSR with jackknife for model:", s, "\n")
    
    # Prepare data for site_eff_uncertainties
    df_unobserved_taxon <- new_data %>%
      mutate(siteID = row_number(), model_id = s) %>%
      select(siteID, all_of(pred_list), model_id)
    
    plsr_result <- tryCatch({
      site_eff_uncertainties(species_dat, df_unobserved_taxon)
    }, error = function(e) {
      cat("  PLSR error for model", s, ":", e$message, "\n")
      return(NULL)
    })
    
    if(is.null(plsr_result)) {
      cat("  Skipping model due to PLSR error:", s, "\n")
      next
    }
    
    # OPTIMIZATION: Limited dredge search space
    cat("  Running dredge for model:", s, "\n")
    
    # Create model formula with limited terms
    formula_str <- paste("TargetVar ~", paste(pred_list, collapse = " + "))
    formula_obj <- as.formula(formula_str)
    
    # Set global options for dredge
    old_options <- options()
    options(na.action = "na.fail")
    
    # Fit base model
    base_model <- lm(formula_obj, data = species_dat)
    
    # OPTIMIZATION: Limit dredge to max 5 predictors to reduce computation
    dredge_result <- tryCatch({
      dredge(base_model, 
             subset = dc(MAT, MAP, latitude_scaled, caNh4d, kNh4d, mgNh4d, naNh4d, cecdNh4, feMjelm, mnMjelm, pMjelm, siMjelm, totalP) <= 5,
             rank = "AICc",
             trace = FALSE)
    }, error = function(e) {
      cat("  Dredge error for model", s, ":", e$message, "\n")
      return(NULL)
    })
    
    # Restore original options
    options(old_options)
    
    if(is.null(dredge_result) || nrow(dredge_result) == 0) {
      cat("  No dredge results for model:", s, "\n")
      next
    }
    
    # OPTIMIZATION: Model averaging with error handling
    if(nrow(dredge_result) > 1) {
      model_avg <- tryCatch({
        model.avg(dredge_result, subset = delta < 2)
      }, error = function(e) {
        cat("  Model averaging error for", s, ":", e$message, "\n")
        return(NULL)
      })
    } else {
      model_avg <- NULL
    }
    
    # Store results
    dredge.out[[s]] <- list(
      dredge_result = dredge_result,
      model_avg = model_avg,
      plsr_result = plsr_result,
      n_sites = nrow(species_dat),
      n_unobserved = nrow(new_data)
    )
    
    cat("  Successfully processed model:", s, "\n")
    
  }, error = function(e) {
    cat("  Error processing model", s, ":", e$message, "\n")
  })
}

cat("Sequential processing completed!\n")
cat("Number of models successfully processed:", length(dredge.out), "\n")

# Save results
if(length(dredge.out) > 0) {
  saveRDS(dredge.out, here("data/summary/site_effects_dredged_CLR.rds"))
  cat("Results saved successfully!\n")
} else {
  cat("No results to save!\n")
}

cat("Performance improvements:\n")
cat("  - Sequential processing (can be switched to parallel)\n")
cat("  - Pre-loaded unobserved data files\n")
cat("  - Full PLSR with jackknife uncertainty estimation\n")
cat("  - Standard dredge search space (max 5 predictors)\n")
cat("  - Streamlined data operations\n")
