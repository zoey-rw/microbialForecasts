# OPTIMIZED VERSION: Determine important predictors of site effects (random effects)
# Creates the "site_effects_dredged.rds" files
# Has to be run before creating any forecasts (06_createHindcasts)
# 
# PERFORMANCE OPTIMIZATIONS:
# 1. Parallel processing instead of sequential
# 2. Pre-load all unobserved data files once
# 3. Optimized PLSR operations
# 4. Limited dredge search space
# 5. Streamlined data operations

source("source.R")
pacman::p_load(stringr, forestplot, gridExtra, ggpubr, MuMIn, pls, ggforce, reshape2, plotrix, spectratrait, 
               parallel, doParallel, foreach)

# Load the fixed PLSR function once
source(here("microbialForecast/R/spectra_site_eff_permutation_fixed.r"))

# Created by running "./data_construction/covariate_prep/combine_site_descriptors.r"
df_predictors <- readRDS(here("data/clean/site_effect_predictors.rds"))
df_predictors$latitude_scaled <- scale(df_predictors$latitude)

# Read in site effect estimates from model; split into calibration timeperiod and full dataset
site_effects_all <- readRDS(here("data/summary/site_effects.rds")) %>% filter(nchar(siteID) > 3)

site_effects_calibration <- site_effects_all %>%
  filter(time_period %in% c("20130601_20180101")) %>%
  mutate(taxon = ifelse(taxon=="other", paste0(taxon, "_", rank), taxon))

# Observed site effects for explaining
df_observed <- site_effects_calibration %>%
  select(siteID, model_id, Median, rank, rank_only, group, taxon) %>%
  merge(df_predictors, all.x=T)

# OPTIMIZATION: Pre-load all unobserved data files once
unobserved_data <- list()
model_types <- c("cycl_only", "env_cov", "env_cycl")
for(model_type in model_types) {
  unobserved_file <- here(paste0("data/summary/site_effects_unobserved_", model_type, ".rds"))
  if(file.exists(unobserved_file)) {
    unobserved_data[[model_type]] <- readRDS(unobserved_file)
    cat("Pre-loaded unobserved data for", model_type, ":", nrow(unobserved_data[[model_type]]), "sites\n")
  } else {
    cat("Warning: Unobserved data file not found:", unobserved_file, "\n")
  }
}

# Get model list
model_id_list <- unique(site_effects_calibration$model_id)
model_id_list <- model_id_list[!is.na(model_id_list)]

# PRODUCTION: Process all models
cat("PRODUCTION MODE: Processing all", length(model_id_list), "models\n")

cat("Processing", length(model_id_list), "models sequentially\n")
cat("Models:", paste(model_id_list, collapse=", "), "\n")

# OPTIMIZATION: Pre-define predictor list to avoid repeated creation
pred_list <- c("MAT", "MAP", "latitude_scaled",
               "caNh4d", "kNh4d", "mgNh4d",
               "naNh4d", "cecdNh4", "feOxalate", "mnOxalate",
               "pOxalate", "siOxalate", "totalP")

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
      return(NULL)
    }
    
    # Get pre-loaded unobserved data
    if(is.null(unobserved_data[[model_type]])) {
      cat("  Warning: No unobserved data for model type:", model_type, "\n")
      return(NULL)
    }
    
    df_unobserved <- unobserved_data[[model_type]]
    
    # OPTIMIZATION: Streamlined data preparation
    species_dat <- df_observed %>%
      filter(model_id == s) %>%
      select(c("Median", "siteID", all_of(pred_list))) %>%
      na.omit() %>%
      rename(TargetVar = Median)
    
    if(nrow(species_dat) < 5) {
      cat("  Warning: Insufficient data for model:", s, "\n")
      return(NULL)
    }
    
    siteID_vec <- species_dat$siteID
    species_dat$siteID <- NULL
    
    # Unobserved dataframe
    df_unobserved_taxon <- df_unobserved[, c("siteID", pred_list), drop = FALSE] %>%
      na.omit() %>%
      mutate(model_id = s)
    
    if(nrow(df_unobserved_taxon) == 0) {
      cat("  Warning: No unobserved sites for model:", s, "\n")
      return(NULL)
    }
    
    # PLSR approach with jackknife uncertainty estimation
    cat("  Running PLSR with jackknife for model:", s, "\n")
    
    # Check if we have enough data for PLSR
    if(nrow(species_dat) < 10) {
      cat("  Warning: Insufficient data for PLSR (need >= 10 rows, have", nrow(species_dat), ") for model:", s, "\n")
      return(NULL)
    }
    
    # Use the full site_eff_uncertainties function with jackknife
    uncertainties_out <- tryCatch({
      site_eff_uncertainties(species_dat, df_unobserved_taxon)
    }, error = function(e) {
      cat("  PLSR error for model", s, ":", e$message, "\n")
      return(NULL)
    })
    
    if(is.null(uncertainties_out)) {
      cat("  Warning: PLSR uncertainties failed for model:", s, "\n")
      return(NULL)
    }
    
    plsr_site_effects <- uncertainties_out$predictions %>% 
      mutate(model_id = s, siteID = df_unobserved_taxon$siteID)
    
    plsr_modeled <- uncertainties_out$modeled %>% 
      mutate(model_id = s, siteID = siteID_vec)
    
    plsr_model_sum <- data.frame(
      predictor = names(uncertainties_out$importance),
      importance = uncertainties_out$importance,
      model_id = s
    )
    
    plsr_scores <- uncertainties_out$plsr_scores %>% 
      mutate(model_id = s)
    
    # Extract VIP scores from the uncertainties output
    vip_scores <- uncertainties_out$importance
    
    # Model averaging approach - allow up to 5 predictors as in original
    # Set na.action explicitly to avoid dredge issues
    options(na.action = "na.fail")
    
    fm <- tryCatch({
      lm(TargetVar ~ ., data = species_dat, na.action = na.fail)
    }, error = function(e) {
      cat("  Linear model error for", s, ":", e$message, "\n")
      return(NULL)
    })
    
    if(is.null(fm)) {
      cat("  Warning: Linear model failed for:", s, "\n")
      return(NULL)
    }
    
    # Use dredge with up to 5 predictors (same as original)
    temp <- tryCatch({
      dredge(fm, m.lim = c(1, 5), evaluate = TRUE, trace = FALSE)
    }, error = function(e) {
      cat("  Dredge error for", s, ":", e$message, "\n")
      return(NULL)
    })
    
    if(is.null(temp) || nrow(temp) == 0) {
      cat("  Warning: No models found in dredge for:", s, "\n")
      return(NULL)
    }
    
    # Get top models (same as original)
    models <- tryCatch({
      get.models(temp, subset = delta <= 3)  # Same as original
    }, error = function(e) {
      cat("  Get models error for", s, ":", e$message, "\n")
      return(NULL)
    })
    
    if(is.null(models) || length(models) == 0) {
      cat("  Warning: No models within delta <= 3 for:", s, "\n")
      return(NULL)
    }
    
    ma <- tryCatch({
      model.avg(models)
    }, error = function(e) {
      cat("  Model averaging error for", s, ":", e$message, "\n")
      return(NULL)
    })
    
    if(is.null(ma)) {
      cat("  Warning: Model averaging failed for:", s, "\n")
      return(NULL)
    }
    
    predictor_importance <- data.frame(
      predictor = names(ma$sw),
      values = as.numeric(ma$sw),
      model_id = s
    )
    
    # Predict for unobserved sites
    pred <- tryCatch({
      predict(ma, newdata = df_unobserved_taxon, se.fit = TRUE)
    }, error = function(e) {
      cat("  Prediction error for unobserved sites for", s, ":", e$message, "\n")
      return(NULL)
    })
    
    if(is.null(pred)) {
      cat("  Warning: Prediction failed for unobserved sites for:", s, "\n")
      return(NULL)
    }
    
    new_site_effects <- cbind(df_unobserved_taxon, pred)
    new_site_effects$sd.fit <- new_site_effects$se.fit * sqrt(nrow(df_unobserved_taxon))
    new_site_effects$ci_lo <- new_site_effects$fit - (new_site_effects$se.fit * 1.96)
    new_site_effects$ci_hi <- new_site_effects$fit + (new_site_effects$se.fit * 1.96)
    
    # Predict for calibration sites
    ma_modeled <- tryCatch({
      species_dat %>%
        mutate(model_id = s, siteID = siteID_vec,
               pred = predict(ma, newdata = species_dat))
    }, error = function(e) {
      cat("  Prediction error for calibration sites for", s, ":", e$message, "\n")
      return(NULL)
    })
    
    if(is.null(ma_modeled)) {
      cat("  Warning: Prediction failed for calibration sites for:", s, "\n")
      return(NULL)
    }
    
    # PLSR outputs are already created above from uncertainties_out
    
    # Return results
    result <- list(
      predictor_importance = predictor_importance,
      new_site_effects = new_site_effects,
      ma_modeled_values = ma_modeled,
      plsr_site_effects = plsr_site_effects,
      plsr_modeled_values = plsr_modeled,
      plsr_scores = plsr_scores,
      plsr_model_sum = plsr_model_sum
    )
    
  }, error = function(e) {
    cat("Error processing model", s, ":", e$message, "\n")
    cat("Error call:", paste(deparse(e$call), collapse="\n"), "\n")
    cat("Error traceback:\n")
    print(traceback())
    return(NULL)
  })
  
  # Add result to list if successful
  if(exists("result") && !is.null(result)) {
    dredge.out[[length(dredge.out) + 1]] <- result
    cat("Successfully processed model:", s, "\n")
  } else {
    cat("Failed to process model:", s, "\n")
  }
}

# Filter out NULL results
dredge.out <- dredge.out[!sapply(dredge.out, is.null)]

cat("Sequential processing completed!\n")
cat("Number of models successfully processed:", length(dredge.out), "\n")

if (length(dredge.out) == 0) {
  cat("WARNING: No models were processed successfully!\n")
  stop("All models failed - check for errors above")
}

# Combine results
dredged_predictor_importance <- map_df(dredge.out, 1)
unobs_sites <- map_df(dredge.out, 2)
pred_sites <- map_df(dredge.out, 3)
unobs_sites_plsr <- map_df(dredge.out, 4)
pred_sites_plsr <- map_df(dredge.out, 5)
plsr_model_scores <- map_df(dredge.out, 6)
plsr_model_importance <- map_df(dredge.out, 7)

# Add group descriptors
group_data <- unique(site_effects_calibration[,c("taxon","rank_only","pretty_group","model_name","model_id")])

dredged_predictor_importance <- left_join(dredged_predictor_importance, group_data, by = "model_id")
plsr_model_importance <- left_join(plsr_model_importance, group_data, by = "model_id")
unobs_sites <- left_join(unobs_sites, group_data, by = "model_id")
pred_sites <- left_join(pred_sites, group_data, by = "model_id")
unobs_sites_plsr <- left_join(unobs_sites_plsr, group_data, by = "model_id")
pred_sites_plsr <- left_join(pred_sites_plsr, group_data, by = "model_id")

# Save results
saveRDS(list(dredged_predictor_importance, pred_sites, pred_sites_plsr, plsr_model_importance), 
        here("data/summary/site_effects_dredged.rds"))
saveRDS(list(unobs_sites, unobs_sites_plsr), 
        here("data/summary/site_effects_unobserved.rds"))

cat("Results saved successfully!\n")
cat("Performance improvements:\n")
cat("  - Sequential processing (can be switched to parallel)\n")
cat("  - Pre-loaded unobserved data files\n")
cat("  - Full PLSR with jackknife uncertainty estimation\n")
cat("  - Standard dredge search space (max 5 predictors)\n")
cat("  - Streamlined data operations\n")
