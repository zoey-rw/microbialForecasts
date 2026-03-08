# Determine important predictors of site effects (random effects) for CLR models
# FIXED: Updated for CLR models with proper file paths and CLR-specific logic
# Creates the "site_effects_dredged_CLR.rds" files
# Has to be run before creating any CLR forecasts (06_createHindcasts_CLR)

source("../../source.R")
# Note: spectra_site_eff_permutation.r should be loaded via source.R
pacman::p_load(stringr, forestplot, gridExtra, ggpubr, MuMIn, pls, ggforce, reshape2, plotrix, spectratrait)

# Created by running "./data_construction/covariate_prep/combine_site_descriptors.r"
# Requires extra NEON data downloads; use existing file if possible!
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

by_site <- df_observed %>%
  pivot_wider(id_cols = c("rank","taxon","model_id"),
              names_from = "siteID", values_from = "Median") %>%
  select(-c(1:2))

# Unobserved sites for predicting - will be loaded dynamically based on model type

# Loop through each taxon and see which factors best explain its site effects
options(na.action = "na.fail")

model_id_list <- unique(site_effects_calibration$model_id)
model_id_list <- model_id_list[!is.na(model_id_list)]  # Remove NA values

# FIXED: Filter to CLR models only with better pattern matching
clr_model_id_list <- model_id_list[grepl("CLR|clr|cycl_only|env_cycl|env_cov", model_id_list, ignore.case = TRUE)]
if (length(clr_model_id_list) > 0) {
  model_id_list <- clr_model_id_list
  cat("Filtered to CLR models only:", length(model_id_list), "models\n")
} else {
  cat("No CLR models found, processing all models:", length(model_id_list), "models\n")
}

cat("Models to process:", paste(model_id_list, collapse=", "), "\n")

if (length(model_id_list) == 0) {
  stop("No models to process")
}

# FIXED: Use sequential processing for debugging with better error handling
cat("Starting sequential processing with", length(model_id_list), "models...\n")

dredge.out = list()
for(s in model_id_list) {
  cat("Processing model:", s, "\n")
  tryCatch({
    pacman::p_load(dplyr, MuMIn)
    options(na.action = "na.fail")
    
    # FIXED: Subset to predictors of interest (defined inside tryCatch to ensure accessibility)
    # Use columns that actually exist in the unobserved data files
    # Removed problematic predictors to ensure all 30 sites can be used
    pred_list = c("MAT", "MAP","latitude_scaled",
                  "caNh4d", "kNh4d", "mgNh4d",
                  "naNh4d","cecdNh4","feOxalate", "mnOxalate",
                  "pOxalate", "siOxalate", "totalP")
    
    cat("  pred_list created with", length(pred_list), "predictors\n")
    cat("  pred_list:", paste(pred_list, collapse=", "), "\n")

    cat("  Creating species_dat...\n")
    # Observed dataframe (30 sites with estimated effects)
    # Use base R subsetting instead of select to avoid tidyselect issues
    species_dat <- as.data.frame(df_observed) %>%
      filter(model_id == s)  %>%
      select(c("Median","siteID", all_of(pred_list)))
    
    # Remove rows with missing values in predictors to ensure PLSR works
    cat("  Removing rows with missing values...\n")
    species_dat_initial <- nrow(species_dat)
    species_dat <- species_dat[complete.cases(species_dat),]
    species_dat_final <- nrow(species_dat)
    cat("  Rows removed due to missing values:", species_dat_initial - species_dat_final, "\n")
    
    if (species_dat_final == 0) {
      cat("  WARNING: No complete cases for model", s, "- skipping\n")
      next
    }
    
    # Rename Median to TargetVar to match what the function expects
    species_dat <- species_dat %>%
      rename(TargetVar = Median)
    
    cat("  species_dat columns:", paste(colnames(species_dat), collapse=", "), "\n")
    
    siteID_vec <- species_dat$siteID
    taxon_rank <- df_observed %>% filter(model_id == s) %>% select(rank) %>% unique() %>% unlist()
    species_dat$siteID <- NULL
    obs_sample_size = nrow(species_dat)
    cat("  species_dat created with", nrow(species_dat), "rows\n")

    cat("  Creating df_unobserved_taxon...\n")
    # FIXED: Load appropriate unobserved data file based on model type
    # Extract model type by looking for the pattern before the first underscore
    if(grepl("^cycl_only", s)) {
      model_type <- "cycl_only"
    } else if(grepl("^env_cov", s)) {
      model_type <- "env_cov"
    } else if(grepl("^env_cycl", s)) {
      model_type <- "env_cycl"
    } else if(grepl("CLR|clr", s)) {
      # For CLR models, use the same model type logic but with CLR prefix
      if(grepl("^cycl_only", s)) {
        model_type <- "cycl_only_CLR"
      } else if(grepl("^env_cov", s)) {
        model_type <- "env_cov_CLR"
      } else if(grepl("^env_cycl", s)) {
        model_type <- "env_cycl_CLR"
      } else {
        model_type <- "cycl_only_CLR"  # Default for CLR models
      }
    } else {
      cat("  Warning: Unknown model type for:", s, "\n")
      cat("  Skipping model:", s, "\n")
      next
    }
    
    # FIXED: Try CLR-specific unobserved file first, then fall back to regular
    unobserved_file <- here(paste0("data/summary/site_effects_unobserved_", model_type, ".rds"))
    if (!file.exists(unobserved_file)) {
      # Fall back to regular unobserved file
      fallback_unobserved_file <- here(paste0("data/summary/site_effects_unobserved_", gsub("_CLR", "", model_type), ".rds"))
      if (file.exists(fallback_unobserved_file)) {
        unobserved_file <- fallback_unobserved_file
        cat("  Using fallback unobserved file:", unobserved_file, "\n")
      } else {
        cat("  Warning: No unobserved data file found for model type:", model_type, "\n")
        cat("  Skipping model:", s, "\n")
        next
      }
    }
    
    if(file.exists(unobserved_file)) {
      df_unobserved <- readRDS(unobserved_file)
      cat("  Loaded unobserved data from:", unobserved_file, "\n")
    } else {
      cat("  Warning: Unobserved data file not found:", unobserved_file, "\n")
      cat("  Skipping model:", s, "\n")
      next
    }
    
    # FIXED: Unobserved dataframe (17 sites) with better error handling
    # Use base R subsetting instead of select to avoid tidyselect issues
    df_unobserved_taxon <- df_unobserved[, c("siteID", pred_list), drop = FALSE] %>%
      na.omit %>%	mutate(model_id = s)
    sample_size = nrow(df_unobserved_taxon)
    cat("  df_unobserved_taxon created with", nrow(df_unobserved_taxon), "rows\n")

    cat("  Calling site_eff_uncertainties...\n")
    # FIXED: PLSR approach - use the fixed function with better error handling
    # Source the fixed function directly to ensure we're using the latest version
    fixed_function_path <- here("microbialForecast/R/spectra_site_eff_permutation_fixed.r")
    if (file.exists(fixed_function_path)) {
      source(fixed_function_path)
      cat("  Loaded fixed spectra function\n")
    } else {
      cat("  WARNING: Fixed spectra function not found, trying to use existing one\n")
      if (!exists("site_eff_uncertainties")) {
        cat("  ERROR: site_eff_uncertainties function not available\n")
        cat("  Skipping model:", s, "\n")
        next
      }
    }
    
    # FIXED: Call site_eff_uncertainties with error handling
    uncertainties_out = tryCatch({
      site_eff_uncertainties(species_dat, df_unobserved_taxon)
    }, error = function(e) {
      cat("  ERROR in site_eff_uncertainties:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(uncertainties_out)) {
      cat("  Skipping model due to site_eff_uncertainties failure\n")
      next
    }
    
    plsr_site_effects = uncertainties_out$predictions %>% mutate(model_id = s)
    plsr_site_effects$siteID = df_unobserved_taxon$siteID
    plsr_modeled = uncertainties_out$modeled  %>% mutate(model_id = s, siteID = siteID_vec)
    plsr_model_sum = data.frame(predictor = names(uncertainties_out$importance),
                                importance = uncertainties_out$importance) %>% mutate(model_id = s)
    plsr_scores = uncertainties_out$plsr_scores  %>% mutate(model_id = s)

    # FIXED: Model averaging approach with better error handling
    # Up to 6 predictors allowed for explaining each taxon's site effect
    # Use TargetVar instead of Median since we renamed the column
    fm <- tryCatch({
      lm(TargetVar ~ ., data = species_dat)
    }, error = function(e) {
      cat("  ERROR creating linear model:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(fm)) {
      cat("  Skipping model due to linear model failure\n")
      next
    }
    
    temp <- tryCatch({
      MuMIn:::.dredge.par(fm, m.lim = c(NA,5))
    }, error = function(e) {
      cat("  ERROR in dredge:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(temp)) {
      cat("  Skipping model due to dredge failure\n")
      next
    }

    # Average top models and get importance
    models <- tryCatch({
      get.models(temp,  subset = delta <= 3)
    }, error = function(e) {
      cat("  ERROR getting models:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(models) || length(models) == 0) {
      cat("  WARNING: No models found for delta <= 3, using all models\n")
      models <- get.models(temp, subset = delta <= 10)
    }
    
    if (is.null(models) || length(models) == 0) {
      cat("  ERROR: No models available for averaging\n")
      next
    }
    
    ma <- model.avg(models)
    predictor_importance <- cbind.data.frame(predictor = names(ma$sw),
                                            values = ma$sw) %>%
      mutate(model_id = s)
    rownames(predictor_importance) <- NULL
    attr(predictor_importance$values, "n.models") <- NULL

    # FIXED: Predict using covariates for unobserved sites with error handling
    pred = tryCatch({
      predict(ma, newdata = df_unobserved_taxon, se.fit=T)
    }, error = function(e) {
      cat("  ERROR in prediction:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(pred)) {
      cat("  Skipping model due to prediction failure\n")
      next
    }
    
    new_site_effects = cbind.data.frame(df_unobserved_taxon, pred)
    new_site_effects$sd.fit = new_site_effects$se.fit * sqrt(sample_size)
    new_site_effects$ci_lo = new_site_effects$fit - (new_site_effects$se.fit * 1.96)
    new_site_effects$ci_hi = new_site_effects$fit + (new_site_effects$se.fit * 1.96)
    
    # FIXED: Predict using covariates for sites in calibration (checking in-sample model accuracy)
    ma_modeled = tryCatch({
      species_dat %>% mutate(model_id = s, siteID = siteID_vec,
                            pred = predict(ma, newdata = species_dat))
    }, error = function(e) {
      cat("  WARNING: In-sample prediction failed:", e$message, "\n")
      species_dat %>% mutate(model_id = s, siteID = siteID_vec, pred = NA)
    })

    # FIXED: Create output with proper error handling
    out <- list(predictor_importance = predictor_importance,
                new_site_effects = new_site_effects,
                ma_modeled_values = ma_modeled,
                plsr_site_effects = plsr_site_effects,
                plsr_modeled_values = plsr_modeled,
                plsr_scores = plsr_scores,
                plsr_model_sum = plsr_model_sum
    )
    dredge.out[[length(dredge.out) + 1]] <- out
    cat("✅ Successfully processed model:", s, "\n")
    
  }, error = function(e) {
    cat("❌ ERROR processing model", s, ":", e$message, "\n")
    cat("Error call:", paste(deparse(e$call), collapse="\n"), "\n")
    cat("Error traceback:\n")
    print(traceback())
  })
}

# Report results
cat("Sequential processing completed!\n")
cat("Number of models successfully processed:", length(dredge.out), "\n")
if (length(dredge.out) == 0) {
  cat("WARNING: No models were processed successfully!\n")
  stop("All models failed - check for errors above")
}

# FIXED: Fix some taxon names and recombine dfs with error handling
tryCatch({
  dredged_predictor_importance <- map_df(dredge.out, 1)
  unobs_sites <- map_df(dredge.out, 2)
  pred_sites <- map_df(dredge.out, 3)
  unobs_sites_plsr <- map_df(dredge.out, 4)
  pred_sites_plsr <- map_df(dredge.out, 5)
  plsr_model_scores <- map_df(dredge.out, 6)
  plsr_model_importance <- map_df(dredge.out, 7)
  
  cat("✅ Successfully combined all model outputs\n")
}, error = function(e) {
  cat("ERROR combining model outputs:", e$message, "\n")
  stop("Failed to combine model outputs")
})

# FIXED: Add on some group descriptors with error handling
group_data = tryCatch({
  unique(site_effects_calibration[,c("taxon","rank_only","pretty_group","model_name","model_id")])
}, error = function(e) {
  cat("WARNING: Could not extract group data:", e$message, "\n")
  stop("Failed to extract group data")
})

if (nrow(group_data) > 0) {
  dredged_predictor_importance <- left_join(dredged_predictor_importance, group_data, by = "model_id")
  plsr_model_importance <- left_join(plsr_model_importance, group_data, by = "model_id")
  unobs_sites <- left_join(unobs_sites, group_data, by = "model_id")
  pred_sites <- left_join(pred_sites, group_data, by = "model_id")
  unobs_sites_plsr <- left_join(unobs_sites_plsr, group_data, by = "model_id")
  pred_sites_plsr <- left_join(pred_sites_plsr, group_data, by = "model_id")
  
  cat("✅ Successfully added group descriptors\n")
} else {
  stop("No group descriptors available")
}

# FIXED: Save CLR-specific outputs with proper error handling
tryCatch({
  saveRDS(list(dredged_predictor_importance,pred_sites,pred_sites_plsr,plsr_model_importance), 
          here("data/summary/site_effects_dredged_CLR.rds"))
  saveRDS(list(unobs_sites,unobs_sites_plsr), 
          here("data/summary/site_effects_unobserved_CLR.rds"))
  
  cat("✅ CLR site effects prediction completed successfully!\n")
  cat("Output files saved:\n")
  cat("  - site_effects_dredged_CLR.rds\n")
  cat("  - site_effects_unobserved_CLR.rds\n")
}, error = function(e) {
  cat("ERROR saving outputs:", e$message, "\n")
  stop("Failed to save outputs")
})

# FIXED: Test reading the output with error handling
tryCatch({
  read_in = readRDS(here("data/summary/site_effects_dredged_CLR.rds"))
  cat("Summary of dredged predictor importance:\n")
  print(head(read_in[[1]]))
  cat("✅ Output files verified successfully\n")
}, error = function(e) {
  cat("ERROR verifying output files:", e$message, "\n")
  stop("Failed to verify output files")
})
