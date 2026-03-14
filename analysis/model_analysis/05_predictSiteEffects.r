# Determine important predictors of site effects (random effects)
# Creates the "site_effects_dredged.rds" files
# 06_createHindcasts must be run after this.

source("source.R")
pacman::p_load(stringr, forestplot, gridExtra, ggpubr, MuMIn, pls, ggforce, 
               reshape2, plotrix, spectratrait, parallel, doParallel, foreach, dplyr, tidyr, here)

# --- 1. DATA PREPARATION (DONE ONCE) ---

# Load predictors
df_predictors <- readRDS(here("data/clean/site_effect_predictors.rds"))
df_predictors$latitude_scaled <- scale(df_predictors$latitude)

# Predictor list
# Mehlich III extractable metals replace former oxalate columns (megapit update)
pred_list <- c("MAT", "MAP", "latitude_scaled",
               "caNh4d", "kNh4d", "mgNh4d",
               "naNh4d", "cecdNh4", "feMjelm", "mnMjelm",
               "pMjelm", "siMjelm", "totalP")

# Load and Filter Site Effects
site_effects_all <- readRDS(here("data/summary/site_effects.rds")) %>% 
  filter(nchar(siteID) > 3)

driver_pattern <- "^((cycl_only|env_cov|env_cycl)_.*_20130601_20180101_with_legacy_covariate_beta_regression)$"
site_effects_calibration <- site_effects_all %>%
  filter(grepl(driver_pattern, model_id)) %>%
  filter(time_period %in% c("20130601_20180101", "2013-2018", "2013-06_2018-01")) %>%
  mutate(taxon = ifelse(taxon=="other", paste0(taxon, "_", rank), taxon))

if (nrow(site_effects_calibration) == 0) stop("No valid models found.")

# Prepare Observed Data
df_observed <- site_effects_calibration %>%
  select(siteID, model_id, Median, rank, rank_only, group, taxon) %>%
  left_join(df_predictors, by = "siteID")

# Prepare Unobserved Data
observed_sites <- unique(site_effects_calibration$siteID)
all_sites <- unique(df_predictors$siteID)
unobserved_sites <- setdiff(all_sites, observed_sites)

unobserved_base <- df_predictors %>% 
    filter(siteID %in% unobserved_sites) %>%
  select(siteID, all_of(pred_list)) %>%
  na.omit()

# Store this in a list to export to workers
unobserved_data_list <- list(
  "cycl_only" = unobserved_base,
  "env_cov"   = unobserved_base,
  "env_cycl"  = unobserved_base
)

# --- 2. RESUME LOGIC & PRE-SPLITTING ---

source_file <- here("data/summary/site_effects.rds")
output_file <- here("data/summary/site_effects_dredged.rds")

all_models <- unique(site_effects_calibration$model_id)
processed_models <- character(0)

force_regenerate <- FALSE

# Check timestamps
if(file.exists(source_file) && file.exists(output_file)) {
  if(file.mtime(source_file) > file.mtime(output_file)) force_regenerate <- TRUE
}

# Check existing data
if(!force_regenerate && file.exists(output_file)) {
  tryCatch({
    existing <- readRDS(output_file)
    if(length(existing) >= 4) processed_models <- unique(existing[[4]]$model_id)
  }, error = function(e) warning("Could not read existing results to resume."))
}

models_to_run <- if(force_regenerate) all_models else setdiff(all_models, processed_models)

# TEST MODE: Limit to 20 models, sampling across all model types
if(exists("TEST_MODE") && TEST_MODE) {
  # Sample 7 from each type to get ~20 total (ensures representation)
  cycl_models <- grep("^cycl_only", models_to_run, value = TRUE)
  env_cov_models <- grep("^env_cov", models_to_run, value = TRUE)
  env_cycl_models <- grep("^env_cycl", models_to_run, value = TRUE)
  
  n_per_type <- min(7, min(length(cycl_models), length(env_cov_models), length(env_cycl_models)))
  if(n_per_type > 0) {
    models_to_run <- c(
      sample(cycl_models, min(n_per_type, length(cycl_models))),
      sample(env_cov_models, min(n_per_type, length(env_cov_models))),
      sample(env_cycl_models, min(n_per_type, length(env_cycl_models)))
    )
    cat("TEST MODE: Limited to", length(models_to_run), "models (sampled across all types)\n")
  }
}

cat("Models to process:", length(models_to_run), "\n")

if(length(models_to_run) == 0) stop("All models processed.")

# Split data for parallel processing
data_split <- df_observed %>%
  filter(model_id %in% models_to_run) %>%
  group_by(model_id) %>%
  group_split()

names(data_split) <- sapply(data_split, function(x) x$model_id[1])

# --- 3. PARALLEL CONFIGURATION ---

# Detect cores and leave one free (or use specified number for testing)
n_cores <- if(exists("TEST_CORES")) TEST_CORES else (parallel::detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

cat("Starting parallel processing on", n_cores, "cores...\n")
cat("Processing", length(data_split), "models...\n")

# --- 4. MAIN COMPUTATIONAL LOOP ---

results_list <- foreach(species_dat_raw = data_split, 
                        .packages = c("dplyr", "MuMIn", "pls", "stats", "stringr", "here", "spectratrait"),
                        .export = c("pred_list", "unobserved_data_list")) %dopar% {
  
  # CRITICAL FIX: Source the function INSIDE the worker
  # The here() context might shift in workers, so we ensure the path is absolute
  tryCatch({
    source(here("microbialForecast/R/spectra_site_eff_permutation_fixed.r"))
  }, error = function(e) {
    # If source fails, we can't proceed. Return NULL.
    return(NULL)
  })
  
  tryCatch({
    # --- A. Setup ---
    s <- species_dat_raw$model_id[1]
    
    model_type <- str_extract(s, "^(cycl_only|env_cov|env_cycl)")
    if(is.na(model_type)) return(NULL)
    
    df_unobserved <- unobserved_data_list[[model_type]]
    if(nrow(df_unobserved) == 0) return(NULL)
    
    # --- B. Aggregation ---
    species_dat <- species_dat_raw %>%
      group_by(siteID) %>%
      summarize(
        TargetVar = mean(Median, na.rm = TRUE),
        across(all_of(pred_list), ~ mean(.x, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      na.omit()
    
    # ROBUSTNESS 1: Min Sample Size
    if(nrow(species_dat) < 5) return(NULL)
    
    siteID_vec <- species_dat$siteID
    species_dat$siteID <- NULL
    
    #  Convert latitude_scaled from matrix to numeric
    if("latitude_scaled" %in% names(species_dat)) {
      species_dat$latitude_scaled <- as.numeric(species_dat$latitude_scaled)
    } 
    
    # Prepare unobserved for PLSR (needs siteID and predictors)
    df_unobserved_taxon <- df_unobserved %>% 
      select(siteID, all_of(pred_list)) %>%
      na.omit()
    
    if(nrow(df_unobserved_taxon) == 0) return(NULL)
    
    # --- C. PLSR Modeling ---
    uncertainties_out <- tryCatch({
      site_eff_uncertainties(observed_sites = species_dat, new_sites = df_unobserved_taxon)
    }, error = function(e) NULL)
    
    if(is.null(uncertainties_out)) return(NULL)
    
    # Extract PLSR results
    plsr_site_effects <- uncertainties_out$predictions %>% 
      mutate(model_id = s, siteID = df_unobserved_taxon$siteID)
    
    plsr_modeled <- uncertainties_out$modeled %>% 
      mutate(model_id = s, siteID = siteID_vec)
    
    plsr_model_sum <- data.frame(
      predictor = names(uncertainties_out$importance),
      importance = uncertainties_out$importance,
      model_id = s
    )
    
    plsr_scores <- as.data.frame(uncertainties_out$plsr_scores)
    plsr_scores$predictor <- rownames(plsr_scores)
    plsr_scores$model_id <- s

    # Per-component CV stats
    plsr_cv_stats <- uncertainties_out$cv_stats_per_ncomp
    plsr_cv_stats$model_id <- s

    # --- D. Model Averaging (Robust Version) ---
    
    curr_pred_list <- pred_list
    
    # ROBUSTNESS 2: Remove Zero Variance Predictors
    # If a predictor has 0 variance (constant), lm will fail or produce NAs
    nzv <- sapply(species_dat[, curr_pred_list], var) > 1e-6
    curr_pred_list <- curr_pred_list[nzv]
    
    # If no predictors left, return PLSR only
    if(length(curr_pred_list) == 0) {
      return(list(NULL, NULL, NULL, plsr_site_effects, plsr_modeled, plsr_scores, plsr_model_sum))
    }
    
    # ROBUSTNESS 3: Sample Size vs Predictors
    # Ensure N >= k + 2
    max_k <- max(1, nrow(species_dat) - 2)
    
    if(length(curr_pred_list) > max_k) {
      # Select top K correlated predictors
      cor_vals <- abs(cor(species_dat$TargetVar, species_dat[, curr_pred_list]))
      # Handle case where cor returns NA (perfect fit or flat line)
      cor_vals[is.na(cor_vals)] <- 0 
      
      # Select top 5 or max_k, whichever is smaller
      n_select <- min(5, max_k)
      curr_pred_list <- names(sort(cor_vals[1,], decreasing = TRUE)[seq_len(n_select)])
    }
    
    # Fit Global Model
    form <- as.formula(paste("TargetVar ~", paste(curr_pred_list, collapse = "+")))
    fm <- tryCatch({
      lm(form, data = species_dat, na.action = na.fail)
    }, error = function(e) NULL)
    
    # ROBUSTNESS 4: Alias Check (Multicollinearity)
    # If fm is not NULL but has NA coefficients
    if(!is.null(fm) && any(is.na(coef(fm)))) {
      # Identify aliased terms
      aliased <- tryCatch({ attributes(alias(fm)$Complete)$dimnames[[1]] }, error = function(e) NULL)
      if(!is.null(aliased)) {
        curr_pred_list <- setdiff(curr_pred_list, aliased)
        if(length(curr_pred_list) > 0) {
           form <- as.formula(paste("TargetVar ~", paste(curr_pred_list, collapse = "+")))
           fm <- tryCatch({ lm(form, data = species_dat, na.action = na.fail) }, error = function(e) NULL)
    } else {
           fm <- NULL
        }
      }
    }
    
    ma_results_list <- list(imp=NULL, eff=NULL, mod=NULL)
    
    if(!is.null(fm)) {
      # Run Dredge
      # m.lim must be strictly less than sample size
      dredge_res <- tryCatch({
        dredge(fm, m.lim = c(1, min(5, length(curr_pred_list))), evaluate = TRUE, trace = FALSE)
      }, error = function(e) NULL)
      
      if(!is.null(dredge_res) && nrow(dredge_res) > 0) {
        models <- tryCatch({ get.models(dredge_res, subset = delta <= 3) }, error = function(e) NULL)
        
        if(!is.null(models) && length(models) > 0) {
          
          # Model Averaging
          ma <- NULL
          sw_vec <- NULL
          
          if(length(models) == 1) {
            ma <- models[[1]]
            cfs <- coef(ma)
            sw_vec <- rep(1.0, length(cfs)-1)
            names(sw_vec) <- names(cfs)[names(cfs) != "(Intercept)"]
    } else {
            ma <- tryCatch({ model.avg(models) }, error = function(e) NULL)
            if(!is.null(ma)) sw_vec <- ma$sw
          }
          
          if(!is.null(ma) && !is.null(sw_vec)) {
            ma_results_list$imp <- data.frame(
        predictor = names(sw_vec),
        values = as.numeric(sw_vec),
        model_id = s
      )
            
            # Predictions for unobserved sites
            # Need to add model_id back for predictions
            # Also ensure latitude_scaled is numeric (not matrix)
            df_unobserved_ma <- df_unobserved %>% 
              mutate(model_id = s)
            
            # CRITICAL FIX: Convert latitude_scaled to numeric if present
            if("latitude_scaled" %in% names(df_unobserved_ma)) {
              df_unobserved_ma$latitude_scaled <- as.numeric(df_unobserved_ma$latitude_scaled)
            }
            
            # Ensure only predictors used in model are present
            model_vars <- names(coef(ma))
            model_vars <- model_vars[model_vars != "(Intercept)"]
            df_unobserved_ma <- df_unobserved_ma[, c("siteID", "model_id", model_vars), drop = FALSE]
            
    pred <- tryCatch({
              predict(ma, newdata = df_unobserved_ma, se.fit = TRUE)
            }, error = function(e) NULL)
            
            if(!is.null(pred)) {
              ma_results_list$eff <- cbind(df_unobserved_ma, pred) %>%
                mutate(
                  sd.fit = se.fit * sqrt(nrow(df_unobserved_ma)),
                  ci_lo = fit - (se.fit * 1.96),
                  ci_hi = fit + (se.fit * 1.96)
                )
              
              # For calibration predictions, ensure latitude_scaled is numeric
              species_dat_pred <- species_dat
              if("latitude_scaled" %in% names(species_dat_pred)) {
                species_dat_pred$latitude_scaled <- as.numeric(species_dat_pred$latitude_scaled)
              }
              
              ma_results_list$mod <- species_dat_pred %>%
                mutate(
                  model_id = s, 
                  siteID = siteID_vec,
                  pred = predict(ma, newdata = species_dat_pred)
                )
            }
          }
        }
      }
    }
    
    # Return formatted list (order matters for extraction later)
    return(list(
      ma_results_list$imp,
      ma_results_list$eff,
      ma_results_list$mod,
      plsr_site_effects,
      plsr_modeled,
      plsr_scores,
      plsr_model_sum,
      plsr_cv_stats
    ))
    
  }, error = function(e) {
    return(NULL) 
  })
}

stopCluster(cl)
cat("Parallel processing complete.\n")
cat("Raw results collected:", length(results_list), "items\n")

# --- 5. RESULT AGGREGATION ---

results_list <- results_list[!sapply(results_list, is.null)]
cat("Non-null results:", length(results_list), "models\n")

if(length(results_list) == 0 && length(processed_models) == 0) {
  stop("All models failed.")
}

extract_df <- function(lst, idx) {
  res <- lapply(lst, `[[`, idx)
  res <- res[!sapply(res, is.null)]
  if(length(res) > 0) dplyr::bind_rows(res) else data.frame()
}

# Extract results
new_dredged_imp <- extract_df(results_list, 1)
new_unobs_sites <- extract_df(results_list, 2)
new_pred_sites  <- extract_df(results_list, 3)
new_unobs_plsr  <- extract_df(results_list, 4)
new_pred_plsr   <- extract_df(results_list, 5)
new_plsr_scores <- extract_df(results_list, 6)
new_plsr_imp    <- extract_df(results_list, 7)
new_plsr_cv_stats <- extract_df(results_list, 8)

# Add group metadata
group_data <- site_effects_calibration %>%
  select(taxon, rank_only, pretty_group, model_name, model_id) %>%
  distinct()

join_groups <- function(df) {
  if(nrow(df) == 0) return(df)
  # Remove duplicates if they exist before joining
  if("taxon" %in% names(df)) return(df)
  if(!"model_id" %in% names(df)) return(df)
  dplyr::left_join(df, group_data, by = "model_id")
}

new_dredged_imp   <- join_groups(new_dredged_imp)
new_unobs_sites   <- join_groups(new_unobs_sites)
new_pred_sites    <- join_groups(new_pred_sites)
new_unobs_plsr    <- join_groups(new_unobs_plsr)
new_pred_plsr     <- join_groups(new_pred_plsr)
new_plsr_imp      <- join_groups(new_plsr_imp)
new_plsr_cv_stats <- join_groups(new_plsr_cv_stats)

# --- 6. MERGE AND SAVE ---

if(length(processed_models) > 0 && !force_regenerate) {
  ex_d <- readRDS(here("data/summary/site_effects_dredged.rds"))
  ex_u <- readRDS(here("data/summary/site_effects_unobserved.rds"))

  final_dredged_imp <- dplyr::bind_rows(ex_d[[1]], new_dredged_imp)
  final_pred_sites  <- dplyr::bind_rows(ex_d[[2]], new_pred_sites)
  final_pred_plsr   <- dplyr::bind_rows(ex_d[[3]], new_pred_plsr)
  final_plsr_imp    <- dplyr::bind_rows(ex_d[[4]], new_plsr_imp)
  final_plsr_scores <- if(length(ex_d) >= 5) dplyr::bind_rows(ex_d[[5]], new_plsr_scores) else new_plsr_scores
  final_plsr_cv_stats <- if(length(ex_d) >= 6) dplyr::bind_rows(ex_d[[6]], new_plsr_cv_stats) else new_plsr_cv_stats

  final_unobs_sites <- dplyr::bind_rows(ex_u[[1]], new_unobs_sites)
  final_unobs_plsr  <- dplyr::bind_rows(ex_u[[2]], new_unobs_plsr)
} else {
  final_dredged_imp   <- new_dredged_imp
  final_pred_sites    <- new_pred_sites
  final_pred_plsr     <- new_pred_plsr
  final_plsr_imp      <- new_plsr_imp
  final_plsr_scores   <- new_plsr_scores
  final_plsr_cv_stats <- new_plsr_cv_stats
  final_unobs_sites   <- new_unobs_sites
  final_unobs_plsr    <- new_unobs_plsr
}

saveRDS(list(final_dredged_imp, final_pred_sites, final_pred_plsr, final_plsr_imp, final_plsr_scores, final_plsr_cv_stats),
        here("data/summary/site_effects_dredged.rds"))

saveRDS(list(final_unobs_sites, final_unobs_plsr),
        here("data/summary/site_effects_unobserved.rds"))

# --- 7. SAVE HINDCAST FILES ---

if(nrow(final_unobs_plsr) > 0) {
  hindcast_data <- final_unobs_plsr
  if("Median" %in% names(hindcast_data)) hindcast_data$fit <- hindcast_data$Median
  
  if(!"model_name" %in% names(hindcast_data)) {
    hindcast_data$model_name <- str_extract(hindcast_data$model_id, "^(env_cycl|cycl_only|env_cov)")
  }
  
  for(mt in c("env_cycl", "env_cov", "cycl_only")) {
    sub_dat <- hindcast_data %>% 
      filter(model_name == mt) %>%
      select(siteID, fit, any_of(c("se_fit", "model_id", "taxon", "rank_only", "pretty_group")))
    
    if(nrow(sub_dat) > 0) {
      saveRDS(sub_dat, here("data/summary", paste0("site_effects_unobserved_", mt, ".rds")))
    }
  }
}

cat("Process Completed Successfully.\n")
