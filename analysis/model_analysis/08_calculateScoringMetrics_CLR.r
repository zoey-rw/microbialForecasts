# Calculate scoring metrics for CLR hindcasts
library(here)
source(here("source.R"))
source(here("microbialForecast/R/statsFunctions.r"))
source(here("analysis/model_analysis/robust_add_scoring_metrics.R"))
library(data.table)

cat("Loading CLR hindcast data...\n")

# Load CLR hindcast data only
clr_hindcasts_file <- here("data/summary/all_hindcasts_CLR.rds")
clr_hindcasts_parquet <- here("data/summary/parquet/all_hindcasts_CLR.parquet")

if (file.exists(clr_hindcasts_file)) {
  hindcast_data <- readRDS(clr_hindcasts_file)
  cat("✅ CLR hindcasts loaded from RDS\n")
} else if (file.exists(clr_hindcasts_parquet)) {
  if (require(arrow, quietly = TRUE)) {
    hindcast_data <- arrow::read_parquet(clr_hindcasts_parquet)
    cat("✅ CLR hindcasts loaded from Parquet\n")
  } else {
    stop("Arrow package required for Parquet files")
  }
} else {
  stop("CLR hindcast file not found. Please run 07_tidyHindcasts_CLR.r first to generate CLR hindcast data.")
}

cat(sprintf("Loaded %d rows\n", nrow(hindcast_data)))

# Create subsets
hindcast_only <- hindcast_data %>% 
  filter(fcast_period == "hindcast", !is.na(truth), !is.na(mean))

calibration_only <- hindcast_data %>% 
  filter(fcast_period == "calibration", !is.na(truth), !is.na(mean))

# Filter calibration to exclude first dates if plot_start_date exists
if ("plot_start_date" %in% colnames(calibration_only)) {
  calibration_only_not_first <- calibration_only %>% 
    filter(date_num > plot_start_date)
  cat("Filtered calibration by plot_start_date\n")
} else {
  cat("Warning: plot_start_date not found, using all calibration data\n")
  calibration_only_not_first <- calibration_only
}

cat(sprintf("Hindcast: %d, Calibration: %d, Calibration (not first): %d\n", 
            nrow(hindcast_only), nrow(calibration_only), nrow(calibration_only_not_first)))

# Calculate scoring metrics for hindcast data (by site_prediction)
scoring_metrics <- as.data.table(hindcast_only) %>%
  .[!is.na(site_prediction)] %>%
  .[, {
    if (.N == 0) {
      list(RMSE = NA_real_, BIAS = NA_real_, MAE = NA_real_, 
           CRPS = NA_real_, RSQ = NA_real_, RSQ.1 = NA_real_,
           RMSE.iqr = NA_real_, RMSE.norm = NA_real_, 
           CRPS_truncated = NA_real_, residual_variance = NA_real_,
           predictive_variance = NA_real_, total_PL = NA_real_)
    } else {
      tryCatch({
        result <- robust_add_scoring_metrics(
          observed = truth,
          mean_predicted = mean,
          median_predicted = med,
          sd_predicted = sd
        )
        as.list(result)
      }, error = function(e) {
        cat("Warning in group:", paste(.BY, collapse = ", "), "-", e$message, "\n")
        list(RMSE = NA_real_, BIAS = NA_real_, MAE = NA_real_, 
             CRPS = NA_real_, RSQ = NA_real_, RSQ.1 = NA_real_,
             RMSE.iqr = NA_real_, RMSE.norm = NA_real_,
             CRPS_truncated = NA_real_, residual_variance = NA_real_,
             predictive_variance = NA_real_, total_PL = NA_real_)
      })
    }
  }, by = .(model_id, fcast_type, pretty_group, model_name, 
            pretty_name, rank_name, taxon, site_prediction)] %>%
  as.data.frame()

# Add mean_crps_sample for skill score analysis
scoring_metrics <- scoring_metrics %>%
  mutate(mean_crps_sample = CRPS)

# Calibration metrics (overall)
calibration_metrics <- as.data.table(calibration_only_not_first) %>%
  .[, {
    if (.N == 0) {
      list(RMSE = NA_real_, BIAS = NA_real_, MAE = NA_real_, 
           CRPS = NA_real_, RSQ = NA_real_, RSQ.1 = NA_real_,
           RMSE.iqr = NA_real_, RMSE.norm = NA_real_,
           CRPS_truncated = NA_real_, residual_variance = NA_real_,
           predictive_variance = NA_real_, total_PL = NA_real_)
    } else {
      tryCatch({
        result <- robust_add_scoring_metrics(
          observed = truth,
          median_predicted = med,
          mean_predicted = mean,
          sd_predicted = sd
        )
        as.list(result)
      }, error = function(e) {
        list(RMSE = NA_real_, BIAS = NA_real_, MAE = NA_real_, 
             CRPS = NA_real_, RSQ = NA_real_, RSQ.1 = NA_real_,
             RMSE.iqr = NA_real_, RMSE.norm = NA_real_,
             CRPS_truncated = NA_real_, residual_variance = NA_real_,
             predictive_variance = NA_real_, total_PL = NA_real_)
      })
    }
  }, by = .(model_id, fcast_type, pretty_group, model_name, 
            rank_name, pretty_name, taxon)] %>%
  as.data.frame()

# Calibration metrics by site
calibration_metrics_site <- as.data.table(calibration_only_not_first) %>%
  .[, {
    if (.N == 0) {
      list(RMSE = NA_real_, BIAS = NA_real_, MAE = NA_real_, 
           CRPS = NA_real_, RSQ = NA_real_, RSQ.1 = NA_real_,
           RMSE.iqr = NA_real_, RMSE.norm = NA_real_,
           CRPS_truncated = NA_real_, residual_variance = NA_real_,
           predictive_variance = NA_real_, total_PL = NA_real_)
    } else {
      tryCatch({
        result <- robust_add_scoring_metrics(
          observed = truth,
          median_predicted = med,
          mean_predicted = mean,
          sd_predicted = sd
        )
        as.list(result)
      }, error = function(e) {
        list(RMSE = NA_real_, BIAS = NA_real_, MAE = NA_real_, 
             CRPS = NA_real_, RSQ = NA_real_, RSQ.1 = NA_real_,
             RMSE.iqr = NA_real_, RMSE.norm = NA_real_,
             CRPS_truncated = NA_real_, residual_variance = NA_real_,
             predictive_variance = NA_real_, total_PL = NA_real_)
      })
    }
  }, by = .(model_id, fcast_type, pretty_group, model_name, 
            rank_name, pretty_name, taxon, siteID)] %>%
  as.data.frame()

# Hindcast metrics by site
scoring_metrics_site <- as.data.table(hindcast_only) %>%
  .[!is.na(site_prediction)] %>%
  .[, {
    if (.N <= 1) {
      list(RMSE = NA_real_, BIAS = NA_real_, MAE = NA_real_, 
           CRPS = NA_real_, RSQ = NA_real_, RSQ.1 = NA_real_,
           RMSE.iqr = NA_real_, RMSE.norm = NA_real_,
           CRPS_truncated = NA_real_, residual_variance = NA_real_,
           predictive_variance = NA_real_, total_PL = NA_real_)
    } else {
      tryCatch({
        result <- robust_add_scoring_metrics(
          observed = truth,
          mean_predicted = mean,
          median_predicted = med,
          sd_predicted = sd
        )
        as.list(result)
      }, error = function(e) {
        list(RMSE = NA_real_, BIAS = NA_real_, MAE = NA_real_, 
             CRPS = NA_real_, RSQ = NA_real_, RSQ.1 = NA_real_,
             RMSE.iqr = NA_real_, RMSE.norm = NA_real_,
             CRPS_truncated = NA_real_, residual_variance = NA_real_,
             predictive_variance = NA_real_, total_PL = NA_real_)
      })
    }
  }, by = .(model_id, siteID, site_prediction)] %>%
  as.data.frame()

# Pivot to long format
scoring_metrics_long <- scoring_metrics %>% pivot_metrics()
calibration_metrics_long <- calibration_metrics %>% pivot_metrics()
calibration_metrics_site_long <- calibration_metrics_site %>% 
  pivot_metrics() %>% 
  filter(!is.infinite(score), !is.nan(score))
scoring_metrics_site_long <- scoring_metrics_site %>% 
  pivot_metrics() %>% 
  filter(!is.infinite(score), !is.nan(score))

# Calculate coefficient of variation (full granularity)
truth_vals <- calibration_only %>%
  filter(model_name == "env_cycl")

# Per-plot CV, then average per site
cv_tax_per_plot <- truth_vals %>%
  group_by(pretty_group, rank_name, taxon, siteID, plotID) %>%
  summarize(per_plot_cv = calc_cv(truth), .groups = "drop")

cv_tax_per_plot_site <- cv_tax_per_plot %>%
  group_by(pretty_group, rank_name, taxon, siteID) %>%
  summarize(mean_per_plot_site_cv = mean(per_plot_cv, na.rm = TRUE), .groups = "drop")

cv_tax_per_plot_taxon <- cv_tax_per_plot %>%
  group_by(pretty_group, rank_name, taxon) %>%
  summarize(mean_per_plot_cv = mean(per_plot_cv, na.rm = TRUE), .groups = "drop")

# Per-site CV, then average per taxon
cv_tax_per_site <- truth_vals %>%
  group_by(pretty_group, rank_name, taxon, siteID) %>%
  summarize(per_site_cv = calc_cv(truth), .groups = "drop")

cv_tax_per_site_taxon <- cv_tax_per_site %>%
  group_by(pretty_group, rank_name, taxon) %>%
  summarize(mean_per_site_cv = mean(per_site_cv, na.rm = TRUE), .groups = "drop")

# Overall CV
cv_tax_overall <- truth_vals %>%
  group_by(pretty_group, pretty_name, rank_name, taxon) %>%
  summarize(overall_cv = calc_cv(truth), .groups = "drop")

# Combine all CV metrics
cv_tax <- cv_tax_overall %>%
  left_join(cv_tax_per_site_taxon, by = c("pretty_group", "rank_name", "taxon")) %>%
  left_join(cv_tax_per_site, by = c("pretty_group", "rank_name", "taxon"))

# Merge CV with scoring metrics
scoring_metrics_cv <- merge(scoring_metrics_long, cv_tax, all = TRUE)
scoring_metrics_cv_site <- merge(scoring_metrics_site_long, cv_tax_per_site)

# Pivot CV to long format for scaled analysis
scoring_metrics_cv_long <- scoring_metrics_cv %>% 
  pivot_longer(cols = c(overall_cv, per_site_cv, mean_per_site_cv),
               names_to = "cv_type", values_to = "cv")

# CV scaled by rank
cv_metric_scaled <- scoring_metrics_cv_long %>%
  group_by(pretty_group, pretty_name, rank_name, cv_type, metric) %>%
  mutate(CV_scale = scale(cv)[, 1],
         metric_scale = scale(score)[, 1]) %>%
  ungroup()

# Calculate skill scores if possible
available_types <- unique(scoring_metrics_long$site_prediction)
cat("Site prediction types:", paste(available_types, collapse = ", "), "\n")

if (length(available_types) > 1 && 
    "New time (observed site)" %in% available_types) {
  
  # CRPS-based skill score
  skill_score_taxon <- scoring_metrics_long %>%
    filter(metric == "CRPS_truncated") %>%
    pivot_wider(
      id_cols = c("model_id", "fcast_type", "pretty_group", "model_name", 
                  "pretty_name", "rank_name", "taxon"),
      values_from = "score", 
      names_from = "site_prediction"
    ) %>%
    mutate(
      skill_score = if ("New time x site (modeled effect)" %in% colnames(.)) {
        (1 - (`New time x site (modeled effect)` / `New time (observed site)`))
      } else NA_real_,
      skill_score_random = if ("New time x site (random effect)" %in% colnames(.)) {
        (1 - (`New time x site (random effect)` / `New time (observed site)`))
      } else NA_real_
    )
  
  # RMSE-based skill score
  skill_score_taxon_RMSE <- scoring_metrics_long %>%
    filter(metric == "RMSE.norm") %>%
    pivot_wider(
      id_cols = c("model_id", "fcast_type", "pretty_group", "model_name", 
                  "pretty_name", "rank_name", "taxon"),
      values_from = "score", 
      names_from = "site_prediction"
    ) %>%
    mutate(
      skill_score = if ("New time x site (modeled effect)" %in% colnames(.)) {
        (1 - (`New time x site (modeled effect)` / `New time (observed site)`))
      } else NA_real_,
      skill_score_random = if ("New time x site (random effect)" %in% colnames(.)) {
        (1 - (`New time x site (random effect)` / `New time (observed site)`))
      } else NA_real_
    )
  
  skill_score_rank <- skill_score_taxon %>%
    group_by(model_id, pretty_group, pretty_name, rank_name) %>%
    summarize(mean_skill_score = mean(skill_score, na.rm = TRUE), .groups = "drop")
} else {
  cat("Insufficient site_prediction types for skill score calculation\n")
  skill_score_taxon <- data.frame()
  skill_score_taxon_RMSE <- data.frame()
  skill_score_rank <- data.frame()
}

# Load convergence lists if available
unconverged <- if (file.exists(here("data/summary/unconverged_taxa_list.rds"))) {
  readRDS(here("data/summary/unconverged_taxa_list.rds"))
} else NULL

converged <- if (file.exists(here("data/summary/weak_converged_taxa_list.rds"))) {
  readRDS(here("data/summary/weak_converged_taxa_list.rds"))
} else NULL

converged_strict <- if (file.exists(here("data/summary/converged_taxa_list.rds"))) {
  readRDS(here("data/summary/converged_taxa_list.rds"))
} else NULL

# Save all outputs
to_save <- list(
  cv_metric_scaled = cv_metric_scaled,
  scoring_metrics_cv = scoring_metrics_cv,
  scoring_metrics_cv_site = scoring_metrics_cv_site,
  scoring_metrics_cv_long = scoring_metrics_cv_long,
  calibration_truth_vals = truth_vals,
  
  scoring_metrics = scoring_metrics,
  scoring_metrics_long = scoring_metrics_long,
  scoring_metrics_site = scoring_metrics_site,
  scoring_metrics_site_long = scoring_metrics_site_long,
  
  calibration_metrics = calibration_metrics,
  calibration_metrics_long = calibration_metrics_long,
  calibration_metrics_site = calibration_metrics_site,
  calibration_metrics_site_long = calibration_metrics_site_long,
  
  skill_score_taxon = skill_score_taxon,
  skill_score_taxon_RMSE = skill_score_taxon_RMSE,
  skill_score_rank = skill_score_rank,
  unconverged_list = unconverged,
  converged_list = converged,
  converged_strict_list = converged_strict
)

saveRDS(to_save, here("data/summary/scoring_metrics_CLR.rds"))
cat("Saved: scoring_metrics_CLR.rds with", length(to_save), "components\n")
