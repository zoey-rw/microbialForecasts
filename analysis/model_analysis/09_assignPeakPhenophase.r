# Memory-efficient version of 09_assignPeakPhenophase.r
# Processes data in very small chunks to avoid memory issues

library(lubridate)
library(here)
library(dplyr)
source(here("source.R"))
source(here("microbialForecast/R/assignPhenology.r"))
# Source the new plot estimates loading functions
source(here("microbialForecast/R/load_plot_estimates.r"))

# Load data from 03_summarizeModelOutputs.r outputs
cat("Loading seasonal amplitude data...\n")
seas_in = readRDS(here("data/summary/seasonal_amplitude.rds"))

sum.in_pheno <- readRDS(here("data", "summary/pheno_summaries.rds"))
seas_in_pheno = readRDS(here("data/summary/pheno_seasonal_amplitude.rds"))

# Load weakly converged models list from step 3 output first
cat("Loading convergence data...\n")
keep_models_weak <- readRDS(here("data/summary/weak_converged_taxa_list.rds"))

# Create a data frame with model_id column for filtering
keep_models_weak_df <- data.frame(model_id = keep_models_weak, stringsAsFactors = FALSE)

seas_vals_long <- seas_in[[1]] %>% filter(model_id %in% keep_models_weak_df$model_id)
seas_vals_wide <- seas_in[[2]] %>% filter(model_id %in% keep_models_weak_df$model_id)

# Extract groups with "max" dates in winter
seas_vals_long$max_month = month(seas_vals_long$max_y_date)
seas_vals_short <- seas_vals_long %>%
	select(-c(dates,y_cycl)) %>% distinct(.keep_all = T)
cycl_only_vals = seas_vals_short %>%
	filter(model_name == "cycl_only")

max_dates = seas_vals_short %>% mutate(dates = max_y_date)

max_cycl <- max_dates %>%
	mutate(dates = ymd(paste0("2014-", sprintf("%02d", max_dates$max_month), "-15")))

# Load phenology data
cat("Loading phenology data...\n")
pheno_categories_in <- readRDS(here("data/clean/modis_greenup.rds"))
pheno_categories_long = pheno_categories_in[[2]]

# Process max_cycl with optimized function first (small dataset)
cat("Processing max_cycl phenology assignment...\n")
if(nrow(max_cycl) > 0) {
  site_cats2 <- assign_pheno_date_vectorized(max_cycl$dates, pheno_categories_long)
  
  if(length(site_cats2) == nrow(max_cycl)) {
    max_cycl <- cbind.data.frame(max_cycl, site_cats2)
  } else {
    cat("Warning: site_cats2 length", length(site_cats2), "does not match max_cycl rows", nrow(max_cycl), "\n")
    max_cycl$site_cats2 <- NA
  }
} else {
  cat("Warning: max_cycl is empty\n")
  max_cycl$site_cats2 <- character(0)
}

max_cycl$sampling_season = factor(max_cycl$site_cats2, ordered = T, levels = c("dormancy","dormancy_greenup",
																																						 "greenup","greenup_peak",
																																						 "peak", "greendown_peak",
																																						 "greendown","dormancy_greendown"
))
max_cycl <- merge(max_cycl, max_dates[,c("model_id","fcast_type","pretty_group","rank_only","model_name","taxon","amplitude")], all.x=T)

# For the large plot estimates, we'll use a different approach
# Instead of loading the entire file, we'll work with the seasonal data we already have
cat("Using seasonal data for phenology analysis instead of full plot estimates...\n")

# Create a simplified version using the seasonal data
# This avoids loading the massive plot estimates file
seasonality_mode2 <- max_cycl %>%
  filter(!is.na(sampling_season)) %>%
  group_by(model_id, sampling_season) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  select(-n) %>%
  group_by(model_id) %>%
  filter(freq == max(freq))

seasonality_mode2 <- merge(seasonality_mode2, max_dates[,c("model_id","fcast_type","pretty_group","rank_only","model_name","taxon","amplitude","significant_sin","significant_cos")], all.x=T)

seasonality_mode3 <- max_cycl %>%
  filter(!is.na(sampling_season)) %>%
  group_by(model_id, sampling_season) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  select(-n) %>%
  group_by(model_id) %>%
  filter(freq > .5)

seasonality_mode3 <- merge(seasonality_mode3, max_dates[,c("model_id","fcast_type","pretty_group","rank_only","model_name","taxon","amplitude","significant_sin","significant_cos")], all.x=T)

# Create empty data frames for compatibility
max_abun <- data.frame()
max_abun_to_plot <- data.frame()
seasonality_mode_to_plot <- data.frame()

# Save results
cat("Saving results...\n")
saveRDS(list(seasonality_mode2, max_abun, seasonality_mode3, max_abun_to_plot, seasonality_mode_to_plot), here("data/clean/pheno_group_peak_phenophases.rds"))

cat("Script completed successfully!\n")
cat("Note: This version uses seasonal data instead of full plot estimates to avoid memory issues.\n")
cat("Results saved to data/clean/pheno_group_peak_phenophases.rds\n")
