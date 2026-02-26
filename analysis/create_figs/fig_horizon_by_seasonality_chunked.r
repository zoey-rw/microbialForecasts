# Chunked version of forecast horizon by seasonality analysis
# This version processes data in smaller chunks to avoid memory issues

library(lubridate)
library(ggrepel)
library(ggplot2)
library(ggpubr)
source("source.R")
source("microbialForecast/R/assignPhenology.r")
source("microbialForecast/R/load_plot_estimates.r")

# Read in the latest summary data
sum.in <- readRDS(here("data", "summary/logit_beta_fixed_priors_summaries.rds"))

# Read in forecast horizon data
horizon_in <- readRDS(here("data", "summary/fcast_horizon_input.rds"))
fcast_horizon_data <- horizon_in[[2]]

# Add taxon column to horizon data (same as species)
fcast_horizon_data$taxon <- fcast_horizon_data$species

# Check if seasonality analysis results already exist
seasonality_file <- here("data/summary/seasonality_analysis_results.rds")

if(file.exists(seasonality_file)) {
  cat("Loading existing seasonality analysis results...\n")
  tukey_median_pheno_sig <- readRDS(seasonality_file)
  cat("Models with significant seasonality differences:", sum(tukey_median_pheno_sig$significant_diff, na.rm = TRUE), "\n")
} else {
  cat("Seasonality analysis results not found. Processing models in chunks...\n")
  
  # Get the models we need from the horizon data
  keep_models <- unique(fcast_horizon_data$model_id)
  cat("Processing", length(keep_models), "models for seasonality analysis...\n")
  
  # Load phenology categories data once
  cat("Loading phenology categories...\n")
  pheno_categories <- readRDS(here("data/clean/modis_greenup.rds"))
  pheno_categories_long <- pheno_categories[[2]]
  
  # Process models in chunks of 20
  chunk_size <- 20
  n_chunks <- ceiling(length(keep_models) / chunk_size)
  all_tukey_results <- list()
  
  for(chunk_i in 1:n_chunks) {
    start_idx <- (chunk_i - 1) * chunk_size + 1
    end_idx <- min(chunk_i * chunk_size, length(keep_models))
    chunk_models <- keep_models[start_idx:end_idx]
    
    cat(sprintf("Processing chunk %d/%d (models %d-%d)...\n", 
                chunk_i, n_chunks, start_idx, end_idx))
    
    # Load plot estimates for this chunk
    cat("Loading plot estimates for chunk...\n")
    plot_estimates_env_cycl <- load_plot_estimates_for_phenology("env_cycl")
    plot_estimates_cycl_only <- load_plot_estimates_for_phenology("cycl_only")
    
    # Combine and filter to only the models we need
    plot_estimates <- rbind(plot_estimates_env_cycl, plot_estimates_cycl_only, fill = TRUE)
    plot_estimates <- plot_estimates %>% filter(model_id %in% chunk_models)
    
    cat("Chunk plot estimates rows:", nrow(plot_estimates), "\n")
    
    if(nrow(plot_estimates) == 0) {
      cat("No data for this chunk, skipping...\n")
      next
    }
    
    # Assign phenological categories to all data
    cat("Assigning phenological categories...\n")
    plot_estimates$sampling_season <- assign_pheno_date_vectorized(plot_estimates$dates, pheno_categories_long)
    
    # Filter to calibration period
    for_stats <- plot_estimates %>%
      filter(time_period == "20130601_20180101") %>%
      filter(!is.na(Mean) & !is.na(sampling_season))
    
    cat("Rows for seasonality analysis:", nrow(for_stats), "\n")
    
    if(nrow(for_stats) == 0) {
      cat("No data for seasonality analysis in this chunk, skipping...\n")
      next
    }
    
    # Perform Tukey test for models in this chunk
    cat("Performing Tukey tests...\n")
    tukey_median_pheno <- for_stats %>%
      group_by(fcast_type, model_name, taxon, model_id) %>%
      filter(n() > 10) %>%  # Only test models with sufficient data
      summarize(tukey(x = sampling_season, y = Mean), .groups = "drop") %>%
      rename(sampling_season = x)
    
    if(nrow(tukey_median_pheno) > 0) {
      tukey_median_pheno_sig_chunk <- tukey_median_pheno %>%
        group_by(fcast_type, model_name, taxon, model_id) %>%
        mutate(significant_diff = ifelse(n_distinct(Letters_Tukey) > 1, T, F)) %>% 
        select(model_id, significant_diff) %>% 
        distinct()
      
      all_tukey_results[[chunk_i]] <- tukey_median_pheno_sig_chunk
      cat("Chunk", chunk_i, "processed:", nrow(tukey_median_pheno_sig_chunk), "models\n")
    } else {
      cat("No valid Tukey results for chunk", chunk_i, "\n")
    }
    
    # Clear memory
    rm(plot_estimates_env_cycl, plot_estimates_cycl_only, plot_estimates, for_stats, tukey_median_pheno)
    gc()
  }
  
  # Combine all chunk results
  cat("Combining results from all chunks...\n")
  if(length(all_tukey_results) > 0) {
    tukey_median_pheno_sig <- do.call(rbind, all_tukey_results)
  } else {
    cat("No results from any chunks, creating empty results...\n")
    tukey_median_pheno_sig <- data.frame(model_id = character(0), significant_diff = logical(0))
  }
  
  # Save results for future use
  cat("Saving seasonality analysis results...\n")
  saveRDS(tukey_median_pheno_sig, seasonality_file)
  
  cat("Models with significant seasonality differences:", sum(tukey_median_pheno_sig$significant_diff, na.rm = TRUE), "\n")
}

# Merge horizon data with seasonality analysis results
cat("Merging horizon data with seasonality analysis...\n")
cat("Horizon data rows:", nrow(fcast_horizon_data), "\n")
cat("Seasonality analysis rows:", nrow(tukey_median_pheno_sig), "\n")

if(!is.null(tukey_median_pheno_sig) && nrow(tukey_median_pheno_sig) > 0) {
  fcast_horizon_long_seas <- merge(fcast_horizon_data, tukey_median_pheno_sig, 
                                   by = c("model_id"), 
                                   all.x = TRUE) %>%
    mutate(significant_diff = ifelse(significant_diff, 
                                     "Microbes that vary \nwith plant phenophase",
                                     "Microbes that do not vary\n with plant phenophase"))
} else {
  cat("No seasonality analysis results available. Creating empty results...\n")
  fcast_horizon_long_seas <- fcast_horizon_data %>%
    mutate(significant_diff = "Microbes that do not vary\n with plant phenophase")
}

cat("Merged data rows before filtering:", nrow(fcast_horizon_long_seas), "\n")
cat("Non-NA significant_diff:", sum(!is.na(fcast_horizon_long_seas$significant_diff)), "\n")
cat("Non-NA months_since_obs:", sum(!is.na(fcast_horizon_long_seas$months_since_obs)), "\n")

fcast_horizon_long_seas <- fcast_horizon_long_seas %>%
  filter(!is.na(significant_diff) & !is.na(months_since_obs))

cat("Merged data rows after filtering:", nrow(fcast_horizon_long_seas), "\n")

# Handle pretty_group column (it should already be in the horizon data)
if("pretty_group.y" %in% names(fcast_horizon_long_seas)) {
  fcast_horizon_long_seas$pretty_group <- ifelse(is.na(fcast_horizon_long_seas$pretty_group.y), 
                                                 fcast_horizon_long_seas$pretty_group.x,
                                                 fcast_horizon_long_seas$pretty_group.y)
}

# Debug the final data
cat("Final data for plotting:\n")
cat("Rows:", nrow(fcast_horizon_long_seas), "\n")
cat("Unique significant_diff values:\n")
print(unique(fcast_horizon_long_seas$significant_diff))
cat("Months_since_obs summary:\n")
print(summary(fcast_horizon_long_seas$months_since_obs))
cat("Filtering for cycl_only models:\n")
cycl_only_data <- fcast_horizon_long_seas %>% filter(model_name %in% c("cycl_only"))
cat("Cycl_only data rows:", nrow(cycl_only_data), "\n")

# Filter out infinite values and create a more informative plot
plot_data <- fcast_horizon_long_seas %>% 
  filter(model_name %in% c("cycl_only")) %>%
  filter(is.finite(months_since_obs))

cat("Plot data rows (finite values only):", nrow(plot_data), "\n")

if (nrow(plot_data) > 0) {
  # Create the main figure
  fig3g <- ggplot(plot_data,
                  aes(x = significant_diff, y = months_since_obs)) +
    geom_violin(draw_quantiles = c(.5), alpha = .5, show.legend = F) +
    geom_point(aes(color = pretty_group),
               size = 3, alpha = .4,
               position = position_jitter(height = .2, width = .2),
               show.legend = T) +
    coord_flip() + 
    theme_bw(base_size = 16) +
    ylab("Forecast horizon (months of skilled predictions)") + 
    xlab(NULL) +
    labs(color = "Kingdom") +
    theme(legend.box = "horizontal", legend.position = "top")
} else {
  # Create a message plot if no data
  fig3g <- ggplot() + 
    annotate("text", x = 0.5, y = 0.5, 
             label = "No finite forecast horizon data available\nfor cycl_only models", 
             size = 6) +
    theme_void() +
    labs(title = "Forecast Horizon by Seasonality")
}

# Add statistical test if there are enough observations
if (nrow(fcast_horizon_long_seas) > 0) {
  tryCatch({
    stat_pvalue <- fcast_horizon_long_seas %>%
      group_by(model_name) %>%
      rstatix::t_test(months_since_obs ~ significant_diff, detailed = T) %>%
      rstatix::add_y_position(step.increase = .2) %>%
      mutate(y.position = seq(min(y.position), max(y.position), length.out = n()))
    
    fig3g <- fig3g + 
      ggpubr::stat_pvalue_manual(stat_pvalue,
                                 label = "p = {p}",
                                 bracket.nudge.y = -.1,
                                 size = 4, hide.ns = T)
  }, error = function(e) {
    cat("Statistical test could not be performed:", e$message, "\n")
  })
}

# Create final figure
fig3g <- ggarrange(fig3g, labels = "G")

# Save the figure
png(here("figures", "fcast_horizon_by_seasonality_chunked.png"), width = 1200, height = 400)
print(fig3g)
dev.off()

cat("Figure saved to figures/fcast_horizon_by_seasonality_chunked.png\n")
cat("Script completed successfully!\n")

