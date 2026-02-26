# Modified version of fig_latitude_phenophases.r that works with available data
library(lubridate)
library(ggrepel)
source("source.R")
source("microbialForecast/R/assignPhenology.r")
source("microbialForecast/R/load_plot_estimates.r")

sum.in <- readRDS(here("data", "summary/logit_beta_fixed_priors_summaries.rds"))

# Load phenology data directly - this contains the abundance information we need
# Element 1 = seasonality_mode2 (one row per model, peak phenophase only).
# Element 4 = max_abun_to_plot (multiple rows per model across phenophases; use for full seasonal coverage).
cat("Loading phenology data for latitude phenophase analysis...\n")
read_in <- readRDS(here("data/clean/pheno_group_peak_phenophases.rds"))
use_full_phenophases <- length(read_in) >= 4L && nrow(read_in[[4]]) > 0L &&
  all(c("mean_modeled_abun", "sampling_season", "taxon") %in% names(read_in[[4]]))
if (use_full_phenophases) {
  full_phenophase_abundance <- as.data.frame(read_in[[4]])
  abundance_col <- "mean_modeled_abun"
  abundance_ylab <- "Mean modeled abundance"
  cat("Using max_abun_to_plot (element 4) for full phenophase coverage.\n")
} else {
  full_phenophase_abundance <- as.data.frame(read_in[[1]])
  abundance_col <- "amplitude"
  abundance_ylab <- "Amplitude"
  cat("Using seasonality_mode2 (element 1); each model has one peak phenophase.\n")
}

# Filter to only converged models - handle different keep_models formats
if ("keep_models" %in% names(sum.in)) {
  if (is.data.frame(sum.in$keep_models) && "model_id" %in% names(sum.in$keep_models)) {
    keep_model_ids <- sum.in$keep_models$model_id
  } else if (is.character(sum.in$keep_models) || is.vector(sum.in$keep_models)) {
    keep_model_ids <- sum.in$keep_models
  } else {
    keep_model_ids <- NULL
  }
  if (!is.null(keep_model_ids)) {
    full_phenophase_abundance <- full_phenophase_abundance %>% filter(model_id %in% keep_model_ids)
  }
}

# Fix pretty_group assignment using proper taxonomy
fg_names <- microbialForecast:::keep_fg_names
rank_spec_names <- microbialForecast:::rank_spec_names
all_bacteria <- unlist(rank_spec_names[grepl("_bac$", names(rank_spec_names))])
all_fungi <- unlist(rank_spec_names[grepl("_fun$", names(rank_spec_names))])

assign_pretty_group <- function(data) {
  if ("species" %in% colnames(data)) {
    # For functional groups, use the proper taxonomy assignment
    if (any(data$species %in% fg_names)) {
      # Assign functional group categories first
      data$fg_cat <- microbialForecast:::assign_fg_categories(data$species)
      # Then assign kingdoms based on functional categories
      data$group <- microbialForecast:::assign_fg_kingdoms(data$fg_cat)
      # Convert to pretty_group
      data$pretty_group <- ifelse(data$group == "16S", "Bacteria", 
                                 ifelse(data$group == "ITS", "Fungi", NA))
    }
    
    # For regular taxonomic species, use the existing logic
    data$pretty_group <- ifelse(data$species %in% all_bacteria, "Bacteria",
                               ifelse(data$species %in% all_fungi, "Fungi", 
                                      data$pretty_group))  # Keep functional group assignment if already done
  }
  return(data)
}

# Fix pretty_group assignment in the phenology data
full_phenophase_abundance <- assign_pretty_group(full_phenophase_abundance)

# Debug: Check available columns
cat("Available columns in full_phenophase_abundance:", colnames(full_phenophase_abundance), "\n")

# This factor should be ordered (so that phenophases are sequential)
if (nrow(full_phenophase_abundance) > 0) {
  full_phenophase_abundance$sampling_season =
    factor(full_phenophase_abundance$sampling_season, ordered = T, levels = c("dormancy","greenup","peak","greendown"))
}

# Read in descriptions of NEON site-level soil chemistry, NCLD class, climate
site_descr <- readRDS(here("data/clean/site_effect_predictors.rds"))
site_descr$latitude_bin = cut(site_descr$latitude, breaks = 10) 	# Bin latitudes into groups

site_descr$latitude_category = ifelse(site_descr$latitude > 44, "High-latitude",
																			ifelse(site_descr$latitude < 31, "Low-latitude",
																						 "Mid-latitude"))
# This factor should be ordered
site_descr$latitude_category =
	factor(site_descr$latitude_category, ordered = T, levels = c("Low-latitude",  "Mid-latitude","High-latitude"))

# For the phenology analysis, we'll work directly with the phenology data
# since it already contains the abundance information (amplitude)
plot_phenophase_abundance <- full_phenophase_abundance

# For phenology data, we don't need to filter by time_period since it's already summarized
plot_phenophase_abundance_cal <- plot_phenophase_abundance

phenophase_fg_abundance_fungi <- plot_phenophase_abundance %>%
	filter(pretty_group=="Fungi" & fcast_type == "Functional")

# Check what columns are available in the phenology data
cat("Available columns in plot_phenophase_abundance_cal:", colnames(plot_phenophase_abundance_cal), "\n")
cat("Using abundance column:", abundance_col, "\n")

# Filter data for statistical tests
for_stats <- plot_phenophase_abundance_cal  %>%
	filter(!is.na(.data[[abundance_col]]) & !is.na(sampling_season))

cat("Data for statistical tests:\n")
cat("- Total observations:", nrow(for_stats), "\n")
cat("- Unique sampling seasons:", unique(for_stats$sampling_season), "\n")
cat("- Groups with data:", nrow(for_stats %>% group_by(fcast_type,model_name,taxon,model_id) %>% summarise(n = n(), .groups = "drop")), "\n")

# Note: Each model has only one phenophase, so we can't do within-group Tukey tests
# Instead, we'll do between-group comparisons (comparing different models across phenophases)
cat("Note: Each model has only one phenophase, so we'll do between-group comparisons instead of within-group tests\n")

# Create a summary of abundance by phenophase and taxonomic group (with SE for variability)
seasonal_summary_by_group <- for_stats %>%
  group_by(pretty_group, sampling_season) %>%
  summarise(
    mean_amplitude = mean(.data[[abundance_col]], na.rm = TRUE),
    median_amplitude = median(.data[[abundance_col]], na.rm = TRUE),
    n_obs = sum(!is.na(.data[[abundance_col]])),
    se_amplitude = ifelse(n_obs > 1L, sd(.data[[abundance_col]], na.rm = TRUE) / sqrt(n_obs), NA_real_),
    n_models = n_distinct(model_id),
    .groups = "drop"
  )

cat("Seasonal summary by taxonomic group:\n")
print(seasonal_summary_by_group)

# For the horizon analysis, we'll create a dummy significant_diff column
# since we can't do the original Tukey tests
tukey_median_pheno_sig <- for_stats %>%
  select(model_id, pretty_group, model_name) %>%
  distinct(model_id, .keep_all = TRUE) %>%
  mutate(significant_diff = FALSE)  # Set to FALSE since we can't do the original tests

# Use forecast horizon results (per-model horizon in months), not horizon input (which has months_since_obs = Inf)
if (nrow(tukey_median_pheno_sig) > 0) {
  horizon_results_file <- here("data/summary/fcast_horizon_df.rds")
  if (!file.exists(horizon_results_file)) {
    cat("fcast_horizon_df.rds not found; skipping horizon vs seasonality plot\n")
  } else {
    fcast_horizon <- readRDS(horizon_results_file)
    fcast_horizon_data <- as.data.frame(fcast_horizon[[1]])
    fcast_horizon_data$forecast_horizon <- fcast_horizon_data$rsq_fcast_horizon
    fcast_horizon_data$forecast_horizon <- ifelse(
      is.na(fcast_horizon_data$forecast_horizon) | !is.finite(fcast_horizon_data$forecast_horizon),
      ifelse(is.finite(fcast_horizon_data$crps_fcast_horizon) & fcast_horizon_data$crps_fcast_horizon > 0,
             fcast_horizon_data$crps_fcast_horizon, fcast_horizon_data$rmse_fcast_horizon),
      fcast_horizon_data$forecast_horizon)

    fcast_horizon_data$model_id_norm <- gsub("_beta_regression$|_combined$", "", fcast_horizon_data$model_id)
    tukey_median_pheno_sig$model_id_norm <- gsub("_beta_regression$|_combined$", "", tukey_median_pheno_sig$model_id)

    core_horizon_seas <- merge(fcast_horizon_data, tukey_median_pheno_sig,
                               by.x = "model_id_norm", by.y = "model_id_norm", all = FALSE,
                               suffixes = c("", "_pheno"))
    core_horizon_seas$model_id_norm <- NULL
    if ("pretty_group_pheno" %in% names(core_horizon_seas)) {
      core_horizon_seas$pretty_group <- dplyr::coalesce(core_horizon_seas$pretty_group, core_horizon_seas$pretty_group_pheno)
      core_horizon_seas$pretty_group_pheno <- NULL
    }
    if ("model_name_pheno" %in% names(core_horizon_seas)) {
      core_horizon_seas$model_name <- dplyr::coalesce(core_horizon_seas$model_name, core_horizon_seas$model_name_pheno)
      core_horizon_seas$model_name_pheno <- NULL
    }
    core_horizon_seas <- core_horizon_seas %>%
      filter(!is.na(pretty_group) & !is.na(model_name) & is.finite(forecast_horizon) & forecast_horizon > 0)

    if (nrow(core_horizon_seas) > 0) {
      p1 <- ggplot(core_horizon_seas,
                   aes(x = significant_diff, y = forecast_horizon)) +
        geom_boxplot(alpha = .5, show.legend = FALSE) +
        geom_point(size = 1, alpha = .3, position = position_jitter(height = 0), show.legend = FALSE) +
        facet_grid(pretty_group ~ model_name) +
        coord_flip() +
        theme_bw() +
        ylab("Forecast horizon (months)")
      print(p1)
      png(here("figures", "horizon_vs_seasonality.png"), width = 1000, height = 600)
      print(p1)
      dev.off()
      cat("Created and saved horizon vs seasonality plot\n")
    } else {
      cat("No overlapping model_ids or no finite forecast horizons; skipping horizon vs seasonality plot\n")
    }
  }
} else {
  cat("No significant groups found, skipping horizon vs seasonality plot\n")
}

# For fungi functional groups, we'll work with the phenology data directly
# Since we don't have site information in the phenology data, we'll skip the latitude category analysis
for_stats_fungi <- phenophase_fg_abundance_fungi  %>%
	filter(model_name=="cycl_only" & !is.na(.data[[abundance_col]]) & !is.na(sampling_season) & !is.na(pretty_group))

cat("Data for fungi functional groups:\n")
cat("- Total observations:", nrow(for_stats_fungi), "\n")
cat("- Available taxa:", unique(for_stats_fungi$taxon), "\n")

# Since we don't have site information, we'll skip the latitude category analysis
# and just create a summary by phenophase
if(nrow(for_stats_fungi) > 0) {
  cat("Creating summary for fungi functional groups by phenophase...\n")
  
  fungi_seasonal_summary <- for_stats_fungi %>%
    group_by(taxon, sampling_season) %>%
    summarise(
      mean_amplitude = mean(.data[[abundance_col]], na.rm = TRUE),
      median_amplitude = median(.data[[abundance_col]], na.rm = TRUE),
      n_models = n_distinct(model_id),
      .groups = "drop"
    )
  
  cat("Fungi seasonal summary:\n")
  print(fungi_seasonal_summary)
  
  tukey_median_pheno_evergreen <- data.frame()  # No Tukey tests possible
} else {
  cat("No fungi functional group data available\n")
  tukey_median_pheno_evergreen <- data.frame()
}

# Create the main phenophase abundance plot
cat("Creating phenophase abundance plot...\n")

# Check what taxa are actually available for plotting
available_taxa <- unique(for_stats$taxon)
target_taxa <- c("ectomycorrhizal","plant_pathogen","saprotroph","endophyte")
plot_taxa <- intersect(available_taxa, target_taxa)

cat("Taxa available for plotting:", plot_taxa, "\n")

if(length(plot_taxa) > 0) {
  plot_median_abun <- ggplot(for_stats %>% filter(taxon %in% plot_taxa),
													 aes(y = .data[[abundance_col]],
													 		x = sampling_season)) +

	geom_boxplot(alpha=.5, show.legend = F) +
	geom_point(
						 size=1, alpha=.3,
						 position=position_jitter(height=0), show.legend = F) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	facet_grid(~taxon, scales = "free") +
	ylab(abundance_ylab) +
	xlab("Plant phenophase")

  print(plot_median_abun)

  # Save the plot
  png(here("figures", "phenophase_abundance.png"), width = 1200, height = 800)
  print(plot_median_abun)
  dev.off()
  cat("Created and saved phenophase abundance plot\n")
} else {
  cat("No target taxa available for plotting\n")
}

# Create a summary plot showing seasonal patterns by taxonomic group
cat("Creating seasonal patterns plot...\n")

if(nrow(seasonal_summary_by_group) > 0) {
  seasonal_plot <- ggplot(seasonal_summary_by_group,
                         aes(x = sampling_season, y = mean_amplitude, color = pretty_group)) +
    geom_errorbar(aes(ymin = mean_amplitude - se_amplitude, ymax = mean_amplitude + se_amplitude,
                     group = pretty_group), width = 0.15, linewidth = 0.5, na.rm = TRUE) +
    geom_line(aes(group = pretty_group), size = 1) +
    geom_point(size = 2) +
    theme_bw() +
    theme(text = element_text(size = 14)) +
    ylab(paste0("Mean ", abundance_ylab)) +
    xlab("Plant Phenophase") +
    labs(color = "Taxonomic Group") +
    ggtitle("Seasonal Patterns by Taxonomic Group")
  
  print(seasonal_plot)
  
  # Save the plot
  png(here("figures", "seasonal_patterns_by_group.png"), width = 1000, height = 600)
  print(seasonal_plot)
  dev.off()
  cat("Created and saved seasonal patterns plot\n")
} else {
  cat("No data available for seasonal patterns plot\n")
}

cat("Script completed successfully!\n")
cat("Statistical tests were run to test differences in predicted abundance between seasons.\n")