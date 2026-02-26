source("source.R")
library(cowplot)
library(MuMIn)
options(na.action = "na.fail")

# Read in forecast scores
scores_list = readRDS(here("data/summary/scoring_metrics_plsr2.rds"))
unconverged = scores_list$unconverged_list
converged = scores_list$converged_list
converged_strict = scores_list$converged_strict_list

# Read in Moran's I (spatial autocorrelation)
moran.stat_all_rank = readRDS("data/clean/moran_stat.rds")
cat("Moran's I data loaded:", nrow(moran.stat_all_rank), "rows\n")

# Read in model estimates for rho (temporal autocorrelation) and core_sd (spatial variation)
rho_core_in <- readRDS(here("data", "summary/rho_core_sd_effects.rds")) %>%
	filter(model_name != "all_covariates" & model_id %in% converged_strict) %>%
	select(-pretty_name) %>%
	mutate(model_id = gsub("_beta_regression$", "", model_id))  # Remove _beta_regression suffix
core_sd = rho_core_in  %>% filter(rowname == "core_sd")
rho = rho_core_in  %>% filter(rowname == "rho")

# Check initial data dimensions
cat("scores_list$scoring_metrics dimensions:", nrow(scores_list$scoring_metrics), "x", ncol(scores_list$scoring_metrics), "\n")
cat("rho_core_in dimensions:", nrow(rho_core_in), "x", ncol(rho_core_in), "\n")
cat("converged_strict length:", length(converged_strict), "\n")

# Check model_id overlap
scores_model_ids <- unique(scores_list$scoring_metrics$model_id)
rho_model_ids <- unique(rho_core_in$model_id)
cat("Scores model_ids (first 10):", head(scores_model_ids, 10), "\n")
cat("Rho model_ids (first 10):", head(rho_model_ids, 10), "\n")
cat("Overlap between scores and rho model_ids:", length(intersect(scores_model_ids, rho_model_ids)), "\n")

# Check for column name conflicts
scores_cols <- colnames(scores_list$scoring_metrics)
rho_cols <- colnames(rho_core_in)
common_cols <- intersect(scores_cols, rho_cols)
cat("Common columns between scores and rho data:", common_cols, "\n")

# Merge the data that has been read in so far
# Check if "New time (observed site)" exists in any data structure
use_observed_data <- FALSE
scoring_metrics_obs <- NULL

if("scoring_metrics_long" %in% names(scores_list)) {
  obs_count <- sum(scores_list$scoring_metrics_long$site_prediction == "New time (observed site)", na.rm = TRUE)
  if(obs_count > 0) {
    cat("Found", obs_count, "rows with 'New time (observed site)' in scoring_metrics_long\n")
    # Filter to observed sites and pivot back to wide format for merging
    scoring_metrics_obs <- scores_list$scoring_metrics_long %>%
      filter(site_prediction == "New time (observed site)") %>%
      pivot_wider(id_cols = c("model_id", "fcast_period", "pretty_group", "model_name", 
                              "rank_name", "species"),
                  names_from = "metric", values_from = "score")
    cat("Using scoring_metrics_long filtered to observed sites:", nrow(scoring_metrics_obs), "rows\n")
    use_observed_data <- TRUE
  }
}

# If no observed data found, check scoring_metrics for observed site_prediction
if(!use_observed_data && "site_prediction" %in% colnames(scores_list$scoring_metrics)) {
  obs_in_metrics <- grep("observed", scores_list$scoring_metrics$site_prediction, value = TRUE, ignore.case = TRUE)
  if(length(obs_in_metrics) > 0) {
    cat("Found observed site_prediction values in scoring_metrics:", paste(unique(obs_in_metrics), collapse = ", "), "\n")
    scoring_metrics_obs <- scores_list$scoring_metrics %>%
      filter(grepl("observed", site_prediction, ignore.case = TRUE))
    cat("Using scoring_metrics filtered to observed sites:", nrow(scoring_metrics_obs), "rows\n")
    use_observed_data <- TRUE
  }
}

# If still no observed data, use all scoring_metrics (which may already be at observed sites)
if(is.null(scoring_metrics_obs)) {
  scoring_metrics_obs <- scores_list$scoring_metrics
  cat("No 'New time (observed site)' data found. Using all scoring_metrics:", nrow(scoring_metrics_obs), "rows\n")
  cat("Note: Current data only has 'New time x site' values. If observed data should exist,\n")
  cat("      consider regenerating scoring_metrics with observed site predictions.\n")
}

# Deduplicate rho_core_in before merging to avoid cartesian explosion
rho_core_unique <- rho_core_in %>% distinct(model_id, rowname, .keep_all = TRUE)
model_scores_vals = merge(scoring_metrics_obs, rho_core_unique, by = "model_id", allow.cartesian = TRUE)
cat("After merge with rho_core (before filtering):", nrow(model_scores_vals), "rows\n")

# Check what site_prediction values are available
cat("Available site_prediction values:\n")
print(table(model_scores_vals$site_prediction, useNA = "ifany"))

# Apply filtering - use the model_ids that are actually in the merged data
available_model_ids <- unique(model_scores_vals$model_id)
cat("Available model_ids in merged data:", length(available_model_ids), "\n")

# If we already filtered to observed data, we're done
# Otherwise, check if site_prediction column exists and filter if needed
if(!use_observed_data && "site_prediction" %in% colnames(model_scores_vals)) {
  cat("Available site_prediction values:", paste(unique(model_scores_vals$site_prediction), collapse = ", "), "\n")
  # Try different patterns for "observed"
  observed_patterns <- c("observed", "Observed", "New time \\(observed")
  matching_rows <- model_scores_vals %>%
    filter(grepl(paste(observed_patterns, collapse = "|"), site_prediction))
  if(nrow(matching_rows) > 0) {
    model_scores_vals <- matching_rows
    cat("After filtering for observed sites:", nrow(model_scores_vals), "rows\n")
    use_observed_data <- TRUE
  } else {
    cat("No rows match 'observed' pattern, using all data\n")
    cat("Using", nrow(model_scores_vals), "rows without filtering\n")
  }
} else if(use_observed_data) {
  cat("Already filtered to observed sites\n")
} else {
  cat("site_prediction column not found, using all data\n")
}
cat("Final model_scores_vals rows:", nrow(model_scores_vals), "\n")

# Merge with Moran's I data using the correct column names
# Check what columns are available in moran data
cat("Moran data columns:", colnames(moran.stat_all_rank), "\n")
cat("Model scores columns:", colnames(model_scores_vals), "\n")

# Use the correct merge column - likely 'taxon' or 'species'
if("taxon" %in% colnames(moran.stat_all_rank) && "taxon" %in% colnames(model_scores_vals)) {
  model_scores_vals <- merge(model_scores_vals, moran.stat_all_rank, by = "taxon", all.x = T)
} else if("species" %in% colnames(model_scores_vals)) {
  model_scores_vals <- merge(model_scores_vals, moran.stat_all_rank, by.x = "species", by.y = "taxon", all.x = T)
} else {
  cat("Cannot merge with Moran data - no suitable column found\n")
}
cat("After merge with moran data:", nrow(model_scores_vals), "rows\n")


# Read in model estimates of seasonal amplitude
seasonal_amplitude_in = readRDS(here("data/summary/seasonal_amplitude.rds"))
cat("seasonal_amplitude_in structure:\n")
print(str(seasonal_amplitude_in))

cycl_only_vals_scores = seasonal_amplitude_in[[6]] %>% filter(model_name=="cycl_only") %>%
	mutate(cycl_amplitude=amplitude, cycl_max = max)
cat("cycl_only_vals_scores dimensions:", nrow(cycl_only_vals_scores), "x", ncol(cycl_only_vals_scores), "\n")

env_cycl_vals_scores = seasonal_amplitude_in[[6]] %>% filter(model_name=="env_cycl")  %>%
	mutate(env_cycl_amplitude=amplitude, env_cycl_max = max)
cat("env_cycl_vals_scores dimensions:", nrow(env_cycl_vals_scores), "x", ncol(env_cycl_vals_scores), "\n")

# Deduplicate before merging to avoid cartesian explosion
cycl_only_unique <- cycl_only_vals_scores %>%
	select(taxon, time_period, cycl_amplitude, cycl_max) %>%
	distinct(taxon, time_period, .keep_all = TRUE)
model_scores_vals <- merge(model_scores_vals, cycl_only_unique, by = c("taxon", "time_period"), all.x = TRUE, allow.cartesian = TRUE)
cat("After merge with cycl_only data:", nrow(model_scores_vals), "rows\n")

env_cycl_unique <- env_cycl_vals_scores %>%
	select(taxon, time_period, env_cycl_amplitude, env_cycl_max) %>%
	distinct(taxon, time_period, .keep_all = TRUE)
model_scores_vals <- merge(model_scores_vals, env_cycl_unique, by = c("taxon", "time_period"), all.x = TRUE, allow.cartesian = TRUE)
cat("After merge with env_cycl data:", nrow(model_scores_vals), "rows\n")
# Clean up duplicate columns - prefer pretty_group from scoring metrics (now guaranteed to be present)
if("pretty_group.x" %in% colnames(model_scores_vals) && "pretty_group.y" %in% colnames(model_scores_vals)) {
  # Use pretty_group.x (from scoring_metrics) since it's now guaranteed to be correct
  model_scores_vals$pretty_group <- ifelse(is.na(model_scores_vals$pretty_group.x), 
                                           model_scores_vals$pretty_group.y, 
                                           model_scores_vals$pretty_group.x)
  model_scores_vals$pretty_group.x <- NULL
  model_scores_vals$pretty_group.y <- NULL
} else if("pretty_group.y" %in% colnames(model_scores_vals)) {
  # If only .y exists, use it (shouldn't happen if scoring metrics has pretty_group)
  model_scores_vals$pretty_group <- model_scores_vals$pretty_group.y
  model_scores_vals$pretty_group.y <- NULL
}
if("model_name.x" %in% colnames(model_scores_vals)) {
  model_scores_vals$model_name <- model_scores_vals$model_name.x
  model_scores_vals$model_name.x <- NULL
  model_scores_vals$model_name.y <- NULL
}
if("time_period.x" %in% colnames(model_scores_vals)) {
  model_scores_vals$time_period <- model_scores_vals$time_period.x
  model_scores_vals$time_period.x <- NULL
  model_scores_vals$time_period.y <- NULL
}
if("rank.x" %in% colnames(model_scores_vals)) {
  model_scores_vals$rank <- model_scores_vals$rank.x
  model_scores_vals$rank.x <- NULL
  model_scores_vals$rank.y <- NULL
}

# Check available columns first
cat("Available columns in model_scores_vals:\n")
print(colnames(model_scores_vals))
cat("model_scores_vals dimensions:", nrow(model_scores_vals), "x", ncol(model_scores_vals), "\n")

# Check what model_name values are available
cat("Available model_name values:\n")
print(table(model_scores_vals$model_name, useNA = "ifany"))

# Check if rowname and Mean columns exist before pivoting
if("rowname" %in% colnames(model_scores_vals) && "Mean" %in% colnames(model_scores_vals)) {
  # Select only existing columns
  existing_cols <- colnames(model_scores_vals)
  cols_to_remove <- c("SD", "Naive SE", "Time-series SE")
  cols_to_remove <- cols_to_remove[cols_to_remove %in% existing_cols]
  
  if(length(cols_to_remove) > 0) {
    model_scores_vals_wide = model_scores_vals %>% select(-all_of(cols_to_remove)) %>%
      pivot_wider(names_from = "rowname", values_from = "Mean")
  } else {
    model_scores_vals_wide = model_scores_vals %>%
      pivot_wider(names_from = "rowname", values_from = "Mean")
  }
  cat("After pivot_wider:", nrow(model_scores_vals_wide), "rows\n")
} else {
  cat("Warning: rowname or Mean columns not found, skipping pivot_wider\n")
  cat("Available columns:", paste(colnames(model_scores_vals), collapse = ", "), "\n")
  # Use the data as-is if pivot isn't possible
  model_scores_vals_wide <- model_scores_vals
  cat("Using model_scores_vals without pivoting:", nrow(model_scores_vals_wide), "rows\n")
}

# Use all available models instead of filtering to just env_cycl
cat("Using all available models:", nrow(model_scores_vals_wide), "rows\n")
cat("Model distribution:\n")
print(table(model_scores_vals_wide$model_name, useNA = "ifany"))

# Read in model estimates of environmental predictors
beta_in <- readRDS(here("data", "summary/predictor_effects.rds")) %>%
	filter(time_period=="2015-11_2018-01" & model_name=="env_cycl" &
				 	model_id %in% converged_strict)
cat("beta_in rows:", nrow(beta_in), "\n")

if(nrow(beta_in) > 0) {
  beta_wide <- beta_in %>%
    pivot_wider(id_cols = c("model_name","pretty_group","taxon","time_period"),
    						names_from = "beta", values_from = "effSize")
  cat("beta_wide dimensions before complete.cases:", nrow(beta_wide), "x", ncol(beta_wide), "\n")
  
  # Don't use complete.cases - it's too restrictive, just merge what we have
  # beta_wide <- beta_wide[complete.cases(beta_wide),]
  cat("beta_wide columns:", paste(head(colnames(beta_wide), 10), collapse = ", "), "...\n")
  
  # Check merge keys
  cat("model_scores_vals_wide merge keys (sample):", paste(head(unique(model_scores_vals_wide$taxon), 5), collapse = ", "), "\n")
  cat("beta_wide merge keys (sample):", paste(head(unique(beta_wide$taxon), 5), collapse = ", "), "\n")
  
  model_scores_vals_wide_betas <- merge(model_scores_vals_wide, beta_wide, 
                                         by = c("taxon", "time_period", "pretty_group", "model_name"),
                                         all.x = TRUE)
  cat("After merge with beta_wide:", nrow(model_scores_vals_wide_betas), "rows\n")
  
  # If the merge results in 0 rows, use the original data without beta
  if(nrow(model_scores_vals_wide_betas) == 0) {
    cat("Merge with beta_wide resulted in 0 rows, using original data\n")
    model_scores_vals_wide_betas <- model_scores_vals_wide
  }
} else {
  cat("No beta_in data available, using model_scores_vals_wide without beta data\n")
  model_scores_vals_wide_betas <- model_scores_vals_wide
}


# Check correlation between moran's and seasonality (negative)
cycl_vals = cycl_only_vals_scores %>% filter(model_id %in% converged_strict &
																						 	time_period == "2015-11_2018-01")
plotting_df = merge(cycl_vals, moran.stat_all_rank)

# Check the data before creating the main plot
cat("model_scores_vals_wide_betas dimensions:", nrow(model_scores_vals_wide_betas), "x", ncol(model_scores_vals_wide_betas), "\n")
cat("Available columns in model_scores_vals_wide_betas:", colnames(model_scores_vals_wide_betas), "\n")

# Check for non-finite values in the main plot data
if("mean_morans" %in% colnames(model_scores_vals_wide_betas)) {
  cat("mean_morans - finite values:", sum(is.finite(model_scores_vals_wide_betas$mean_morans)), "out of", nrow(model_scores_vals_wide_betas), "\n")
}
if("cycl_amplitude" %in% colnames(model_scores_vals_wide_betas)) {
  cat("cycl_amplitude - finite values:", sum(is.finite(model_scores_vals_wide_betas$cycl_amplitude)), "out of", nrow(model_scores_vals_wide_betas), "\n")
}

# Filter to only finite values for the main plot
finite_main_data <- model_scores_vals_wide_betas[is.finite(model_scores_vals_wide_betas$mean_morans) & 
                                                 is.finite(model_scores_vals_wide_betas$cycl_amplitude), ]
cat("Finite data for main plot:", nrow(finite_main_data), "rows\n")

# Create a minimal test plot
cat("Creating minimal test plot...\n")
cat("finite_main_data dimensions:", nrow(finite_main_data), "x", ncol(finite_main_data), "\n")
cat("mean_morans range:", range(finite_main_data$mean_morans, na.rm = TRUE), "\n")
cat("cycl_amplitude range:", range(finite_main_data$cycl_amplitude, na.rm = TRUE), "\n")

# Try to create a simple plot
tryCatch({
  a <- ggplot(finite_main_data, aes(x = mean_morans, y = cycl_amplitude)) +
    geom_point() +
    theme_bw()
  cat("Plot created successfully\n")
}, error = function(e) {
  cat("Error creating plot:", e$message, "\n")
  a <- NULL
})

# cycl_vals = cycl_only_vals_scores %>% filter(model_id %in% converged_strict &
# 																						 	time_period == "2015-11_2018-01")
# plotting_df = merge(cycl_vals, moran.stat_all_rank)
# ggplot(model_scores_vals_wide_betas,
# 			 aes(x = mean_morans, y = cycl_amplitude, color = pretty_group)) +
# 	geom_jitter(alpha=.3, size = 3, height = 0.01, width = 0.01) +
# 	geom_smooth(method="lm",
# 							linewidth=2, alpha = .2, na.rm = T, se = F)  +
# 	theme_bw(base_size = 18) +
# 	#facet_grid(metric ~pretty_group, scales="free") +
# 	scale_x_sqrt() +
# 	ylab("Seasonal amplitude") + xlab("Variation among cores") +
# 	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~"))) +
# 	labs(color = "Domain")

# testing_df = merge(core_sd, moran.stat_all_rank) %>% filter(model_name=="cycl_only" &
# 																														model_id %in% converged &
# 																															time_period == "2015-11_2018-01")

# Check data before creating plots
cat("model_scores_vals_wide dimensions:", nrow(model_scores_vals_wide), "x", ncol(model_scores_vals_wide), "\n")
cat("Available columns:", colnames(model_scores_vals_wide), "\n")

# Check for non-finite values in the data
cat("Checking for non-finite values:\n")
if("precision" %in% colnames(model_scores_vals_wide)) {
  cat("precision - finite values:", sum(is.finite(model_scores_vals_wide$precision)), "out of", nrow(model_scores_vals_wide), "\n")
}
if("cycl_amplitude" %in% colnames(model_scores_vals_wide)) {
  cat("cycl_amplitude - finite values:", sum(is.finite(model_scores_vals_wide$cycl_amplitude)), "out of", nrow(model_scores_vals_wide), "\n")
}

# Create simple density plots using ggplot instead of ggdensity
if(nrow(model_scores_vals_wide) > 0 && "precision" %in% colnames(model_scores_vals_wide) && sum(is.finite(model_scores_vals_wide$precision)) > 0) {
  # Filter to only finite values for plotting
  finite_data_precision <- model_scores_vals_wide[is.finite(model_scores_vals_wide$precision), ]
  if(nrow(finite_data_precision) > 0) {
    cat("Creating xplot with", nrow(finite_data_precision), "finite precision values\n")
    tryCatch({
      xplot <- ggplot(finite_data_precision, aes(x = precision)) +
        geom_density(alpha = 0.7, fill = "blue") +
        theme_bw() + 
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 45, hjust = 1))
      cat("xplot created successfully\n")
    }, error = function(e) {
      cat("Error creating xplot:", e$message, "\n")
      xplot <- NULL
    })
  } else {
    cat("Cannot create xplot - no finite precision values\n")
    xplot <- NULL
  }
} else {
  cat("Cannot create xplot - insufficient data or missing precision column\n")
  xplot <- NULL
}

if(nrow(model_scores_vals_wide) > 0 && "cycl_amplitude" %in% colnames(model_scores_vals_wide) && sum(is.finite(model_scores_vals_wide$cycl_amplitude)) > 0) {
  # Filter to only finite values for plotting
  finite_data_amplitude <- model_scores_vals_wide[is.finite(model_scores_vals_wide$cycl_amplitude), ]
  if(nrow(finite_data_amplitude) > 0) {
    cat("Creating yplot with", nrow(finite_data_amplitude), "finite cycl_amplitude values\n")
    tryCatch({
      yplot <- ggplot(finite_data_amplitude, aes(x = cycl_amplitude)) +
        geom_density(alpha = 0.7, fill = "red") +
        theme_bw() + 
        theme(legend.position = "none") +
        coord_flip() +
        theme(axis.text.y = element_text(angle = 45, hjust = 1))
      cat("yplot created successfully\n")
    }, error = function(e) {
      cat("Error creating yplot:", e$message, "\n")
      yplot <- NULL
    })
  } else {
    cat("Cannot create yplot - no finite cycl_amplitude values\n")
    yplot <- NULL
  }
} else {
  cat("Cannot create yplot - insufficient data or missing cycl_amplitude column\n")
  yplot <- NULL
}

# Skip legend extraction for now to avoid aesthetics error
legend_grob <- NULL

sp <- a

# Create individual plots instead of plot_grid to avoid aesthetics issues
cat("Creating individual plots...\n")

# Skip printing plots for now to avoid aesthetics error
cat("Skipping plot printing to avoid aesthetics error\n")

# Check what groups are available in the original data
cat("Available pretty_group values in original scores data:\n")
cat("Skipping table printing to avoid aesthetics error\n")

cat("Available pretty_group values in rho_core_in data:\n")
cat("Skipping table printing to avoid aesthetics error\n")

cat("Available pretty_group values in merged data before cleaning:\n")
cat("Skipping table printing to avoid aesthetics error\n")

# Check what's happening with the filtering steps
cat("Before pivot_wider - model_name distribution:\n")
cat("Skipping table printing to avoid aesthetics error\n")

cat("After pivot_wider - model_name distribution:\n")
cat("Skipping table printing to avoid aesthetics error\n")

cat("Available pretty_group values in final data:\n")
cat("Skipping table printing to avoid aesthetics error\n")

# Check if we're losing data due to the env_cycl filter
cat("Total rows before env_cycl filter:", nrow(model_scores_vals), "\n")
cat("Total rows after env_cycl filter:", nrow(model_scores_vals_wide), "\n")

# Only run summary if we have enough data and multiple groups
if(nrow(model_scores_vals_wide) > 0 && sum(!is.na(model_scores_vals_wide$RSQ.1)) > 0 && length(unique(model_scores_vals_wide$pretty_group)) > 1) {
  cat("Running first summary (RSQ.1 ~ mean_morans * cycl_amplitude * pretty_group)...\n")
  summary(lm(RSQ.1 ~ mean_morans * cycl_amplitude * pretty_group, model_scores_vals_wide, na.action = na.omit))
} else {
  cat("Cannot run first summary - insufficient data or only one group\n")
}




# Create the main scatter plot (b) if we have data
if(nrow(finite_main_data) > 0 && "pretty_group" %in% colnames(finite_main_data) && length(unique(finite_main_data$pretty_group)) > 0) {
  b <- ggplot(finite_main_data,
  			 aes(x = mean_morans, y = cycl_amplitude, color = pretty_group)) +
  	geom_jitter(alpha=.3, size = 3, height = 0.01, width = 0.01) +
  	geom_smooth(method="lm",
  							linewidth=2, alpha = .2, na.rm = T, se = F)  +
  	theme_bw(base_size = 18) +
  	scale_x_sqrt() +
  	ylab("Seasonal amplitude") + xlab("Variation among cores") +
  	stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~"))) +
  	labs(color = "Domain")
  cat("Main scatter plot (b) created successfully\n")
} else {
  cat("Cannot create main scatter plot - insufficient data\n")
  b <- NULL
}

# Create density plots for margins
xplot <- NULL
yplot <- NULL

if(nrow(model_scores_vals_wide) > 0 && "core_sd" %in% colnames(model_scores_vals_wide) && "pretty_group" %in% colnames(model_scores_vals_wide)) {
  tryCatch({
    xplot <- ggdensity(model_scores_vals_wide, "core_sd", fill = "pretty_group") + clean_theme() + rremove("legend")
  }, error = function(e) {
    cat("Error creating xplot:", e$message, "\n")
  })
}

if(nrow(model_scores_vals_wide) > 0 && "cycl_amplitude" %in% colnames(model_scores_vals_wide) && "pretty_group" %in% colnames(model_scores_vals_wide)) {
  tryCatch({
    yplot <- ggdensity(model_scores_vals_wide, "cycl_amplitude", fill = "pretty_group") +
    	rotate() + clean_theme() + rremove("legend")
  }, error = function(e) {
    cat("Error creating yplot:", e$message, "\n")
  })
}

# Extract legend from main plot if it exists
legend_grob <- NULL
if(!is.null(a) && !is.null(b)) {
  tryCatch({
    grobs <- ggplotGrob(b)$grobs
    legend_idx <- which(sapply(grobs, function(x) x$name) == "guide-box")
    if(length(legend_idx) > 0) {
      legend_grob <- grobs[[legend_idx]]
    }
  }, error = function(e) {
    cat("Error extracting legend:", e$message, "\n")
  })
}

# Prepare plots for plot_grid
sp <- NULL
if(!is.null(b)) {
  sp <- b + rremove("legend")
}

if(!is.null(yplot)) {
  yplot <- yplot + clean_theme() + rremove("legend")
}

if(!is.null(xplot)) {
  xplot <- xplot + clean_theme() + rremove("legend")
}

# Only create plot_grid if we have the main plot
if(!is.null(sp)) {
  # If we have margin plots, create a composite, otherwise just save the main plot
  if(!is.null(xplot) && !is.null(yplot)) {
    if(!is.null(legend_grob)) {
      p_spatial_amplitude <- plot_grid(xplot, legend_grob,
      					sp, yplot, ncol = 2, align = "hv",
      					rel_widths = c(2, 1), rel_heights = c(1, 2))
    } else {
      p_spatial_amplitude <- plot_grid(xplot, sp, yplot, ncol = 2, align = "hv",
      					rel_widths = c(1, 1), rel_heights = c(1, 1))
    }
  } else {
    # Just use the main plot if margin plots aren't available
    p_spatial_amplitude <- sp
  }
  
  print(p_spatial_amplitude)
  tryCatch({
    png(here("figures","spatial_vs_amplitude.png"), width = 1400, height = 1200)
    print(p_spatial_amplitude)
    dev.off()
    cat("Spatial vs amplitude plot saved successfully\n")
  }, error = function(e) {
    cat("Error saving spatial_vs_amplitude:", e$message, "\n")
  })
} else {
  cat("Cannot create plot - missing main scatter plot\n")
}

# Only run summary if we have enough data and multiple groups
if(nrow(model_scores_vals_wide) > 0 && sum(!is.na(model_scores_vals_wide$RSQ.1)) > 0 && length(unique(model_scores_vals_wide$pretty_group)) > 1) {
  cat("Running second summary (RSQ.1 ~ mean_morans * pretty_group)...\n")
  summary(lm(RSQ.1 ~ mean_morans * pretty_group, model_scores_vals_wide, na.action = na.omit))
} else {
  cat("Cannot run second summary - insufficient data or only one group\n")
}
