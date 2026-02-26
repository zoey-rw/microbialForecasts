# Fig Patchiness Predictability by Group (Fungi vs Bacteria)
# This script analyzes how trends vary between fungi and bacteria

source("source.R")
library(cowplot)
library(MuMIn)
library(sjPlot)
library(ggplot2)
options(na.action = "na.omit")

cat("=== ANALYZING PATCHINESS PREDICTABILITY BY GROUP ===\n")

# Step 1: Read in forecast scores
cat("Step 1: Reading forecast scores...\n")
scores_list = readRDS(here("data/summary/scoring_metrics_plsr2.rds"))
converged = unique(scores_list$scoring_metrics$model_id)
cat("Found", length(converged), "models with scoring data\n")

# Step 2: Read in Moran's I (spatial autocorrelation)
cat("Step 2: Reading Moran's I data...\n")
moran.stat_all_rank = readRDS("data/clean/moran_stat.rds")
cat("Moran's I data loaded\n")

# Step 3: Read in model estimates for rho (temporal autocorrelation) and precision (core heterogeneity)
cat("Step 3: Reading rho and precision data...\n")
rho_core_in <- readRDS(here("data", "summary/rho_core_sd_effects.rds")) %>%
  filter(time_period=="20130601_20180101") %>%
  filter(model_name != "all_covariates") %>%
  select(model_id, taxon, rowname, Mean) %>%
  pivot_wider(names_from = "rowname", values_from = "Mean", values_fn = mean) %>%
  mutate(model_id_base = gsub("_(combined|beta_regression)$", "", model_id))

# Filter to only include models that have scoring data
converged_base = gsub('_(combined|beta_regression)$', '', converged)
rho_core_in <- rho_core_in %>% filter(model_id_base %in% converged_base)
cat("Rho and precision data loaded:", nrow(rho_core_in), "rows\n")

# Step 4: Get scores data
cat("Step 4: Filtering scores data...\n")
scores_df = scores_list$scoring_metrics %>%
  filter(model_id %in% converged)

# Add group information based on proper taxonomy structure from microbialForecast package
cat("Step 4a: Adding group information based on proper taxonomy...\n")

# Get the correct taxonomy lists
fg_names <- microbialForecast:::keep_fg_names
rank_spec_names <- microbialForecast:::rank_spec_names

# Create vectors of all bacteria and fungi species
all_bacteria <- unlist(rank_spec_names[grepl("_bac$", names(rank_spec_names))])
all_fungi <- unlist(rank_spec_names[grepl("_fun$", names(rank_spec_names))])

# Assign pretty_group based on species membership using proper taxonomy
assign_pretty_group <- function(data) {
  if ("species" %in% colnames(data)) {
    # For functional groups, use the proper taxonomy assignment
    if (any(data$species %in% fg_names)) {
      # Assign functional group categories first
      data$fg_cat <- assign_fg_categories(data$species)
      # Then assign kingdoms based on functional categories
      data$group <- assign_fg_kingdoms(data$fg_cat)
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

scores_df <- assign_pretty_group(scores_df)

cat("Scores data filtered:", nrow(scores_df), "rows\n")
cat("Group distribution:\n")
print(table(scores_df$pretty_group))

# Step 4c: Also get predictor effects data which has better group information
cat("Step 4c: Reading predictor effects data for better group information...\n")
predictor_effects <- readRDS(here("data", "summary/predictor_effects.rds"))
cat("Predictor effects data loaded:", nrow(predictor_effects), "rows\n")
cat("Predictor effects group distribution:\n")
print(table(predictor_effects$pretty_group))

# Step 4b: Get CV data for spatial heterogeneity
cat("Step 4b: Getting CV data for spatial heterogeneity...\n")
cv_data = scores_list$cv_metric_scaled %>%
  filter(model_id %in% converged &
         cv_type == "mean_per_site_cv" &
         !is.na(cv)) %>%
  select(model_id, species, pretty_group, cv) %>%
  distinct()
cat("CV data filtered:", nrow(cv_data), "rows\n")

# Step 5: Read in model estimates of seasonal amplitude
cat("Step 5: Reading seasonal amplitude data...\n")
seasonal_amplitude_in = readRDS(here("data/summary/seasonal_amplitude.rds"))
seas_amplitude_long <- seasonal_amplitude_in[[6]] %>% filter(time_period=="20130601_20180101")  %>%
  select(model_id, taxon, model_name, amplitude) %>%
  mutate(model_id = gsub("_(combined|beta_regression)$", "", model_id))
cat("Seasonal amplitude data loaded:", nrow(seas_amplitude_long), "rows\n")

# Step 6: Read in model estimates of environmental predictors
cat("Step 6: Reading environmental predictors data...\n")
beta_in <- readRDS(here("data", "summary/predictor_effects.rds")) %>%
  filter(time_period=="20130601_20180101") %>%
  select(model_id, beta, effSize) %>%
  group_by(model_id, beta) %>%
  summarise(effSize = mean(effSize, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(id_cols = c("model_id"),
              names_from = "beta", values_from = "effSize") %>%
  mutate(model_id_base = gsub("_(combined|beta_regression)$", "", model_id))

# Filter to only include models that have scoring data
beta_in <- beta_in %>% filter(model_id_base %in% converged_base)
cat("Environmental predictors data loaded:", nrow(beta_in), "rows\n")

# Step 7: Merge all data
cat("Step 7: Merging all data...\n")
master_df <- merge(scores_df, rho_core_in, by.x = "model_id", by.y = "model_id_base")
cat("After merge with rho_core:", nrow(master_df), "rows\n")

master_df <- merge(master_df, moran.stat_all_rank, by.x = "species", by.y = "taxon", all.x = T)
cat("After merge with Moran:", nrow(master_df), "rows\n")

master_df <- as.data.frame(master_df)
seas_amplitude_long <- as.data.frame(seas_amplitude_long)
master_df <- merge(master_df, seas_amplitude_long, all.x = T)
cat("After merge with seasonal amplitude:", nrow(master_df), "rows\n")

# Clean up duplicate columns after each merge
if("model_id.y" %in% colnames(master_df)) {
  master_df$model_id <- master_df$model_id.y
  master_df$model_id.y <- NULL
}
# Remove any other duplicate columns
master_df <- master_df[, !duplicated(colnames(master_df))]

beta_in <- as.data.frame(beta_in)
master_df <- merge(master_df, beta_in, by.x = "model_id", by.y = "model_id_base", all.x = T)
cat("After merge with environmental predictors:", nrow(master_df), "rows\n")

# Clean up duplicate columns
if("model_id.y" %in% colnames(master_df)) {
  master_df$model_id <- master_df$model_id.y
  master_df$model_id.y <- NULL
}
# Remove any other duplicate columns
master_df <- master_df[, !duplicated(colnames(master_df))]

cv_data <- as.data.frame(cv_data)
master_df <- merge(master_df, cv_data, all.x = T)
cat("After merge with CV data:", nrow(master_df), "rows\n")

# Clean up duplicate columns
if("model_id.y" %in% colnames(master_df)) {
  master_df$model_id <- master_df$model_id.y
  master_df$model_id.y <- NULL
}
# Remove any other duplicate columns
master_df <- master_df[, !duplicated(colnames(master_df))]
cat("Cleaned up duplicate columns\n")

# Step 8: Create the analysis dataset
cat("Step 8: Creating analysis dataset...\n")
to_model <- master_df %>% mutate(
  `temporal memory` = rho,
  seasonality = amplitude,
  `variation among soil cores` = precision,
  `spatial heterogeneity` = cv,
  `spatial autocorrelation` = mean_morans
) %>%
filter(!is.na(RSQ))

cat("Analysis dataset created:", nrow(to_model), "rows\n")

# Step 9: Add additional columns
cat("Step 9: Adding additional columns...\n")
to_model$seasonal_predictors = ifelse(to_model$model_name %in% c("env_cycl","cycl_only"), 1, 0)
to_model$environmental_predictors = ifelse(to_model$model_name %in% c("env_cycl","env_cov"), 1, 0)
to_model$new_site = ifelse(to_model$site_prediction =="New time x site (random effect)", 1, 0)
to_model$Fungi = ifelse(to_model$pretty_group =="Fungi", 1, 0)
# Handle missing fcast_type column
if("fcast_type" %in% colnames(to_model)) {
  to_model$Taxonomic = ifelse(to_model$fcast_type =="Taxonomic", 1, 0)
} else {
  to_model$Taxonomic = 0  # Default value when fcast_type is missing
}

# Add missing columns with default values
to_model$latitude <- NA
to_model$`spatial autocorrelation` <- NA

cat("Additional columns added\n")

# Step 9b: Analyze spatial metrics separately
cat("Step 9b: Analyzing spatial metrics separately...\n")

# Check data availability for spatial metrics
spatial_summary <- to_model %>% 
  group_by(pretty_group) %>%
  summarise(
    n_total = n(),
    n_cv = sum(!is.na(cv)),
    n_moran = sum(!is.na(mean_morans)),
    n_precision = sum(!is.na(precision)),
    mean_cv = mean(cv, na.rm = TRUE),
    mean_moran = mean(mean_morans, na.rm = TRUE),
    mean_precision = mean(precision, na.rm = TRUE)
  )
print("Spatial metrics summary by group:")
print(spatial_summary)

# Debug: Check if mean_morans column exists
cat("Checking mean_morans column:\n")
cat("mean_morans column exists:", "mean_morans" %in% colnames(to_model), "\n")
if("mean_morans" %in% colnames(to_model)) {
  cat("mean_morans non-NA values:", sum(!is.na(to_model$mean_morans)), "out of", nrow(to_model), "\n")
  cat("Sample mean_morans values:", head(to_model$mean_morans[!is.na(to_model$mean_morans)], 5), "\n")
} else {
  cat("mean_morans column not found in to_model\n")
}

# Step 10: Create scaled dataset
cat("Step 10: Creating scaled dataset...\n")
to_model_scaled = to_model %>% group_by(pretty_group, model_name) %>%
  filter(site_prediction == "New time (observed site)") %>%
  filter(!is.na(RMSE.norm) & !is.na(RSQ)) %>%  # Remove rows with missing key variables
  mutate(
    RMSE.norm = scale(RMSE.norm),
    mean_crps_sample = scale(mean_crps_sample),
    RSQ = scale(RSQ),
    `seasonal amplitude` = scale(seasonality),
    `temporal memory` = scale(`temporal memory`),
    amplitude = scale(amplitude),
    temperature = scale(Temperature),
    moisture = scale(Moisture),
    pH = scale(pH),
    `percent carbon` = scale(pC),
    `ectomycorrhizal trees` = scale(`Ectomycorrhizal\ntrees`),
    LAI = scale(LAI),
    `spatial heterogeneity` = scale(`spatial heterogeneity`),
    `spatial autocorrelation` = scale(mean_morans),
    `variation among soil cores` = scale(`variation among soil cores`)
  )

cat("Scaled dataset created:", nrow(to_model_scaled), "rows\n")

# Step 11: Check data distribution by group
cat("Step 11: Checking data distribution by group...\n")
group_summary <- to_model_scaled %>% 
  group_by(pretty_group) %>% 
  summarise(
    n_models = n(),
    n_taxa = length(unique(taxon)),
    mean_rmse = mean(RMSE.norm, na.rm = TRUE),
    mean_rsq = mean(RSQ, na.rm = TRUE),
    mean_temporal_memory = mean(`temporal memory`, na.rm = TRUE),
    mean_seasonal_amp = mean(`seasonal amplitude`, na.rm = TRUE)
  )
print(group_summary)

# Step 12: Create separate models for fungi and bacteria
cat("Step 12: Creating separate models for fungi and bacteria...\n")

# Check data availability for each group
fungi_data <- to_model_scaled %>% filter(pretty_group == "Fungi")
bacteria_data <- to_model_scaled %>% filter(pretty_group == "Bacteria")

cat("Fungi data:", nrow(fungi_data), "rows\n")
cat("Bacteria data:", nrow(bacteria_data), "rows\n")

# Check which spatial metrics are available
cat("Checking spatial metrics availability:\n")
cat("Fungi - CV available:", sum(!is.na(fungi_data$`spatial heterogeneity`)), "out of", nrow(fungi_data), "\n")
cat("Fungi - Moran available:", sum(!is.na(fungi_data$`spatial autocorrelation`)), "out of", nrow(fungi_data), "\n")
cat("Fungi - Precision available:", sum(!is.na(fungi_data$`variation among soil cores`)), "out of", nrow(fungi_data), "\n")
cat("Bacteria - CV available:", sum(!is.na(bacteria_data$`spatial heterogeneity`)), "out of", nrow(bacteria_data), "\n")
cat("Bacteria - Moran available:", sum(!is.na(bacteria_data$`spatial autocorrelation`)), "out of", nrow(bacteria_data), "\n")
cat("Bacteria - Precision available:", sum(!is.na(bacteria_data$`variation among soil cores`)), "out of", nrow(bacteria_data), "\n")

# Create models with available data
if (nrow(fungi_data) > 0 && sum(!is.na(fungi_data$RMSE.norm)) > 0) {
  # Check if we have enough non-NA values for all predictors
  predictor_cols <- c('temperature', 'moisture', 'pH', 'percent carbon', 'ectomycorrhizal trees', 'LAI', 'temporal memory', 'seasonal amplitude')
  complete_cases <- complete.cases(fungi_data[, predictor_cols, drop = FALSE])
  
  if(sum(complete_cases) >= 2) {  # Need at least 2 complete cases
    tryCatch({
      fungi_lm <- lm(RMSE.norm ~
        temperature +
        moisture + pH + `percent carbon` +
        `ectomycorrhizal trees` + LAI +
        `temporal memory` + `seasonal amplitude`,
        data = fungi_data)
      
      cat("Fungi basic model created successfully\n")
      cat("Fungi basic model summary:\n")
      print(summary(fungi_lm))
    }, error = function(e) {
      cat("Error creating fungi model:", e$message, "\n")
      fungi_lm <<- NULL
    })
  } else {
    cat("Insufficient complete cases for fungi model (need at least 2, have", sum(complete_cases), ")\n")
    fungi_lm <- NULL
  }
} else {
  cat("Insufficient data for fungi model\n")
  fungi_lm <- NULL
}

# Bacteria model
if (nrow(bacteria_data) > 0 && sum(!is.na(bacteria_data$RMSE.norm)) > 0) {
  # Check if we have enough non-NA values for all predictors
  predictor_cols <- c('temperature', 'moisture', 'pH', 'percent carbon', 'ectomycorrhizal trees', 'LAI', 'temporal memory', 'seasonal amplitude')
  complete_cases <- complete.cases(bacteria_data[, predictor_cols, drop = FALSE])
  
  if(sum(complete_cases) >= 2) {  # Need at least 2 complete cases
    tryCatch({
      bacteria_lm <- lm(RMSE.norm ~
        temperature +
        moisture + pH + `percent carbon` +
        `ectomycorrhizal trees` + LAI +
        `temporal memory` + `seasonal amplitude`,
        data = bacteria_data)
      
      cat("Bacteria basic model created successfully\n")
      cat("Bacteria basic model summary:\n")
      print(summary(bacteria_lm))
    }, error = function(e) {
      cat("Error creating bacteria model:", e$message, "\n")
      bacteria_lm <<- NULL
    })
  } else {
    cat("Insufficient complete cases for bacteria model (need at least 2, have", sum(complete_cases), ")\n")
    bacteria_lm <- NULL
  }
} else {
  cat("Insufficient data for bacteria model\n")
  bacteria_lm <- NULL
}

# Step 12b: Analyze spatial metrics separately
cat("Step 12b: Analyzing spatial metrics separately...\n")

# Create models with spatial metrics where available
if (!is.null(fungi_lm)) {
  # Fungi model with spatial metrics
  fungi_spatial_data <- fungi_data %>% filter(!is.na(`spatial heterogeneity`) & !is.na(`variation among soil cores`))
  cat("Fungi spatial data available:", nrow(fungi_spatial_data), "rows\n")
  
  if (nrow(fungi_spatial_data) > 10) {
    fungi_spatial_lm <- lm(RMSE.norm ~
      temperature + moisture + pH + `percent carbon` +
      `ectomycorrhizal trees` + LAI +
      `temporal memory` + `seasonal amplitude` +
      `spatial heterogeneity` + `spatial autocorrelation` + `variation among soil cores`,
      data = fungi_spatial_data)
    
    cat("Fungi spatial model created successfully\n")
    cat("Fungi spatial model summary:\n")
    print(summary(fungi_spatial_lm))
  } else {
    cat("Insufficient spatial data for fungi spatial model\n")
    fungi_spatial_lm <- NULL
  }
}

if (!is.null(bacteria_lm)) {
  # Bacteria model with spatial metrics
  bacteria_spatial_data <- bacteria_data %>% filter(!is.na(`spatial heterogeneity`) & !is.na(`variation among soil cores`))
  cat("Bacteria spatial data available:", nrow(bacteria_spatial_data), "rows\n")
  
  if (nrow(bacteria_spatial_data) > 10) {
    bacteria_spatial_lm <- lm(RMSE.norm ~
      temperature + moisture + pH + `percent carbon` +
      `ectomycorrhizal trees` + LAI +
      `temporal memory` + `seasonal amplitude` +
      `spatial heterogeneity` + `spatial autocorrelation` + `variation among soil cores`,
      data = bacteria_spatial_data)
    
    cat("Bacteria spatial model created successfully\n")
    cat("Bacteria spatial model summary:\n")
    print(summary(bacteria_spatial_lm))
  } else {
    cat("Insufficient spatial data for bacteria spatial model\n")
    bacteria_spatial_lm <- NULL
  }
}

# Step 13: Create comparison plot
cat("Step 13: Creating comparison plot...\n")

# Debug: Check model diagnostics
if (!is.null(fungi_lm)) {
  cat("Fungi model diagnostics:\n")
  cat("Degrees of freedom:", fungi_lm$df.residual, "\n")
  cat("Number of observations used:", length(fungi_lm$fitted.values), "\n")
  cat("Number of observations deleted:", sum(is.na(fungi_lm$model)), "\n")
}

if (!is.null(bacteria_lm)) {
  cat("Bacteria model diagnostics:\n")
  cat("Degrees of freedom:", bacteria_lm$df.residual, "\n")
  cat("Number of observations used:", length(bacteria_lm$fitted.values), "\n")
  cat("Number of observations deleted:", sum(is.na(bacteria_lm$model)), "\n")
}

if (!is.null(fungi_lm) && !is.null(bacteria_lm)) {
  # Create side-by-side plots with better error bar and significance handling
  fungi_plot <- sjPlot::plot_models(model = list(fungi_lm),
    vline.color = "black", show.legend = F, colors = "gs",
    show.p = TRUE, p.shape = TRUE, p.threshold = c(0.05, 0.01, 0.001),
    ci.lvl = 0.95, show.values = FALSE) +
    theme_bw(base_size = 16) + 
    ylab("Parameter effect on forecast error (nRMSE)") +
    ggtitle("Fungi") +
    ylim(c(-1.5,1.5)) +
    theme(plot.margin = unit(c(1,1,1,0), "cm")) +
    annotate("rect", xmin = 8.6, xmax = 9.7, ymin = -1.5, ymax = 1.5,
      color="black", fill = "lightgrey")+
    annotate(label = "Environmental predictors\n& seasonality", x= 9.2, y=0, geom="text", size=4) +  
    coord_flip(expand = FALSE)
  
  bacteria_plot <- sjPlot::plot_models(model = list(bacteria_lm),
    vline.color = "black", show.legend = F, colors = "gs",
    show.p = TRUE, p.shape = TRUE, p.threshold = c(0.05, 0.01, 0.001),
    ci.lvl = 0.95, show.values = FALSE) +
    theme_bw(base_size = 16) + 
    ylab("Parameter effect on forecast error (nRMSE)") +
    ggtitle("Bacteria") +
    ylim(c(-1.5,1.5)) +
    theme(plot.margin = unit(c(1,1,1,0), "cm")) +
    annotate("rect", xmin = 8.6, xmax = 9.7, ymin = -1.5, ymax = 1.5,
      color="black", fill = "lightgrey")+
    annotate(label = "Environmental predictors\n& seasonality", x= 9.2, y=0, geom="text", size=4) +  
    coord_flip(expand = FALSE)
  
  # Combine plots
  combined_plot <- cowplot::plot_grid(fungi_plot, bacteria_plot, 
                            ncol = 2, 
                            labels = c("A", "B"),
                            label_size = 20)
  
  cat("Saving comparison plot...\n")
  png(here("figures","explain_predictability_by_group.png"), width = 1600, height=800)
  print(combined_plot)
  dev.off()
  cat("Plot saved to figures/explain_predictability_by_group.png\n")
  
} else {
  cat("Cannot create comparison plot - missing model data\n")
  
  # Create a simple plot using the predictor effects data instead
  cat("Creating alternative plot using predictor effects data...\n")
  
  # Create a simple comparison plot using predictor effects
  library(ggplot2)
  
  # Filter predictor effects data for the main environmental predictors
  beta_names <- c("Ectomycorrhizal\ntrees", "LAI", "pC", "pH", "Temperature", "Moisture")
  plot_data <- predictor_effects %>%
    filter(beta %in% beta_names & 
           model_name == "env_cycl" &
           time_period == "20130601_20180101") %>%
    filter(!is.na(effSize) & !is.na(pretty_group))
  
  if(nrow(plot_data) > 0) {
    # Create a simple comparison plot
    simple_plot <- ggplot(plot_data, aes(x = beta, y = effSize, fill = pretty_group)) +
      geom_boxplot(alpha = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      theme_bw(base_size = 14) +
      labs(title = "Effect Sizes by Group (Fungi vs Bacteria)",
           x = "Environmental Predictor",
           y = "Effect Size",
           fill = "Group") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom") +
      scale_fill_manual(values = c("Bacteria" = "lightblue", "Fungi" = "lightcoral"))
    
    # Save the simple plot
    ggsave(here("figures", "predictability_by_group_simple.png"), simple_plot, 
           width = 12, height = 8, dpi = 300)
    cat("Simple plot saved to figures/predictability_by_group_simple.png\n")
  } else {
    cat("No data available for simple plot either\n")
  }
}

# Step 14: Statistical comparison
cat("Step 14: Statistical comparison between groups...\n")

if (!is.null(fungi_lm) && !is.null(bacteria_lm)) {
  # Extract coefficients for comparison
  fungi_coef <- summary(fungi_lm)$coefficients
  bacteria_coef <- summary(bacteria_lm)$coefficients
  
  # Create comparison table
  comparison_df <- data.frame(
    Parameter = rownames(fungi_coef)[-1], # Exclude intercept
    Fungi_Estimate = fungi_coef[-1, "Estimate"],
    Fungi_Pvalue = fungi_coef[-1, "Pr(>|t|)"],
    Bacteria_Estimate = bacteria_coef[-1, "Estimate"],
    Bacteria_Pvalue = bacteria_coef[-1, "Pr(>|t|)"],
    Difference = fungi_coef[-1, "Estimate"] - bacteria_coef[-1, "Estimate"]
  )
  
  cat("Coefficient comparison:\n")
  print(comparison_df)
  
  # Identify significant differences (where one group is significant and other is not)
  significant_diffs <- comparison_df %>%
    filter((Fungi_Pvalue < 0.05 & Bacteria_Pvalue >= 0.05) | 
           (Fungi_Pvalue >= 0.05 & Bacteria_Pvalue < 0.05))
  
  if (nrow(significant_diffs) > 0) {
    cat("\nSignificant differences between groups:\n")
    print(significant_diffs)
  } else {
    cat("\nNo significant differences in significance patterns between groups\n")
  }
}

cat("=== ANALYSIS COMPLETED ===\n")