# Clean Saprophytic Fungi Phenology Visualization
# Uses existing hindcast data with proper confidence intervals (no zigzag)
# Follows patterns from existing visualization scripts

cat("=== CLEAN SAPROPHYTIC FUNGI PHENOLOGY VISUALIZATION ===\n\n")

# Load required packages and environment
source("../../source.R")
library(patchwork)

# Function to get core-level observed data for any species
get_core_data <- function(species, min_date, max_date) {
  tryCatch({
    # Load fungal functional group data
    fungi <- readRDS(here("data/clean/groupAbundances_ITS_2023.rds"))
    
    if (!species %in% names(fungi)) {
      cat("Species", species, "not found in data\n")
      return(NULL)
    }
    
    rank.df <- fungi[[species]]
    
    if (is.null(rank.df)) {
      cat("No data for species", species, "\n")
      return(NULL)
    }
    
    # Filter and prepare data
    core_data <- rank.df %>%
      select(siteID, plotID, dateID, sampleID, dates, plot_date, !!species) %>%
      filter(dates >= as.Date(min_date, format = "%Y%m%d") & 
             dates <= as.Date(max_date, format = "%Y%m%d")) %>%
      filter(!is.na(!!sym(species))) %>%
      mutate(
        month = month(dates),
        month_name = month.abb[month],
        year = year(dates),
        season = case_when(
          month >= 6 & month < 9 ~ "Summer",
          month >= 9 & month < 12 ~ "Fall",
          month >= 12 | month < 3 ~ "Winter",
          TRUE ~ "Spring"
        ),
        abundance = !!sym(species)
      )
    
    return(core_data)
    
  }, error = function(e) {
    cat("Error processing", species, ":", e$message, "\n")
    return(NULL)
  })
}

cat("✅ Environment loaded successfully\n")

# Load environmental parameters from both logit and CLR models
cat("📊 Loading environmental parameters...\n")
pred_effects <- readRDS(here("data/summary/predictor_effects.rds"))
sapr_effects <- pred_effects[pred_effects$taxon == "saprotroph" & !is.na(pred_effects$taxon), ]

# Also load CLR model results if available
clr_effects <- NULL
tryCatch({
  clr_summaries <- readRDS(here("data/summary/CLR_regression_summaries.rds"))
  if ("saprotroph" %in% names(clr_summaries)) {
    clr_effects <- clr_summaries[["saprotroph"]]
    cat("✅ Loaded CLR model results for saprotroph\n")
  }
}, error = function(e) {
  cat("Note: CLR summaries not available, proceeding with logit models only\n")
})

cat("Found", nrow(sapr_effects), "parameter estimates for saprotroph\n")

# Extract seasonal parameters for different model types (logit and CLR)
seasonal_params <- sapr_effects[sapr_effects$beta %in% c("sin", "cos"), ]

# Also extract CLR seasonal parameters if available
if (!is.null(clr_effects)) {
  for (model_name in names(clr_effects)) {
    if (model_name %in% c("cycl_only", "env_cycl")) {
      model_data <- clr_effects[[model_name]]
      if ("sin" %in% names(model_data) && "cos" %in% names(model_data)) {
        sin_coef <- model_data$sin$Mean
        cos_coef <- model_data$cos$Mean
        
        if (!is.na(sin_coef) && !is.na(cos_coef)) {
          amplitude <- sqrt(sin_coef^2 + cos_coef^2)
          peak_timing <- atan2(cos_coef, sin_coef) * 12/(2*pi)
          if (peak_timing < 0) peak_timing <- peak_timing + 12
          
          seasonal_params <- rbind(seasonal_params, data.frame(
            beta = c("sin", "cos"),
            Mean = c(sin_coef, cos_coef),
            SD = c(model_data$sin$SD, model_data$cos$SD),
            effSize = c(abs(sin_coef), abs(cos_coef)),
            significant = c(0, 0),  # Will determine based on CI
            model_name = c(paste0("CLR_", model_name), paste0("CLR_", model_name)),
            stringsAsFactors = FALSE
          ))
          
          cat("CLR Model", model_name, ": Amplitude =", round(amplitude, 3), 
              ", Peak month =", round(peak_timing, 1), "\n")
        }
      }
    }
  }
}

models_with_seasonality <- list()
for (model in unique(seasonal_params$model_name)) {
  sin_coef <- seasonal_params[seasonal_params$beta == "sin" & seasonal_params$model_name == model, "Mean"]
  cos_coef <- seasonal_params[seasonal_params$beta == "cos" & seasonal_params$model_name == model, "Mean"]
  
  if (length(sin_coef) > 0 && length(cos_coef) > 0) {
    amplitude <- sqrt(sin_coef^2 + cos_coef^2)
    peak_timing <- atan2(cos_coef, sin_coef) * 12/(2*pi)
    if (peak_timing < 0) peak_timing <- peak_timing + 12
    
    # Get SD for uncertainty
    sin_sd <- seasonal_params[seasonal_params$beta == "sin" & seasonal_params$model_name == model, "SD"]
    cos_sd <- seasonal_params[seasonal_params$beta == "cos" & seasonal_params$model_name == model, "SD"]
    
    models_with_seasonality[[model]] <- list(
      sin_coef = sin_coef,
      cos_coef = cos_coef,
      sin_sd = sin_sd,
      cos_sd = cos_sd,
      amplitude = amplitude,
      peak_timing = peak_timing,
      peak_month = round(peak_timing, 1)
    )
    
    cat("Model", model, ": Amplitude =", round(amplitude, 3), 
        ", Peak month =", round(peak_timing, 1), "\n")
  }
}

# Extract most significant environmental parameters (logit and CLR)
env_params <- sapr_effects[!sapr_effects$beta %in% c("sin", "cos"), ]

# Also extract CLR environmental parameters if available
if (!is.null(clr_effects)) {
  for (model_name in names(clr_effects)) {
    if (model_name %in% c("env_cov", "env_cycl")) {
      model_data <- clr_effects[[model_name]]
      # Extract environmental parameters (excluding sin/cos)
      env_vars <- names(model_data)[!names(model_data) %in% c("sin", "cos")]
      
      for (var in env_vars) {
        if (var %in% names(model_data) && "Mean" %in% names(model_data[[var]])) {
          mean_val <- model_data[[var]]$Mean
          sd_val <- model_data[[var]]$SD
          
          if (!is.na(mean_val) && !is.na(sd_val)) {
            # Add to env_params
            new_row <- data.frame(
              rowname = "beta",
              beta_num = "1",
              Mean = mean_val,
              SD = sd_val,
              Naive_SE = NA,
              Time_series_SE = NA,
              beta = var,
              taxon = "saprotroph",
              significant = ifelse(abs(mean_val) > 2 * sd_val, 1, 0),  # Rough significance test
              effSize = abs(mean_val),
              Median = NA,
              site_num = NA,
              siteID = NA,
              rank = var,
              model_name = paste0("CLR_", model_name),
              group = "ITS",
              rank_only = "functional",
              time_period = "2015-11_2018-01",
              fcast_type = "Functional",
              pretty_group = "Fungi",
              model_id = paste0("CLR_", model_name, "_saprotroph"),
              tax_rank = var,
              pretty_name = NA,
              only_rank = NA,
              stringsAsFactors = FALSE
            )
            
            env_params <- rbind(env_params, new_row)
          }
        }
      }
    }
  }
}

env_params <- env_params[order(-abs(env_params$effSize)), ]

cat("\n🌍 Most significant environmental parameters:\n")
top_env_params <- head(env_params, 6)
for (i in 1:nrow(top_env_params)) {
  param <- top_env_params[i, ]
  cat(sprintf("%d. %s (%s): Effect = %.3f%s\n", 
              i, param$beta, param$model_name, param$effSize,
              ifelse(param$significant == 1, " *significant*", "")))
}

# Load existing hindcast data for visualization
cat("\n📊 Loading existing hindcast data for visualization...\n")

# Check if we have the main hindcast data
hindcast_file <- here("data/summary/hindcast_saprotroph.rds")
if (!file.exists(hindcast_file)) {
  stop("hindcast_saprotroph.rds not found. Please ensure the hindcast data is available.")
}

hindcast_data <- readRDS(hindcast_file)
cat("Loaded hindcast data:", nrow(hindcast_data), "rows\n")

# Check what columns are available for confidence intervals
cat("Available columns:", names(hindcast_data), "\n")

# Look for confidence interval columns
ci_cols <- names(hindcast_data)[grepl("lo|hi|2\\.5|97\\.5|25|75", names(hindcast_data))]
cat("Confidence interval columns found:", ci_cols, "\n")

# Use existing hindcast data for visualization
plot_data <- hindcast_data %>%
  filter(!is.na(med)) %>%  # Remove rows without predictions
  mutate(
    # Create model type labels
    model_type = case_when(
      grepl("CLR_", model_name) ~ "CLR Model",
      model_name == "cycl_only" ~ "Seasonal Only",
      model_name == "env_cov" ~ "Environmental Only",
      model_name == "env_cycl" ~ "Environmental + Seasonal",
      model_name == "all_covariates" ~ "All Covariates",
      TRUE ~ model_name
    )
  )

cat("Plot data prepared:", nrow(plot_data), "rows\n")

# Also load core observation data for comparison
cat("📊 Loading core-level observation data...\n")
core_data <- get_core_data("saprotroph", "20150101", "20200101")
if (!is.null(core_data) && nrow(core_data) > 0) {
  combined_core_data <- core_data
  cat("Loaded", nrow(core_data), "core observations\n")
} else {
  cat("No core data found, proceeding with hindcast data only\n")
  combined_core_data <- data.frame()
}

cat("\n🎨 Creating phenology-style visualization...\n")

# Create the main phenology plot using existing hindcast data
main_plot <- ggplot() +
  # Add seasonal background (subtle)
  annotate("rect", xmin = 5.5, xmax = 8.5, ymin = -Inf, ymax = Inf, 
           alpha = 0.1, fill = "gold") +
  annotate("rect", xmin = 8.5, xmax = 11.5, ymin = -Inf, ymax = Inf, 
           alpha = 0.1, fill = "orange") +
  annotate("rect", xmin = 11.5, xmax = 2.5, ymin = -Inf, ymax = Inf, 
           alpha = 0.1, fill = "lightblue") +
  annotate("rect", xmin = 2.5, xmax = 5.5, ymin = -Inf, ymax = Inf, 
           alpha = 0.1, fill = "lightgreen") +
  
  # Add seasonal labels
  annotate("text", x = 7, y = max(plot_data$med, na.rm = TRUE) * 0.9, 
           label = "Summer", color = "darkgoldenrod", fontface = "bold", size = 3) +
  annotate("text", x = 10, y = max(plot_data$med, na.rm = TRUE) * 0.9, 
           label = "Fall", color = "darkorange", fontface = "bold", size = 3) +
  annotate("text", x = 1, y = max(plot_data$med, na.rm = TRUE) * 0.9, 
           label = "Winter", color = "darkblue", fontface = "bold", size = 3) +
  annotate("text", x = 4, y = max(plot_data$med, na.rm = TRUE) * 0.9, 
           label = "Spring", color = "darkgreen", fontface = "bold", size = 3) +
  
  # Add uncertainty bands using existing hindcast data (no zigzag!)
  geom_ribbon(data = plot_data, 
              aes(x = month(dates), ymin = lo, ymax = hi, fill = model_type), 
              alpha = 0.15) +
  
  # Add predicted trend lines using existing hindcast data
  geom_line(data = plot_data, 
            aes(x = month(dates), y = med, color = model_type, linetype = model_type), 
            linewidth = 1.2, alpha = 0.8) +
  
  # Add individual core observations with jitter
  geom_point(data = combined_core_data, 
             aes(x = month, y = abundance), 
             size = 1.5, alpha = 0.6, color = "darkgreen",
             position = position_jitter(width = 0.2)) +
  
  # Enhanced styling
  scale_color_manual(values = c(
    "CLR Model" = "#d62728",
    "Seasonal Only" = "#1f77b4",
    "Environmental Only" = "#ff7f0e", 
    "Environmental + Seasonal" = "#2ca02c",
    "All Covariates" = "#9467bd"
  )) +
  scale_fill_manual(values = c(
    "CLR Model" = "#d62728",
    "Seasonal Only" = "#1f77b4",
    "Environmental Only" = "#ff7f0e", 
    "Environmental + Seasonal" = "#2ca02c",
    "All Covariates" = "#9467bd"
  )) +
  scale_linetype_manual(values = c(
    "CLR Model" = "solid",
    "Seasonal Only" = "solid",
    "Environmental Only" = "dashed",
    "Environmental + Seasonal" = "dotdash",
    "All Covariates" = "longdash"
  )) +
  scale_x_continuous(breaks = 1:12, 
                    labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"),
                    limits = c(1, 12)) +
  
  # Labels and theme
  labs(title = "Saprophytic Fungi: Seasonal Phenology with Core-Level Observations",
       subtitle = paste("Peak timing:", paste(sapply(names(models_with_seasonality), function(name) {
         params <- models_with_seasonality[[name]]
         paste(month.abb[round(params$peak_timing)], "(", name, ")")
       }), collapse = ", "),
         "| Green points = observed values, lines = model predictions with uncertainty bands"),
       x = "Month", 
       y = "Relative Abundance",
       color = "Model Type",
       linetype = "Model Type") +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, color = "darkgray", hjust = 0.5),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  )

# Create environmental parameters forest plot
cat("📊 Creating environmental parameters forest plot...\n")

# Prepare data for forest plot with confidence intervals
forest_data <- top_env_params %>%
  mutate(
    # Calculate confidence intervals (95% CI = mean ± 1.96 * SD)
    ci_lower = Mean - 1.96 * SD,
    ci_upper = Mean + 1.96 * SD,
    # Determine significance based on CI not crossing zero
    significant = ifelse(sign(ci_lower) == sign(ci_upper), 1, 0),
    # Create model type labels
    model_type = case_when(
      grepl("CLR_", model_name) ~ "CLR Model",
      model_name == "env_cov" ~ "Logit Environmental",
      model_name == "env_cycl" ~ "Logit Environmental + Seasonal",
      TRUE ~ model_name
    ),
    # Create parameter labels
    param_label = case_when(
      beta == "pC" ~ "Soil Carbon",
      beta == "Ectomycorrhizal\ntrees" ~ "Ectomycorrhizal Trees",
      beta == "Temperature" ~ "Temperature",
      beta == "Moisture" ~ "Soil Moisture",
      beta == "pH" ~ "Soil pH",
      beta == "LAI" ~ "Leaf Area Index",
      TRUE ~ as.character(beta)
    )
  )

# Create forest plot
env_forest_plot <- ggplot(forest_data, 
                          aes(x = reorder(param_label, abs(effSize)), y = Mean, 
                              color = model_type, shape = factor(significant))) +
  # Add reference line at zero
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7) +
  
  # Add confidence intervals
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper, color = model_type), 
                width = 0.3, linewidth = 0.8, alpha = 0.8) +
  
  # Add point estimates
  geom_point(aes(fill = model_type), size = 3, alpha = 0.9) +
  
  # Add significance indicators
  geom_text(aes(label = ifelse(significant == 1, "*", ""), 
                y = ifelse(Mean >= 0, ci_upper + max(abs(ci_upper)) * 0.05, 
                           ci_lower - max(abs(ci_lower)) * 0.05)), 
            size = 4, color = "red", fontface = "bold") +
  
  # Color and shape scales
  scale_color_manual(values = c(
    "CLR Model" = "#d62728",
    "Logit Environmental" = "#ff7f0e", 
    "Logit Environmental + Seasonal" = "#2ca02c"
  )) +
  scale_fill_manual(values = c(
    "CLR Model" = "#d62728",
    "Logit Environmental" = "#ff7f0e", 
    "Logit Environmental + Seasonal" = "#2ca02c"
  )) +
  scale_shape_manual(values = c("0" = 21, "1" = 19), 
                     labels = c("0" = "Non-significant", "1" = "Significant")) +
  
  # Flip coordinates for horizontal forest plot
  coord_flip() +
  
  # Labels and theme
  labs(title = "Environmental Parameter Effects on Saprophytic Fungi",
       subtitle = "Forest plot showing effect sizes with 95% confidence intervals\n* indicates significant effects (CI does not cross zero)",
       x = "Environmental Parameter",
       y = "Effect Size (logit scale)",
       color = "Model Type",
       fill = "Model Type",
       shape = "Significance") +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, color = "darkgray", hjust = 0.9),
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  ) +
  guides(color = guide_legend(nrow = 1), 
         fill = guide_legend(nrow = 1),
         shape = guide_legend(nrow = 1))

# Create combined plot
cat("🎨 Creating combined visualization...\n")

combined_plot <- (main_plot / env_forest_plot) +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(
    title = "Saprophytic Fungi: Seasonal Patterns and Environmental Drivers",
    subtitle = paste("Analysis based on", nrow(plot_data), "hindcast predictions across NEON sites"),
    caption = "Top: Seasonal phenology with uncertainty bands | Bottom: Environmental parameter effects (forest plot)",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      plot.caption = element_text(hjust = 0.5, size = 10, color = "gray50")
    )
  )

# Save all plots
cat("💾 Saving improved visualizations...\n")

# Create output directory
output_dir <- here("figures", "saprotroph_phenology_clean")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Save main phenology plot
ggsave(
  filename = file.path(output_dir, "saprotroph_phenology_main_clean.png"),
  plot = main_plot,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

# Save environmental parameters forest plot
ggsave(
  filename = file.path(output_dir, "saprotroph_environmental_effects_forest_clean.png"),
  plot = env_forest_plot,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)

# Save combined plot
ggsave(
  filename = file.path(output_dir, "saprotroph_phenology_combined_clean.png"),
  plot = combined_plot,
  width = 12,
  height = 12,
  dpi = 300,
  bg = "white"
)

cat("✅ All clean visualizations saved to:", output_dir, "\n")

# Create summary report
cat("\n📊 SAPROTROPH PHENOLOGY ANALYSIS SUMMARY\n")
cat("=========================================\n")

cat("\n🍃 SEASONAL PATTERNS:\n")
for (model_name in names(models_with_seasonality)) {
  params <- models_with_seasonality[[model_name]]
  month_name <- month.abb[round(params$peak_timing)]
  cat(sprintf("  %s: Peak in %s (month %.1f), Amplitude = %.3f\n",
              model_name, month_name, params$peak_timing, params$amplitude))
}

cat("\n🌍 TOP ENVIRONMENTAL DRIVERS:\n")
for (i in 1:min(3, nrow(top_env_params))) {
  param <- top_env_params[i, ]
  direction <- ifelse(param$Mean > 0, "positive", "negative")
  significance <- ifelse(param$significant == 1, " (significant)", "")
  cat(sprintf("  %d. %s: %s effect (%.3f)%s\n", 
              i, param$beta, direction, abs(param$effSize), significance))
}

cat("\n📊 DATA COVERAGE:\n")
cat("  Total hindcast predictions:", nrow(plot_data), "\n")
if (nrow(combined_core_data) > 0) {
  cat("  Core observations:", nrow(combined_core_data), "\n")
  cat("  Sites represented:", length(unique(combined_core_data$siteID)), "\n")
  cat("  Seasonal distribution:\n")
  season_counts <- table(combined_core_data$season)
  for (season in names(season_counts)) {
    cat("    ", season, ":", season_counts[season], "\n")
  }
}

cat("\n=== CLEAN SAPROPHYTIC FUNGI VISUALIZATION COMPLETE ===\n")
cat("✅ No more zigzag patterns - using existing confidence intervals!\n")
