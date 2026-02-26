source("source.R")

# Check if required data files exist
required_files <- c(
  "data/clean/pheno_group_peak_phenophases.rds",
  "data/summary/pheno_predictor_effects.rds",
  "data/summary/predictor_effects.rds",
  "data/clean/all_predictor_data.rds"
)

optional_files <- c(
  "data/clean/group_peak_phenophases.rds"
)

missing_files <- c()
for (file in required_files) {
  if (!file.exists(here(file))) {
    missing_files <- c(missing_files, file)
  }
}

if (length(missing_files) > 0) {
  cat("Missing required files:", paste(missing_files, collapse = ", "), "\n")
  stop("Cannot proceed without required data files")
}

# Check optional files
for (file in optional_files) {
  if (!file.exists(here(file))) {
    cat("Optional file not found:", file, "\n")
  }
}

# Load data
phenophase_in_legacy = readRDS(here("data/clean/pheno_group_peak_phenophases.rds"))
beta_effects_legacy = readRDS(here("data/summary/pheno_predictor_effects.rds"))
beta_effects_calibration = readRDS(here("data/summary/predictor_effects.rds"))

# Try to load optional phenophase file
if (file.exists(here("data/clean/group_peak_phenophases.rds"))) {
  phenophase_in = readRDS(here("data/clean/group_peak_phenophases.rds"))
} else {
  cat("Using only legacy phenophase data\n")
  phenophase_in <- NULL
}

# Check data structure
if (!is.null(phenophase_in) && length(phenophase_in) >= 4) {
  max_abun = phenophase_in[[4]] %>% filter(siteID %in% c("CPER","DSNY","HARV","OSBS","STER"))
} else {
  cat("phenophase_in not available or does not have expected structure\n")
  max_abun <- data.frame()
}

if (length(phenophase_in_legacy) >= 4) {
  # Check if siteID column exists
  if ("siteID" %in% names(phenophase_in_legacy[[4]])) {
    max_abun_legacy = phenophase_in_legacy[[4]] %>% filter(siteID %in% c("CPER","DSNY","HARV","OSBS","STER"))
  } else {
    cat("siteID column not found in phenophase_in_legacy[[4]]\n")
    max_abun_legacy <- data.frame()
  }
} else {
  cat("phenophase_in_legacy does not have expected structure\n")
  max_abun_legacy <- data.frame()
}

# Merge beta effects
all_betas_merged = rbind(beta_effects_legacy, beta_effects_calibration)

# Check available columns for merging
cat("Legacy columns:", paste(names(beta_effects_legacy), collapse = ", "), "\n")
cat("Calibration columns:", paste(names(beta_effects_calibration), collapse = ", "), "\n")

# Find common columns for merging
common_cols <- intersect(names(beta_effects_legacy), names(beta_effects_calibration))
cat("Common columns for merging:", paste(common_cols, collapse = ", "), "\n")

# Create legacy vs calibration comparison with available columns
if (length(common_cols) > 0 && nrow(beta_effects_legacy) > 0) {
  legacy_vs_calibration = merge(beta_effects_legacy, 
                               beta_effects_calibration %>% filter(time_period == "2015-11_2018-01"), 
                               by = common_cols, suffixes = c("legacy","calibration"))
  
  # Only filter if the column exists
  if ("model_name" %in% names(legacy_vs_calibration)) {
    legacy_vs_calibration <- legacy_vs_calibration %>% filter(model_name=="env_cycl")
  }
} else {
  cat("No common columns found for merging or legacy data is empty\n")
  legacy_vs_calibration <- data.frame()
}

# Create calibration vs refit comparison
if (length(common_cols) > 0) {
  calibration_vs_refit = merge(beta_effects_calibration %>% filter(time_period == "2015-11_2018-01"), 
                               beta_effects_calibration %>% filter(time_period == "2015-11_2020-01"), 
                               by = common_cols, suffixes = c("calibration","refit"))
  
  # Only filter if the column exists
  if ("model_name" %in% names(calibration_vs_refit)) {
    calibration_vs_refit <- calibration_vs_refit %>% filter(model_name=="env_cycl")
  }
} else {
  cat("No common columns found for calibration vs refit comparison\n")
  calibration_vs_refit <- data.frame()
}

# Check if we have data for plotting
if (nrow(legacy_vs_calibration) > 0) {
  # Create legacy vs calibration plot
  p1 <- ggplot(legacy_vs_calibration %>% filter(significantcalibration==1 & significantlegacy==1), 
               aes(x = Meancalibration, y = Meanlegacy, color=beta)) + 
    geom_point() + 
    geom_smooth(method="lm", se=F) +
    geom_abline(slope = 1, intercept = 0, linetype=2) +
    theme_bw() +
    ggtitle("Legacy vs Calibration Effects")
  
  print(p1)
  png(here("figures","effect_consistency_legacy_vs_calibration.png"), width = 1000, height = 800)
  print(p1)
  dev.off()
  cat("Legacy vs calibration plot created\n")
} else {
  cat("No data available for legacy vs calibration comparison\n")
}

if (nrow(calibration_vs_refit) > 0) {
  # Create calibration vs refit plot
  p2 <- ggplot(calibration_vs_refit %>% filter(significantcalibration==1 & significantrefit==1), 
               aes(x = Meancalibration, y = Meanrefit, color=beta)) + 
    geom_point() + 
    geom_smooth(method="lm", se=F) +
    geom_abline(slope = 1, intercept = 0, linetype=2) +
    theme_bw() +
    ggtitle("Calibration vs Refit Effects")
  
  print(p2)
  png(here("figures","effect_consistency_calibration_vs_refit.png"), width = 1000, height = 800)
  print(p2)
  dev.off()
  cat("Calibration vs refit plot created\n")
} else {
  cat("No data available for calibration vs refit comparison\n")
}

# Create env_cov comparisons
if (length(common_cols) > 0 && nrow(beta_effects_legacy) > 0) {
  legacy_vs_calibration_env_cov = merge(beta_effects_legacy, 
                                       beta_effects_calibration %>% filter(time_period == "2015-11_2018-01"), 
                                       by = common_cols, suffixes = c("legacy","calibration"))
  
  # Only filter if the column exists
  if ("model_name" %in% names(legacy_vs_calibration_env_cov)) {
    legacy_vs_calibration_env_cov <- legacy_vs_calibration_env_cov %>% filter(model_name=="env_cov")
  }
} else {
  cat("No common columns found for env_cov legacy vs calibration comparison\n")
  legacy_vs_calibration_env_cov <- data.frame()
}

if (length(common_cols) > 0) {
  calibration_vs_refit_env_cov = merge(beta_effects_calibration %>% filter(time_period == "2015-11_2018-01"), 
                                      beta_effects_calibration %>% filter(time_period == "2015-11_2020-01"), 
                                      by = common_cols, suffixes = c("calibration","refit"))
  
  # Only filter if the column exists
  if ("model_name" %in% names(calibration_vs_refit_env_cov)) {
    calibration_vs_refit_env_cov <- calibration_vs_refit_env_cov %>% filter(model_name=="env_cov")
  }
} else {
  cat("No common columns found for env_cov calibration vs refit comparison\n")
  calibration_vs_refit_env_cov <- data.frame()
}

# Create env_cov plots
if (nrow(legacy_vs_calibration_env_cov) > 0) {
  p3 <- ggplot(legacy_vs_calibration_env_cov %>% filter(significantcalibration==1 & significantlegacy==1), 
               aes(x = Meancalibration, y = Meanlegacy, color=beta)) + 
    geom_point() + 
    geom_smooth(method="lm", se=F) +
    geom_abline(slope = 1, intercept = 0, linetype=2) +
    theme_bw() +
    ggtitle("Legacy vs Calibration Effects (env_cov)")
  
  print(p3)
  png(here("figures","effect_consistency_legacy_vs_calibration_env_cov.png"), width = 1000, height = 800)
  print(p3)
  dev.off()
  cat("env_cov legacy vs calibration plot created\n")
}

if (nrow(calibration_vs_refit_env_cov) > 0) {
  p4 <- ggplot(calibration_vs_refit_env_cov %>% filter(significantcalibration==1 & significantrefit==1), 
               aes(x = Meancalibration, y = Meanrefit, color=beta)) + 
    geom_point() + 
    geom_smooth(method="lm", se=F) +
    geom_abline(slope = 1, intercept = 0, linetype=2) +
    theme_bw() +
    ggtitle("Calibration vs Refit Effects (env_cov)")
  
  print(p4)
  png(here("figures","effect_consistency_calibration_vs_refit_env_cov.png"), width = 1000, height = 800)
  print(p4)
  dev.off()
  cat("env_cov calibration vs refit plot created\n")
}

# Load predictor data for nitrifier analysis
predictor_data <- readRDS(here("data/clean/all_predictor_data.rds"))

# Check if pH data exists
if ("pH" %in% names(predictor_data)) {
  pH_df <- predictor_data$pH %>% as.data.frame() %>% 
    select("pH" = 88) %>% 
    rownames_to_column("plotID")
  
  # Load plot estimates from separate file (plot_est is intentionally empty in summaries)
  if (file.exists(here("data/summary/plot_estimates.rds"))) {
    cat("Loading plot_estimates from separate file for nitrifier analysis...\n")
    calibration_df = readRDS(here("data/summary/plot_estimates.rds"))
    cat("Loaded", nrow(calibration_df), "rows\n")
    
    # Check if truth column exists
    if ("truth" %in% names(calibration_df)) {
      calibration_df = calibration_df %>% filter(!is.na(truth))
      cat("Filtered to", nrow(calibration_df), "rows with truth values\n")
    } else {
      cat("Warning: truth column not found in plot_estimates\n")
      cat("Available columns:", paste(names(calibration_df), collapse = ", "), "\n")
    }
  } else {
    # Fallback to summaries file
    sum.in <- readRDS(here("data/summary/logit_beta_fixed_priors_summaries.rds"))
    if ("plot_est" %in% names(sum.in) && nrow(sum.in$plot_est) > 0) {
      if ("truth" %in% names(sum.in$plot_est)) {
        calibration_df = sum.in$plot_est %>% filter(!is.na(truth))
      } else {
        cat("Warning: truth column not found in plot_est\n")
        calibration_df = sum.in$plot_est
      }
    } else {
      cat("Warning: plot_est not available\n")
      calibration_df = data.frame()
    }
  }
  
  # Merge with pH data (only if we have calibration_df)
  if (exists("calibration_df") && nrow(calibration_df) > 0) {
    calibration_df <- merge(calibration_df, pH_df, by = "plotID", all.x = TRUE)
    
    # Filter for nitrification data
    nitr <- calibration_df %>% filter(taxon=="nitrification")
    
    if (nrow(nitr) > 0) {
      # Create nitrifier pH plot
      p5 <- ggplot(nitr %>% filter(model_name=="env_cycl"), 
                   aes(x = pH, y = truth, color=time_period)) + 
        geom_point(show.legend = FALSE) + 
        geom_smooth(method="lm", se=F, color=1) +
        facet_grid(~time_period) +
        theme_minimal() +
        theme(text = element_text(size = 16),
              axis.title=element_text(size=18),
              strip.text.y = element_text(size=12)) + 
        ylab("Observed nitrifier plot-level abundance") + 
        xlab("pH, mean-scaled") +
        ggtitle("Nitrifier Abundance vs pH")
      
      print(p5)
      png(here("figures","effect_consistency_nitrifier_pH.png"), width = 1200, height = 800)
      print(p5)
      dev.off()
      cat("Nitrifier pH plot created\n")
      
      # Create model comparison plot - use Mean instead of 50% if 50% doesn't exist
      y_col <- if ("50%" %in% names(nitr)) "50%" else if ("Mean" %in% names(nitr)) "Mean" else stop("No y column found")
      p6 <- ggplot(nitr, aes(x = pH, y = .data[[y_col]], color=model_name)) + 
        geom_point() + 
        geom_smooth(se=F) +
        facet_wrap(model_name~time_period, scales="free") +
        theme_bw() +
        ggtitle("Nitrifier Model Predictions vs pH")
      
      print(p6)
      png(here("figures","effect_consistency_nitrifier_model_comparison.png"), width = 1400, height = 1000)
      print(p6)
      dev.off()
      cat("Nitrifier model comparison plot created\n")
    } else {
      cat("No nitrification data available\n")
    }
  } else {
    cat("No calibration data available for nitrifier analysis\n")
  }
} else {
  cat("No pH data available in predictor data\n")
}

