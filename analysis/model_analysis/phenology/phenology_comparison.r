# Phenology comparisons: how do conclusions differ between CLR or logit beta models?
# Options: comparison_types, model_types, output_formats

source("../../source.R")


# Load required packages
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(viridis)
library(patchwork)
library(corrplot)

# Configuration options
phenology_comparison <- function(
  # Comparison options
  comparison_types = "both",                # "clr_vs_logit", "time_periods", "both"
  model_types = "both",                     # "logit", "clr", "both"
  output_formats = "both",                  # "rds", "csv", "both"
  
  # Input/output options
  input_dir = NULL,                         # NULL = auto-detect
  output_dir = NULL,                        # NULL = auto-detect
  species_list = NULL,                      # NULL = auto-detect
  
  # Analysis options
  include_convergence = TRUE,               # Include convergence metrics in comparison
  create_correlation_plots = TRUE,          # Create correlation plots
  create_heatmaps = TRUE,                   # Create parameter heatmaps
  
  # Visualization options
  plot_theme = "minimal",                   # "minimal", "classic", "bw"
  color_palette = "viridis"                 # "viridis", "plasma", "magma"
) {
  
  cat("🎯 Master Phenology Comparison\n")
  cat("Configuration:\n")
  cat("  Comparison types:", comparison_types, "\n")
  cat("  Model types:", model_types, "\n")
  cat("  Output formats:", output_formats, "\n")
  cat("  Include convergence:", include_convergence, "\n\n")
  
  # Auto-detect input directories
  if (is.null(input_dir)) {
    input_dirs <- list()
    if (model_types %in% c("logit", "both")) {
      # Check for improved models first
      improved_dir <- here("data/model_outputs/logit_phenology_improved_2013_2015")
      if (dir.exists(improved_dir)) {
        input_dirs$logit_improved <- improved_dir
      }
      
      # Check for standard models
      standard_dir <- here("data/model_outputs/functional_groups/cycl_only")
      if (dir.exists(standard_dir)) {
        input_dirs$logit_standard <- standard_dir
      }
    }
    if (model_types %in% c("clr", "both")) {
      clr_dir <- here("data/model_outputs/CLR_phenology_2013_2015")
      if (dir.exists(clr_dir)) {
        input_dirs$clr <- clr_dir
      }
    }
  } else {
    input_dirs <- list(logit = input_dir, clr = input_dir)
  }
  
  # Auto-detect output directory
  if (is.null(output_dir)) {
    output_dir <- here("data/model_outputs/phenology_comparison")
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Auto-detect species list
  if (is.null(species_list)) {
    species_list <- c("saprotroph", "ectomycorrhizal", "endophyte", "plant_pathogen", 
                      "animal_pathogen", "lichenized", "wood_saprotroph", 
                      "soil_saprotroph", "litter_saprotroph", "cellulolytic", 
                      "assim_nitrite_reduction")
  }
  
  # Function to extract seasonal parameters
  extract_seasonal_params <- function(file_path, model_type) {
    tryCatch({
      model_data <- readRDS(file_path)
      
      if (!"samples" %in% names(model_data)) return(NULL)
      
      samples <- model_data$samples
      if (length(samples) < 1) return(NULL)
      
      # Extract species name
      if (model_type == "clr") {
        species <- str_extract(basename(file_path), "CLR_phenology_(.+)_\\d{8}_\\d{8}")
        species <- str_remove(species, "CLR_phenology_")
        species <- str_remove(species, "_\\d{8}_\\d{8}")
        sin_param <- "beta[1, 1]"
        cos_param <- "beta[2, 1]"
      } else {
        species <- str_extract(basename(file_path), "samples_cycl_only_(.+)_\\d{8}_\\d{8}")
        species <- str_remove(species, "samples_cycl_only_")
        species <- str_remove(species, "_\\d{8}_\\d{8}")
        sin_param <- "beta[1]"
        cos_param <- "beta[2]"
      }
      
      param_names <- colnames(samples[[1]])
      if (!all(c(sin_param, cos_param) %in% param_names)) return(NULL)
      
      # Calculate seasonal parameters
      sin_samples <- lapply(samples, function(x) x[, sin_param])
      cos_samples <- lapply(samples, function(x) x[, cos_param])
      
      sin_coef <- mean(unlist(sin_samples))
      cos_coef <- mean(unlist(cos_samples))
      amplitude <- sqrt(sin_coef^2 + cos_coef^2)
      peak_timing <- atan2(cos_coef, sin_coef) * 12/(2*pi)
      if (peak_timing < 0) peak_timing <- peak_timing + 12
      
      # Calculate convergence diagnostics if requested
      convergence_metrics <- NULL
      if (include_convergence && length(samples) >= 2) {
        tryCatch({
          # Gelman-Rubin diagnostic
          mcmc_objects <- lapply(samples, function(x) {
            if (is.matrix(x)) coda::mcmc(x) else NULL
          })
          mcmc_objects <- mcmc_objects[!sapply(mcmc_objects, is.null)]
          
          if (length(mcmc_objects) >= 2) {
            gelman_result <- coda::gelman.diag(mcmc_objects)
            ess_result <- coda::effectiveSize(mcmc_objects)
            
            convergence_metrics <- list(
              gelman_rubin = max(gelman_result$psrf[, "Upper C.I."], na.rm = TRUE),
              min_ess = min(ess_result, na.rm = TRUE),
              mean_ess = mean(ess_result, na.rm = TRUE)
            )
          }
        }, error = function(e) {
          convergence_metrics <- list(
            gelman_rubin = NA,
            min_ess = NA,
            mean_ess = NA
          )
        })
      }
      
      return(list(
        species = species,
        model_type = model_type,
        sin_coef = sin_coef,
        cos_coef = cos_coef,
        amplitude = amplitude,
        peak_timing = peak_timing,
        peak_month = round(peak_timing, 1),
        convergence_metrics = convergence_metrics,
        n_chains = length(samples),
        n_samples = nrow(samples[[1]])
      ))
      
    }, error = function(e) {
      return(NULL)
    })
  }
  
  # Function to compare CLR vs Logit Beta approaches
  compare_clr_vs_logit <- function() {
    cat("🔄 Comparing CLR vs Logit Beta approaches...\n")
    
    # Extract parameters from both approaches
    clr_params <- list()
    logit_params <- list()
    
    # Process CLR models
    if ("clr" %in% names(input_dirs)) {
      clr_files <- list.files(input_dirs$clr, 
                             pattern = "CLR_phenology_.*\\.rds$",
                             full.names = TRUE)
      
      for (file in clr_files) {
        params <- extract_seasonal_params(file, "clr")
        if (!is.null(params)) {
          clr_params[[params$species]] <- params
        }
      }
    }
    
    # Process Logit Beta models
    if (any(str_detect(names(input_dirs), "logit"))) {
      for (logit_type in names(input_dirs)[str_detect(names(input_dirs), "logit")]) {
        logit_files <- list.files(input_dirs[[logit_type]], 
                                 pattern = "samples_cycl_only_.*\\.rds$",
                                 full.names = TRUE)
        
        for (file in logit_files) {
          params <- extract_seasonal_params(file, "logit")
          if (!is.null(params)) {
            logit_params[[params$species]] <- params
          }
        }
      }
    }
    
    # Find common species
    common_species <- intersect(names(clr_params), names(logit_params))
    
    if (length(common_species) == 0) {
      cat("  ❌ No common species found for comparison\n")
      return(NULL)
    }
    
    cat("  📊 Found", length(common_species), "common species for comparison\n")
    
    # Create comparison dataframe
    comparison_data <- data.frame()
    
    for (species in common_species) {
      clr_data <- clr_params[[species]]
      logit_data <- logit_params[[species]]
      
      # CLR data
      clr_row <- data.frame(
        species = species,
        approach = "CLR",
        group_type = "bacterial_phylum",
        sin_coef = clr_data$sin_coef,
        cos_coef = clr_data$cos_coef,
        amplitude = clr_data$amplitude,
        peak_timing = clr_data$peak_timing,
        peak_month = clr_data$peak_month,
        n_chains = clr_data$n_chains,
        n_samples = clr_data$n_samples,
        stringsAsFactors = FALSE
      )
      
      # Logit data
      logit_row <- data.frame(
        species = species,
        approach = "Logit Beta",
        group_type = "fungal_functional_group",
        sin_coef = logit_data$sin_coef,
        cos_coef = logit_data$cos_coef,
        amplitude = logit_data$amplitude,
        peak_timing = logit_data$peak_timing,
        peak_month = logit_data$peak_month,
        n_chains = logit_data$n_chains,
        n_samples = logit_data$n_samples,
        stringsAsFactors = FALSE
      )
      
      # Add convergence metrics if available
      if (include_convergence) {
        if (!is.null(clr_data$convergence_metrics)) {
          clr_row$gelman_rubin <- clr_data$convergence_metrics$gelman_rubin
          clr_row$min_ess <- clr_data$convergence_metrics$min_ess
          clr_row$mean_ess <- clr_data$convergence_metrics$mean_ess
        } else {
          clr_row$gelman_rubin <- NA
          clr_row$min_ess <- NA
          clr_row$mean_ess <- NA
        }
        
        if (!is.null(logit_data$convergence_metrics)) {
          logit_row$gelman_rubin <- logit_data$convergence_metrics$gelman_rubin
          logit_row$min_ess <- logit_data$convergence_metrics$min_ess
          logit_row$mean_ess <- logit_data$convergence_metrics$mean_ess
        } else {
          logit_row$gelman_rubin <- NA
          logit_row$min_ess <- NA
          logit_row$mean_ess <- NA
        }
      }
      
      comparison_data <- rbind(comparison_data, clr_row, logit_row)
    }
    
    return(comparison_data)
  }
  
  # Function to compare time periods
  compare_time_periods <- function() {
    cat("🔄 Comparing time periods...\n")
    
    # This would require data from multiple time periods
    # For now, we'll focus on the 2013-2015 period
    cat("  ⚠️ Time period comparison requires multiple time period data\n")
    cat("  📊 Currently only 2013-2015 data available\n")
    
    return(NULL)
  }
  
  # Function to create correlation plots
  create_correlation_plots <- function(comparison_data) {
    if (!create_correlation_plots || is.null(comparison_data)) return(NULL)
    
    cat("  🎨 Creating correlation plots...\n")
    
    # Separate CLR and Logit data
    clr_data <- comparison_data[comparison_data$approach == "CLR", ]
    logit_data <- comparison_data[comparison_data$approach == "Logit Beta", ]
    
    if (nrow(clr_data) == 0 || nrow(logit_data) == 0) {
      cat("    ⚠️ Insufficient data for correlation plots\n")
      return(NULL)
    }
    
    # Create correlation matrix
    numeric_cols <- c("sin_coef", "cos_coef", "amplitude", "peak_timing")
    correlation_matrix <- cor(clr_data[, numeric_cols], logit_data[, numeric_cols], 
                            use = "complete.obs")
    
    # Plot 1: Correlation heatmap
    p1 <- ggplot(data = reshape2::melt(correlation_matrix), 
                 aes(x = Var1, y = Var2, fill = value)) +
      geom_tile() +
      scale_fill_gradient2(low = "red", mid = "white", high = "blue", 
                          midpoint = 0, limits = c(-1, 1)) +
      labs(title = "CLR vs Logit Beta Parameter Correlations",
           x = "CLR Parameters",
           y = "Logit Beta Parameters",
           fill = "Correlation") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5, face = "bold"))
    
    # Plot 2: Amplitude comparison
    p2 <- ggplot() +
      geom_point(data = clr_data, aes(x = reorder(species, amplitude), y = amplitude, 
                                     color = "CLR"), size = 3, alpha = 0.7) +
      geom_point(data = logit_data, aes(x = reorder(species, amplitude), y = amplitude, 
                                       color = "Logit Beta"), size = 3, alpha = 0.7) +
      scale_color_manual(values = c("CLR" = "blue", "Logit Beta" = "red")) +
      labs(title = "Seasonal Amplitude Comparison",
           x = "Species",
           y = "Amplitude",
           color = "Approach") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5, face = "bold"))
    
    # Plot 3: Peak timing comparison
    p3 <- ggplot() +
      geom_point(data = clr_data, aes(x = reorder(species, peak_month), y = peak_month, 
                                     color = "CLR"), size = 3, alpha = 0.7) +
      geom_point(data = logit_data, aes(x = reorder(species, peak_month), y = peak_month, 
                                       color = "Logit Beta"), size = 3, alpha = 0.7) +
      scale_color_manual(values = c("CLR" = "blue", "Logit Beta" = "red")) +
      scale_y_continuous(breaks = 1:12, 
                        labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
      labs(title = "Peak Timing Comparison",
           x = "Species",
           y = "Peak Month",
           color = "Approach") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5, face = "bold"))
    
    # Save plots
    ggsave(file.path(output_dir, "clr_vs_logit_correlations.png"), p1, 
           width = 10, height = 8, dpi = 300)
    ggsave(file.path(output_dir, "amplitude_comparison.png"), p2, 
           width = 12, height = 8, dpi = 300)
    ggsave(file.path(output_dir, "peak_timing_comparison.png"), p3, 
           width = 12, height = 8, dpi = 300)
    
    # Combined plot
    combined_plot <- (p1 + p2) / p3 +
      plot_annotation(
        title = "CLR vs Logit Beta Phenology Comparison",
        subtitle = paste("Species:", length(unique(comparison_data$species))),
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
      )
    
    ggsave(file.path(output_dir, "combined_clr_vs_logit_comparison.png"), combined_plot, 
           width = 20, height = 16, dpi = 300)
    
    cat("    ✅ Saved correlation plots\n")
    
    return(list(p1, p2, p3, combined_plot))
  }
  
  # Function to create parameter heatmaps
  create_parameter_heatmaps <- function(comparison_data) {
    if (!create_heatmaps || is.null(comparison_data)) return(NULL)
    
    cat("  🎨 Creating parameter heatmaps...\n")
    
    # Create heatmap for sin coefficients
    p1 <- ggplot(comparison_data, aes(x = approach, y = reorder(species, sin_coef), fill = sin_coef)) +
      geom_tile() +
      scale_fill_gradient2(low = "red", mid = "white", high = "blue", 
                          midpoint = 0, name = "Sin Coefficient") +
      labs(title = "Sin Coefficient Comparison",
           x = "Approach",
           y = "Species") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    # Create heatmap for cos coefficients
    p2 <- ggplot(comparison_data, aes(x = approach, y = reorder(species, cos_coef), fill = cos_coef)) +
      geom_tile() +
      scale_fill_gradient2(low = "red", mid = "white", high = "blue", 
                          midpoint = 0, name = "Cos Coefficient") +
      labs(title = "Cos Coefficient Comparison",
           x = "Approach",
           y = "Species") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    # Create heatmap for amplitudes
    p3 <- ggplot(comparison_data, aes(x = approach, y = reorder(species, amplitude), fill = amplitude)) +
      geom_tile() +
      scale_fill_viridis(option = "plasma", name = "Amplitude") +
      labs(title = "Amplitude Comparison",
           x = "Approach",
           y = "Species") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    # Save plots
    ggsave(file.path(output_dir, "sin_coefficient_heatmap.png"), p1, 
           width = 10, height = 8, dpi = 300)
    ggsave(file.path(output_dir, "cos_coefficient_heatmap.png"), p2, 
           width = 10, height = 8, dpi = 300)
    ggsave(file.path(output_dir, "amplitude_heatmap.png"), p3, 
           width = 10, height = 8, dpi = 300)
    
    # Combined heatmap
    combined_heatmap <- (p1 + p2) / p3 +
      plot_annotation(
        title = "Parameter Comparison Heatmaps",
        subtitle = paste("Species:", length(unique(comparison_data$species))),
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
      )
    
    ggsave(file.path(output_dir, "combined_parameter_heatmaps.png"), combined_heatmap, 
           width = 20, height = 16, dpi = 300)
    
    cat("    ✅ Saved parameter heatmaps\n")
    
    return(list(p1, p2, p3, combined_heatmap))
  }
  
  # Main processing
  all_results <- list()
  
  # CLR vs Logit Beta comparison
  if (comparison_types %in% c("clr_vs_logit", "both")) {
    clr_vs_logit_data <- compare_clr_vs_logit()
    if (!is.null(clr_vs_logit_data)) {
      all_results$clr_vs_logit <- clr_vs_logit_data
      
      # Save comparison data
      if (output_formats %in% c("rds", "both")) {
        rds_file <- file.path(output_dir, "clr_vs_logit_comparison.rds")
        saveRDS(clr_vs_logit_data, rds_file)
        cat("  💾 Saved CLR vs Logit comparison to:", basename(rds_file), "\n")
      }
      
      if (output_formats %in% c("csv", "both")) {
        csv_file <- file.path(output_dir, "clr_vs_logit_comparison.csv")
        write.csv(clr_vs_logit_data, csv_file, row.names = FALSE)
        cat("  💾 Saved CLR vs Logit comparison CSV to:", basename(csv_file), "\n")
      }
      
      # Create visualizations
      correlation_plots <- create_correlation_plots(clr_vs_logit_data)
      parameter_heatmaps <- create_parameter_heatmaps(clr_vs_logit_data)
      
      all_results$correlation_plots <- correlation_plots
      all_results$parameter_heatmaps <- parameter_heatmaps
    }
  }
  
  # Time period comparison
  if (comparison_types %in% c("time_periods", "both")) {
    time_period_data <- compare_time_periods()
    if (!is.null(time_period_data)) {
      all_results$time_periods <- time_period_data
    }
  }
  
  # Summary
  cat("\n🎯 PHENOLOGY COMPARISON COMPLETE\n")
  cat("=", paste(rep("=", 60), collapse = ""), "\n")
  
  cat("Configuration used:\n")
  cat("  Comparison types:", comparison_types, "\n")
  cat("  Model types:", model_types, "\n")
  cat("  Output formats:", output_formats, "\n")
  cat("  Include convergence:", include_convergence, "\n")
  
  cat("\nResults:\n")
  cat("  Input directories:", length(input_dirs), "\n")
  cat("  Output directory:", output_dir, "\n")
  cat("  Comparisons completed:", length(all_results), "\n")
  
  if ("clr_vs_logit" %in% names(all_results)) {
    cat("  CLR vs Logit species:", length(unique(all_results$clr_vs_logit$species)), "\n")
  }
  
  cat("\n=== MASTER PHENOLOGY COMPARISON COMPLETE ===\n")
  
  return(all_results)
}

# Example usage:
# phenology_comparison(
#   comparison_types = "both",
#   model_types = "both",
#   output_formats = "both"
# )
