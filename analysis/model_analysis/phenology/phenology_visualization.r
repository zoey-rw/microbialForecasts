# Multiple visualization approaches 
# Options: data_type, aesthetic_level, model_integration, output_format

source("../../source.R")

# Load visualization packages
library(ggplot2)
library(viridis)
library(dplyr)
library(patchwork)
library(ggrepel)
library(tidyr)
library(lubridate)

# Configuration options
phenology_visualization <- function(
  # Data options
  data_type = "individual_cores",           # "individual_cores", "aggregated", "both"
  aesthetic_level = "enhanced",             # "basic", "enhanced", "publication"
  model_integration = "improved",           # "improved", "standard", "none"
  output_format = "both",                   # "individual", "combined", "both"
  
  # Species and time period options
  species_list = NULL,                      # NULL = all available species
  time_periods = NULL,                      # NULL = all available periods
  
  # Output options
  output_dir = NULL,                        # NULL = default figures/phenology_comparison
  plot_width = 12,
  plot_height = 8,
  dpi = 300
) {
  
  cat("Configuration:\n")
  cat("  Data type:", data_type, "\n")
  cat("  Aesthetic level:", aesthetic_level, "\n")
  cat("  Model integration:", model_integration, "\n")
  cat("  Output format:", output_format, "\n\n")
  
  # Set default species list if not provided
  if (is.null(species_list)) {
    species_list <- c("saprotroph", "ectomycorrhizal", "endophyte", "plant_pathogen", 
                      "animal_pathogen", "lichenized", "wood_saprotroph", 
                      "soil_saprotroph", "litter_saprotroph", "cellulolytic", 
                      "assim_nitrite_reduction")
  }
  
  # Set default time periods if not provided
  if (is.null(time_periods)) {
    time_periods <- list(
      "2013-2015" = list(min = "20130601", max = "20151101"),
      "2015-2018" = list(min = "20150101", max = "20180101"),
      "2015-2020" = list(min = "20150101", max = "20200101")
    )
  }
  
  # Set default output directory
  if (is.null(output_dir)) {
    output_dir <- here("figures/phenology_comparison")
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Function to get core-level observed data
  get_core_data <- function(species, min_date, max_date) {
    tryCatch({
      # Load appropriate data
      if (species %in% c("saprotroph", "ectomycorrhizal", "endophyte", "plant_pathogen", 
                         "animal_pathogen", "lichenized", "wood_saprotroph", 
                         "soil_saprotroph", "litter_saprotroph", "cellulolytic", 
                         "assim_nitrite_reduction")) {
        fungi <- readRDS(here("data/clean/groupAbundances_ITS_2023.rds"))
        rank.df <- fungi[[species]]
      } else {
        bacteria <- readRDS(here("data/clean/groupAbundances_16S_2023.rds"))
        rank.df <- bacteria[[species]]
      }
      
      if (is.null(rank.df)) return(NULL)
      
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
      cat("  Error processing", species, ":", e$message, "\n")
      return(NULL)
    })
  }
  
  # Function to extract seasonal parameters
  extract_seasonal_params <- function(file_path, model_type) {
    tryCatch({
      model_data <- readRDS(file_path)
      
      if (!"samples" %in% names(model_data)) return(NULL)
      
      samples <- model_data$samples
      if (length(samples) < 2) return(NULL)
      
      # Extract species name
      if (model_type == "CLR") {
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
      
      return(list(
        species = species,
        sin_coef = sin_coef,
        cos_coef = cos_coef,
        amplitude = amplitude,
        peak_timing = peak_timing,
        peak_month = round(peak_timing, 1)
      ))
      
    }, error = function(e) {
      return(NULL)
    })
  }
  
  # Function to generate predicted trends
  generate_predicted_trends <- function(sin_coef, cos_coef, months = 1:12) {
    seasonal_effect <- sin_coef * sin(2 * pi * months / 12) + 
                      cos_coef * cos(2 * pi * months / 12)
    return(seasonal_effect)
  }
  
  # Function to create phenology plot
  create_phenology_plot <- function(species, time_periods) {
    cat("  Creating plot for", species, "\n")
    
    all_data <- list()
    
    for (period in names(time_periods)) {
      min_date <- time_periods[[period]]$min
      max_date <- time_periods[[period]]$max
      
      # Get core data
      core_data <- get_core_data(species, min_date, max_date)
      
      if (!is.null(core_data) && nrow(core_data) > 0) {
        # Look for model files based on integration type
        model_file <- NULL
        model_type <- "Logit"
        
        if (model_integration != "none") {
          if (period == "2013-2015") {
            if (model_integration == "improved") {
              improved_file <- here("data/model_outputs/logit_phenology_improved_2013_2015", 
                                   paste0("samples_cycl_only_", species, "_20130601_20151101.rds"))
              if (file.exists(improved_file)) {
                model_file <- improved_file
                cat("      Found improved model file\n")
              }
            } else {
              standard_file <- here("data/model_outputs/functional_groups/cycl_only", 
                                   paste0("samples_cycl_only_", species, "_20130601_20151101.rds"))
              if (file.exists(standard_file)) {
                model_file <- standard_file
                cat("      Found standard model file\n")
              }
            }
          } else if (period %in% c("2015-2018", "2015-2020")) {
            model_file <- here("data/model_outputs/functional_groups/cycl_only", 
                              paste0("samples_cycl_only_", species, "_20150101_", 
                                    ifelse(period == "2015-2018", "20180101", "20200101"), ".rds"))
            if (file.exists(model_file)) {
              cat("      Found", period, "model file\n")
            }
          }
        }
        
        # Prepare data based on data_type
        if (data_type == "individual_cores") {
          # Keep individual core data
          period_data <- core_data %>%
            mutate(
              species = species,
              time_period = period,
              data_type = "Real Core Data"
            )
        } else if (data_type == "aggregated") {
          # Aggregate by month
          period_data <- core_data %>%
            group_by(month) %>%
            summarise(
              mean_abundance = mean(abundance, na.rm = TRUE),
              se_abundance = sd(abundance, na.rm = TRUE) / sqrt(n()),
              n_obs = n(),
              .groups = 'drop'
            ) %>%
            mutate(
              month_name = month.abb[month],
              species = species,
              time_period = period,
              data_type = "Aggregated Data"
            )
        } else {
          # Both individual and aggregated
          individual_data <- core_data %>%
            mutate(
              species = species,
              time_period = period,
              data_type = "Real Core Data"
            )
          
          aggregated_data <- core_data %>%
            group_by(month) %>%
            summarise(
              mean_abundance = mean(abundance, na.rm = TRUE),
              se_abundance = sd(abundance, na.rm = TRUE) / sqrt(n()),
              n_obs = n(),
              .groups = 'drop'
            ) %>%
            mutate(
              month_name = month.abb[month],
              species = species,
              time_period = period,
              data_type = "Aggregated Data"
            )
          
          period_data <- list(
            individual = individual_data,
            aggregated = aggregated_data
          )
        }
        
        # Add predicted trends if model file exists
        if (!is.null(model_file) && file.exists(model_file)) {
          params <- extract_seasonal_params(model_file, model_type)
          
          if (!is.null(params)) {
            months <- 1:12
            predicted_trends <- generate_predicted_trends(params$sin_coef, params$cos_coef, months)
            
            predicted_data <- data.frame(
              month = months,
              month_name = month.abb[months],
              abundance = predicted_trends,
              species = species,
              time_period = period,
              data_type = "Model Prediction",
              amplitude = params$amplitude,
              peak_month = params$peak_month
            )
            
            if (data_type == "both") {
              period_data$predicted <- predicted_data
            } else {
              period_data <- list(
                data = period_data,
                predicted = predicted_data
              )
            }
          }
        }
        
        all_data[[period]] <- period_data
      }
    }
    
    if (length(all_data) == 0) {
      cat("    No data found for", species, "\n")
      return(NULL)
    }
    
    # Combine data based on data_type
    if (data_type == "both") {
      all_individual <- do.call(rbind, lapply(all_data, function(x) x$individual))
      all_aggregated <- do.call(rbind, lapply(all_data, function(x) x$aggregated))
      all_predicted <- do.call(rbind, lapply(all_data, function(x) x$predicted))
    } else {
      all_data_combined <- do.call(rbind, lapply(all_data, function(x) {
        if (is.list(x) && "data" %in% names(x)) x$data else x
      }))
      all_predicted <- do.call(rbind, lapply(all_data, function(x) {
        if (is.list(x) && "predicted" %in% names(x)) x$predicted else NULL
      }))
    }
    
    # Create the plot
    p <- ggplot() +
      # Add seasonal background based on aesthetic level
      if (aesthetic_level %in% c("enhanced", "publication")) {
        list(
          annotate("rect", xmin = 5.5, xmax = 8.5, ymin = -Inf, ymax = Inf, 
                   alpha = 0.1, fill = "gold"),
          annotate("rect", xmin = 8.5, xmax = 11.5, ymin = -Inf, ymax = Inf, 
                   alpha = 0.1, fill = "orange"),
          annotate("rect", xmin = 11.5, xmax = 2.5, ymin = -Inf, ymax = Inf, 
                   alpha = 0.1, fill = "lightblue"),
          annotate("rect", xmin = 2.5, xmax = 5.5, ymin = -Inf, ymax = Inf, 
                   alpha = 0.1, fill = "lightgreen"),
          
          # Seasonal labels
          annotate("text", x = 7, y = max(all_data_combined$abundance, na.rm = TRUE) * 0.8, 
                   label = "Summer", color = "darkgoldenrod", fontface = "bold"),
          annotate("text", x = 10, y = max(all_data_combined$abundance, na.rm = TRUE) * 0.8, 
                   label = "Fall", color = "darkorange", fontface = "bold"),
          annotate("text", x = 1, y = max(all_data_combined$abundance, na.rm = TRUE) * 0.8, 
                   label = "Winter", color = "darkblue", fontface = "bold"),
          annotate("text", x = 4, y = max(all_data_combined$abundance, na.rm = TRUE) * 0.8, 
                   label = "Spring", color = "darkgreen", fontface = "bold")
        )
      } else {
        NULL
      } +
      
      # Add predicted trends
      if (!is.null(all_predicted) && nrow(all_predicted) > 0) {
        geom_line(data = all_predicted, 
                  aes(x = month, y = abundance, color = time_period), 
                  linewidth = 1.2, alpha = 0.8)
      } +
      
      # Add data points based on data_type
      if (data_type == "individual_cores") {
        geom_point(data = all_data_combined, 
                   aes(x = month, y = abundance, color = time_period), 
                   size = 1.5, alpha = 0.6, position = position_jitter(width = 0.2))
      } else if (data_type == "aggregated") {
        geom_point(data = all_data_combined, 
                   aes(x = month, y = mean_abundance, color = time_period), 
                   size = 2.5, alpha = 0.8) +
        geom_errorbar(data = all_data_combined, 
                      aes(x = month, y = mean_abundance, 
                          ymin = mean_abundance - se_abundance, 
                          ymax = mean_abundance + se_abundance, 
                          color = time_period), 
                      width = 0.3, alpha = 0.6)
      } else {
        # Both individual and aggregated
        geom_point(data = all_individual, 
                   aes(x = month, y = abundance, color = time_period), 
                   size = 1, alpha = 0.4, position = position_jitter(width = 0.2)) +
        geom_point(data = all_aggregated, 
                   aes(x = month, y = mean_abundance, color = time_period), 
                   size = 3, alpha = 0.8) +
        geom_errorbar(data = all_aggregated, 
                      aes(x = month, y = mean_abundance, 
                          ymin = mean_abundance - se_abundance, 
                          ymax = mean_abundance + se_abundance, 
                          color = time_period), 
                      width = 0.3, alpha = 0.6)
      } +
      
      # Enhanced styling
      scale_color_viridis(discrete = TRUE, option = "plasma") +
      scale_x_continuous(breaks = 1:12, 
                        labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
      
      # Labels and theme
      labs(title = paste("Seasonal Phenology:", species),
           subtitle = paste("Data type:", data_type, "| Aesthetics:", aesthetic_level, "| Models:", model_integration),
           x = "Month",
           y = "Abundance (proportion)",
           color = "Time Period") +
      
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            plot.subtitle = element_text(hjust = 0.5, size = 11),
            legend.position = "bottom",
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5))
    
    return(p)
  }
  
  # Create plots for each species
  all_plots <- list()
  for (species in species_list) {
    cat("\n--- Processing", species, "---\n")
    plot_obj <- create_phenology_plot(species, time_periods)
    if (!is.null(plot_obj)) {
      all_plots[[species]] <- plot_obj
    }
  }
  
  # Save plots based on output_format
  if (output_format %in% c("individual", "both")) {
    cat("\n💾 Saving individual species plots...\n")
    for (species in names(all_plots)) {
      plot_file <- file.path(output_dir, paste0("phenology_", species, ".png"))
      ggsave(plot_file, all_plots[[species]], width = plot_width, height = plot_height, dpi = dpi)
      cat("  ✅ Saved", species, "plot\n")
    }
  }
  
  if (output_format %in% c("combined", "both") && length(all_plots) > 1) {
    cat("\n🎯 Creating combined phenology plot...\n")
    
    combined_plot <- wrap_plots(all_plots, ncol = 2) +
      plot_annotation(
        title = "Complete Microbial Phenology Analysis",
        subtitle = paste("Data:", data_type, "| Aesthetics:", aesthetic_level, "| Models:", model_integration),
        caption = paste("Configuration: Data type =", data_type, "| Aesthetics =", aesthetic_level, "| Models =", model_integration),
        theme = theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
          plot.subtitle = element_text(hjust = 0.5, size = 14),
          plot.caption = element_text(hjust = 0.5, size = 10, color = "gray50")
        )
      )
    
    combined_file <- file.path(output_dir, "combined_phenology_all_species.png")
    ggsave(combined_file, combined_plot, width = 20, height = 16, dpi = dpi)
    cat("  ✅ Saved combined phenology plot\n")
  }
  
  # Summary
  cat("\n🎯 PHENOLOGY VISUALIZATION COMPLETE\n")
  cat("=", paste(rep("=", 60), collapse = ""), "\n")
  
  cat("Configuration used:\n")
  cat("  Data type:", data_type, "\n")
  cat("  Aesthetic level:", aesthetic_level, "\n")
  cat("  Model integration:", model_integration, "\n")
  cat("  Output format:", output_format, "\n")
  
  cat("\nResults:\n")
  cat("  Species processed:", length(all_plots), "\n")
  cat("  Time periods:", length(time_periods), "\n")
  cat("  Output directory:", output_dir, "\n")
  
  cat("\n=== MASTER PHENOLOGY VISUALIZATION COMPLETE ===\n")
  
  return(all_plots)
}

# Example usage:
# phenology_visualization(
#   data_type = "individual_cores",
#   aesthetic_level = "enhanced", 
#   model_integration = "improved",
#   output_format = "both"
# )
