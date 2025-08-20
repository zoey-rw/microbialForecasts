# Seasonality consistency

# Evaluates consistency of seasonality trends across different models and approaches
# Shows observed core data and estimated trends from logit beta regression vs CLR models

source("../../source.R")


# Load required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(viridis)
library(lubridate)

# Function to extract seasonal parameters from logit beta regression models
extract_logit_seasonal_params <- function(file_path) {
  tryCatch({
    # Check if this is a summary file or MCMC samples file
    if (grepl("summary_", basename(file_path))) {
      # Summary file structure (2015-2018 models)
      summary_data <- readRDS(file_path)
      
      # The summary file has 3 elements, first element contains the beta coefficients
      beta_data <- summary_data[[1]]
      
      # Find rows with sin and cos coefficients
      sin_row <- which(beta_data$beta == "sin")
      cos_row <- which(beta_data$beta == "cos")
      
      if (length(sin_row) > 0 && length(cos_row) > 0) {
        sin_coef <- beta_data$Mean[sin_row[1]]
        cos_coef <- beta_data$Mean[cos_row[1]]
        
        # Calculate amplitude and peak timing
        amplitude <- sqrt(sin_coef^2 + cos_coef^2)
        peak_timing <- atan2(cos_coef, sin_coef) * 12 / (2 * pi)
        if (peak_timing < 0) peak_timing <- peak_timing + 12
        
        # Extract group name and model type from filename
        filename <- basename(file_path)
        
        # Parse filename to get model type and group name
        if (grepl("summary_cycl_only_", filename)) {
          model_type <- "cycl_only"
          group_name <- str_extract(filename, "summary_cycl_only_(.+?)_2015")
          group_name <- str_remove(group_name, "summary_cycl_only_")
          group_name <- str_remove(group_name, "_2015.*")
        } else if (grepl("summary_env_cov_", filename)) {
          model_type <- "env_cov"
          group_name <- str_extract(filename, "summary_env_cov_(.+?)_2015")
          group_name <- str_remove(group_name, "summary_env_cov_")
          group_name <- str_remove(group_name, "_2015.*")
        } else if (grepl("summary_env_cycl_", filename)) {
          model_type <- "env_cycl"
          group_name <- str_extract(filename, "summary_env_cycl_(.+?)_2015")
          group_name <- str_remove(group_name, "summary_env_cycl_")
          group_name <- str_remove(group_name, "_2015.*")
        } else {
          model_type <- "unknown"
          group_name <- "unknown"
        }
        
        return(list(
          sin_coef = sin_coef,
          cos_coef = cos_coef,
          amplitude = amplitude,
          peak_timing = peak_timing,
          model_type = model_type,
          group_name = group_name,
          time_period = "2015-2018"
        ))
      } else {
        cat("    No sin/cos coefficients found in", basename(file_path), "\n")
        return(NULL)
      }
    } else if (grepl("samples_", basename(file_path))) {
      # MCMC samples file structure (2013-2015 models)
      samples_data <- readRDS(file_path)
      
      # The samples are in samples_data$samples which is a list of chains
      if ("samples" %in% names(samples_data) && is.list(samples_data$samples)) {
        # Combine all chains and extract beta parameters
        all_chains <- do.call(rbind, samples_data$samples)
        
        # Extract sin and cos coefficients from MCMC samples
        sin_coef <- mean(all_chains[, "beta[1]"], na.rm = TRUE)
        cos_coef <- mean(all_chains[, "beta[2]"], na.rm = TRUE)
        
        # Calculate amplitude and peak timing
        amplitude <- sqrt(sin_coef^2 + cos_coef^2)
        peak_timing <- atan2(cos_coef, sin_coef) * 12 / (2 * pi)
        if (peak_timing < 0) peak_timing <- peak_timing + 12
        
        # Extract group name from filename
        filename <- basename(file_path)
        group_name <- str_extract(filename, "samples_cycl_only_(.+?)_2013")
        group_name <- str_remove(group_name, "samples_cycl_only_")
        group_name <- str_remove(group_name, "_2013.*")
        
        return(list(
          sin_coef = sin_coef,
          cos_coef = cos_coef,
          amplitude = amplitude,
          peak_timing = peak_timing,
          model_type = "cycl_only",
          group_name = group_name,
          time_period = "2013-2015"
        ))
      } else {
        cat("    Invalid MCMC samples structure in", basename(file_path), "\n")
        return(NULL)
      }
    } else {
      cat("    Unknown file type:", basename(file_path), "\n")
      return(NULL)
    }
  }, error = function(e) {
    cat("  Error extracting parameters from", file_path, ":", e$message, "\n")
    return(NULL)
  })
}

# Function to extract seasonal parameters from CLR models
extract_clr_seasonal_params <- function(file_path) {
  tryCatch({
    samples_data <- readRDS(file_path)
    
    # Extract sin and cos coefficients from MCMC samples
    sin_coef <- mean(samples_data$samples[, "beta[1]"])
    cos_coef <- mean(samples_data$samples[, "beta[2]"])
    
    # Calculate amplitude and peak timing
    amplitude <- sqrt(sin_coef^2 + cos_coef^2)
    peak_timing <- atan2(cos_coef, sin_coef) * 12 / (2 * pi)
    if (peak_timing < 0) peak_timing <- peak_timing + 12
    
    # Extract group name from filename
    filename <- basename(file_path)
    group_name <- str_extract(filename, "samples_CLR_final_cycl_only_(.+?)_2015")
    group_name <- str_remove(group_name, "samples_CLR_final_cycl_only_")
    group_name <- str_remove(group_name, "_2015.*")
    
    # Check if group name is valid
    if (is.na(group_name) || group_name == "") {
      return(NULL)
    }
    
    return(list(
      sin_coef = sin_coef,
      cos_coef = cos_coef,
      amplitude = amplitude,
      peak_timing = peak_timing,
      model_type = "CLR_cycl_only",
      group_name = group_name
    ))
  }, error = function(e) {
    cat("  Error extracting parameters from", file_path, ":", e$message, "\n")
    return(NULL)
  })
}

# Function to get core-level observed data
get_core_data <- function(species, min_date, max_date) {
  tryCatch({
    # Load appropriate data based on species type
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

# Function to generate predicted seasonal trends
generate_predicted_trends <- function(sin_coef, cos_coef, months = 1:12) {
  seasonal_effect <- sin_coef * sin(2 * pi * months / 12) + 
                    cos_coef * cos(2 * pi * months / 12)
  return(seasonal_effect)
}

# Function to create phenology plot with observed data and trend lines
create_phenology_plot <- function(species, time_periods, logit_params_2013_2015, logit_params_2015_2018 = NULL, clr_params = NULL) {
  cat("  Creating phenology plot for", species, "\n")
  
  all_data <- list()
  all_predicted <- list()
  
  # Process each time period
  for (period in names(time_periods)) {
    min_date <- time_periods[[period]]$min
    max_date <- time_periods[[period]]$max
    
    # Get core data
    core_data <- get_core_data(species, min_date, max_date)
    
    if (!is.null(core_data) && nrow(core_data) > 0) {
      # Add time period info
      period_data <- core_data %>%
        mutate(
          species = species,
          time_period = period,
          data_type = "Observed Core Data"
        )
      
      all_data[[period]] <- period_data
      
      # Generate predicted trends for this period
      if (period == "2013-2015" && !is.null(logit_params_2013_2015)) {
        # Use 2013-2015 Logit Beta Regression parameters
        months <- 1:12
        predicted_trends <- generate_predicted_trends(logit_params_2013_2015$sin_coef, logit_params_2013_2015$cos_coef, months)
        
        predicted_data <- data.frame(
          month = months,
          month_name = month.abb[months],
          abundance = predicted_trends,
          species = species,
          time_period = period,
          data_type = "Logit Beta Regression (2013-2015)",
          model_approach = "Logit Beta 2013-2015"
        )
        
        all_predicted[[paste0(period, "_logit_2013_2015")]] <- predicted_data
      }
      
      if (period == "2015-2018" && !is.null(logit_params_2015_2018)) {
        # Use 2015-2018 Logit Beta Regression parameters
        months <- 1:12
        predicted_trends <- generate_predicted_trends(logit_params_2015_2018$sin_coef, logit_params_2015_2018$cos_coef, months)
        
        predicted_data <- data.frame(
          month = months,
          month_name = month.abb[months],
          abundance = predicted_trends,
          species = species,
          time_period = period,
          data_type = "Logit Beta Regression (2015-2018)",
          model_approach = "Logit Beta 2015-2018"
        )
        
        all_predicted[[paste0(period, "_logit_2015_2018")]] <- predicted_data
      }
      
      if (period == "2015-2018" && !is.null(clr_params)) {
        # Use CLR parameters for 2015-2018
        months <- 1:12
        predicted_trends <- generate_predicted_trends(clr_params$sin_coef, clr_params$cos_coef, months)
        
        predicted_data <- data.frame(
          month = months,
          month_name = month.abb[months],
          abundance = predicted_trends,
          species = species,
          time_period = period,
          data_type = "CLR Regression",
          model_approach = "CLR"
        )
        
        all_predicted[[paste0(period, "_clr")]] <- predicted_data
      }
    }
  }
  
  if (length(all_data) == 0) {
    cat("    No data found for", species, "\n")
    return(NULL)
  }
  
  # Combine all data
  all_cores <- do.call(rbind, all_data)
  all_predicted_combined <- do.call(rbind, all_predicted)
  
  # Create the phenology plot
  p <- ggplot() +
    # Add seasonal background
    annotate("rect", xmin = 5.5, xmax = 8.5, ymin = -Inf, ymax = Inf, 
             alpha = 0.1, fill = "gold") +
    annotate("rect", xmin = 8.5, xmax = 11.5, ymin = -Inf, ymax = Inf, 
             alpha = 0.1, fill = "orange") +
    annotate("rect", xmin = 11.5, xmax = 2.5, ymin = -Inf, ymax = Inf, 
             alpha = 0.1, fill = "lightblue") +
    annotate("rect", xmin = 2.5, xmax = 5.5, ymin = -Inf, ymax = Inf, 
             alpha = 0.1, fill = "lightgreen") +
    
    # Add seasonal labels
    annotate("text", x = 7, y = max(all_cores$abundance, na.rm = TRUE) * 0.8, 
             label = "Summer", color = "darkgoldenrod", fontface = "bold") +
    annotate("text", x = 10, y = max(all_cores$abundance, na.rm = TRUE) * 0.8, 
             label = "Fall", color = "darkorange", fontface = "bold") +
    annotate("text", x = 1, y = max(all_cores$abundance, na.rm = TRUE) * 0.8, 
             label = "Winter", color = "darkblue", fontface = "bold") +
    annotate("text", x = 4, y = max(all_cores$abundance, na.rm = TRUE) * 0.8, 
             label = "Spring", color = "darkgreen", fontface = "bold") +
    
    # Add predicted trend lines
    {
      if (!is.null(all_predicted_combined) && nrow(all_predicted_combined) > 0) {
        geom_line(data = all_predicted_combined, 
                  aes(x = month, y = abundance, color = model_approach, linetype = time_period), 
                  linewidth = 1.2, alpha = 0.8)
      }
    } +
    
    # Add individual core observations
    geom_point(data = all_cores, 
               aes(x = month, y = abundance, color = time_period), 
               size = 1.5, alpha = 0.6, position = position_jitter(width = 0.2)) +
    
    # Enhanced styling
    scale_color_viridis(discrete = TRUE, option = "plasma") +
    scale_x_continuous(breaks = 1:12, 
                      labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
    
    # Labels and theme
    labs(title = paste("Seasonal Phenology:", species),
         subtitle = "Observed Core Data + Model Predictions (2013-2015 vs 2015-2018)",
         x = "Month", 
         y = "Abundance",
         color = "Data Source",
         linetype = "Time Period") +
    
    theme_minimal() +
    theme(
      text = element_text(size = 12),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "darkgray"),
      legend.position = "bottom",
      legend.box = "horizontal",
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# Function to analyze seasonality consistency
analyze_seasonality_consistency <- function(all_parameters) {
  cat("📊 Analyzing seasonality consistency...\n")
  
  # Convert peak_month to season categories
  consistency_analysis <- all_parameters %>%
    mutate(
      season = case_when(
        peak_timing >= 12 | peak_timing < 3 ~ "Winter",
        peak_timing >= 3 & peak_timing < 6 ~ "Spring", 
        peak_timing >= 6 & peak_timing < 9 ~ "Summer",
        peak_timing >= 9 & peak_timing < 12 ~ "Fall",
        TRUE ~ "Unknown"
      ),
      season_order = case_when(
        season == "Winter" ~ 1,
        season == "Spring" ~ 2,
        season == "Summer" ~ 3,
        season == "Fall" ~ 4,
        TRUE ~ 5
      )
    ) %>%
    group_by(group_name) %>%
    summarise(
      n_models = n(),
      n_approaches = n_distinct(model_type),
      season_consistency = n_distinct(season),
      month_consistency = n_distinct(round(peak_timing)),
      consistency_score = case_when(
        season_consistency == 1 ~ "High",
        season_consistency == 2 ~ "Medium", 
        TRUE ~ "Low"
      ),
      .groups = 'drop'
    ) %>%
    arrange(desc(n_models), desc(consistency_score))
  
  return(consistency_analysis)
}

# Function to create phenology visualization
create_phenology_visualization <- function(all_parameters, phenology_plots) {
  cat("🎨 Creating phenology visualization...\n")
  
  # Create peak timing vs amplitude plot
  p1 <- ggplot(all_parameters, aes(x = peak_timing, y = amplitude, color = model_type, shape = model_type)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_viridis(discrete = TRUE, option = "plasma") +
    scale_shape_manual(values = c(16, 17, 18, 19)) +
    labs(title = "Peak Timing vs Amplitude by Model Type",
         x = "Peak Month", y = "Seasonal Amplitude",
         color = "Model Type", shape = "Model Type") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Create season distribution plot
  season_data <- all_parameters %>%
    mutate(
      season = case_when(
        peak_timing >= 12 | peak_timing < 3 ~ "Winter",
        peak_timing >= 3 & peak_timing < 6 ~ "Spring", 
        peak_timing >= 6 & peak_timing < 9 ~ "Summer",
        peak_timing >= 9 & peak_timing < 12 ~ "Fall",
        TRUE ~ "Unknown"
      )
    ) %>%
    filter(season != "Unknown") %>%
    group_by(season, model_type) %>%
    summarise(count = n(), .groups = 'drop')
  
  p2 <- ggplot(season_data, aes(x = season, y = count, fill = model_type)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    scale_fill_viridis(discrete = TRUE, option = "plasma") +
    labs(title = "Seasonal Peak Distribution by Model Type",
         x = "Season", y = "Number of Groups",
         fill = "Model Type") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Create consistency by group type plot
  consistency_data <- all_parameters %>%
    group_by(group_name) %>%
    summarise(
      n_models = n(),
      n_approaches = n_distinct(model_type),
      season_consistency = n_distinct(case_when(
        peak_timing >= 12 | peak_timing < 3 ~ "Winter",
        peak_timing >= 3 & peak_timing < 6 ~ "Spring", 
        peak_timing >= 6 & peak_timing < 9 ~ "Summer",
        peak_timing >= 9 & peak_timing < 12 ~ "Fall",
        TRUE ~ "Unknown"
      )),
      .groups = 'drop'
    ) %>%
    mutate(
      consistency_score = case_when(
        season_consistency == 1 ~ "High",
        season_consistency == 2 ~ "Medium", 
        TRUE ~ "Low"
      )
    )
  
  p3 <- ggplot(consistency_data, aes(x = consistency_score, fill = consistency_score)) +
    geom_bar(alpha = 0.8) +
    scale_fill_manual(values = c("High" = "darkgreen", "Medium" = "orange", "Low" = "red")) +
    labs(title = "Consistency Scores Across Groups",
         x = "Consistency Score", y = "Number of Groups",
         fill = "Score") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # Save individual plots
  ggsave("../../figures/seasonality_peak_timing_amplitude.png", p1, width = 10, height = 8, dpi = 300)
  cat("✅ Saved peak timing plot to figures/seasonality_peak_timing_amplitude.png\n")
  
  ggsave("../../figures/seasonality_distribution_by_season.png", p2, width = 10, height = 8, dpi = 300)
  cat("✅ Saved season distribution plot to figures/seasonality_distribution_by_season.png\n")
  
  ggsave("../../figures/seasonality_consistency_by_group_type.png", p3, width = 10, height = 8, dpi = 300)
  cat("✅ Saved consistency plot to figures/seasonality_consistency_by_group_type.png\n")
  
  # Save phenology plots
  if (length(phenology_plots) > 0) {
    # Save first few phenology plots as examples
    for (i in 1:min(4, length(phenology_plots))) {
      species_name <- names(phenology_plots)[i]
      plot_obj <- phenology_plots[[i]]
      
      if (!is.null(plot_obj)) {
        filename <- paste0("../../figures/phenology_", gsub("[^a-zA-Z0-9]", "_", species_name), ".png")
        ggsave(filename, plot_obj, width = 12, height = 8, dpi = 300)
        cat("✅ Saved phenology plot for", species_name, "to", filename, "\n")
      }
    }
  }
  
  cat("✅ All phenology visualizations saved to figures/ directory\n")
}

# Function to create consistency report
create_consistency_report <- function(all_parameters, consistency_analysis) {
  cat("📊 Creating detailed consistency report...\n")
  
  # Create comprehensive report
  report <- list(
    all_parameters = all_parameters,
    consistency_analysis = consistency_analysis,
    summary_stats = list(
      total_groups = n_distinct(all_parameters$group_name),
      total_models = nrow(all_parameters),
      groups_with_multiple_models = sum(consistency_analysis$n_models > 1),
      high_consistency_groups = sum(consistency_analysis$consistency_score == "High"),
      medium_consistency_groups = sum(consistency_analysis$consistency_score == "Medium"),
      low_consistency_groups = sum(consistency_analysis$consistency_score == "Low")
    ),
    model_coverage = all_parameters %>%
      group_by(model_type) %>%
      summarise(count = n(), .groups = 'drop'),
    seasonal_distribution = all_parameters %>%
      mutate(
        season = case_when(
          peak_timing >= 12 | peak_timing < 3 ~ "Winter",
          peak_timing >= 3 & peak_timing < 6 ~ "Spring", 
          peak_timing >= 6 & peak_timing < 9 ~ "Summer",
          peak_timing >= 9 & peak_timing < 12 ~ "Fall",
          TRUE ~ "Unknown"
        )
      ) %>%
      group_by(season, model_type) %>%
      summarise(count = n(), .groups = 'drop')
  )
  
  # Save report
  saveRDS(report, "../../data/summary/seasonality_consistency_report.rds")
  cat("✅ Saved detailed report to data/summary/seasonality_consistency_report.rds\n")
  
  return(report)
}

# Main function
main <- function() {
  cat("🚀 Starting Seasonality Consistency Analysis\n\n")
  
  # Define time periods
  time_periods <- list(
    "2013-2015" = list(min = "20130601", max = "20151101"),
    "2015-2018" = list(min = "20150101", max = "20180101")
  )
  
  # Extract parameters from Logit Beta Regression models
  cat("📊 Extracting Logit Beta Regression parameters...\n")
  
  # 2013-2015 models (improved phenology models)
  logit_dirs_2013_2015 <- c(
    "../../data/model_outputs/logit_phenology_improved_2013_2015"
  )
  
  # 2015-2018 models (standard models)
  logit_dirs_2015_2018 <- c(
    "../../data/model_outputs/logit_beta_regression/cycl_only",
    "../../data/model_outputs/logit_beta_regression/env_cov",
    "../../data/model_outputs/logit_beta_regression/env_cycl"
  )
  
  logit_params <- list()
  
  # Process 2013-2015 models
  for (logit_dir in logit_dirs_2013_2015) {
    if (dir.exists(logit_dir)) {
      # Look for combined samples files (without chain numbers)
      logit_files <- list.files(logit_dir, pattern = "samples_.*_2013.*\\.rds", full.names = TRUE)
      # Filter out individual chain files - look for files that don't end with _chainN.rds
      logit_files <- logit_files[!grepl("_chain[0-9]+\\.rds$", logit_files)]
      
      cat("    Found", length(logit_files), "combined files in", logit_dir, "\n")
      
      params_2013_2015_count <- 0
      for (file in logit_files) {
        params <- extract_logit_seasonal_params(file)
        if (!is.null(params) && !is.na(params$group_name) && params$group_name != "unknown") {
          logit_params[[paste0(params$group_name, "_2013_2015")]] <- params
          params_2013_2015_count <- params_2013_2015_count + 1
          cat("      ✅ Extracted params from", basename(file), "->", params$group_name, "\n")
        }
      }
      cat("    Total 2013-2015 params extracted:", params_2013_2015_count, "\n")
    }
  }
  
  # Process 2015-2018 models
  params_2015_2018_count <- 0
  for (logit_dir in logit_dirs_2015_2018) {
    if (dir.exists(logit_dir)) {
      logit_files <- list.files(logit_dir, pattern = "summary_.*_2015.*\\.rds", full.names = TRUE)
      
      for (file in logit_files) {
        params <- extract_logit_seasonal_params(file)
        if (!is.null(params) && !is.na(params$group_name) && params$group_name != "unknown") {
          logit_params[[paste0(params$group_name, "_2015_2018")]] <- params
          params_2015_2018_count <- params_2015_2018_count + 1
        }
      }
    }
  }
  
  cat("    Total 2015-2018 params extracted:", params_2015_2018_count, "\n")
  cat("  Extracted parameters from", length(logit_params), "Logit Beta Regression models total\n")
  
  # Extract parameters from CLR models
  cat("📊 Extracting CLR parameters...\n")
  clr_dirs <- c("../../data/model_outputs/CLR_regression_final/cycl_only")
  
  clr_params <- list()
  for (clr_dir in clr_dirs) {
    if (dir.exists(clr_dir)) {
      clr_files <- list.files(clr_dir, pattern = "samples_.*_2015.*\\.rds", full.names = TRUE)
      
      for (file in clr_files) {
        params <- extract_clr_seasonal_params(file)
        if (!is.null(params) && !is.na(params$group_name) && params$group_name != "") {
          clr_params[[params$group_name]] <- params
        }
      }
    }
  }
  
  cat("  Extracted parameters from", length(clr_params), "CLR models\n")
  
  # Combine all parameters
  all_parameters <- rbind(
    do.call(rbind, lapply(logit_params, function(x) {
      data.frame(
        group_name = x$group_name,
        model_type = x$model_type,
        sin_coef = x$sin_coef,
        cos_coef = x$cos_coef,
        amplitude = x$amplitude,
        peak_timing = x$peak_timing,
        approach = "Logit Beta Regression"
      )
    })),
    do.call(rbind, lapply(clr_params, function(x) {
      data.frame(
        group_name = x$group_name,
        model_type = x$model_type,
        sin_coef = x$sin_coef,
        cos_coef = x$cos_coef,
        amplitude = x$amplitude,
        peak_timing = x$peak_timing,
        approach = "CLR Regression"
      )
    }))
  )
  
  # Create phenology plots showing both Logit Beta and CLR approaches
  cat("🎨 Creating phenology plots with both approaches...\n")
  phenology_plots <- list()
  
  # Get groups that have 2013-2015 Logit Beta models
  logit_groups_2013_2015 <- gsub("_2013_2015$", "", names(logit_params)[grepl("_2013_2015$", names(logit_params))])
  
  # Get CLR fungal phyla for comparison
  clr_fungal_groups <- c("ascomycota", "basidiomycota", "chytridiomycota", "mortierellomycota")
  available_clr_fungal <- intersect(clr_fungal_groups, names(clr_params))
  
  cat("  Logit Beta 2013-2015 groups:", paste(logit_groups_2013_2015, collapse = ", "), "\n")
  cat("  Available CLR fungal phyla:", paste(available_clr_fungal, collapse = ", "), "\n")
  
  # Check which groups have abundance data
  fungi <- readRDS("../../data/clean/groupAbundances_ITS_2023.rds")
  available_fungal_groups <- names(fungi)
  cat("  Available fungal groups in data:", paste(available_fungal_groups, collapse = ", "), "\n")
  
  # Find groups that have both models AND data
  groups_with_data <- intersect(logit_groups_2013_2015, available_fungal_groups)
  cat("  Groups with both models and data:", paste(groups_with_data, collapse = ", "), "\n")
  
  # Create phenology plots for each group that has both models and data
  for (group in groups_with_data[1:min(4, length(groups_with_data))]) {  # Limit to 4 plots
    # Get parameters for 2013-2015 Logit Beta
    params_2013_2015 <- logit_params[[paste0(group, "_2013_2015")]]
    
    # Use ascomycota as representative CLR fungal phylum for comparison
    clr_params_for_comparison <- NULL
    if ("ascomycota" %in% names(clr_params)) {
      clr_params_for_comparison <- clr_params[["ascomycota"]]
    }
    
    plot_obj <- create_phenology_plot(
      species = group,  # Use the actual group name from data
      time_periods = time_periods,
      logit_params_2013_2015 = params_2013_2015,  # Use 2013-2015 for trend lines
      logit_params_2015_2018 = NULL,  # No 2015-2018 Logit Beta for these groups
      clr_params = clr_params_for_comparison  # Include CLR for comparison
    )
    
    if (!is.null(plot_obj)) {
      phenology_plots[[group]] <- plot_obj
    }
  }
  
  # Also create a dedicated comparison plot
  if ("ascomycota" %in% names(clr_params) && "saprotroph_2013_2015" %in% names(logit_params)) {
    cat("  Creating comprehensive comparison plot...\n")
    
    comparison_plot <- create_phenology_plot(
      species = "saprotroph",  # Use actual group name
      time_periods = time_periods,
      logit_params_2013_2015 = logit_params[["saprotroph_2013_2015"]],
      logit_params_2015_2018 = NULL,
      clr_params = clr_params[["ascomycota"]]
    )
    
    if (!is.null(comparison_plot)) {
      phenology_plots[["Comprehensive_Comparison"]] <- comparison_plot
    }
  }
  
  cat("  Created", length(phenology_plots), "phenology plots\n")
  
  # Analyze consistency
  consistency_analysis <- analyze_seasonality_consistency(all_parameters)
  
  # Create visualizations
  create_phenology_visualization(all_parameters, phenology_plots)
  
  # Create report
  report <- create_consistency_report(all_parameters, consistency_analysis)
  
  cat("🎉 Seasonality Consistency Analysis Complete!\n")
  cat("📊 Summary:\n")
  cat("  Total groups analyzed:", report$summary_stats$total_groups, "\n")
  cat("  Groups with multiple models:", report$summary_stats$groups_with_multiple_models, "\n")
  cat("  High consistency groups:", report$summary_stats$high_consistency_groups, "\n")
  cat("  Medium consistency groups:", report$summary_stats$medium_consistency_groups, "\n")
  cat("  Low consistency groups:", report$summary_stats$low_consistency_groups, "\n")
}

# Run the analysis
main()
