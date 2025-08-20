#!/usr/bin/env Rscript

# Script to create traceplots for MCMC chains to assess convergence
# Focus on key parameters: beta, beta_sd, rho, sigma, intercept

library(here)
library(coda)
library(ggplot2)
library(dplyr)
library(tidyr)

# Set up environment
source(here("source_local.R"))

# Function to create traceplots for a specific model
create_traceplots <- function(model_id, output_dir = "figures/traceplots") {
  
  cat("Creating traceplots for:", model_id, "\n")
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Find chain files
  model_type <- strsplit(model_id, "_")[[1]][1:2]
  model_type <- paste(model_type, collapse = "_")
  
  chain_files <- list.files(
    here("data", "model_outputs", "logit_beta_regression", model_type),
    pattern = paste0(model_id, "_chain[1-4]\\.rds"),
    full.names = TRUE
  )
  
  if (length(chain_files) == 0) {
    cat("No chain files found for:", model_id, "\n")
    return(NULL)
  }
  
  cat("Found", length(chain_files), "chain files\n")
  
  # Load chains
  chains <- list()
  for (i in seq_along(chain_files)) {
    chain_data <- readRDS(chain_files[i])
    if (is.list(chain_data) && "samples" %in% names(chain_data)) {
      chains[[i]] <- as.mcmc(chain_data$samples)
    } else {
      chains[[i]] <- as.mcmc(chain_data)
    }
  }
  
  # Combine into mcmc.list
  mcmc_list <- as.mcmc.list(chains)
  
  # Key parameters to plot
  key_params <- c("beta", "beta_sd", "rho", "sigma", "intercept", "sig", "core_sd")
  
  # Create traceplots for each key parameter
  for (param in key_params) {
    # Find all parameters that match this pattern
    param_names <- grep(paste0("^", param, "\\[?[0-9]*\\]?$"), colnames(mcmc_list[[1]]), value = TRUE)
    
    if (length(param_names) == 0) next
    
    for (param_name in param_names) {
      cat("  Plotting:", param_name, "\n")
      
      # Extract data for plotting
      plot_data <- data.frame()
      for (chain_idx in seq_along(mcmc_list)) {
        chain_data <- data.frame(
          iteration = 1:nrow(mcmc_list[[chain_idx]]),
          value = mcmc_list[[chain_idx]][, param_name],
          chain = paste0("Chain ", chain_idx)
        )
        plot_data <- rbind(plot_data, chain_data)
      }
      
      # Create traceplot
      p <- ggplot(plot_data, aes(x = iteration, y = var1, color = chain)) +
        geom_line(alpha = 0.7) +
        #facet_wrap(~chain, ncol = 1, scales = "free_y") +
        labs(
          title = paste("Traceplot:", param_name, "-", model_id),
          x = "Iteration",
          y = "Parameter Value",
          color = "Chain"
        ) +
        theme_minimal() +
        theme(
          legend.position = "none",
          plot.title = element_text(size = 12, face = "bold"),
          strip.text = element_text(size = 10, face = "bold")
        )
      
      # Save plot
      filename <- paste0(gsub("\\[|\\]", "_", param_name), "_", model_id, "_traceplot.png")
      filepath <- file.path(output_dir, filename)
      ggsave(filepath, p, width = 10, height = 8, dpi = 150)
      
      cat("    Saved:", filepath, "\n")
    }
  }
  
  # Create density plots for key parameters
  for (param in key_params) {
    param_names <- grep(paste0("^", param, "\\[?[0-9]*\\]?$"), colnames(mcmc_list[[1]]), value = TRUE)
    
    if (length(param_names) == 0) next
    
    for (param_name in param_names) {
      cat("  Creating density plot for:", param_name, "\n")
      
      # Extract data for plotting
      plot_data <- data.frame()
      for (chain_idx in seq_along(mcmc_list)) {
        chain_data <- data.frame(
          value = mcmc_list[[chain_idx]][, param_name],
          chain = paste0("Chain ", chain_idx)
        )
        plot_data <- rbind(plot_data, chain_data)
      }
      
      # Create density plot
      p <- ggplot(plot_data, aes(x = value, fill = chain)) +
        geom_density(alpha = 0.6) +
        labs(
          title = paste("Density Plot:", param_name, "-", model_id),
          x = "Parameter Value",
          y = "Density",
          fill = "Chain"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 12, face = "bold"),
          legend.position = "bottom"
        )
      
      # Save plot
      filename <- paste0(gsub("\\[|\\]", "_", param_name), "_", model_id, "_density.png")
      filepath <- file.path(output_dir, filename)
      ggsave(filepath, p, width = 8, height = 6, dpi = 150)
      
      cat("    Saved:", filepath, "\n")
    }
  }
  
  # Create summary statistics table
  cat("  Creating summary statistics\n")
  summary_stats <- summary(mcmc_list)
  
  # Save summary to file
  summary_file <- file.path(output_dir, paste0(model_id, "_summary.txt"))
  sink(summary_file)
  print(summary_stats)
  sink()
  
  cat("    Saved summary to:", summary_file, "\n")
  
  return(mcmc_list)
}

# Main execution
if (!interactive()) {
  # Test with the ectomycorrhizal model
  model_id <- "env_cov_ectomycorrhizal_20151101_20180101"
  
  cat("Creating traceplots for:", model_id, "\n")
  mcmc_list <- create_traceplots(model_id)
  
  if (!is.null(mcmc_list)) {
    cat("\nTraceplot creation completed successfully!\n")
    cat("Check the 'figures/traceplots' directory for output files.\n")
  } else {
    cat("\nFailed to create traceplots.\n")
  }
}
