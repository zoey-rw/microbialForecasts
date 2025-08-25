#!/usr/bin/env Rscript

# Diagnose convergence issues in cycl_only models
# Examine traceplots and parameter behavior

library(coda)
library(ggplot2)
library(gridExtra)

# Function to load and analyze MCMC output
analyze_mcmc_output <- function(file_path, model_name) {
  cat("Analyzing:", file_path, "\n")
  
  # Load output
  output <- readRDS(file_path)
  samples <- output$samples
  
  cat("Sample dimensions:", dim(samples), "\n")
  cat("Parameters monitored:", length(colnames(samples)), "\n")
  
  # Key parameters to examine
  key_params <- c("beta[1]", "beta[2]", "intercept", "rho", "sigma")
  available_params <- key_params[key_params %in% colnames(samples)]
  
  if (length(available_params) == 0) {
    cat("No key parameters found in this output\n")
    return(NULL)
  }
  
  # Parameter summaries
  cat("\n=== Parameter Summaries ===\n")
  print(summary(samples[, available_params]))
  
  # Effective sample sizes
  cat("\n=== Effective Sample Sizes ===\n")
  ess <- effectiveSize(samples[, available_params])
  print(ess)
  
  # Convergence diagnostics
  cat("\n=== Convergence Diagnostics ===\n")
  cat("Parameters with ESS < 100:", names(ess[ess < 100]), "\n")
  cat("Parameters with ESS < 50:", names(ess[ess < 50]), "\n")
  cat("Parameters with ESS < 10:", names(ess[ess < 10]), "\n")
  
  # Create traceplots for key parameters
  traceplots <- list()
  for (param in available_params) {
    if (param %in% colnames(samples)) {
      traceplots[[param]] <- ggplot(data.frame(
        iteration = 1:nrow(samples),
        value = samples[, param]
      ), aes(x = iteration, y = value)) +
        geom_line(alpha = 0.7) +
        labs(title = paste("Traceplot:", param),
             x = "Iteration", y = "Value") +
        theme_minimal()
    }
  }
  
  # Create density plots
  density_plots <- list()
  for (param in available_params) {
    if (param %in% colnames(samples)) {
      density_plots[[param]] <- ggplot(data.frame(
        value = samples[, param]
      ), aes(x = value)) +
        geom_density(fill = "steelblue", alpha = 0.7) +
        labs(title = paste("Density:", param),
             x = "Value", y = "Density") +
        theme_minimal()
    }
  }
  
  # Save plots
  output_dir <- "figures/convergence_diagnostics"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Combine traceplots
  if (length(traceplots) > 0) {
    trace_combined <- do.call(grid.arrange, c(traceplots, ncol = 2))
    ggsave(file.path(output_dir, paste0(model_name, "_traceplots.png")), 
           trace_combined, width = 12, height = 8)
  }
  
  # Combine density plots
  if (length(density_plots) > 0) {
    density_combined <- do.call(grid.arrange, c(density_plots, ncol = 2))
    ggsave(file.path(output_dir, paste0(model_name, "_densities.png")), 
           density_combined, width = 12, height = 8)
  }
  
  return(list(
    samples = samples,
    ess = ess,
    traceplots = traceplots,
    density_plots = density_plots
  ))
}

# Analyze recent cycl_only outputs
cat("=== Diagnosing cycl_only Model Convergence ===\n")

# Find recent output files
output_dir <- "data/model_outputs/old_mcmc/logit_beta_regression/cycl_only"
files <- list.files(output_dir, pattern = "samples_.*chain1.rds", full.names = TRUE)
recent_files <- files[order(file.info(files)$mtime, decreasing = TRUE)]

cat("Found", length(recent_files), "recent output files\n")

# Analyze the most recent file
if (length(recent_files) > 0) {
  most_recent <- recent_files[1]
  model_name <- gsub(".*samples_cycl_only_(.+)_chain1.rds", "\\1", basename(most_recent))
  
  cat("\nAnalyzing most recent file:", basename(most_recent), "\n")
  results <- analyze_mcmc_output(most_recent, model_name)
  
  if (!is.null(results)) {
    cat("\n=== Analysis Complete ===\n")
    cat("Check figures/convergence_diagnostics/ for plots\n")
  }
} else {
  cat("No recent output files found\n")
}

