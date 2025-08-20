#!/usr/bin/env Rscript

# Evaluate MCMC Convergence Status
# This script examines all MCMC outputs to assess convergence quality

# Load required packages
library(here)
library(coda)
library(dplyr)
library(purrr)

# Source environment
source("source.R")

# Function to evaluate convergence for a single model
evaluate_convergence <- function(file_path) {
  tryCatch({
    # Load the MCMC samples
    samples <- readRDS(file_path)
    
    # Extract model info from filename
    filename <- basename(file_path)
    model_info <- strsplit(filename, "_")[[1]]
    
    # Determine if this is a combined file or individual chain
    is_combined <- !grepl("chain", filename)
    
    if (is_combined) {
      # Combined file - check if it's a list or matrix
      if (is.list(samples)) {
        # Convert list to mcmc.list
        mcmc_list <- mcmc.list(lapply(samples, mcmc))
      } else {
        # Single matrix - convert to mcmc.list
        mcmc_list <- mcmc.list(mcmc(samples))
      }
    } else {
      # Individual chain - convert to mcmc
      mcmc_list <- mcmc.list(mcmc(samples))
    }
    
    # Calculate convergence diagnostics
    gelman_rubin <- gelman.diag(mcmc_list, multivariate = FALSE)
    effective_size <- effectiveSize(mcmc_list)
    
    # Extract key statistics
    max_psrf <- max(gelman_rubin$psrf[, "Point est."], na.rm = TRUE)
    min_eff_size <- min(effective_size, na.rm = TRUE)
    n_chains <- length(mcmc_list)
    n_iter <- nrow(mcmc_list[[1]])
    
    # Determine convergence status
    if (max_psrf <= 1.1 && min_eff_size >= 100) {
      status <- "CONVERGED"
    } else if (max_psrf <= 1.2 && min_eff_size >= 50) {
      status <- "PARTIALLY_CONVERGED"
    } else {
      status <- "NOT_CONVERGED"
    }
    
    # Return results
    list(
      filename = filename,
      model_type = paste(model_info[1:2], collapse = "_"),
      species = model_info[3],
      time_period = paste(model_info[4:5], collapse = "_"),
      n_chains = n_chains,
      n_iterations = n_iter,
      max_psrf = max_psrf,
      min_effective_size = min_eff_size,
      convergence_status = status,
      file_size_mb = file.size(file_path) / 1024^2
    )
    
  }, error = function(e) {
    # Return error info
    list(
      filename = basename(file_path),
      model_type = "ERROR",
      species = "ERROR",
      time_period = "ERROR",
      n_chains = NA,
      n_iterations = NA,
      max_psrf = NA,
      min_effective_size = NA,
      convergence_status = "ERROR",
      file_size_mb = file.size(file_path) / 1024^2,
      error_message = e$message
    )
  })
}

# Function to find all MCMC output files
find_mcmc_files <- function() {
  base_dir <- here("data", "model_outputs")
  
  # Find all .rds files in model output directories
  files <- list.files(
    path = base_dir,
    pattern = "\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  # Filter for MCMC output files (exclude summaries)
  mcmc_files <- files[!grepl("summary_", basename(files))]
  
  return(mcmc_files)
}

# Main execution
cat("Evaluating MCMC convergence status...\n")

# Find all MCMC files
mcmc_files <- find_mcmc_files()
cat("Found", length(mcmc_files), "MCMC output files\n")

# Evaluate convergence for each file
cat("Evaluating convergence...\n")
results <- map_dfr(mcmc_files, evaluate_convergence, .id = "file_id")

# Save results
output_file <- here("data", "summary", "mcmc_convergence_status.rds")
saveRDS(results, output_file)

# Display summary
cat("\n=== MCMC CONVERGENCE SUMMARY ===\n")
cat("Total files evaluated:", nrow(results), "\n")

# Convergence status summary
status_summary <- results %>%
  filter(!grepl("ERROR", convergence_status)) %>%
  group_by(convergence_status) %>%
  summarise(
    count = n(),
    avg_iterations = mean(n_iterations, na.rm = TRUE),
    avg_chains = mean(n_chains, na.rm = TRUE),
    avg_psrf = mean(max_psrf, na.rm = TRUE),
    avg_eff_size = mean(min_effective_size, na.rm = TRUE)
  )

print(status_summary)

# Model type summary
model_summary <- results %>%
  filter(!grepl("ERROR", convergence_status)) %>%
  group_by(model_type) %>%
  summarise(
    count = n(),
    converged = sum(convergence_status == "CONVERGED"),
    partially_converged = sum(convergence_status == "PARTIALLY_CONVERGED"),
    not_converged = sum(convergence_status == "NOT_CONVERGED")
  )

cat("\n=== MODEL TYPE SUMMARY ===\n")
print(model_summary)

# Detailed results for non-converged models
cat("\n=== NON-CONVERGED MODELS ===\n")
not_converged <- results %>%
  filter(convergence_status == "NOT_CONVERGED") %>%
  select(filename, model_type, species, max_psrf, min_effective_size, n_iterations)

if (nrow(not_converged) > 0) {
  print(not_converged)
} else {
  cat("All models are converged or partially converged!\n")
}

# Save detailed results as CSV for easy viewing
csv_file <- here("data", "summary", "mcmc_convergence_status.csv")
write.csv(results, csv_file, row.names = FALSE)

cat("\nResults saved to:\n")
cat("RDS:", output_file, "\n")
cat("CSV:", csv_file, "\n")

cat("\n=== RECOMMENDATIONS ===\n")

# Identify models that need more iterations
needs_more_iterations <- results %>%
  filter(!grepl("ERROR", convergence_status)) %>%
  filter(n_iterations < 1000 | min_effective_size < 100)

if (nrow(needs_more_iterations) > 0) {
  cat("Models needing more iterations:\n")
  print(needs_more_iterations %>% 
          select(filename, n_iterations, min_effective_size, convergence_status))
} else {
  cat("All models have sufficient iterations!\n")
}

cat("\nEvaluation complete!\n")
