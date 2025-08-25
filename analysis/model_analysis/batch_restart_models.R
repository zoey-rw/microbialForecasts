#!/usr/bin/env Rscript

# Batch Restart Script for MCMC Models
# Automatically finds models with existing chains and restarts them with improved initial values
# Designed to "rescue" models that were fit with loose priors

# Load required libraries
library(here)
library(tidyverse)
library(parallel)
library(doParallel)

# Load restart functions
source("model_restart_functions.R")

#' Batch restart models that have existing chain files
#'
#' @param model_input_df Path to model input dataframe CSV
#' @param output_base_dir Base directory for model outputs
#' @param max_models Maximum number of models to process (for testing)
#' @param restart_params List of restart parameters
#' @param parallel_cores Number of cores for parallel processing
#' @return Summary of batch restart results
batch_restart_models <- function(model_input_df = "/Users/zoeywerbin/Documents/microbialForecasts/data/clean/model_input_df.csv",
                                output_base_dir = "/Users/zoeywerbin/Documents/microbialForecasts/data/model_outputs/logit_beta_regression",
                                max_models = NULL,
                                restart_params = list(),
                                parallel_cores = 4) {

  cat("=== BATCH MODEL RESTART ===\n")
  cat("Input file:", model_input_df, "\n")
  cat("Output base:", output_base_dir, "\n")
  cat("Parallel cores:", parallel_cores, "\n")
  if (!is.null(max_models)) {
    cat("Max models to process:", max_models, "\n")
  }

  # Set default restart parameters
  default_params <- list(
    min_ess = 50,
    max_rhat = 1.2,
    burnin_proportion = 0.4,
    use_fallback_for_extreme = TRUE,
    fallback_strategy = "conservative",
    niter = 2000,
    nburnin = 500,
    thin = 5,
    nchains = 2
  )

  restart_params <- modifyList(default_params, restart_params)

  # Load model input data
  if (!file.exists(model_input_df)) {
    stop("Model input file not found: ", model_input_df)
  }

  cat("Loading model input data...\n")
  params_in <- read.csv(model_input_df,
                       colClasses = c(rep("character", 4),
                                    rep("logical", 2),
                                    rep("character", 4)))

  cat("Found", nrow(params_in), "models in input data\n")

  # Find models that have existing chain files
  cat("Scanning for existing chain files...\n")

  models_with_chains <- list()
  total_chain_files <- 0

  for (i in 1:nrow(params_in)) {
    if (!is.null(max_models) && max_models > 0 && length(models_with_chains) >= max_models) {
      break
    }

    row <- params_in[i, ]

    # Create model ID for file naming
    legacy_indicator <- ifelse(grepl("Legacy with covariate", row$scenario),
                             "with_legacy_covariate",
                             "without_legacy_covariate")
    model_id <- paste(row$model_name, row$species, row$min.date, row$max.date, legacy_indicator, sep = "_")

    # Look for existing chain files
    output_dir <- file.path(output_base_dir, row$model_name)
    pattern <- paste0("samples_", model_id, "_chain[0-9]+\\.rds$")
    chain_files <- list.files(output_dir, pattern = pattern, full.names = TRUE)

    if (length(chain_files) > 0) {
      models_with_chains[[length(models_with_chains) + 1]] <- list(
        model_info = row,
        model_id = model_id,
        chain_files = chain_files,
        n_chains = length(chain_files)
      )
      total_chain_files <- total_chain_files + length(chain_files)

      cat("  ✓ Model", i, ":", row$species, "(", row$model_name, ") -", length(chain_files), "chains\n")
    }
  }

  cat("\nFound", length(models_with_chains), "models with existing chains\n")
  cat("Total chain files:", total_chain_files, "\n")

  if (length(models_with_chains) == 0) {
    cat("❌ No models with existing chains found. Nothing to restart.\n")
    return(NULL)
  }

  # Create restart tasks
  restart_tasks <- list()
  task_counter <- 1

  for (model_idx in seq_along(models_with_chains)) {
    model_data <- models_with_chains[[model_idx]]

    # Create a task for each chain
    for (chain_no in 1:model_data$n_chains) {
      restart_tasks[[task_counter]] <- list(
        task_id = task_counter,
        model_idx = model_idx,
        chain_no = chain_no,
        model_info = model_data$model_info,
        model_id = model_data$model_id,
        chain_files = model_data$chain_files,
        restart_params = restart_params
      )
      task_counter <- task_counter + 1
    }
  }

  cat("Created", length(restart_tasks), "restart tasks\n")

  # For now, use sequential processing to test the core functionality
  # TODO: Fix parallel processing in future version
  cat("Using sequential processing for testing (parallel processing needs refinement)...\n")

  # Define processing function inside the main function scope
  process_restart_task <- function(task_idx) {
    # Access the variables from the parent environment
    task <- restart_tasks[[task_idx]]
    model_info <- task$model_info

    cat("Worker", task_idx, ": Processing model", task$model_idx, "chain", task$chain_no, "\n")

    tryCatch({
      # Extract final values from all chains for this model
      cat("  Extracting final values from", length(task$chain_files), "chain files...\n")
      extraction_result <- extract_final_values_from_chains(
        task$chain_files,
        min_ess = task$restart_params$min_ess,
        max_rhat = task$restart_params$max_rhat,
        burnin_proportion = task$restart_params$burnin_proportion
      )

      # Create restart initial values
      restart_inits <- create_restart_inits(
        extraction_result,
        use_fallback_for_extreme = task$restart_params$use_fallback_for_extreme,
        fallback_strategy = task$restart_params$fallback_strategy
      )

      # Return successful result
      return(list(
        task_id = task_idx,
        model_idx = task$model_idx,
        chain_no = task$chain_no,
        status = "SUCCESS",
        model_info = model_info,
        extraction_result = extraction_result,
        restart_inits = restart_inits,
        extreme_values_detected = sum(restart_inits$extreme_flags),
        fallback_values_used = sum(restart_inits$fallback_used)
      ))

    }, error = function(e) {
      cat("Worker", task_idx, ": Error processing model", task$model_idx, "chain", task$chain_no, "\n")
      cat("  Error:", e$message, "\n")

      return(list(
        task_id = task_idx,
        model_idx = task$model_idx,
        chain_no = task$chain_no,
        status = "ERROR",
        model_info = model_info,
        error = e$message
      ))
    })
  }

  # Run restart tasks sequentially
  cat("Starting sequential restart processing...\n")
  start_time <- Sys.time()

  restart_results <- list()
  for (task_idx in 1:length(restart_tasks)) {
    cat("Processing task", task_idx, "of", length(restart_tasks), "...\n")
    restart_results[[task_idx]] <- process_restart_task(task_idx)
  }

  end_time <- Sys.time()
  processing_time <- difftime(end_time, start_time, units = "mins")

  # Summarize results
  successful_restarts <- 0
  failed_restarts <- 0
  total_extreme_values <- 0
  total_fallback_values <- 0

  cat("\n=== RESTART PROCESSING COMPLETE ===\n")
  cat("Processing time:", round(processing_time, 2), "minutes\n")
  cat("Tasks processed:", length(restart_results), "\n")

  # Process results by model
  model_results <- list()

  for (result in restart_results) {
    model_idx <- result$model_idx
    model_key <- paste0("model_", model_idx)

    if (is.null(model_results[[model_key]])) {
      model_results[[model_key]] <- list(
        model_info = result$model_info,
        chains = list()
      )
    }

    model_results[[model_key]]$chains[[length(model_results[[model_key]]$chains) + 1]] <- result

    if (result$status == "SUCCESS") {
      successful_restarts <- successful_restarts + 1
      total_extreme_values <- total_extreme_values + result$extreme_values_detected
      total_fallback_values <- total_fallback_values + result$fallback_values_used

      cat("✅ Model", model_idx, "Chain", result$chain_no, ": SUCCESS\n")
      cat("   - Extreme values:", result$extreme_values_detected, "\n")
      cat("   - Fallback values:", result$fallback_values_used, "\n")
    } else {
      failed_restarts <- failed_restarts + 1
      cat("❌ Model", model_idx, "Chain", result$chain_no, ": FAILED -", result$error, "\n")
    }
  }

  # Create summary report
  summary_report <- list(
    total_models = length(models_with_chains),
    total_tasks = length(restart_tasks),
    successful_restarts = successful_restarts,
    failed_restarts = failed_restarts,
    total_extreme_values = total_extreme_values,
    total_fallback_values = total_fallback_values,
    processing_time_minutes = as.numeric(processing_time),
    model_results = model_results,
    restart_params = restart_params
  )

  cat("\n=== SUMMARY REPORT ===\n")
  cat("Models processed:", length(models_with_chains), "\n")
  cat("Restart tasks:", length(restart_tasks), "\n")
  cat("Successful restarts:", successful_restarts, "\n")
  cat("Failed restarts:", failed_restarts, "\n")
  cat("Success rate:", round(successful_restarts / length(restart_tasks) * 100, 1), "%\n")
  cat("Total extreme values detected:", total_extreme_values, "\n")
  cat("Total fallback values used:", total_fallback_values, "\n")
  cat("Processing time:", round(processing_time, 2), "minutes\n")

  return(summary_report)
}

#' Create a report of models that can be restarted
#'
#' @param model_input_df Path to model input dataframe CSV
#' @param output_base_dir Base directory for model outputs
#' @return Dataframe with restart candidates
create_restart_report <- function(model_input_df = "/Users/zoeywerbin/Documents/microbialForecasts/data/clean/model_input_df.csv",
                                 output_base_dir = "/Users/zoeywerbin/Documents/microbialForecasts/data/model_outputs/logit_beta_regression") {

  cat("=== CREATING RESTART REPORT ===\n")

  # Load model input data
  if (!file.exists(model_input_df)) {
    stop("Model input file not found: ", model_input_df)
  }

  params_in <- read.csv(model_input_df,
                       colClasses = c(rep("character", 4),
                                    rep("logical", 2),
                                    rep("character", 4)))

  cat("Found", nrow(params_in), "models in input data\n")

  # Check each model for existing chains
  restart_candidates <- list()

  for (i in 1:nrow(params_in)) {
    row <- params_in[i, ]

    # Create model ID for file naming
    legacy_indicator <- ifelse(grepl("Legacy with covariate", row$scenario),
                             "with_legacy_covariate",
                             "without_legacy_covariate")
    model_id <- paste(row$model_name, row$species, row$min.date, row$max.date, legacy_indicator, sep = "_")

    # Look for existing chain files
    output_dir <- file.path(output_base_dir, row$model_name)
    pattern <- paste0("samples_", model_id, "_chain[0-9]+\\.rds$")
    chain_files <- list.files(output_dir, pattern = pattern, full.names = TRUE)

    if (length(chain_files) > 0) {
      # Try to get basic info from first chain file
      chain_info <- tryCatch({
        chain_data <- readRDS(chain_files[1])
        if (is.list(chain_data) && "metadata" %in% names(chain_data)) {
          metadata <- chain_data$metadata
          niter <- if ("niter" %in% names(metadata)) metadata$niter else NA
          convergence_status <- "unknown"
        } else {
          niter <- NA
          convergence_status <- "unknown"
        }
        list(
          niter = niter,
          convergence_status = convergence_status,
          samples_dim = if (is.list(chain_data) && "samples" %in% names(chain_data)) dim(chain_data$samples) else c(NA, NA)
        )
      }, error = function(e) {
        list(
          niter = NA,
          convergence_status = "error_loading",
          samples_dim = c(NA, NA)
        )
      })

      restart_candidates[[length(restart_candidates) + 1]] <- c(
        row,
        list(
          model_id = model_id,
          n_chains = length(chain_files),
          chain_files = paste(basename(chain_files), collapse = "; "),
          n_iterations = chain_info$niter,
          convergence_status = chain_info$convergence_status,
          samples_dim = paste(chain_info$samples_dim, collapse = " x ")
        )
      )
    }
  }

  # Convert to dataframe
  if (length(restart_candidates) > 0) {
    result_df <- do.call(rbind, lapply(restart_candidates, function(x) {
      as.data.frame(x, stringsAsFactors = FALSE)
    }))

    cat("Found", nrow(result_df), "models that can be restarted\n")

    # Summary by model type
    cat("\nModels by type:\n")
    print(table(result_df$model_name))

    # Summary by rank
    cat("\nModels by rank:\n")
    print(table(result_df$rank.name))

    return(result_df)
  } else {
    cat("❌ No models found that can be restarted\n")
    return(NULL)
  }
}

# Main execution
if (sys.nframe() == 0) {
  cat("Batch Restart Models Script\n")
  cat("Usage: Rscript batch_restart_models.R [max_models] [cores]\n")

  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  max_models <- if (length(args) >= 1) as.numeric(args[1]) else NULL
  parallel_cores <- if (length(args) >= 2) as.numeric(args[2]) else 4

  # Create restart report first
  cat("Creating restart report...\n")
  restart_report <- create_restart_report()

  if (!is.null(restart_report)) {
    # Save report
    report_file <- paste0("restart_candidates_report_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    write.csv(restart_report, report_file, row.names = FALSE)
    cat("Report saved to:", report_file, "\n")

    # Check if we should proceed automatically or ask user
    if (length(args) >= 1) {
      # Non-interactive mode - proceed automatically
      cat("\nFound", nrow(restart_report), "models that can be restarted.\n")
      cat("Proceeding with batch restart in non-interactive mode...\n")
      proceed <- TRUE
    } else {
      # Interactive mode - ask user
      cat("\nFound", nrow(restart_report), "models that can be restarted.\n")
      cat("Proceed with batch restart? (y/N): ")
      response <- tolower(readLines(con = stdin(), n = 1))
      proceed <- (response == "y" || response == "yes")
    }

    if (proceed) {
      cat("Starting batch restart...\n")

      # Run batch restart
      restart_results <- batch_restart_models(
        max_models = max_models,
        parallel_cores = parallel_cores,
        restart_params = list(
          min_ess = 50,
          max_rhat = 1.2,
          fallback_strategy = "conservative",
          niter = 2000,
          nchains = 2
        )
      )

      # Save results
      if (!is.null(restart_results)) {
        results_file <- paste0("batch_restart_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
        saveRDS(restart_results, results_file)
        cat("Results saved to:", results_file, "\n")
      }

    } else {
      cat("Batch restart cancelled by user.\n")
    }

  } else {
    cat("No models available for restart.\n")
  }
}
