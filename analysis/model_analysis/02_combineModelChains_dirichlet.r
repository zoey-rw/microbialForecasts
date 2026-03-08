# Combine chains from each Dirichlet model, and create basic summary stats
source("source.R")
source(here::here("analysis/model_analysis/dirichlet_covariance/dirichlet_helper_functions.r"))

# Function to combine Dirichlet chains and create summary structure
combine_dirichlet_chains <- function(chain_paths, model_id) {
  message("Combining ", length(chain_paths), " chains for ", model_id)
  
  # All chain files contain both samples and samples2 in a list
  param_chains <- chain_paths

  # Read parameter chains and plot chains from the unified files
  param_samples <- list()
  plot_samples <- list()
  metadata_list <- list()

  for (i in 1:length(param_chains)) {
    chain_data <- readRDS(param_chains[[i]])
    message("Parameter chain ", i, " dimensions: ", nrow(chain_data$samples), " x ", ncol(chain_data$samples))

    if (is.list(chain_data)) {
      param_samples[[i]] <- chain_data$samples

      # FIX: Extract samples2 natively from the chain list
      if ("samples2" %in% names(chain_data)) {
        plot_samples[[i]] <- chain_data$samples2
      }
      if ("metadata" %in% names(chain_data)) {
        metadata_list[[i]] <- chain_data$metadata
      }
    } else {
      # Fallback
      param_samples[[i]] <- chain_data
      metadata_list[[i]] <- list(model_id = model_id)
    }
  }

  # Combine parameter chains
  library(coda)
  param_mcmc_chains <- lapply(param_samples, mcmc)
  combined_param_chains <- mcmc.list(param_mcmc_chains)

  # Combine plot chains if they exist
  if (length(plot_samples) > 0 && !is.null(plot_samples[[1]])) {
    plot_mcmc_chains <- lapply(plot_samples, mcmc)
    combined_plot_chains <- mcmc.list(plot_mcmc_chains)
  } else {
    combined_plot_chains <- NULL
  }
  
  # Calculate parameter summaries (use fast.summary.mcmc to tolerate NA/NaN in chains)
  param_sum <- fast.summary.mcmc(combined_param_chains)
  param_summary <- list(param_sum$statistics, param_sum$quantiles)

  # Calculate plot summaries if available
  plot_summary <- list()
  if (!is.null(combined_plot_chains)) {
    plot_sum <- fast.summary.mcmc(combined_plot_chains)
    plot_summary[[1]] <- plot_sum$statistics
    plot_summary[[2]] <- plot_sum$quantiles
  } else {
    # Create minimal plot summaries for Dirichlet models
    plot_summary[[1]] <- data.frame()  # Placeholder for means
    plot_summary[[2]] <- data.frame()  # Placeholder for quantiles
  }
  
  # Extract real metadata from the first chain
  metadata <- metadata_list[[1]]
  if (is.null(metadata)) {
    metadata <- list(model_id = model_id)
  }
  
  # Ensure metadata has required fields
  if (is.null(metadata$niteration)) {
    metadata$niteration <- nrow(param_samples[[1]])
  }
  
  # Create the combined output structure (samples2 needed for post-combine plot_summary/gelman)
  out <- list(
    samples = combined_param_chains,
    samples2 = combined_plot_chains,
    param_summary = param_summary,
    plot_summary = plot_summary,
    metadata = metadata
  )

  return(out)
}

# For Dirichlet models: focus on driver uncertainty outputs
target_dir <- here("data/model_outputs/dirichlet_driver_uncertainty/")

cat("Processing Dirichlet regression models in:", target_dir, "\n")

# Dynamically find all sample chain files (exclude checkpoints and progress files)
file.list <- list.files(path = target_dir,
                        pattern = "_chain[0-9]+\\.rds$",
                        recursive = TRUE,
                        full.names = TRUE)
file.list <- file.list[grepl("^samples_", basename(file.list))]

# Extract unique model IDs directly from the file names
model_id_list <- unique(gsub("_chain[0-9]+.rds$", "", gsub("^samples_", "", basename(file.list))))

cat("Total models to process:", length(model_id_list), "\n\n")
cat("Model IDs:", paste(model_id_list, collapse = ", "), "\n")

# Use sequential processing for debugging
cat("Starting sequential processing with", length(model_id_list), "models...\n")

for(model_id in model_id_list) {
  # Do we want to keep all the chain files separately? Deleting them will save space
  delete_samples_files = F

  # Subset to chain files for this model_id
  chain_paths <- file.list[grepl(model_id, file.list, fixed = T)]

  message("Searching model outputs for ", model_id)
  message("Total files in file.list: ", length(file.list))
  message("Found ", length(chain_paths), " chain files for ", model_id)
  if (length(chain_paths) > 0) {
    message("Chain files: ", paste(basename(chain_paths), collapse = ", "))
  }
  
  if (length(chain_paths) == 0) {
    message("Skipping ", model_id, "; no chains found")
    next
  }
  if (length(chain_paths) == 1) {
    message("Skipping ", model_id, "; only one chain found (need multiple chains to combine)")
    next
  }
  savepath <- gsub("_chain[1234567]", "", chain_paths[[1]])
  
  # Don't run loop if the samples file is already newer than the chain files
  if (samples_exists(chain_paths[[1]])) {
    message("Summary samples file already exists")
  } else {
    # Calculate summary on each output subset, using custom summary function
    message("Combining chains for ", model_id, " with ", length(chain_paths), " chains...")
    out <- tryCatch({
      combine_dirichlet_chains(chain_paths, model_id)
    }, error = function(e) {
      message("Error combining chains for ", model_id, ": ", e$message)
      return(NULL)
    })
    
    if (is.null(out)) {
      message("Failed to combine chains for ", model_id, ", skipping...")
      next
    }
    
    # Remove chains if one has very different values than the others
    chains <- out$samples
    if (length(chains) == 0) {
      message("No chains found in output for ", model_id, ", skipping...")
      next
    }
    
    # Check if intercept parameter exists (NIMBLE uses bracketed names like intercept[1])
    if (!any(grepl("^intercept", colnames(chains[[1]])))) {
      message("No intercept parameter found in chains for ", model_id, ", skipping...")
      next
    }
    
    # Ensure chains have reasonable sample sizes before outlier detection
    min_chain_size <- 100  # Minimum reasonable chain size
    if (min(sapply(chains, nrow)) < min_chain_size) {
      message("Chains too small for reliable outlier detection in ", model_id, ", skipping outlier removal")
    } else {
      intercept_col <- grep("^intercept", colnames(chains[[1]]), value = TRUE)[1]
      means <- lapply(chains, function(x) mean(x[, intercept_col], na.rm=T))
      scaled_means = scale(unlist(means))
      potential_outlier <- which(abs(scaled_means) > 1.3)
      if (length(potential_outlier) %in% c(1,2)){
        chains_without_outlier <- chains[-c(potential_outlier)]
        new_gelman= gelman.diag(chains_without_outlier, multivariate = F)[[1]][,1] %>%  mean(na.rm=T)
        old_gelman= gelman.diag(chains, multivariate = F)[[1]][,1] %>%  mean(na.rm=T)
        improvement = old_gelman - new_gelman
        remove <- ifelse(improvement > .1, T, F)
        if (remove) {
          message(model_id, " removing outlier chain: ", potential_outlier,
                  "\nGelman diagnostic improves from ", round(old_gelman, 3), " to ", round(new_gelman, 3))
          out$samples = chains_without_outlier
          out$samples2 = out$samples2[-c(potential_outlier)]
        }
      }
    }
    
    param_summary <- fast.summary.mcmc(out$samples)
    plot_summary <- fast.summary.mcmc(out$samples2)
    es <- effectiveSize(out$samples)
    if (length(out$samples) > 1) {
      gelman_out <- cbind(gelman.diag(out$samples, multivariate = F)[[1]], es)
    } else gelman_out = cbind(`Point est.`=NA, `Upper C.I.`=NA, es)
    
    # Combine and save output
    out_summary <- list(
      samples = out$samples,
      param_summary = param_summary,
      plot_summary = plot_summary,
      metadata = out$metadata,
      gelman = gelman_out
    )
    
    saveRDS(out_summary, savepath, compress = F)
    message("Saved combined samples output for ", model_id, " to: ", savepath)
    
    if (min(es) < 500) {
      message("Low effective sample sizes - check for unconverged parameters in model: ", model_id)
      print(head(gelman_out, 50))
    }
    
    # If the summary now exists, delete the chains
    if (delete_samples_files){
      if (samples_exists(chain_paths[[1]])) {
        unlink(chain_paths)
        message("Deleting samples files, e.g.: ", chain_paths[[1]])
      }
    }
  }
}

cat("\n=== All Dirichlet models completed ===\n")
cat("Total models processed:", length(model_id_list), "\n")
