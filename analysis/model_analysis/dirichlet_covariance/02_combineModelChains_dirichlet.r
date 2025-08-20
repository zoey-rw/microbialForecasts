# Combine chains from each Dirichlet model, and create basic summary stats
source("source.R")
source("dirichlet_helper_functions.r")

# Function to combine Dirichlet chains and create summary structure
combine_dirichlet_chains <- function(chain_paths, model_id) {
  message("Combining ", length(chain_paths), " chains for ", model_id)
  
  # Separate parameter samples and plot predictions
  param_chains <- chain_paths[grepl("samples_dirichlet", chain_paths)]
  plot_chains <- chain_paths[grepl("samples2_dirichlet", chain_paths)]
  
  message("Parameter chains: ", length(param_chains))
  message("Plot prediction chains: ", length(plot_chains))
  
  # Read parameter chains
  param_samples <- list()
  metadata_list <- list()
  
  for (i in 1:length(param_chains)) {
    chain_data <- readRDS(param_chains[[i]])
    message("Parameter chain ", i, " dimensions: ", nrow(chain_data), " x ", ncol(chain_data))
    
    # Extract metadata from the first chain
    if (i == 1 && is.list(chain_data) && "metadata" %in% names(chain_data)) {
      metadata_list[[i]] <- chain_data$metadata
      param_samples[[i]] <- chain_data$samples
    } else if (is.list(chain_data) && "metadata" %in% names(chain_data)) {
      metadata_list[[i]] <- chain_data$metadata
      param_samples[[i]] <- chain_data$samples
    } else {
      # Fallback: treat as raw matrix
      param_samples[[i]] <- chain_data
      metadata_list[[i]] <- list(model_id = model_id)
    }
  }
  
  # Read plot prediction chains if they exist
  plot_samples <- list()
  if (length(plot_chains) > 0) {
    for (i in 1:length(plot_chains)) {
      plot_data <- readRDS(plot_chains[[i]])
      if (is.list(plot_data) && "samples2" %in% names(plot_data)) {
        plot_samples[[i]] <- plot_data$samples2
      } else {
        plot_samples[[i]] <- plot_data
      }
    }
  }
  
  # Combine parameter chains using mcmc.list
  library(coda)
  param_mcmc_chains <- lapply(param_samples, mcmc)
  combined_param_chains <- mcmc.list(param_mcmc_chains)
  
  # Combine plot chains if they exist
  if (length(plot_samples) > 0) {
    plot_mcmc_chains <- lapply(plot_samples, mcmc)
    combined_plot_chains <- mcmc.list(plot_mcmc_chains)
  } else {
    combined_plot_chains <- NULL
  }
  
  # Calculate parameter summaries
  param_summary <- list()
  param_summary[[1]] <- summary(combined_param_chains)$statistics  # Means
  param_summary[[2]] <- summary(combined_param_chains)$quantiles   # Quantiles
  
  # Calculate plot summaries if available
  plot_summary <- list()
  if (!is.null(combined_plot_chains)) {
    plot_summary[[1]] <- summary(combined_plot_chains)$statistics  # Means
    plot_summary[[2]] <- summary(combined_plot_chains)$quantiles   # Quantiles
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
  
  # Create the combined output structure
  out <- list(
    samples = combined_param_chains,
    param_summary = param_summary,
    plot_summary = plot_summary,
    metadata = metadata
  )
  
  return(out)
}

# For Dirichlet models: focus on fungal phylum models for 2015-2018
model_name <- "env_cycl"  # Can be "env_cov", "env_cycl", or "cycl_only"

# For Dirichlet models: use phylum composition models
model_id_list = c(
  "env_cycl_phylum_fun_20151101_20180101",
  "env_cov_phylum_fun_20151101_20180101"
)

cat("Processing Dirichlet regression models for 2015-2018 data\n")
cat("Model type:", model_name, "\n")
cat("Total models to process:", length(model_id_list), "\n\n")

cat("Model ID list length:", length(model_id_list), "\n")
cat("Model IDs:", paste(model_id_list, collapse = ", "), "\n")

# Use sequential processing for debugging
cat("Starting sequential processing with", length(model_id_list), "models...\n")

for(model_id in model_id_list) {
  # Do we want to keep all the chain files separately? Deleting them will save space
  delete_samples_files = F
  
  # Find chain files for Dirichlet models
  file.list <- list.files(path = here("data/model_outputs/dirichlet_regression/"),
                          pattern = "_chain",
                          recursive = T,
                          full.names = T)
  
  # Subset to newest output files
  info <- file.info(file.list)
  
  # Subset more than 100KB (for our test files, reasonable threshold)
  large_enough <- rownames(info[which(info$size > 100000), ])
  
  # Remove date filtering for now - all files are recent
  newer <- rownames(info)  # Include all files
  
  # Don't want files still being written - at least 1 min old (for testing)
  older <- rownames(info[which(info$mtime < (Sys.time()-60)), ])
  
  file.list <- file.list[file.list %in% newer & file.list %in% large_enough]
  
  message("Searching model outputs for ", model_id)
  message("Total files in file.list: ", length(file.list))
  message("First few files: ", paste(head(file.list, 3), collapse = ", "))
  
  # Subset to files of interest
  chain_paths <- file.list[grepl(model_id, file.list, fixed = T)]
  
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
    
    # Check if intercept parameter exists (for Dirichlet models, this might be different)
    if (!"intercept" %in% colnames(chains[[1]])) {
      message("No intercept parameter found in chains for ", model_id, ", skipping...")
      next
    }
    
    # Ensure chains have reasonable sample sizes before outlier detection
    min_chain_size <- 100  # Minimum reasonable chain size
    if (min(sapply(chains, nrow)) < min_chain_size) {
      message("Chains too small for reliable outlier detection in ", model_id, ", skipping outlier removal")
    } else {
      means <- lapply(chains, function(x) mean(x[,"intercept"], na.rm=T))
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
