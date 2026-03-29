# Combine chains from each Dirichlet model, and create basic summary stats
source("../../source.R")
source(here::here("analysis/model_analysis/dirichlet_covariance/dirichlet_helper_functions.r"))

# Memory-efficient chain combination: loads one chain at a time,
# keeps only the second half (post-warmup), thins, then computes summary stats.
# This avoids holding all raw samples in memory at once.
combine_dirichlet_chains_lowmem <- function(chain_paths, model_id,
                                            keep_frac = 0.5, thin_by = 5) {
  library(coda)
  message("Combining ", length(chain_paths), " chains for ", model_id,
          " (keep_frac=", keep_frac, ", thin_by=", thin_by, ")")

  param_chains <- list()
  plot_chains <- list()
  metadata <- NULL

  for (i in seq_along(chain_paths)) {
    message("  Loading chain ", i, ": ", basename(chain_paths[i]))
    chain_data <- readRDS(chain_paths[i])

    samp <- chain_data$samples
    n_total <- nrow(samp)
    start_row <- floor(n_total * (1 - keep_frac)) + 1
    keep_idx <- seq(start_row, n_total, by = thin_by)
    message("    samples: ", n_total, " rows -> keeping ", length(keep_idx),
            " (rows ", start_row, "-", n_total, ", thin=", thin_by, ")")
    param_chains[[i]] <- mcmc(samp[keep_idx, , drop = FALSE])
    rm(samp)

    if (!is.null(chain_data$samples2)) {
      s2 <- chain_data$samples2
      n2 <- nrow(s2)
      # samples2 may have different thinning; use same fraction
      start2 <- floor(n2 * (1 - keep_frac)) + 1
      thin2 <- max(1, round(thin_by * n2 / n_total))
      keep2 <- seq(start2, n2, by = thin2)
      message("    samples2: ", n2, " rows -> keeping ", length(keep2))
      plot_chains[[i]] <- mcmc(s2[keep2, , drop = FALSE])
      rm(s2)
    }

    if (is.null(metadata) && !is.null(chain_data$metadata)) {
      metadata <- chain_data$metadata
    }

    rm(chain_data); gc()
  }

  # Trim all chains to the length of the shortest so mcmc.list() works
  # (required for Gelman-Rubin Rhat across chains)
  min_param <- min(sapply(param_chains, nrow))
  param_chains <- lapply(param_chains, function(ch) mcmc(ch[(nrow(ch) - min_param + 1):nrow(ch), , drop = FALSE]))
  message("  Trimmed param chains to ", min_param, " rows each")
  combined_params <- mcmc.list(param_chains)

  if (length(plot_chains) > 0) {
    min_plot <- min(sapply(plot_chains, nrow))
    plot_chains <- lapply(plot_chains, function(ch) mcmc(ch[(nrow(ch) - min_plot + 1):nrow(ch), , drop = FALSE]))
    message("  Trimmed plot chains to ", min_plot, " rows each")
    combined_plots <- mcmc.list(plot_chains)
  } else {
    combined_plots <- NULL
  }

  param_sum <- fast.summary.mcmc(combined_params)
  param_summary <- list(param_sum$statistics, param_sum$quantiles)

  if (!is.null(combined_plots)) {
    plot_sum <- fast.summary.mcmc(combined_plots)
    plot_summary <- list(plot_sum$statistics, plot_sum$quantiles)
  } else {
    plot_summary <- list(data.frame(), data.frame())
  }

  if (is.null(metadata)) metadata <- list(model_id = model_id)

  return(list(
    samples = combined_params,
    samples2 = combined_plots,
    param_summary = param_summary,
    plot_summary = plot_summary,
    metadata = metadata
  ))
}

# For Dirichlet models: use the reparameterized 75k warm-started outputs
target_dir <- here("data/model_outputs/dirichlet_driver_uncertainty_reparam_75k/")

cat("Processing Dirichlet regression models in:", target_dir, "\n")

# Find all chain files: final samples_* AND checkpoint files
all_rds <- list.files(path = target_dir,
                      pattern = "\\.rds$",
                      recursive = TRUE,
                      full.names = TRUE)
all_rds <- all_rds[!grepl("progress_|warmstart_", basename(all_rds))]

# For each chain, prefer the final samples file; fall back to the latest checkpoint.
# Extract (model_id, chain_num) from each file, then pick the best per chain.
parse_chain_info <- function(f) {
  bn <- basename(f)
  is_final <- grepl("^samples_", bn)
  # Extract chain number
  chain_num <- as.integer(sub(".*chain(\\d+).*", "\\1", bn))
  # Extract model_id: strip prefix, chain suffix, and checkpoint suffix
  model_id <- bn
  model_id <- sub("^(samples|checkpoint)_", "", model_id)
  model_id <- sub("_chain\\d+.*\\.rds$", "", model_id)
  # For checkpoints, extract loop number for sorting
  loop_num <- if (grepl("_loop(\\d+)", bn)) as.integer(sub(".*_loop(\\d+).*", "\\1", bn)) else NA
  data.frame(file = f, model_id = model_id, chain_num = chain_num,
             is_final = is_final, loop_num = loop_num, stringsAsFactors = FALSE)
}

file_info <- do.call(rbind, lapply(all_rds, parse_chain_info))

# For each (model_id, chain_num), pick: final samples if available, else latest checkpoint
best_per_chain <- file_info %>%
  group_by(model_id, chain_num) %>%
  arrange(desc(is_final), desc(loop_num)) %>%
  slice(1) %>%
  ungroup()

file.list <- best_per_chain$file
model_id_list <- unique(best_per_chain$model_id)

cat("Total models to process:", length(model_id_list), "\n\n")
cat("Model IDs:", paste(model_id_list, collapse = ", "), "\n")

# Use sequential processing for debugging
cat("Starting sequential processing with", length(model_id_list), "models...\n")

for(model_id in model_id_list) {
  # Do we want to keep all the chain files separately? Deleting them will save space
  delete_samples_files = F

  # Subset to chain files for this model_id (already filtered to best per chain)
  chain_paths <- best_per_chain$file[best_per_chain$model_id == model_id]

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
  # Build savepath from model_id in the parent directory of the chain files
  savepath <- file.path(dirname(chain_paths[[1]]),
                        paste0("samples_", model_id, ".rds"))

  # Don't run loop if the combined samples file already exists
  if (file.exists(savepath)) {
    message("Summary samples file already exists")
  } else {
    # Calculate summary on each output subset, using custom summary function
    message("Combining chains for ", model_id, " with ", length(chain_paths), " chains...")
    out <- tryCatch({
      combine_dirichlet_chains_lowmem(chain_paths, model_id)
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
