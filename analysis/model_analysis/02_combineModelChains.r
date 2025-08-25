# Combine chains from each taxon model, and create basic summary stats
source("../../source.R")

# Define model types and time periods for logit beta regression
model_types <- c("cycl_only", "env_cycl", "env_cov")
time_periods <- c("20130601_20151101", "20151101_20180101", "20130601_20180101", "20130601_20200101")

cat("Processing logit beta regression models with legacy effects\n")
cat("Model types:", paste(model_types, collapse = ", "), "\n")
cat("Time periods:", paste(time_periods, collapse = ", "), "\n\n")

# Get all available model outputs from logit_beta_regression directory
base_path <- here("data/model_outputs/logit_beta_regression")
file.list <- list.files(path = base_path,
                       pattern = "_chain",
                       recursive = TRUE,
                       full.names = TRUE)

# Subset to files larger than 100KB
info <- file.info(file.list)
large_enough <- rownames(info[which(info$size > 100000), ])
file.list <- file.list[file.list %in% large_enough]

# Don't want files still being written - at least 1 min old
older <- rownames(info[which(info$mtime < (Sys.time()-60)), ])
file.list <- file.list[file.list %in% older]

cat("Total chain files found:", length(file.list), "\n")
cat("First few files:", paste(basename(head(file.list, 3)), collapse = ", "), "\n\n")

# Group files by model ID (removing chain suffix)
model_files <- list()
for (file in file.list) {
  # Extract model ID by removing chain suffix
  model_id <- gsub("_chain[1234567]", "", file)
  if (!(model_id %in% names(model_files))) {
    model_files[[model_id]] <- list()
  }
  model_files[[model_id]] <- c(model_files[[model_id]], file)
}

# Process each model that has multiple chains
for (model_id in names(model_files)) {
  chain_paths <- model_files[[model_id]]
  
  if (length(chain_paths) < 2) {
    message("Skipping ", basename(model_id), "; need at least 2 chains to combine")
    next
  }
  
  message("Processing ", basename(model_id), " with ", length(chain_paths), " chains...")
  
  # Check if combined file already exists
  # The combined file should have the same name as the model_id but without _chain suffix
  # and should be different from the original individual model files
  combined_file_path <- gsub("_chain[1234567]", "", chain_paths[[1]])
  
  # Check if this is actually a combined file (not the original individual model file)
  # The combined file should have been created by this script, so we'll check if it exists
  # and has the expected structure
  if (file.exists(combined_file_path)) {
    # Check if it's actually a combined file by examining its structure
    tryCatch({
      test_read <- readRDS(combined_file_path)
      if (is.list(test_read) && "samples" %in% names(test_read) && 
          is.list(test_read$samples) && length(test_read$samples) > 1) {
        message("Combined file already exists for ", basename(model_id))
        next
      } else {
        message("File exists but is not a combined file, will overwrite: ", basename(model_id))
      }
    }, error = function(e) {
      message("Error reading existing file, will overwrite: ", basename(model_id))
    })
  }
  
  # Combine chains
  message("Combining chains for ", basename(model_id), "...")
  out <- tryCatch({
    # Use the new combine_chains function that handles different sample sizes robustly
    combine_chains(chain_paths)
  }, error = function(e) {
    message("Error combining chains for ", basename(model_id), ": ", e$message)
    return(NULL)
  })
  
  if (is.null(out)) {
    message("Failed to combine chains for ", basename(model_id), ", skipping...")
    next
  }
  
  # Verify metadata preservation including Nimble model code
  if (!is.null(out$metadata)) {
    message("✅ Metadata preserved during chain combination")
    
    # Check for key metadata components
    metadata_components <- names(out$metadata)
    message("   Metadata components found: ", paste(metadata_components, collapse = ", "))
    
    # Verify Nimble model code is preserved
    if ("nimble_code" %in% metadata_components) {
      message("   ✅ Nimble model code preserved")
    } else if ("model_code" %in% metadata_components) {
      message("   ✅ Model code preserved")
    } else {
      message("   ⚠️  Model code component not found in metadata")
    }
    
    # Verify other important metadata components
    if ("model_name" %in% metadata_components) {
      message("   ✅ Model name preserved: ", out$metadata$model_name)
    }
    if ("use_legacy_covariate" %in% metadata_components) {
      message("   ✅ Legacy covariate flag preserved: ", out$metadata$use_legacy_covariate)
    }
  } else {
    message("❌ ERROR: Metadata was lost during chain combination!")
    next
  }
  
  # Remove chains if one has very different values than the others
  chains <- out$samples
  if (length(chains) == 0) {
    message("No chains found in output for ", basename(model_id), ", skipping...")
    next
  }
  
  # Check if intercept parameter exists
  if (!"intercept" %in% colnames(chains[[1]])) {
    message("No intercept parameter found in chains for ", basename(model_id), ", skipping...")
    next
  }
  
  # Ensure chains have reasonable sample sizes before outlier detection
  min_chain_size <- 100
  chain_sizes <- sapply(chains, nrow)
  if (min(chain_sizes) < min_chain_size) {
    message("Chains too small for reliable outlier detection in ", basename(model_id), ", skipping outlier removal")
  } else {
    means <- lapply(chains, function(x) mean(x[,"intercept"], na.rm=TRUE))
    scaled_means <- scale(unlist(means))
    potential_outlier <- which(abs(scaled_means) > 1.3)
    if (length(potential_outlier) %in% c(1,2)) {
      chains_without_outlier <- chains[-c(potential_outlier)]
      new_gelman <- gelman.diag(chains_without_outlier, multivariate = FALSE)[[1]][,1] %>% mean(na.rm=TRUE)
      old_gelman <- gelman.diag(chains, multivariate = FALSE)[[1]][,1] %>% mean(na.rm=TRUE)
      improvement <- old_gelman - new_gelman
      remove <- ifelse(improvement > 0.1, TRUE, FALSE)
      if (remove) {
        message(basename(model_id), " removing outlier chain: ", potential_outlier,
                "\nGelman diagnostic improves from ", round(old_gelman, 3), " to ", round(new_gelman, 3))
        out$samples <- chains_without_outlier
        out$samples2 <- out$samples2[-c(potential_outlier)]
      }
    }
  }
  
  # Calculate summaries
  param_summary <- fast.summary.mcmc(out$samples)
  plot_summary <- fast.summary.mcmc(out$samples2)
  es <- effectiveSize(out$samples)
  
  if (length(out$samples) > 1) {
    gelman_out <- cbind(gelman.diag(out$samples, multivariate = FALSE)[[1]], es)
  } else {
    gelman_out <- cbind(`Point est.`=NA, `Upper C.I.`=NA, es)
  }
  
  # Combine and save output - ensure metadata is preserved
  out_summary <- list(
    samples = out$samples,
    param_summary = param_summary,
    plot_summary = plot_summary,
    metadata = out$metadata,  # This preserves all metadata including Nimble model code
    gelman = gelman_out
  )
  
  # Final verification that metadata is preserved in output
  if (!is.null(out_summary$metadata)) {
    message("✅ Final verification: Metadata preserved in output summary")
  } else {
    message("❌ ERROR: Metadata lost in final output summary!")
  }
  
  # Save to the same directory as the chains
  savepath <- gsub("_chain[1234567]", "", chain_paths[[1]])
  saveRDS(out_summary, savepath, compress = FALSE)
  message("Saved combined samples output for ", basename(model_id), " to: ", basename(savepath))
  
  if (min(es) < 500) {
    message("Low effective sample sizes - check for unconverged parameters in model: ", basename(model_id))
    print(head(gelman_out, 50))
  }
}

cat("\nChain combining complete!\n")


