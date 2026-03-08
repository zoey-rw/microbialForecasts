# Combine chains from each taxon model, and create basic summary stats
# TRUNCATED NORMAL VERSION - Adapted for truncated normal model outputs
source("../../source.R")

# Define model types and time periods for truncated normal models
model_types <- c("cycl_only", "env_cycl", "env_cov")
time_periods <- c("20130601_20151101", "20151101_20180101", "20130601_20180101", "20130601_20200101")

cat("Processing truncated normal models with legacy effects\n")
cat("Model types:", paste(model_types, collapse = ", "), "\n")
cat("Time periods:", paste(time_periods, collapse = ", "), "\n\n")

# Get all available model outputs from truncated_normal directory
base_path <- here("data/model_outputs/truncated_normal")
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
    # Use the combine_chains function that handles different sample sizes robustly
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
    if ("model_structure" %in% metadata_components) {
      message("   ✅ Model structure preserved: ", out$metadata$model_structure)
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
  
  # Check if core_sd parameter exists (truncNorm specific)
  if (!"core_sd" %in% colnames(chains[[1]])) {
    message("No core_sd parameter found in chains for ", basename(model_id), ", skipping...")
    next
  }
  
  # Ensure chains have reasonable sample sizes before outlier detection
  min_chain_size <- 100
  chain_sizes <- sapply(chains, nrow)
  if (min(chain_sizes) < min_chain_size) {
    message("Chains too small for reliable outlier detection in ", basename(model_id), ", skipping outlier removal")
  } else {
    # Outlier detection for truncNorm models
    message("Checking for outlier chains in ", basename(model_id), "...")
    
    # Calculate median values for key parameters
    key_params <- c("core_sd", "sigma", "intercept", "rho")
    available_params <- key_params[key_params %in% colnames(chains[[1]])]
    
    if (length(available_params) > 0) {
      param_medians <- matrix(NA, nrow = length(chains), ncol = length(available_params))
      colnames(param_medians) <- available_params
      
      for (i in 1:length(chains)) {
        for (j in 1:length(available_params)) {
          param <- available_params[j]
          param_medians[i, j] <- median(chains[[i]][, param], na.rm = TRUE)
        }
      }
      
      # Calculate Mahalanobis distance for outlier detection
      tryCatch({
        mahal_dist <- mahalanobis(param_medians, 
                                  colMeans(param_medians, na.rm = TRUE), 
                                  cov(param_medians, use = "complete.obs"))
        
        # Identify outliers (chains with Mahalanobis distance > 3)
        outlier_threshold <- 3
        outliers <- which(mahal_dist > outlier_threshold)
        
        if (length(outliers) > 0) {
          message("Found ", length(outliers), " outlier chain(s) in ", basename(model_id), ": ", paste(outliers, collapse = ", "))
          message("Removing outlier chains...")
          
          # Remove outlier chains
          out$samples <- out$samples[-outliers]
          message("Remaining chains after outlier removal: ", length(out$samples))
        } else {
          message("No outlier chains detected in ", basename(model_id))
        }
      }, error = function(e) {
        message("Outlier detection failed for ", basename(model_id), ": ", e$message)
        message("Proceeding with all chains...")
      })
    }
  }
  
  # Save combined output
  message("Saving combined chains for ", basename(model_id), "...")
  tryCatch({
    saveRDS(out, combined_file_path)
    message("✅ Successfully saved combined chains to: ", basename(combined_file_path))
  }, error = function(e) {
    message("❌ Failed to save combined chains for ", basename(model_id), ": ", e$message)
  })
}

cat("\nTruncated normal chain combination completed!\n")
cat("Combined files saved in: ", base_path, "\n")
