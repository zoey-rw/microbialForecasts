# Combine chains from each taxon model, and create basic summary stats
# FIXED: Updated for CLR models with proper error handling and CLR-specific logic
source("../../source.R")
project_root <- here()

# Function to filter files to only include those with both 'with_legacy_covariate' and 'clr'
# AND exclude checkpoint files (we only want final combined sample files)
filter_standard_files <- function(file_list) {
  if (length(file_list) == 0) return(file_list)
  
  # For CLR models: filter for chain files with legacy covariate
  # Exclude checkpoints and already-combined files
  standard_files <- file_list[
    grepl('_chain[0-9]+\\.rds$', basename(file_list)) &  # Is a chain file
    grepl('with_legacy_covariate', basename(file_list)) &  # Has legacy covariate
    !grepl('checkpoint', basename(file_list))  # Not a checkpoint
  ]
  
  cat('File filtering applied:\n')
  cat('  Original files:', length(file_list), '\n')
  cat('  Standard CLR chain files:', length(standard_files), '\n')
  cat('  Filtered out:', length(file_list) - length(standard_files), '\n\n')
  
  return(standard_files)
}

# Find parameters that have 0 non-NA cases or zero variance (causes lm.fit / effectiveSize / ar to fail)
check_chain_columns_for_summary <- function(chain_list, label = "samples") {
  x <- as.matrix(chain_list[[1]])
  vnames <- colnames(x)
  bad <- character(0)
  for (j in seq_len(ncol(x))) {
    col_j <- x[, j]
    n_ok <- sum(!is.na(col_j))
    if (n_ok == 0) {
      bad <- c(bad, vnames[j])
    } else if (n_ok > 0 && (is.na(var(col_j, na.rm = TRUE)) || var(col_j, na.rm = TRUE) == 0)) {
      bad <- c(bad, vnames[j])
    }
  }
  # Check across all chains (in case one chain is fine but another has all NA for a param)
  if (length(chain_list) > 1) {
    for (i in seq_along(chain_list)) {
      x_i <- as.matrix(chain_list[[i]])
      for (j in seq_len(ncol(x_i))) {
        if (vnames[j] %in% bad) next
        col_j <- x_i[, j]
        n_ok <- sum(!is.na(col_j))
        if (n_ok == 0) bad <- c(bad, vnames[j])
        else if (is.na(var(col_j, na.rm = TRUE)) || var(col_j, na.rm = TRUE) == 0) bad <- c(bad, vnames[j])
      }
    }
  }
  unique(bad)
}

# Safe effective sample size per column (NA for flat/zero-variance columns; avoids lm.fit crash)
safe_effectiveSize <- function(mat) {
  if (is.list(mat) && inherits(mat[[1]], "mcmc")) {
    mat <- do.call(rbind, lapply(mat, as.matrix))
  }
  mat <- as.matrix(mat)
  apply(mat, 2, function(x) {
    if (sum(!is.na(x)) < 2 || is.na(var(x, na.rm = TRUE)) || var(x, na.rm = TRUE) == 0) {
      return(NA_real_)
    }
    tryCatch(coda::effectiveSize(coda::as.mcmc(x)), error = function(e) NA_real_)
  })
}

# Fallback when fast.summary.mcmc fails (e.g. zero-variance columns); builds statistics and quantiles per column
fallback_summary_mcmc <- function(chain_list) {
  x <- if (length(chain_list) == 1) as.matrix(chain_list[[1]]) else as.matrix(do.call(rbind, lapply(chain_list, as.matrix)))
  nvar <- ncol(x)
  vnames <- colnames(x)
  statnames <- c("Mean", "SD", "Naive SE", "Time-series SE")
  varstats <- matrix(NA_real_, nrow = nvar, ncol = length(statnames), dimnames = list(vnames, statnames))
  qu <- matrix(NA_real_, nrow = nvar, ncol = 5, dimnames = list(vnames, c("2.5%", "25%", "50%", "75%", "97.5%")))
  ntot <- nrow(x)
  for (j in seq_len(nvar)) {
    v <- x[, j]
    n_ok <- sum(!is.na(v))
    if (n_ok < 2 || is.na(var(v, na.rm = TRUE)) || var(v, na.rm = TRUE) == 0) next
    varstats[j, "Mean"] <- mean(v, na.rm = TRUE)
    varstats[j, "SD"] <- sd(v, na.rm = TRUE)
    varstats[j, "Naive SE"] <- varstats[j, "SD"] / sqrt(ntot)
    varstats[j, "Time-series SE"] <- NA_real_
    qu[j, ] <- quantile(v, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
  }
  out <- list(statistics = varstats, quantiles = qu)
  class(out) <- "summary.mcmc"
  out
}

# Process all CLR model types (cycl_only, env_cov, env_cycl)
# No need to hardcode model_name as we process all models found

# Use CLR model IDs from CLR_regression directory
# Get all model IDs from the CLR_regression directory by looking at individual chain files
clr_output_dir <- here("data/model_outputs/CLR_regression")
clr_chain_files <- list.files(clr_output_dir, pattern = "samples_.*_chain[0-9]+\\.rds$", 
                              recursive = TRUE, full.names = FALSE)

# Extract model IDs from chain file names (remove _chain[0-9] suffix)
model_id_list <- unique(gsub("samples_(.*)_chain[0-9]+\\.rds", "\\1", clr_chain_files))
model_id_list <- gsub(".*/", "", model_id_list)  # Remove directory path

cat("Found", length(model_id_list), "CLR models in CLR_regression\n")
cat("Total models to process:", length(model_id_list), "\n\n")

# Set up parallel processing (1–2 cores to limit RAM; each worker holds full chain data)
n_cores <- min(2, detectCores() - 1)
cl <- makeCluster(n_cores, outfile = "")
registerDoParallel(cl)
cat("Using", n_cores, "cores for parallel processing (low core count to reduce RAM)\n")

	# Process each model ID in parallel (.errorhandling = "pass" returns error objects; we check them below)
	combine_results <- foreach(model_id=model_id_list, .errorhandling = "pass", .export = c("project_root", "filter_standard_files", "check_chain_columns_for_summary", "safe_effectiveSize", "fallback_summary_mcmc")) %dopar% {
		source(file.path(project_root, "source.R"))
		
		# Do we want to keep all the chain files separately? Deleting them will save space
		delete_samples_files = F
	
	file.list <- list.files(path = here("data/model_outputs/CLR_regression/"),
													pattern = "_chain",
													recursive = T,
													full.names = T)
	
	# Subset to newest output files (revised 01: only chains from last 7 days)
	info <- file.info(file.list)
	
	# Subset more than 1KB (for our test files, lower threshold)
	large_enough <- rownames(info[which(info$size > 1000), ])
	
	# Only use chain files modified in the last 7 days (so env_cycl and other recent 01 runs are included)
	newer <- rownames(info[which(info$mtime > (Sys.time() - as.difftime(7, units = "days"))), ])
	
	# Don't want files still being written - at least 1 min old
	older <- rownames(info[which(info$mtime < (Sys.time()-60)), ])
	
	file.list <- file.list[file.list %in% newer & file.list %in% large_enough & file.list %in% older]
	
	# Apply filtering to only process standard files
	file.list <- filter_standard_files(file.list)
	
	message("Searching model outputs for ",model_id)
	
	# Subset to files of interest
	chain_paths <- file.list[grepl(model_id, file.list, fixed = T)]
	
	if (length(chain_paths) == 0)
		return(message("Skipping ", model_id, "; no chains"))
	
	# FIXED: Handle single chain models properly for CLR
	if (length(chain_paths) == 1) {
		cat("Single chain model, copying to combined samples file...\n")
		single_chain_path <- chain_paths[[1]]
		combined_samples_path <- gsub("_chain[1234567]", "", single_chain_path)
		
		if (file.exists(combined_samples_path)) {
			cat("Combined samples file already exists, skipping...\n")
			return(message("Skipping ", model_id, "; combined samples already exist"))
		}
		
		cat("Copying single chain to combined samples file...\n")
		chain_data <- readRDS(single_chain_path)
		
		# FIXED: Create the same output structure as multi-chain models for CLR
		if (is.list(chain_data) && "samples" %in% names(chain_data)) {
			# New format with metadata (CLR models)
			if ("samples2" %in% names(chain_data)) {
				out <- list(
					samples = list(mcmc(chain_data$samples, start = 1, end = nrow(chain_data$samples), thin = 1)),
					samples2 = list(mcmc(chain_data$samples2, start = 1, end = nrow(chain_data$samples2), thin = 1)),
					metadata = chain_data$metadata
				)
			} else {
				# CLR models use same samples for both parameter and plot predictions
				out <- list(
					samples = list(mcmc(chain_data$samples, start = 1, end = nrow(chain_data$samples), thin = 1)),
					samples2 = list(mcmc(chain_data$samples, start = 1, end = nrow(chain_data$samples), thin = 1)),
					metadata = chain_data$metadata
				)
			}
		} else {
			# Old format (raw matrix) - handle gracefully
			out <- list(
				samples = list(mcmc(chain_data, start = 1, end = nrow(chain_data), thin = 1)),
				samples2 = list(mcmc(chain_data, start = 1, end = nrow(chain_data), thin = 1)),
				metadata = list(model_id = model_id)
			)
		}
		
		# Single chain: warn (do not stop) on flat params; unmodeled eta padding has 0 variance by design
		bad_param <- check_chain_columns_for_summary(out$samples, "param")
		bad_plot <- check_chain_columns_for_summary(out$samples2, "plot")
		if (length(bad_param) > 0 || length(bad_plot) > 0) {
			message("WARNING: Some parameters/latent states have 0 variance (expected for unmodeled eta padding). model_id: ", model_id)
		}
		cat("Creating parameter summaries for single chain...\n")
		param_summary <- tryCatch(
			fast.summary.mcmc(out$samples),
			error = function(e) {
				message("fast.summary.mcmc failed (e.g. flat columns): ", conditionMessage(e), " — using fallback")
				fallback_summary_mcmc(out$samples)
			}
		)
		cat("Creating plot summaries for single chain...\n")
		plot_summary <- tryCatch(
			fast.summary.mcmc(out$samples2),
			error = function(e) {
				message("fast.summary.mcmc failed for plot: ", conditionMessage(e), " — using fallback")
				fallback_summary_mcmc(out$samples2)
			}
		)
		cat("Calculating effective sample sizes...\n")
		es <- safe_effectiveSize(out$samples[[1]])
		
		# For single chain, Gelman diagnostics are not applicable
		cat("Single chain, setting Gelman diagnostics to NA...\n")
		gelman_out <- cbind(`Point est.` = NA, `Upper C.I.` = NA, es)
		
		# Combine and save output
		out_summary <- list(
			samples = out$samples,
			param_summary = param_summary,
			plot_summary = plot_summary,
			metadata = out$metadata,
			gelman = gelman_out
		)
		
		cat("Saving combined samples file...\n")
		saveRDS(out_summary, combined_samples_path)
		cat("Single chain model processed successfully!\n")
		return(message("Completed single chain model: ", model_id))
	}
	
	savepath <- gsub("_chain[1234567]", "", chain_paths[[1]])
	
	# Don't run loop if the samples file is already newer than the chain files
	cat("Checking if samples already exist...\n")
	if (samples_exists(chain_paths[[1]])) {
		message("Summary samples file already exists")
	} else {
		cat("No existing samples file, proceeding with combination...\n")
		
		# FIXED: For CLR outputs: handle the proper structure with samples and metadata
		cat("Loading", length(chain_paths), "chains...\n")
		chains <- list()
		chains2 <- list()
		metadata_list <- list()
		
		for (i in seq_along(chain_paths)) {
			cat("  Loading chain", i, ":", basename(chain_paths[i]), "\n")
			chain_data <- readRDS(chain_paths[i])
			
			if (is.list(chain_data) && "samples" %in% names(chain_data)) {
				# New format with metadata (CLR models)
				cat("    Chain", i, "dimensions:", dim(chain_data$samples), "\n")
				chains[[i]] <- chain_data$samples
				
				# Handle samples2 - CLR models may have same samples for both
				if ("samples2" %in% names(chain_data)) {
					chains2[[i]] <- chain_data$samples2
					cat("    ✓ samples2 found\n")
				} else {
					# CLR models use same samples for both parameter and plot predictions
					chains2[[i]] <- chain_data$samples
					cat("    ✓ Using samples for samples2 (CLR model)\n")
				}
				
				metadata_list[[i]] <- chain_data$metadata
				cat("    ✓ Metadata loaded\n")
			} else {
				# Old format (raw matrix) - handle gracefully
				cat("    Chain", i, "dimensions:", dim(chain_data), "\n")
				chains[[i]] <- chain_data
				chains2[[i]] <- chain_data
				metadata_list[[i]] <- list(model_id = model_id)
				cat("    ✓ Old format handled\n")
			}
		}
		
		cat("Creating output structure...\n")
		
		# Validate that we have valid chains
		if (length(chains) == 0) {
			cat("ERROR: No valid chains loaded for", model_id, "\n")
			return(message("ERROR: No valid chains loaded for ", model_id))
		}
		
		cat("Successfully loaded", length(chains), "chains\n")
		
		# Validate metadata structure matches between CLR and beta models
		if (!is.null(metadata_list[[1]]) && length(metadata_list[[1]]) > 0) {
			required_fields <- c("model_name", "model_id", "use_legacy_covariate")
			existing_fields <- intersect(required_fields, names(metadata_list[[1]]))
			missing <- setdiff(required_fields, existing_fields)
			if (length(missing) > 0) {
				warning("Metadata missing fields: ", paste(missing, collapse=", "))
			} else {
				cat("✓ Metadata validation passed\n")
			}
		} else {
			cat("WARNING: No metadata available for validation\n")
		}
		
		# Validate chain dimensions
		chain_dims <- sapply(chains, dim)
		cat("Chain dimensions:", paste(apply(chain_dims, 2, paste, collapse="x"), collapse=", "), "\n")
		
		# Check if all chains have the same dimensions
		if (length(unique(apply(chain_dims, 2, paste, collapse="x"))) > 1) {
			cat("WARNING: Chains have different dimensions, this may cause issues\n")
		}
		
		# Convert chains to mcmc objects for compatibility with summary functions
		cat("Converting chains to mcmc format...\n")
		mcmc_chains <- lapply(chains, function(x) {
			mcmc(x, start = 1, end = nrow(x), thin = 1)
		})
		mcmc_chains2 <- lapply(chains2, function(x) {
			mcmc(x, start = 1, end = nrow(x), thin = 1)
		})
		
		out <- list(
			samples = mcmc_chains,
			samples2 = mcmc_chains2,
			metadata = metadata_list[[1]]  # Use metadata from first chain
		)
		
		# FIXED: Process chains for outlier detection with better error handling
		cat("Processing chains for outlier detection...\n")
		chains <- out$samples
		cat("Calculating chain means for outlier detection...\n")
		
		# Check if intercept parameter exists before calculating means
		if ("intercept" %in% colnames(chains[[1]])) {
			means <- lapply(chains, function(x) mean(x[,"intercept"], na.rm=T))
			cat("Chain means:", paste(round(unlist(means), 4), collapse = ", "), "\n")
			scaled_means = scale(unlist(means))
			cat("Scaled means:", paste(round(scaled_means, 4), collapse = ", "), "\n")
			potential_outlier <- which(abs(scaled_means) > 1.3)
			cat("Potential outliers:", potential_outlier, "\n")
			
			if (length(potential_outlier) %in% c(1,2)){
				chains_without_outlier <- chains[-c(potential_outlier)]
				new_gelman <- tryCatch(gelman.diag(chains_without_outlier, multivariate = F)[[1]][,1] %>% mean(na.rm = T), error = function(e) NA_real_)
				old_gelman <- tryCatch(gelman.diag(chains, multivariate = F)[[1]][,1] %>% mean(na.rm = T), error = function(e) NA_real_)
				improvement = old_gelman - new_gelman
				remove <- ifelse(improvement > .1, T, F)
				if (remove) {
					message(model_id, " removing outlier chain: ", potential_outlier,
									"\nGelman diagnostic improves from ", round(old_gelman, 3), " to ", round(new_gelman, 3))
					out$samples = chains_without_outlier
					out$samples2 = out$samples2[-c(potential_outlier)]
				}
			}
		} else {
			cat("WARNING: No intercept parameter found, skipping outlier detection\n")
		}
		
		# Warn (do not stop) on flat params; unmodeled eta padding has 0 variance by design
		bad_param <- check_chain_columns_for_summary(out$samples, "param")
		bad_plot <- check_chain_columns_for_summary(out$samples2, "plot")
		if (length(bad_param) > 0 || length(bad_plot) > 0) {
			message("WARNING: Some parameters/latent states have 0 variance (expected for unmodeled eta padding). model_id: ", model_id)
		}
		cat("Creating parameter summaries...\n")
		param_summary <- tryCatch(
			fast.summary.mcmc(out$samples),
			error = function(e) {
				message("fast.summary.mcmc failed (e.g. flat columns): ", conditionMessage(e), " — using fallback")
				fallback_summary_mcmc(out$samples)
			}
		)
		cat("Creating plot summaries...\n")
		plot_summary <- tryCatch(
			fast.summary.mcmc(out$samples2),
			error = function(e) {
				message("fast.summary.mcmc failed for plot: ", conditionMessage(e), " — using fallback")
				fallback_summary_mcmc(out$samples2)
			}
		)
		cat("Calculating effective sample sizes...\n")
		es <- safe_effectiveSize(out$samples)
		
		if (length(out$samples) > 1) {
			cat("Calculating Gelman diagnostics...\n")
			gelman_out <- tryCatch(
				cbind(gelman.diag(out$samples, multivariate = F)[[1]], es),
				error = function(e) {
					message("Gelman diagnostic failed (e.g. 0 variance): ", conditionMessage(e), " — using NA for Point est. / Upper C.I.")
					cbind(`Point est.` = NA_real_, `Upper C.I.` = NA_real_, es)
				}
			)
		} else {
			cat("Single chain, setting Gelman diagnostics to NA...\n")
			gelman_out = cbind(`Point est.`=NA, `Upper C.I.`=NA, es)
		}
		
		# Combine and save output
		out_summary <- list(
			samples = out$samples,
			param_summary = param_summary,
			plot_summary = plot_summary,
			metadata = out$metadata,
			gelman = gelman_out
		)
		
		saveRDS(out_summary, savepath, compress = F)
		message("Saved combined samples output for ", model_id, " to: ",savepath)
		
		min_es <- min(es, na.rm = TRUE)
		if (length(es[!is.na(es)]) > 0 && is.finite(min_es) && min_es < 500) {
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
	return()
} # End model_id loop

stopCluster(cl)

# Treat worker errors as script failure: do not report success if any model failed
failed <- which(vapply(combine_results, function(x) inherits(x, "error"), logical(1)))
if (length(failed) > 0) {
  failed_ids <- model_id_list[failed]
  err_msgs <- vapply(combine_results[failed], function(x) conditionMessage(x), character(1))
  cat("ERROR: ", length(failed), " model(s) failed during combine/summary.\n", sep = "")
  cat("Failed model_id(s): ", paste(failed_ids, collapse = ", "), "\n", sep = "")
  cat("First error: ", err_msgs[1], "\n", sep = "")
  stop("02_combineModelChains_CLR.r: ", length(failed), " model(s) had invalid chains (all NA or zero variance). Fix 01_fitModels_CLR.r or data for these runs.")
}

cat("✅ CLR chain combination completed successfully!\n")
cat("Combined samples saved to:", clr_output_dir, "\n")


