#' @title 			combine_chains
#' @description Combine MCMC chains with robust handling of different sample sizes and metadata preservation
#' @export
combine_chains <- function(chain_paths,
																	save = FALSE,
																	cut_size1 = NULL,
																	cut_size2 = NULL){
	require(coda)
	require(tidyverse)

	# cut_size1 and cut_size2 default to NULL (no limit) if not specified
	# NULL means use minimum size among chains (no truncation beyond matching sizes)
	# Non-NULL values will truncate chains to that size
	
	message("combine_chains called with ", length(chain_paths), " chain paths")
	if (!is.null(cut_size1)) {
		message("  cut_size1 (parameter samples limit): ", cut_size1)
	}
	if (!is.null(cut_size2)) {
		message("  cut_size2 (plot samples limit): ", cut_size2)
	}

	readInputRdsFile = function(input_rds){
		input = tryCatch({
			x <- readRDS(input_rds)
			# If compat shim object with plot_mu[, rename to mu[ for old utils
			if (is.list(x) && !is.null(x$samples2) && is.matrix(x$samples2)) {
				cn <- colnames(x$samples2)
				if (!is.null(cn) && any(grepl("^plot_mu\\[", cn))) {
					colnames(x$samples2) <- sub("^plot_mu\\[", "mu[", cn)
				}
			}
			# Also handle old format where samples2 is second element
			if (is.list(x) && length(x) >= 2 && is.matrix(x[[2]])) {
				cn <- colnames(x[[2]])
				if (!is.null(cn) && any(grepl("^plot_mu\\[", cn))) {
					colnames(x[[2]]) <- sub("^plot_mu\\[", "mu[", cn)
				}
			}
			x
		}, error = function(c) {
			message("The input *rds is invalid")
			return(NA)
		})
	}

	# initialize
	samples <- list()
	metadata <- list()
	for(i in 1:length(chain_paths)){
		print(i)
		# Read chain file
		chain <- readInputRdsFile(chain_paths[[i]])
		if (any(is.na(chain))) {
			message("Chain ", i, " is NA, skipping...")
			next()
		}
		
		message("Chain ", i, " has ", length(chain), " elements")
		if (is.matrix(chain)) {
			message("Chain ", i, " is a matrix with ", nrow(chain), " rows and ", ncol(chain), " columns")
		} else {
			message("Chain ", i, " is not a matrix, class: ", class(chain))
		}
		
		# Handle both old and new chain formats
		if (is.list(chain) && "samples" %in% names(chain)) {
			# New format: list with samples and metadata
			# Convert samples to mcmc object if it's a matrix
			if (is.matrix(chain$samples)) {
				chain_samples <- chain$samples
				# Apply cut_size1 limit immediately after reading to save memory
				if (!is.null(cut_size1) && nrow(chain_samples) > cut_size1) {
					message("Truncating chain ", i, " immediately after reading: ", nrow(chain_samples), " -> ", cut_size1, " rows")
					start_row <- nrow(chain_samples) - cut_size1 + 1
					chain_samples <- chain_samples[start_row:nrow(chain_samples), , drop = FALSE]
				}
				samples[[i]] <- mcmc(chain_samples)
			} else {
				samples[[i]] <- chain$samples
			}
			if ("metadata" %in% names(chain)) {
				metadata[[i]] <- chain$metadata
			}
		} else if (is.list(chain) && length(chain) >= 3) {
			# Old format: list with samples, samples2, metadata
			# Convert samples to mcmc object if it's a matrix
			if (is.matrix(chain[[1]])) {
				chain_samples <- chain[[1]]
				# Apply cut_size1 limit immediately after reading to save memory
				if (!is.null(cut_size1) && nrow(chain_samples) > cut_size1) {
					message("Truncating chain ", i, " immediately after reading: ", nrow(chain_samples), " -> ", cut_size1, " rows")
					start_row <- nrow(chain_samples) - cut_size1 + 1
					chain_samples <- chain_samples[start_row:nrow(chain_samples), , drop = FALSE]
				}
				samples[[i]] <- mcmc(chain_samples)
			} else {
				samples[[i]] <- chain[[1]]
			}
			if (length(chain) >= 3) {
				metadata[[i]] <- chain[[3]]
			}
		} else {
			message("Chain ", i, " has unexpected format, skipping...")
			next()
		}
	}

	# Remove any NULL elements
	samples <- samples[!sapply(samples, is.null)]
	metadata <- metadata[!sapply(metadata, is.null)]
	
	# Debug: Check if we have any valid samples
	message("After filtering, samples list has ", length(samples), " elements")
	if (length(samples) == 0) {
		message("No valid samples found after filtering")
		return(NULL)
	}

	# Now make them all the same size using window_chain (best functionality from combine_chains_simple)
	nrows <- lapply(samples, nrow) %>% unlist()
	message("nrows: ", paste(nrows, collapse = ", "))
	min_nrow <- min(nrows)
	message("min_nrow: ", min_nrow)
	
	# Apply cut_size1 limit if specified (for parameter samples)
	if (!is.null(cut_size1) && cut_size1 < min_nrow) {
		message("Applying cut_size1 limit: truncating all chains to ", cut_size1, " samples")
		min_nrow <- cut_size1
	}
	
	# Ensure we don't truncate chains to unreasonably small sizes
	min_reasonable_size <- 1000  # Minimum reasonable chain size
	if (min_nrow < min_reasonable_size) {
		message("Warning: Chains are already very small (", min_nrow, " rows). Not truncating further.")
		# Even for small chains, we need to make them the same size
		for(i in 1:length(samples)){
			current_nrow <- nrow(samples[[i]])
			if (min_nrow < current_nrow){
				message("Truncating chain ", i, " from ", current_nrow, " to ", min_nrow, " rows")
				samples[[i]] <- window_chain(samples[[i]], max_size = min_nrow)
			}
		}
	} else {
		for(i in 1:length(samples)){
			current_nrow <- nrow(samples[[i]])
			if (min_nrow < current_nrow){
				message("Truncating chain ", i, " from ", current_nrow, " to ", min_nrow, " rows")
				# Use a reasonable minimum size instead of (min_nrow-1)
				target_size <- max(min_nrow, min_reasonable_size)
				samples[[i]] <- window_chain(samples[[i]], max_size = target_size)
			}
		}
	}
	
	# CRITICAL: After truncation, ensure all chains have the same mcpar
	# The window_chain function should set mcpar, but we need to ensure they all match
	# Use the mcpar from the first chain (which should have been set by window_chain)
	if (length(samples) > 0) {
		reference_mcpar <- attr(samples[[1]], "mcpar")
		if (is.null(reference_mcpar)) {
			# If first chain doesn't have mcpar, create one based on its size
			n_rows <- nrow(samples[[1]])
			reference_mcpar <- c(1, n_rows, 1)  # start, end, thin
			attr(samples[[1]], "mcpar") <- reference_mcpar
			message("Created default mcpar for first chain after truncation: ", paste(reference_mcpar, collapse = ", "))
		}
		
		# Set all chains to have the same mcpar
		for (i in 1:length(samples)) {
			attr(samples[[i]], "mcpar") <- reference_mcpar
		}
		message("Set matching mcpar for all chains after truncation: ", paste(reference_mcpar, collapse = ", "))
	}
	
	# Verify all chains are now the same size
	final_nrows <- lapply(samples, nrow) %>% unlist()
	message("Final chain sizes after truncation: ", paste(final_nrows, collapse = ", "))
	if (length(unique(final_nrows)) > 1) {
		message("ERROR: Chains still have different sizes after truncation!")
		stop("Chain size mismatch after truncation")
	}

	# Make the attributes match up (best functionality from combine_chains_simple)
	# CRITICAL: Ensure all chains have matching mcpar attributes before creating mcmc.list
	# Get or create reference mcpar based on actual chain size
	n_rows <- nrow(samples[[1]])
	reference_mcpar <- c(1, n_rows, 1)  # start, end, thin
	
	# Force all chains to have the same mcpar (overwrite any existing values)
	for (i in 1:length(samples)) {
		attr(samples[[i]], "mcpar") <- reference_mcpar
	}
	
	message("Set matching mcpar for all chains: ", paste(reference_mcpar, collapse = ", "))

	# Read metadata from the first chain file (best functionality from combine_chains_simple_new)
	first_chain_path <- chain_paths[[1]]
	first_chain <- readRDS(first_chain_path)
	
	if (!is.list(first_chain) || !"metadata" %in% names(first_chain)) {
		# Try to get metadata from the metadata list we collected
		if (length(metadata) > 0) {
			metadata_final <- metadata[[1]]
			message("Using metadata from collected metadata list")
		} else {
			stop("Chain file does not contain valid metadata. Cannot proceed without real metadata.")
		}
	} else {
		metadata_final <- first_chain$metadata
		message("Found metadata in first chain file, preserving it")
	}
	
	# Verify that all important metadata components are preserved
	metadata_components <- names(metadata_final)
	message("Metadata components found: ", paste(metadata_components, collapse = ", "))
	
	# Check for critical metadata components
	critical_components <- c("model_name", "use_legacy_covariate", "has_driver_uncertainty")
	for (comp in critical_components) {
		if (comp %in% metadata_components) {
			message("✅ ", comp, " preserved: ", metadata_final[[comp]])
		} else {
			message("⚠️  ", comp, " not found in metadata")
		}
	}
	
	# Check for Nimble model code specifically
	if ("nimble_code" %in% metadata_components) {
		message("✅ Nimble model code preserved (length: ", length(metadata_final$nimble_code), ")")
	} else if ("model_code" %in% metadata_components) {
		message("✅ Model code preserved (length: ", length(metadata_final$model_code), ")")
	} else {
		message("⚠️  Model code component not found in metadata")
	}
	
	# Handle samples2 properly - it should contain plot-level predictions, not parameter samples
	samples2_list <- list()
	for(i in 1:length(chain_paths)){
		chain <- readInputRdsFile(chain_paths[[i]])
		if (any(is.na(chain))) {
			next()
		}
		
		# Extract samples2 if it exists and has proper plot prediction structure
		if (is.list(chain) && "samples2" %in% names(chain)) {
			if (is.matrix(chain$samples2)) {
				col_names <- colnames(chain$samples2)
				# CRITICAL FIX: Check for plot prediction columns (plot_mu, mu, Ex, plot_estimates, plot_predictions)
				has_plot_cols <- any(grepl("plot_mu", col_names)) || any(grepl("^Ex\\[", col_names)) || 
				                 any(grepl("^mu\\[", col_names)) || any(grepl("plot_estimates", col_names)) || 
				                 any(grepl("plot_predictions", col_names))
				
				# CRITICAL FIX: Also check that it does NOT contain parameter columns (beta, intercept, etc.)
				# If samples2 contains parameters instead of plot predictions, skip it
				has_param_cols <- any(grepl("^beta\\[", col_names)) || any(grepl("^intercept", col_names)) || 
				                  any(grepl("^site_effect\\[", col_names)) || any(grepl("^precision", col_names)) ||
				                  any(grepl("^rho", col_names)) || any(grepl("^legacy_effect", col_names))
				
				if (has_plot_cols && !has_param_cols) {
					# samples2 has proper plot prediction structure (plot_mu, Ex[i,j], mu[i,j], plot_estimates, or plot_predictions format)
					samples2_list[[i]] <- chain$samples2
					plot_param_type <- if(any(grepl("plot_mu", col_names))) "plot_mu" 
						else if(any(grepl("^Ex\\[", col_names))) "Ex[i,j]"
						else if(any(grepl("^mu\\[", col_names))) "mu[i,j]"
						else if(any(grepl("plot_estimates", col_names))) "plot_estimates"
						else "plot_predictions"
					message("Chain ", i, " has proper samples2 with ", plot_param_type, " parameters")
				} else if (has_param_cols) {
					# samples2 contains parameter columns instead of plot predictions - this is wrong
					message("Chain ", i, " has samples2 with PARAMETER columns (beta, intercept, etc.) instead of plot predictions - skipping samples2")
					message("  Found columns like:", paste(head(col_names[grepl("^beta\\[|^intercept|^site_effect", col_names)], 5), collapse=", "))
					samples2_list[[i]] <- NULL
				} else {
					# samples2 is malformed - don't use it
					message("Chain ", i, " has malformed samples2 (", ncol(chain$samples2), " cols, no plot predictions), skipping samples2")
					samples2_list[[i]] <- NULL
				}
			}
		} else if (is.list(chain) && length(chain) >= 2) {
			# Old format: check if second element is samples2
			if (is.matrix(chain[[2]]) && (any(grepl("plot_mu", colnames(chain[[2]]))) || any(grepl("^Ex\\[", colnames(chain[[2]]))) || any(grepl("^mu\\[", colnames(chain[[2]]))) || any(grepl("plot_estimates", colnames(chain[[2]]))) || any(grepl("plot_predictions", colnames(chain[[2]]))))) {
				samples2_list[[i]] <- chain[[2]]
				plot_param_type <- if(any(grepl("plot_mu", colnames(chain[[2]])))) "plot_mu" 
					else if(any(grepl("^Ex\\[", colnames(chain[[2]])))) "Ex[i,j]"
					else if(any(grepl("^mu\\[", colnames(chain[[2]])))) "mu[i,j]"
					else if(any(grepl("plot_estimates", colnames(chain[[2]])))) "plot_estimates"
					else "plot_predictions"
				message("Chain ", i, " has proper samples2 from old format with ", plot_param_type, " parameters")
			} else {
				message("Chain ", i, " has malformed samples2 in old format, skipping")
				samples2_list[[i]] <- NULL
			}
		} else {
			message("Chain ", i, " has no samples2, skipping")
			samples2_list[[i]] <- NULL
		}
	}
	
	# Remove NULL elements from samples2
	samples2_list <- samples2_list[!sapply(samples2_list, is.null)]
	
	# Ensure all samples2 chains have the same size (match samples chains)
	if (length(samples2_list) > 0) {
		samples2_nrows <- lapply(samples2_list, nrow) %>% unlist()
		min_samples2_nrow <- min(samples2_nrows)
		
		# Apply cut_size2 limit if specified
		if (!is.null(cut_size2) && cut_size2 < min_samples2_nrow) {
			message("Applying cut_size2 limit: truncating all samples2 chains to ", cut_size2, " samples")
			min_samples2_nrow <- cut_size2
		}
		
		# Truncate all samples2 chains to the same size
		for(i in 1:length(samples2_list)){
			current_nrow <- nrow(samples2_list[[i]])
			if (min_samples2_nrow < current_nrow){
				message("Truncating samples2 chain ", i, " from ", current_nrow, " to ", min_samples2_nrow, " rows")
				if (is.matrix(samples2_list[[i]])) {
					samples2_list[[i]] <- samples2_list[[i]][(current_nrow - min_samples2_nrow + 1):current_nrow, , drop = FALSE]
				}
			}
		}
		
		# Verify all samples2 chains are now the same size
		final_samples2_nrows <- lapply(samples2_list, nrow) %>% unlist()
		if (length(unique(final_samples2_nrows)) > 1) {
			message("WARNING: samples2 chains still have different sizes after truncation: ", paste(final_samples2_nrows, collapse = ", "))
		}
		
		# CRITICAL: Ensure all samples2 chains have the same number of rows as samples chains
		# This is important because samples2 mcpar should match samples mcpar
		if (length(samples) > 0) {
			samples_nrows <- nrow(samples[[1]])
			if (min_samples2_nrow != samples_nrows) {
				message("WARNING: samples2 chains (", min_samples2_nrow, " rows) don't match samples chains (", samples_nrows, " rows)")
				message("Truncating samples2 chains to match samples chains: ", samples_nrows, " rows")
				min_samples2_nrow <- samples_nrows
				for(i in 1:length(samples2_list)){
					current_nrow <- nrow(samples2_list[[i]])
					if (min_samples2_nrow < current_nrow){
						if (is.matrix(samples2_list[[i]])) {
							samples2_list[[i]] <- samples2_list[[i]][(current_nrow - min_samples2_nrow + 1):current_nrow, , drop = FALSE]
						}
					}
				}
			}
		}
	}
	
	# Create output with proper samples and samples2
	# CRITICAL: Create samples mcmc.list first to get reference mcpar
	# Debug: Check mcpar attributes before creating mcmc.list
	message("DEBUG: About to create samples mcmc.list...")
	message("DEBUG: Checking mcpar attributes before creating samples mcmc.list:")
	for (i in 1:length(samples)) {
		current_mcpar <- attr(samples[[i]], "mcpar")
		message("DEBUG:   Chain ", i, " mcpar: ", if (is.null(current_mcpar)) "NULL" else paste(current_mcpar, collapse = ", "))
	}
	
	message("DEBUG: Attempting to create samples mcmc.list...")
	samples_mcmc_list <- tryCatch({
		result <- as.mcmc.list(samples)
		message("DEBUG: Successfully created samples mcmc.list")
		result
	}, error = function(e) {
		message("ERROR creating samples mcmc.list: ", e$message)
		message("DEBUG: Chain mcpar values at error:")
		for (i in 1:length(samples)) {
			current_mcpar <- attr(samples[[i]], "mcpar")
			message("DEBUG:   Chain ", i, ": ", if (is.null(current_mcpar)) "NULL" else paste(current_mcpar, collapse = ", "))
		}
		stop("Failed to create samples mcmc.list: ", e$message)
	})
	
	if (length(samples2_list) > 0) {
		# We have valid samples2 data - convert matrices to mcmc objects
		message("Using ", length(samples2_list), " valid samples2 chains")
		
		# Get reference mcpar from samples (already matched)
		reference_mcpar <- attr(samples[[1]], "mcpar")
		
		# Create samples2 mcmc objects with matching mcpar attributes from the start
		samples2_mcmc <- lapply(samples2_list, function(x) {
			if (is.matrix(x)) {
				mcmc_obj <- mcmc(x)
				# Set mcpar to match samples chains
				attr(mcmc_obj, "mcpar") <- reference_mcpar
				return(mcmc_obj)
			} else if (inherits(x, "mcmc")) {
				# If already an mcmc object, just update mcpar
				attr(x, "mcpar") <- reference_mcpar
				return(x)
			} else {
				return(x)
			}
		})
		
		# Verify all samples2 mcmc objects have matching mcpar
		for (i in 1:length(samples2_mcmc)) {
			if (inherits(samples2_mcmc[[i]], "mcmc")) {
				current_mcpar <- attr(samples2_mcmc[[i]], "mcpar")
				if (!identical(current_mcpar, reference_mcpar)) {
					message("WARNING: samples2 chain ", i, " mcpar doesn't match, fixing...")
					attr(samples2_mcmc[[i]], "mcpar") <- reference_mcpar
				}
			}
		}
		message("Set matching mcpar attributes for all samples2 chains")
		
		# Debug: Check samples2 mcpar before creating mcmc.list
		message("DEBUG: Checking samples2 mcpar before creating mcmc.list:")
		for (i in 1:length(samples2_mcmc)) {
			if (inherits(samples2_mcmc[[i]], "mcmc")) {
				current_mcpar <- attr(samples2_mcmc[[i]], "mcpar")
				message("DEBUG:   samples2 chain ", i, " mcpar: ", if (is.null(current_mcpar)) "NULL" else paste(current_mcpar, collapse = ", "))
			}
		}
		
		# Create samples2 mcmc.list
		message("DEBUG: Attempting to create samples2 mcmc.list...")
		samples2_mcmc_list <- tryCatch({
			result <- as.mcmc.list(samples2_mcmc)
			message("DEBUG: Successfully created samples2 mcmc.list")
			result
		}, error = function(e) {
			message("ERROR creating samples2 mcmc.list: ", e$message)
			message("DEBUG: samples2 chain mcpar values at error:")
			for (i in 1:length(samples2_mcmc)) {
				if (inherits(samples2_mcmc[[i]], "mcmc")) {
					current_mcpar <- attr(samples2_mcmc[[i]], "mcpar")
					message("DEBUG:   Chain ", i, ": ", if (is.null(current_mcpar)) "NULL" else paste(current_mcpar, collapse = ", "))
				} else {
					message("DEBUG:   Chain ", i, ": not an mcmc object, class = ", class(samples2_mcmc[[i]]))
				}
			}
			stop("Failed to create samples2 mcmc.list: ", e$message)
		})
		
		out <- list(samples = samples_mcmc_list,
								samples2 = samples2_mcmc_list,
								metadata = metadata_final)
	} else {
		# No valid samples2 data - create empty structure
		message("No valid samples2 data found - creating empty samples2 structure")
		out <- list(samples = samples_mcmc_list,
								samples2 = list(),  # Empty list instead of wrong data
								metadata = metadata_final)
	}

	if(!isFALSE(save)){
		saveRDS(out, file = save)
	}
	return(out)
}

# Remove the old functions since we now have the hybrid version
# combine_chains_simple and combine_chains_simple_new are deprecated


#' @title 			window_chain
#' @description Shorten chain to a specific size, taking the last samples
#' @export
window_chain = function(chain, thin = 1, max_size = 10000) {
	require(coda)
	if (all(class(chain) != "mcmc")){
		chain <- mcmc(as.matrix(chain))
	}
	nrow_samples <- nrow(chain)
	
	# Always truncate to max_size, taking the last samples
	if (nrow_samples <= max_size) {
		message("Chain already at or below target size (", nrow_samples, " <= ", max_size, ")")
		out_chain <- chain
		attr(out_chain, "mcpar")[[1]] <- 1
		attr(out_chain, "mcpar")[[2]] <- nrow_samples
		return(out_chain)
	} else {
		message("Truncating chain from ", nrow_samples, " to ", max_size, " samples")
		# Take the last max_size samples
		start_row <- nrow_samples - max_size + 1
		out_chain <- chain[start_row:nrow_samples, , drop = FALSE]
		# Convert back to mcmc object
		out_chain <- mcmc(out_chain)
		attr(out_chain, "mcpar")[[1]] <- 1
		attr(out_chain, "mcpar")[[2]] <- max_size
		return(out_chain)
	}
}


