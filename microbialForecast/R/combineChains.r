#' @title 			combine_chains
#' @description Combine MCMC chains with robust handling of different sample sizes and metadata preservation
#' @export
combine_chains <- function(chain_paths,
																	save = FALSE,
																	cut_size1 = NULL,
																	cut_size2 = NULL){
	require(coda)
	require(tidyverse)

	if (is.null(cut_size1)) cut_size1 <- 19999
	if (is.null(cut_size2)) cut_size2 <- 9999
	
	message("combine_chains called with ", length(chain_paths), " chain paths")

	readInputRdsFile = function(input_rds){
		input = tryCatch(readRDS(input_rds),
										 error = function(c) {
										 	message("The input *rds is invalid")
										 	return(NA)
										 }
		)
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
				samples[[i]] <- mcmc(chain$samples)
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
				samples[[i]] <- mcmc(chain[[1]])
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
	
	# Verify all chains are now the same size
	final_nrows <- lapply(samples, nrow) %>% unlist()
	message("Final chain sizes after truncation: ", paste(final_nrows, collapse = ", "))
	if (length(unique(final_nrows)) > 1) {
		message("ERROR: Chains still have different sizes after truncation!")
		stop("Chain size mismatch after truncation")
	}

	# Make the attributes match up (best functionality from combine_chains_simple)
	for (i in 1:length(samples)) {
		attr(samples[[i]], "mcpar") = attr(samples[[1]], "mcpar")
	}

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
			if (is.matrix(chain$samples2) && (any(grepl("plot_mu", colnames(chain$samples2))) || any(grepl("^Ex\\[", colnames(chain$samples2))) || any(grepl("plot_estimates", colnames(chain$samples2))) || any(grepl("plot_predictions", colnames(chain$samples2))))) {
				# samples2 has proper plot prediction structure (plot_mu, Ex[i,j], plot_estimates, or plot_predictions format)
				samples2_list[[i]] <- chain$samples2
				plot_param_type <- if(any(grepl("plot_mu", colnames(chain$samples2)))) "plot_mu" 
					else if(any(grepl("^Ex\\[", colnames(chain$samples2)))) "Ex[i,j]"
					else if(any(grepl("plot_estimates", colnames(chain$samples2)))) "plot_estimates"
					else "plot_predictions"
				message("Chain ", i, " has proper samples2 with ", plot_param_type, " parameters")
			} else {
				# samples2 is malformed - don't use it
				message("Chain ", i, " has malformed samples2 (", ncol(chain$samples2), " cols, no plot predictions), skipping samples2")
				samples2_list[[i]] <- NULL
			}
		} else if (is.list(chain) && length(chain) >= 2) {
			# Old format: check if second element is samples2
			if (is.matrix(chain[[2]]) && (any(grepl("plot_mu", colnames(chain[[2]]))) || any(grepl("^Ex\\[", colnames(chain[[2]]))) || any(grepl("plot_estimates", colnames(chain[[2]]))) || any(grepl("plot_predictions", colnames(chain[[2]]))))) {
				samples2_list[[i]] <- chain[[2]]
				plot_param_type <- if(any(grepl("plot_mu", colnames(chain[[2]])))) "plot_mu" 
					else if(any(grepl("^Ex\\[", colnames(chain[[2]])))) "Ex[i,j]"
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
	
	# Create output with proper samples and samples2
	if (length(samples2_list) > 0) {
		# We have valid samples2 data - convert matrices to mcmc objects
		message("Using ", length(samples2_list), " valid samples2 chains")
		samples2_mcmc <- lapply(samples2_list, function(x) {
			if (is.matrix(x)) {
				return(mcmc(x))
			} else {
				return(x)
			}
		})
		out <- list(samples = as.mcmc.list(samples),
								samples2 = as.mcmc.list(samples2_mcmc),  # Use actual plot predictions as mcmc.list
								metadata = metadata_final)
	} else {
		# No valid samples2 data - create empty structure
		message("No valid samples2 data found - creating empty samples2 structure")
		out <- list(samples = as.mcmc.list(samples),
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

