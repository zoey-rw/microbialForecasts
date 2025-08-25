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
	critical_components <- c("model_name", "use_legacy_covariate")
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
	
	# Create output with both samples and samples2 (compatibility with existing code)
	out <- list(samples = as.mcmc.list(samples),
							samples2 = as.mcmc.list(samples),  # Use same samples for samples2
							metadata = metadata_final)

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

