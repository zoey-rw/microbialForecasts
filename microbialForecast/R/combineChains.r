#' @title 			combine_chains_simple
#' @description Combine MCMC chains using paths, shortening each chain due to RAM constraints
#' @export
combine_chains_simple <- function(chain_paths,
																	save = FALSE,
																	cut_size1 = NULL,
																	cut_size2 = NULL){
	require(coda)
	require(tidyverse)

	if (is.null(cut_size1)) cut_size1 <- 19999
	if (is.null(cut_size2)) cut_size2 <- 9999
	
	message("combine_chains_simple called with ", length(chain_paths), " chain paths")

	readInputRdsFile = function(input_rds){
		input = tryCatch(readRDS(input_rds),
										 error = function(c) {
										 	message("The input *rds is invalid")
										 	return(NA)
										 }
		)
	}


	# initialize
	samples <- samples2 <- metadata <- list()
	first_iter <- last_iter <- list()
	for(i in 1:length(chain_paths)){


		print(i)
		# paste model file path to chain number
		chain <- readInputRdsFile(chain_paths[[i]])
		if (any(is.na(chain))) {
			message("Chain ", i, " is NA, skipping...")
			next()
		}
		
		message("Chain ", i, " has ", length(chain), " elements")
		if (length(chain) < 3) {
			message("Chain ", i, " has insufficient elements, skipping...")
			next()
		}
		
		nrow_samples <- nrow(chain[[1]])
		nrow_samples2 <- nrow(chain[[2]])
		message("Chain ", i, " samples: ", nrow_samples, " rows, samples2: ", nrow_samples2, " rows")

		samples[[i]] <- window_chain(chain[[1]])
		samples2[[i]] <- window_chain(chain[[2]])
		metadata[[i]] <- chain[[3]]
	}

	samples <- samples[!sapply(samples,FUN = function(x) is.null(colnames(x)))]
	samples2 <- samples2[!sapply(samples2,FUN = function(x) is.null(colnames(x)))]
	metadata<-metadata[!sapply(metadata,is.null)]
	
	# Debug: Check if we have any valid samples
	message("After filtering, samples list has ", length(samples), " elements")
	if (length(samples) == 0) {
		message("No valid samples found after filtering")
		return(NULL)
	}

	# Now make them all the same size
	nrows <- lapply(samples, nrow) %>% unlist()
	message("nrows: ", paste(nrows, collapse = ", "))
	min_nrow <- min(nrows)
	message("min_nrow: ", min_nrow)
	for(i in 1:length(samples)){
		current_nrow <- nrow(samples[[i]])
		if (min_nrow < current_nrow){
			print(i)
			samples[[i]] <- window_chain(samples[[i]], max_size = (min_nrow-1))

		}
	}

	# Now make them all the same size, chain 2
	nrows <- lapply(samples2, nrow) %>% unlist()
	min_nrow <- min(nrows)
	for(i in 1:length(samples2)){
		current_nrow <- nrow(samples2[[i]])
		if (min_nrow < current_nrow){
			print(i)
			samples2[[i]] <- window_chain(samples2[[i]], max_size = (min_nrow-1))
		}
	}

	# Make the attributes match up (sort of arbitrary)
	for (i in 1:length(samples)) {
		attr(samples[[i]], "mcpar") = attr(samples[[1]], "mcpar")
	}

	for (i in 1:length(samples2)) {
		attr(samples2[[i]], "mcpar") = attr(samples2[[1]], "mcpar")
	}

	out <- list(samples = as.mcmc.list(samples),
							samples2 = as.mcmc.list(samples2),
							metadata = metadata[[1]])

	if(!isFALSE(save)){
		saveRDS(out, file = save)
	}
	return(out)
}


#' @title 			combine_chains_simple_new
#' @description Combine MCMC chains that are saved as matrices (not lists)
#' @export
combine_chains_simple_new <- function(chain_paths,
																	save = FALSE,
																	cut_size1 = NULL,
																	cut_size2 = NULL){
	require(coda)
	require(tidyverse)

	if (is.null(cut_size1)) cut_size1 <- 19999
	if (is.null(cut_size2)) cut_size2 <- 9999
	
	message("combine_chains_simple_new called with ", length(chain_paths), " chain paths")

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
		
		# Convert to mcmc object and window if needed
		# chain should now be a list with samples and metadata
		if (!is.list(chain) || !"samples" %in% names(chain)) {
			stop("Chain file must be a list with 'samples' and 'metadata' components")
		}
		samples[[i]] <- window_chain(chain$samples)
	}

	# Remove any NULL elements
	samples <- samples[!sapply(samples, is.null)]
	
	# Debug: Check if we have any valid samples
	message("After filtering, samples list has ", length(samples), " elements")
	if (length(samples) == 0) {
		message("No valid samples found after filtering")
		return(NULL)
	}

	# Now make them all the same size
	nrows <- lapply(samples, nrow) %>% unlist()
	message("nrows: ", paste(nrows, collapse = ", "))
	min_nrow <- min(nrows)
	message("min_nrow: ", min_nrow)
	
	# Ensure we don't truncate chains to unreasonably small sizes
	min_reasonable_size <- 1000  # Minimum reasonable chain size
	if (min_nrow < min_reasonable_size) {
		message("Warning: Chains are already very small (", min_nrow, " rows). Not truncating further.")
	} else {
		for(i in 1:length(samples)){
			current_nrow <- nrow(samples[[i]])
			if (min_nrow < current_nrow){
				print(i)
				# Use a reasonable minimum size instead of (min_nrow-1)
				target_size <- max(min_nrow, min_reasonable_size)
				samples[[i]] <- window_chain(samples[[i]], max_size = target_size)
			}
		}
	}

	# Make the attributes match up (sort of arbitrary)
	for (i in 1:length(samples)) {
		attr(samples[[i]], "mcpar") = attr(samples[[1]], "mcpar")
	}

	# Read metadata from the first chain file - this is the only valid approach
	# The function requires real metadata to work properly
	first_chain_path <- chain_paths[[1]]
	first_chain <- readRDS(first_chain_path)
	
	if (!is.list(first_chain) || !"metadata" %in% names(first_chain)) {
		stop("Chain file does not contain valid metadata. Cannot proceed without real metadata.")
	}
	
	message("Found metadata in first chain file, preserving it")
	metadata <- first_chain$metadata
	
	# Verify that all important metadata components are preserved
	metadata_components <- names(metadata)
	message("Metadata components found: ", paste(metadata_components, collapse = ", "))
	
	# Check for critical metadata components
	critical_components <- c("model_name", "use_legacy_covariate")
	for (comp in critical_components) {
		if (comp %in% metadata_components) {
			message("✅ ", comp, " preserved: ", metadata[[comp]])
		} else {
			message("⚠️  ", comp, " not found in metadata")
		}
	}
	
	# Check for Nimble model code specifically
	if ("nimble_code" %in% metadata_components) {
		message("✅ Nimble model code preserved (length: ", length(metadata$nimble_code), ")")
	} else if ("model_code" %in% metadata_components) {
		message("✅ Model code preserved (length: ", length(metadata$model_code), ")")
	} else {
		message("⚠️  Model code component not found in metadata")
	}
	
	out <- list(samples = as.mcmc.list(samples),
							samples2 = as.mcmc.list(samples),  # Use same samples for samples2
							metadata = metadata)

	if(!isFALSE(save)){
		saveRDS(out, file = save)
	}
	return(out)
}


#' @title 			window_chain
#' @description Shorten chain if it contains more than a certain number of samples
#' @export
window_chain = function(chain, thin = 1, max_size = 10000) {
	require(coda)
	if (all(class(chain) != "mcmc")){
		chain <- mcmc(as.matrix(chain))
	}
	nrow_samples <- nrow(chain)
	if (nrow_samples < max_size) {
		message("Sample size not reduced; current size is fewer than ", max_size)
		out_chain <- chain
		attr(out_chain, "mcpar")[[1]] <- 1
		attr(out_chain, "mcpar")[[2]] <- nrow_samples

		return(out_chain)
	} else {

		first_iter <- attr(chain, "mcpar")[[1]]
		last_iter <- attr(chain, "mcpar")[[2]]

		window_start = last_iter - max_size
		out_chain <- window(chain, window_start, last_iter, 1)
		attr(out_chain, "mcpar")[[1]] <- 1
		attr(out_chain, "mcpar")[[2]] <- max_size
		return(out_chain)
	}
}

