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
		if (any(is.na(chain))) next()
		nrow_samples <- nrow(chain[[1]])
		nrow_samples2 <- nrow(chain[[2]])

		samples[[i]] <- window_chain(chain[[1]])
		samples2[[i]] <- window_chain(chain[[2]])
		metadata[[i]] <- chain[[3]]
	}

	samples <- samples[!sapply(samples,FUN = function(x) is.null(colnames(x)))]
	samples2 <- samples2[!sapply(samples2,FUN = function(x) is.null(colnames(x)))]
	metadata<-metadata[!sapply(metadata,is.null)]

	# Now make them all the same size
	nrows <- lapply(samples, nrow) %>% unlist()
	min_nrow <- min(nrows)
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

