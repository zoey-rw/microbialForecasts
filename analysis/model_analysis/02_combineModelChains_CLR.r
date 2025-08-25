# Combine chains from each taxon model, and create basic summary stats
source("../../source.R")

# For phenology analysis: use cycl_only models for 2013-2015 data
model_name <- "cycl_only"

# For testing: use our actual CLR model IDs instead of reading from CSV
# params_in = read.csv(here("data/clean/model_input_df.csv"))
# model_id_list = unique(params_in$model_id)

# Use our expanded CLR model IDs from CLR_regression directory
# Get all model IDs from the CLR_regression directory by looking at individual chain files
clr_output_dir <- here("data/model_outputs/CLR_regression")
clr_chain_files <- list.files(clr_output_dir, pattern = "samples_.*_chain[0-9]+\\.rds$", 
                              recursive = TRUE, full.names = FALSE)
# Extract model IDs from chain file names (remove _chain[0-9] suffix)
model_id_list <- unique(gsub("samples_(.*)_chain[0-9]+\\.rds", "\\1", clr_chain_files))
model_id_list <- gsub(".*/", "", model_id_list)  # Remove directory path

cat("Found", length(model_id_list), "CLR models in CLR_regression\n")

cat("Processing CLR phenology models for 2013-2015 data\n")
cat("Model type:", model_name, "\n")
cat("Total models to process:", length(model_id_list), "\n\n")

# For testing: use just 1 core instead of 28
# cl <- makeCluster(28, outfile="")
# registerDoParallel(cl)
cl <- makeCluster(1, outfile="")
registerDoParallel(cl)
#for(model_id in model_id_list){
foreach(model_id=model_id_list, .errorhandling = "pass") %dopar% {
	source("../../source.R")
	
	# Do we want to keep all the chain files separately? Deleting them will save space
	delete_samples_files = F
	
	#for (model_name in c("all_covariates", "cycl_only")) {
	file.list <- list.files(path = here("data/model_outputs/CLR_regression/"),
													pattern = "_chain",
													recursive = T,
													full.names = T)
	
	# Subset to newest output files
	info <- file.info(file.list)
	
	# Subset more than 1KB (for our test files, lower threshold)
	# large_enough <- rownames(info[which(info$size > 100000000), ])
	# large_enough <- rownames(info[which(info$size > 5000000), ])
	large_enough <- rownames(info[which(info$size > 1000), ])
	
	#newer <- rownames(info[which(info$mtime > "2023-03-09 06:00:00 EDT"), ])
	#newer <- rownames(info[which(info$mtime > "2023-06-09 18:00:00 EDT"), ])
	# For testing: accept all recent files
	newer <- rownames(info[which(info$mtime > "2020-01-01 00:00:00 EDT"), ])
	
	# Don't want files still being written - at least 1 min old (for testing)
	# older <- rownames(info[which(info$mtime < (Sys.time()-3600)), ])
	older <- rownames(info[which(info$mtime < (Sys.time()-60)), ])
	
	# If deleting sample files, we don't want models still being run - at least 18 hours old
	# if (delete_samples_files) {
	# 	older <- rownames(info[which(info$mtime < (Sys.time()-64800)), ])
	# }
	file.list <- file.list[file.list %in% newer & file.list %in% large_enough]
	
	message("Searching model outputs for ",model_id)
	
	
	# Subset to files of interest
	chain_paths <-
		file.list[grepl(model_id, file.list, fixed = T)]
	
	if (length(chain_paths) == 0)
		return(message("Skipping ", model_id, "; no chains"))
	if (length(chain_paths) == 1) {
		cat("Single chain model, copying to combined samples file...\n")
		# For single chain models, just copy the chain file to the samples file
		single_chain_path <- chain_paths[[1]]
		combined_samples_path <- gsub("_chain[1234567]", "", single_chain_path)
		
		if (file.exists(combined_samples_path)) {
			cat("Combined samples file already exists, skipping...\n")
			return(message("Skipping ", model_id, "; combined samples already exist"))
		}
		
		cat("Copying single chain to combined samples file...\n")
		chain_data <- readRDS(single_chain_path)
		
		# Create the same output structure as multi-chain models
		if (is.list(chain_data) && "samples" %in% names(chain_data)) {
			# New format with metadata
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
			# Old format (raw matrix)
			out <- list(
				samples = list(mcmc(chain_data, start = 1, end = nrow(chain_data), thin = 1)),
				samples2 = list(mcmc(chain_data, start = 1, end = nrow(chain_data), thin = 1)),
				metadata = list(model_id = model_id)
			)
		}
		
		# Create summaries for single chain
		cat("Creating parameter summaries for single chain...\n")
		param_summary <- summary(out$samples[[1]])
		cat("Creating plot summaries for single chain...\n")
		plot_summary <- summary(out$samples2[[1]])
		cat("Calculating effective sample sizes...\n")
		es <- effectiveSize(out$samples[[1]])
		
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
		
		
		# For CLR outputs: handle the proper structure with samples and metadata
		cat("Loading", length(chain_paths), "chains...\n")
		chains <- list()
		chains2 <- list()
		metadata_list <- list()
		
		for (i in seq_along(chain_paths)) {
			cat("  Loading chain", i, ":", basename(chain_paths[i]), "\n")
			chain_data <- readRDS(chain_paths[i])
			
			if (is.list(chain_data) && "samples" %in% names(chain_data)) {
				# New format with metadata
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
				# Old format (raw matrix)
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
		
		cat("Processing chains for outlier detection...\n")
		# Remove chains if one has very different values than the other 3
		chains <- out$samples
		cat("Calculating chain means for outlier detection...\n")
		means <- lapply(chains, function(x) mean(x[,"intercept"], na.rm=T))
		cat("Chain means:", paste(round(unlist(means), 4), collapse = ", "), "\n")
		scaled_means = scale(unlist(means))
		cat("Scaled means:", paste(round(scaled_means, 4), collapse = ", "), "\n")
		potential_outlier <- which(abs(scaled_means) > 1.3)
		cat("Potential outliers:", potential_outlier, "\n")
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
		
		
		
		cat("Creating parameter summaries...\n")
		param_summary <- fast.summary.mcmc(out$samples)
		cat("Creating plot summaries...\n")
		plot_summary <- fast.summary.mcmc(out$samples2)
		cat("Calculating effective sample sizes...\n")
		es <- effectiveSize(out$samples)
		if (length(out$samples) > 1) {
			cat("Calculating Gelman diagnostics...\n")
			gelman_out <- cbind(gelman.diag(out$samples, multivariate = F)[[1]], es)
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
		
		
		if (min(es) < 500) {
			message("Low effective sample sizes - check for unconverged parameters in model: ", model_id)
			print(head(gelman_out, 50))
		}
		
		
		# If the summary now exists, delete the chains
		if (delete_samples_files){
			if (samples_exists(chain_paths_time_period[[1]])) {
				unlink(chain_paths_time_period)
				message("Deleting samples files, e.g.: ", chain_paths_time_period[[1]])
			}
		} #else message("Not deleting samples files, e.g.: ", chain_paths_time_period[[1]])
	}
	return()
} # End model_id loop


 stopCluster(cl)



# Removed hardcoded paths and test code for local testing 


