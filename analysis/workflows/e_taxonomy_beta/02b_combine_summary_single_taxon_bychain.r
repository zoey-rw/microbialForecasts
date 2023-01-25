# Combine chains from each taxon model, and create basic summary stats
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

# Do we want to keep all the chain files separately? Deleting them will save space
delete_samples_files = F

# For interactive use/testing.
model_name <- "cycl_only"
model_name <- "all_covariates"
rank <- "order_bac"
rank <- "order_fun"
rank <- "phylum_bac"
rank <- "family_fun"
rank <- "family_bac"
rank <- "class_fun"
rank <- "genus_fun"
rank <- "genus_bac"
rank <- "class_bac"
time_period <- "20151101_20200101"
time_period <- "20151101_20180101"
spec = microbialForecast:::rank_spec_names[[rank]][[1]]


for (model_name in c("all_covariates", "cycl_only")) {
	file.list <- list.files(path = paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/single_taxon/", model_name),
			pattern = "_chain",
			full.names = T)

	# Subset to newer output files
	info <- file.info(file.list)
	newer <- rownames(info[which(info$mtime > "2022-08-02 00:00:00 EDT"), ])
	newer <- rownames(info[which(info$mtime > "2022-12-02 00:00:00 EDT"), ])

	# Don't want files still being written - at least 2 min old
	#older <- rownames(info[which(info$mtime < (Sys.time()-120)), ])
	older <- rownames(info[which(info$mtime < (Sys.time()-3000)), ])

	# If deleting sample files, we don't want models still being run - at least 18 hours old
	if (delete_samples_files) {
	older <- rownames(info[which(info$mtime < (Sys.time()-64800)), ])
	}

	file.list <- file.list[file.list %in% newer & file.list %in% older]

	for (rank in microbialForecast:::tax_names[6]) {
		for (spec in microbialForecast:::rank_spec_names[[rank]]) {
			for (time_period in c("20151101_20180101")) {
				#for (time_period in c("20151101_20180101","20151101_20200101")){
				print(model_name)
				print(rank)
				print(spec)
				print(time_period)

				spec_rank = paste0(rank, "_", spec)
				spec_rank_time = paste0(rank, "_", spec, "_", time_period)
				chain_paths <-
					file.list[grepl(paste0(spec_rank, "_20"), file.list, fixed = T)]
				chain_paths_time_period <-
					chain_paths[grepl(time_period, chain_paths)]


				savepath <- gsub("_chain[1234567]", "", chain_paths_time_period[[1]])

			if (any(grepl("chain[567]", chain_paths_time_period))) {

				chain_paths_time_period <- chain_paths_time_period[grepl("chain[567]", chain_paths_time_period)]

			}


				if (length(chain_paths_time_period) == 0) {
					next()
				}

				if (length(chain_paths_time_period) == 1) {
					message("Skipping ", chain_paths_time_period, "; only one chain")
					next()
				}

				# if (samples_exists(chain_paths_time_period[[1]])) {
				# 	message("Summary file already exists")
				#
				#
				# 	if (delete_samples_files) {
				# 		unlink(chain_paths_time_period)
				# 		message("Deleting samples file: ", chain_paths_time_period)
				# 	}
				# } else {
				# 	print(chain_paths_time_period)
				#
				# 	# Read in and combine each chain
				#
				#
				# 	# If new chains are 5-7, combine with previous
				# 	if (grepl("chain[567]", chain_paths_time_period[[1]])) {
				# 		existing_samples = readRDS(savepath)
				# 		out = combine_chains_existing(input_list = c(chain_paths_time_period, existing_samples$samples))
				# 		# Calculate summary on each output subset, using custom summary function
				# 		param_summary <- fast.summary.mcmc(out)
				# 		plot_summary <- existing_samples$plot_summary
				# 		es <- effectiveSize(out)
				# 		if (length(out) > 1) {
				# 		gelman_out <- cbind(gelman.diag(out)[[1]], es)
				# 		} else gelman_out = NULL
				#
				# 			metadata = existing_samples$metadata
				#
				# 		# Combine and save output
				# 		out_summary <- list(
				# 			samples = out,
				# 			param_summary = param_summary,
				# 			plot_summary = plot_summary,
				# 			metadata = metadata,
				# 			gelman = gelman_out
				# 		)
				#
				# 	} else {

						out <- combine_chains_simple(chain_paths_time_period)
						# Calculate summary on each output subset, using custom summary function
						param_summary <- fast.summary.mcmc(out$samples)
						plot_summary <- fast.summary.mcmc(out$samples2)
						es <- effectiveSize(out$samples)
						if (length(out$samples) > 1) {
							gelman_out <- cbind(gelman.diag(out$samples)[[1]], es)
						} else gelman_out = NULL

					# Combine and save output
					out_summary <- list(
						samples = out$samples,
						param_summary = param_summary,
						plot_summary = plot_summary,
						metadata = out$metadata,
						gelman = gelman_out
					)
					}

					saveRDS(out_summary, savepath, compress = F)
					message("Saved samples and summary output for ",model_name,", ",time_period,", ",rank," to: ",savepath)


if (min(es) < 5) {
	print(es[1:50])
}



# If the summary now exists, delete the chains
if (delete_samples_files){
	if (samples_exists(chain_paths_time_period[[1]])) {
		unlink(chain_paths_time_period)
		message("Deleting samples files, e.g.: ", chain_paths_time_period[[1]])
	}
} else message("Not deleting samples files, e.g.: ", chain_paths_time_period[[1]])

				#}
		} # End time-period loop
	} # End species loop
} # End rank loop
	} # End model type loop
