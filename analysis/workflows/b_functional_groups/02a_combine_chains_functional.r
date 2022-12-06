# Combine chains for functional group models
pacman::p_load(nimble, coda, tidyverse)
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

model_name <- "cycl_only"
model_name <- "all_covariates"
for (model_name in c("cycl_only")){
#for (model_name in c("all_covariates")){

	#for (model_name in c("all_covariates", "cycl_only")){
	print(model_name)

	file.list <- list.files(path = paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/functional_groups//", model_name),
													pattern = "_chain",
													full.names = T)

	info <- file.info(file.list)
	newer <- rownames(info[which(info$mtime > "2022-11-07 00:00:00 EDT"),])
	file.list <- file.list[file.list %in% newer]
	rank <- "denitrification"
	#time_period <- "calibration"
	time_period <- "20151101_20200101"

	for (rank in microbialForecast:::fg_names[1:66]){
		print(rank)


			chain_paths <- file.list[grepl(paste0(rank, "_20"), file.list, fixed=T)]
			if(length(chain_paths)==0) next()
			#for (time_period in c("refit")){
			print(time_period)

			chain_paths_time_period <- chain_paths[grepl(time_period, chain_paths)]
			if(length(chain_paths_time_period)==0) next()

			savepath <- gsub("_chain[1234]","_summary",chain_paths_time_period[[1]])
			print(savepath)
			#if (!file.exists(savepath)){
				# TO DO: remove the problematic chain here
				#if (rank == "phylum_fun" & model_name == "all_covariates" & time_period == "calibration") chain_paths_time_period <- chain_paths_time_period[c(1,3:4)]#
				#if (rank == "class_bac" & model_name == "cycl_only") chain_paths_time_period <- chain_paths_time_period[2:4]#
				#if (rank == "family_bac" & model_name == "cycl_only") chain_paths_time_period <- chain_paths_time_period[2:4]#

				print(chain_paths_time_period)


#				if(length(chain_paths_time_period)==1) next()

				try_read_chains = function(chain_paths_time_period){
					out = tryCatch(combine_chains(chain_paths_time_period),
								 error = function(c) {
								 	message("Couldn't read files")
								 	return(NA)
								 })}
				out = try_read_chains(chain_paths_time_period)
				if (is.na(out)) next()

				es <- effectiveSize(out$samples)
				if (min(es) < 10) {
					#plot(out$samples[,1:10])
					#plot(out$samples[,21:30])
					print(es[1:50])
				}

				# Calculate summary and save output.
				param_summary <- fast.summary.mcmc(out$samples)
				plot_summary <- fast.summary.mcmc(out$samples2)

				out_summary <- list(samples = out$samples,
														param_summary = param_summary,
														plot_summary = plot_summary,
														metadata = out$metadata)

				if (length(out$samples) > 1){
					gelman_out <- cbind(gelman.diag(out$samples)[[1]], es)
					out_summary$gelman = gelman_out
				}

				saveRDS(out_summary, savepath, compress = F)
				cat(paste0("\nSaved summary output for ", model_name, ", ", time_period, ", ",
									 rank, " to: \n", savepath, "\n"))


				#If the summary now exists, delete the chains
				if (file.exists(savepath)){
					cat(paste0("\nDeleting chains for ", model_name, ", ", time_period, ", ",
										 rank, " : \n", chain_paths_time_period, "\n"))
					unlink(chain_paths_time_period)
				}


			#} else print("Summary file already exists")
		}
	}
}
