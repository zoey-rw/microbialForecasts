# Combine chains
pacman::p_load(nimble, coda, tidyverse) 


source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

model_name <- "cycl_only"

for (model_name in c("all_covariates", "cycl_only")){
print(model_name)
	
	file.list <- list.files(path = paste0("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/taxa/", model_name),
													pattern = "chain",
													full.names = T)

	rank <- "genus_fun"
	rank <- "phylum_bac"
	time_period <- "calibration"
	
	for (rank in tax_names[3:10]){
		print(rank)
		chain_paths <- file.list[grepl(rank, file.list)]
		for (time_period in c("calibration","refit")){
			chain_paths_time_period <- chain_paths[grepl(time_period, chain_paths)]
			print(chain_paths_time_period)
			savepath <- gsub("_chain[1234]","_summary",chain_paths_time_period[[1]])
			out <- combine_chains(chain_paths_time_period)
			
			effectiveSize(out$samples)
#			plot(out$samples[,1:10])
#			plot(out$samples[,21:30])
			
			# Calculate summary and save output.
			param_summary <- fast.summary.mcmc(out$samples)
			plot_summary <- fast.summary.mcmc(out$samples2)

			out_summary <- list(samples = out$samples,
									param_summary = param_summary,
									plot_summary = plot_summary,
									gelman = gelman.diag(out$samples))
			saveRDS(out_summary, savepath, compress = F)
			cat(paste0("\nSaved summary output for ", model_name, ", ", time_period, ", ",
								 rank, " to: \n", savepath, "\n"))
		}
	}
}
 



