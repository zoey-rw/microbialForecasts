# Combine chains
pacman::p_load(nimble, coda, tidyverse) 
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

model_name <- "cycl_only"

#for (model_name in c("cycl_only")){
for (model_name in c("all_covariates")){
	
#for (model_name in c("all_covariates", "cycl_only")){
		print(model_name)
	
	file.list <- list.files(path = paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/taxa/", model_name),
													pattern = "full_uncertainty_chain",
													full.names = T)
	
	# info <- file.info(file.list)
	# newer <- rownames(info[which(info$mtime > "2022-04-19 13:00:00 EDT"),])
	# file.list <- file.list[file.list %in% newer]

	rank <- "genus_bac"
	rank <- "phylum_bac"
	rank <- "class_fun"
	time_period <- "calibration"
	
	for (rank in tax_names[3]){
		print(rank)
		
		chain_paths <- file.list[grepl(rank, file.list)]
		if(length(chain_paths)==0) next()
		for (time_period in c("calibration","refit")){
		#for (time_period in c("refit")){
			

				print(time_period)
			
			chain_paths_time_period <- chain_paths[grepl(time_period, chain_paths)]
			if(length(chain_paths_time_period)==0) next()
			
			savepath <- gsub("_chain[1234]","_summary",chain_paths_time_period[[1]])
			print(savepath)
			#if (!file.exists(savepath)){

				
				# TO DO: remove the problematic chain here
				#if (rank == "phylum_fun" & model_name == "all_covariates" & time_period == "calibration") chain_paths_time_period <- chain_paths_time_period[c(1,3:4)]#
				if (rank == "class_bac" & model_name == "cycl_only") chain_paths_time_period <- chain_paths_time_period[2:4]#
				if (rank == "family_bac" & model_name == "cycl_only") chain_paths_time_period <- chain_paths_time_period[2:4]#
				
				print(chain_paths_time_period)
				
				
				if(length(chain_paths_time_period)==1) next()
	
			out <- combine_chains(chain_paths_time_period)
			
			es <- effectiveSize(out$samples)
			if (min(es) < 10) {
				#plot(out$samples[,1:10])
				#plot(out$samples[,21:30])
				print(es[1:50])
			}
			
			# Calculate summary and save output.
			param_summary <- fast.summary.mcmc(out$samples)
			plot_summary <- fast.summary.mcmc(out$samples2)

			gelman_out <- cbind(gelman.diag(out$samples)[[1]], es)
			
			out_summary <- list(samples = out$samples,
									param_summary = param_summary,
									plot_summary = plot_summary,
									metadata = out$metadata,
									gelman = gelman_out)
			
			saveRDS(out_summary, savepath, compress = F)
			cat(paste0("\nSaved summary output for ", model_name, ", ", time_period, ", ",
								 rank, " to: \n", savepath, "\n"))
			#} else print("Summary file already exists")
		}
	}
}



chain_paths = chain_paths_time_period
save = FALSE
cut_size1 = NULL
cut_size2 = NULL

combine_chains <- function(chain_paths, 
													 save = FALSE, 
													 cut_size1 = NULL, 
													 cut_size2 = NULL){
	require(coda)
	require(tidyverse)
	
	
	if (is.null(cut_size1)) cut_size1 <- 19999 
	if (is.null(cut_size2)) cut_size2 <- 9999 
	
	
	# initialize
	samples <- samples2 <- metadata <- list()
	for(i in 1:length(chain_paths)){

		print(i)
		# paste model file path to chain number
		chain <- readRDS(chain_paths[[i]])
		nrow_samples <- nrow(chain[[1]])
		nrow_samples2 <- nrow(chain[[2]])
		
		if (nrow_samples < cut_size1) {
			cut_size1 <- nrow_samples
		}
		
		samples[[i]] <- as.mcmc(as.matrix(window(chain[[1]], nrow_samples-cut_size1, nrow_samples, 1)))
		
		
		if (nrow_samples2 < cut_size2) {
			cut_size2 <- nrow_samples2
		}
		samples2[[i]] <- as.mcmc(as.matrix(window(chain[[2]], nrow_samples2-cut_size1, nrow_samples2, 1)))

		#metadata[[i]] <- chain[[3]]
	}
	
	
		nrows <- lapply(samples, nrow) %>% unlist()
		min_nrow <- min(nrows)
		for(i in 1:length(samples)){
			current_nrow <- nrow(samples[[i]])
			if (min_nrow < current_nrow){
				print(i)
				samples[[i]] <- as.mcmc(as.matrix(window(samples[[i]], (current_nrow-min_nrow+1), current_nrow, 1)))
			}
		}
	
		
		nrows <- lapply(samples2, nrow) %>% unlist()
		min_nrow <- min(nrows)
		for(i in 1:length(samples)){
			current_nrow <- nrow(samples2[[i]])
			if (min_nrow < current_nrow){
				print(i)
				samples2[[i]] <- as.mcmc(as.matrix(window(samples2[[i]], (current_nrow-min_nrow+1), current_nrow, 1)))
			}
		}
		
	
	out <- list(samples = as.mcmc.list(samples),
							samples2 = as.mcmc.list(samples2),
							metadata = metadata[[1]])
	
	if(!isFALSE(save)){
		saveRDS(out, file = save)
	}
	return(out)
}
