# Combine chains from each taxon model, and create basic summary stats
source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")


# # For interactive use/testing.
# model_name <- "cycl_only"
#model_id = "env_cycl_animal_pathogen_20151101_20180101"

# Summarize all available models
params_in = read.csv(here("data/clean/model_input_df.csv")) %>% filter(min.date=="20151101" & max.date == "20180101")
model_id_list = unique(params_in$model_id)

# Or, subset to certain model runs
# params <- params_in %>% filter(min.date=="20151101" & max.date == "20180101")
# params <- params_in %>% filter(min.date=="20151101" &
# 															 	max.date == "20180101",
# 															 model_name=="cycl_only")
# model_id_list <- gsub("env_cov","all_covariates",model_id_list)
# model_id_list = unique(params$model_id)

cl <- makeCluster(28, outfile="")
registerDoParallel(cl)
#for(model_id in model_id_list){
foreach(model_id=model_id_list, .errorhandling = "pass") %dopar% {
	source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")

	# Do we want to keep all the chain files separately? Deleting them will save space
	delete_samples_files = F

	#for (model_name in c("all_covariates", "cycl_only")) {
	file.list <- list.files(path = here("data/model_outputs/logit_beta_regression/"),
													pattern = "_chain",
													recursive = T,
													full.names = T)

	# Subset to newest output files
	info <- file.info(file.list)

	# Subset more than 200MB (otherwise, some parameters probably didn't save)
	large_enough <- rownames(info[which(info$size > 100000000), ])

	#newer <- rownames(info[which(info$mtime > "2023-03-09 06:00:00 EDT"), ])
	newer <- rownames(info[which(info$mtime > "2023-04-28 18:00:00 EDT"), ])

	# Don't want files still being written - at least 2 min old
	older <- rownames(info[which(info$mtime < (Sys.time()-120)), ])

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
		return(message("Skipping ", chain_paths, "; no chains"))
	if (length(chain_paths) == 1) {
		return(message("Skipping ", chain_paths, "; only one chain"))
		#next()
	}
	savepath <- gsub("_chain[1234567]", "", chain_paths[[1]])


	# Don't run loop if the samples file is already newer than the chain files
	if (samples_exists(chain_paths[[1]])) message("Summary samples file already exists") else {
		
		
	# Calculate summary on each output subset, using custom summary function
	out <- combine_chains_simple(chain_paths)

	# Remove chains if one has very different values than the other 3
	chains <- out$samples
	means <- lapply(chains, function(x) mean(x[,"intercept"], na.rm=T))
	scaled_means = scale(unlist(means))
	potential_outlier <- which(abs(scaled_means) > 1.3)
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



	param_summary <- fast.summary.mcmc(out$samples)
	plot_summary <- fast.summary.mcmc(out$samples2)
	es <- effectiveSize(out$samples)
	if (length(out$samples) > 1) {
		gelman_out <- cbind(gelman.diag(out$samples, multivariate = F)[[1]], es)
	} else gelman_out = cbind(`Point est.`=NA, `Upper C.I.`=NA, es)

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






plot_est = parse_plot_mu_vars(out_summary$plot_summary$statistics)


eval = merge(out_summary$metadata$model_data %>% filter(species != "other"), plot_est)

# Check the predicted~observed
plot(eval$truth ~ eval$Mean); abline(0,1)


# True values
ggplot(eval) +
	geom_point(aes(x = timepoint, y = as.numeric(truth), color = siteID)) +
	geom_line(aes(x = timepoint, y = as.numeric(truth), color = siteID, group=plotID))
# Estimated values
ggplot(eval) +
	geom_point(aes(x = timepoint, y = Mean, color = siteID)) +
	geom_line(aes(x = timepoint, y = Mean, color = siteID, group=plotID))

# Check the predicted~observed, for low values in particular
ggplot(eval) +
	geom_point(aes( x = Mean, y = as.numeric(truth),color = siteID))  +
	ylab("Observed") + xlab("Predicted") +
	geom_abline(slope = 1,intercept = 0) +
	theme_minimal(base_size = 18) +
	ggtitle("Inflating plot means at lower than 5% abundance")


ggplot(eval %>% filter(plotID %in% "CPER_004")) +
	geom_point(aes( x = date_num, y = as.numeric(truth),
									color = siteID), show.legend = F, alpha=.5)  +
	geom_ribbon(aes( x = date_num, ymin = Mean-(1.96*SD), ymax=Mean+(1.96*SD)), alpha=.2) +
	geom_line(aes( x = date_num, y = Mean,
								 color = siteID, group=plotID), show.legend = F, alpha=.5)  +
	ylab("Observed") + xlab("Predicted") +
	theme_minimal(base_size = 18) +
	ggtitle("Overestimating plot means") +
	facet_wrap(~siteID)


