# Combine chains from each taxon model, and create basic summary stats
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

# Do we want to keep all the chain files separately? Deleting them will save space
delete_samples_files = T

# For interactive use/testing.
model_name <- "cycl_only"
model_name <- "all_covariates"
rank <- "order_bac"
rank <- "order_fun"
rank <- "phylum_bac"
rank <- "phylum_fun"
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
	file.list <- list.files(path = paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/logit_beta_regression/", model_name),
			pattern = "_chain",
			full.names = T)

	# Subset to newer output files
	info <- file.info(file.list)
	newer <- rownames(info[which(info$mtime > "2022-08-02 00:00:00 EDT"), ])
	newer <- rownames(info[which(info$mtime > "2022-12-02 00:00:00 EDT"), ])

	# Don't want files still being written - at least 2 min old
	older <- rownames(info[which(info$mtime < (Sys.time()-120)), ])
	#older <- rownames(info[which(info$mtime < (Sys.time()-3000)), ])

	# If deleting sample files, we don't want models still being run - at least 18 hours old
	# if (delete_samples_files) {
	# 	older <- rownames(info[which(info$mtime < (Sys.time()-64800)), ])
	# }

	file.list <- file.list[file.list %in% newer & file.list %in% older]

	for (rank in microbialForecast:::tax_names[1]) {
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


				if (length(chain_paths_time_period) == 0) {
					next()
				}


				savepath <- gsub("_chain[1234567]", "", chain_paths_time_period[[1]])



				if (length(chain_paths_time_period) == 1) {
					message("Skipping ", chain_paths_time_period, "; only one chain")
					next()
				}

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
					#}

					saveRDS(out_summary, savepath, compress = F)
					message("Saved samples and summary output for ",model_name,", ",time_period,", ",rank," to: ",savepath)


if (min(es) < 500) {
	print(es[1:50])
}



# If the summary now exists, delete the chains
if (delete_samples_files){
	if (samples_exists(chain_paths_time_period[[1]])) {
		unlink(chain_paths_time_period)
		message("Deleting samples files, e.g.: ", chain_paths_time_period[[1]])
	}
} #else message("Not deleting samples files, e.g.: ", chain_paths_time_period[[1]])

				#}
		} # End time-period loop
	} # End species loop
} # End rank loop
	} # End model type loop






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


