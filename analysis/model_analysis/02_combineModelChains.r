# Combine chains from each taxon model, and create basic summary stats
source("source.R")

# For phenology analysis: use cycl_only models for 2013-2015 data
model_name <- "cycl_only"

# For phenology analysis: use our actual model IDs for 2013-2015 data
model_id_list = c(
  "cycl_only_saprotroph_20130601_20151101",
  "cycl_only_ectomycorrhizal_20130601_20151101", 
  "cycl_only_cellulolytic_20130601_20151101",
  "cycl_only_assim_nitrite_reduction_20130601_20151101"
)

cat("Processing logit beta regression phenology models for 2013-2015 data\n")
cat("Model type:", model_name, "\n")
cat("Total models to process:", length(model_id_list), "\n\n")

# Summarize all available models
# params_in = read.csv(here("data/clean/model_input_df.csv")) %>% filter(!min.date=="20151101" | !max.date == "20180101")

# Since we have our specific model IDs, use them directly
# model_id_list <- chain_model_ids

cat("Model ID list length:", length(model_id_list), "\n")
cat("Model IDs:", paste(model_id_list, collapse = ", "), "\n")

# cl <- makeCluster(28, outfile="")
# registerDoParallel(cl)
#for(model_id in model_id_list){
# foreach(model_id=model_id_list, .errorhandling = "pass") %dopar% {
for(model_id in model_id_list) {
	# source("source.R")  # Commented out to avoid potential issues

	# Do we want to keep all the chain files separately? Deleting them will save space
	delete_samples_files = F

	#for (model_name in c("all_covariates", "cycl_only")) {
	file.list <- list.files(path = here("data/model_outputs/logit_phenology_2013_2015/"),
													pattern = "_chain",
													recursive = T,
													full.names = T)

	# Subset to newest output files
	info <- file.info(file.list)

	# Subset more than 100KB (for our test files, reasonable threshold)
	# large_enough <- rownames(info[which(info$size > 100000000), ])
	# large_enough <- rownames(info[which(info$size > 5000000), ])
	large_enough <- rownames(info[which(info$size > 100000), ])
	
	#newer <- rownames(info[which(info$mtime > "2023-03-09 06:00:00 EDT"), ])
	# Remove date filtering for now - all files are recent
	newer <- rownames(info)  # Include all files

	# Don't want files still being written - at least 2 min old
	# older <- rownames(info[which(info$mtime < (Sys.time()-120)), ])
	
	# Don't want files still being written - at least 1 min old (for testing)
	# older <- rownames(info[which(info$mtime < (Sys.time()-3600)), ]
	older <- rownames(info[which(info$mtime < (Sys.time()-60)), ]

	# If deleting sample files, we don't want models still being run - at least 18 hours old
	# if (delete_samples_files) {
	# }
	file.list <- file.list[file.list %in% newer & file.list %in% large_enough]

	message("Searching model outputs for ",model_id)
	message("Total files in file.list: ", length(file.list))
	message("First few files: ", paste(head(file.list, 3), collapse = ", "))

	# Subset to files of interest
	chain_paths <- file.list[grepl(model_id, file.list, fixed = T)]

	message("Found ", length(chain_paths), " chain files for ", model_id)
	if (length(chain_paths) > 0) {
		message("Chain files: ", paste(basename(chain_paths), collapse = ", "))
	}

	if (length(chain_paths) == 0) {
		message("Skipping ", model_id, "; no chains found")
		next
	}
	if (length(chain_paths) == 1) {
		message("Skipping ", model_id, "; only one chain found (need multiple chains to combine)")
		next
	}
	savepath <- gsub("_chain[1234567]", "", chain_paths[[1]])


	# Don't run loop if the samples file is already newer than the chain files
	if (samples_exists(chain_paths[[1]])) {
		message("Summary samples file already exists")
	} else {
		# Calculate summary on each output subset, using custom summary function
		message("Combining chains for ", model_id, " with ", length(chain_paths), " chains...")
		out <- tryCatch({
			combine_chains_simple_new(chain_paths)
		}, error = function(e) {
			message("Error combining chains for ", model_id, ": ", e$message)
			return(NULL)
		})
		
		if (is.null(out)) {
			message("Failed to combine chains for ", model_id, ", skipping...")
			next
		}

		# Remove chains if one has very different values than the other 3
		chains <- out$samples
		if (length(chains) == 0) {
			message("No chains found in output for ", model_id, ", skipping...")
			next
		}
		
		# Check if intercept parameter exists
		if (!"intercept" %in% colnames(chains[[1]])) {
			message("No intercept parameter found in chains for ", model_id, ", skipping...")
			next
		}
		
		# Ensure chains have reasonable sample sizes before outlier detection
		min_chain_size <- 100  # Minimum reasonable chain size
		if (min(nrow(chains[[1]]), nrow(chains[[2]]), nrow(chains[[3]]), nrow(chains[[4]])) < min_chain_size) {
			message("Chains too small for reliable outlier detection in ", model_id, ", skipping outlier removal")
		} else {
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
			if (samples_exists(chain_paths[[1]])) {
				unlink(chain_paths)
				message("Deleting samples files, e.g.: ", chain_paths[[1]])
			}
		} #else message("Not deleting samples files, e.g.: ", chain_paths[[1]])
	} # End else block
} # End model_id loop

# stopCluster(cl)

# Remove the problematic return() call at the end
# The script should end here without any return() statements

# Comment out plotting code since it requires out_summary which may not exist
# if the script skips processing due to existing files

# plot_est = parse_plot_mu_vars(out_summary$plot_summary$statistics)
# 
# 
# eval = merge(out_summary$metadata$model_data %>% filter(species != "other"), plot_est)
# 
# # Check the predicted~observed
# plot(eval$truth ~ eval$Mean); abline(0,1)
# 
# 
# # True values
# ggplot(eval) +
# 	geom_point(aes(x = timepoint, y = as.numeric(truth), color = siteID)) +
# 	geom_line(aes(x = timepoint, y = as.numeric(truth), color = siteID, group=plotID))
# # Estimated values
# ggplot(eval) +
# 	geom_point(aes(x = timepoint, y = Mean, color = siteID)) +
# 	geom_line(aes(x = timepoint, y = Mean, color = siteID, group=plotID))
# 
# # Check the predicted~observed, for low values in particular
# ggplot(eval) +
# 	geom_point(aes( x = Mean, y = as.numeric(truth),color = siteID))  +
# 	ylab("Observed") + xlab("Predicted") +
# 	theme_minimal(base_size = 18) +
# 	ggtitle("Inflating plot means at lower than 5% abundance")
# 
# 
# ggplot(eval %>% filter(plotID %in% "CPER_004")) +
# 	geom_point(aes( x = date_num, y = as.numeric(truth),
# 									color = siteID), show.legend = F, alpha=.5)  +
# 	geom_ribbon(aes( x = date_num, ymin = Mean-(1.96*SD), ymax=Mean+(1.96*SD)), alpha=.2) +
# 	geom_line(aes( x = date_num, y = Mean,
# 								 color = siteID, group=plotID), show.legend = F, alpha=.5)  +
# 	geom_point(aes( x = date_num, y = Mean,
# 								 color = siteID, group=plotID), show.legend = F, alpha=.5)  +
# 	ylab("Observed") + xlab("Predicted") +
# 	theme_minimal(base_size = 18) +
# 	ggtitle("Overestimating plot means") +
# 	facet_wrap(~siteID)


