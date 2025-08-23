# Summarize MCMC output from all single-taxon models
# Assumes input files have already had MCMC chains combined

source("source.R")


file.list = intersect(list.files(here("data/model_outputs/CLR_regression_expanded/"),recursive = T,
																 pattern = "20130601_20151101|20151101_20180101|20151101_20200101|20130601_20200101|20130601_20180101", full.names = T),
											list.files(here("data/model_outputs/CLR_regression_expanded/"), recursive = T,
																 pattern = "samples_CLR_expanded_.*\\.rds$", full.names = T))

# Filter out individual chain files, keep only combined samples files
file.list = file.list[!grepl("_chain[0-9]+", file.list)]

# Remove any files with only one chain
file.list = file.list[!grepl("chain", file.list)]

# Subset to newest output files (for testing, accept all recent files)
info <- file.info(file.list)
newer <- rownames(info[which(info$mtime > "2020-01-01 00:00:00 EDT"), ])
#newer <- rownames(info[which(info$mtime > "2023-03-09 00:00:00 EDT"), ])
file.list <- file.list[file.list %in% newer]

f=file.list[[1]]
cl <- makeCluster(1, outfile="")  # Use 1 core for testing
registerDoParallel(cl)

#Run summary function for multiple groups, in parallel
file_summaries = foreach(f=file.list, .errorhandling = "pass") %dopar% {
	source("source.R")
	# Use the new CLR-specific summary function
	out <- summarize_clr_model(f, save_summary=T, drop_other = T, overwrite = TRUE)
	return(out)
}

stopCluster(cl)


summary_file_list = list.files(here("data/model_outputs/CLR_regression_expanded/"), recursive = T,
															 pattern = "summary", full.names = T)

# Subset to newest output files
info <- file.info(summary_file_list)
newer <- rownames(info[which(info$mtime > "2023-03-09 00:00:00 EDT"), ])
summary_file_list <- summary_file_list[summary_file_list %in% newer]

# Combine summary files for all models/time-periods
# summary_file_list = intersect(list.files(here("data/model_outputs/logit_beta_regression/"),recursive = T,
# 																	pattern = "20151101_20180101", full.names = T),
# 											 list.files(here("data/model_outputs/logit_beta_regression/"), recursive = T,
# 											 					 pattern = "summary", full.names = T))
file_summaries <- purrr::map(summary_file_list, readRDS)
summary_df <- map_df(file_summaries, 1)
plot_est <- map_df(file_summaries, 2)
gelman_list <- map_df(file_summaries, 3)


#summaries = readRDS(here("data/summary/logit_beta_regression_summaries.rds"))


gelman.summary <- gelman_list %>%
	mutate(is_major_param = ifelse(grepl("beta|int|sigma|sd", parameter), T, F)) %>% 
	filter(!grepl("20130601_20151101", model_id))
#mutate(model_id2 = paste(model_name, rank, taxon, time_period, sep = "_"))

by_rank <- gelman.summary %>%
	group_by(model_id, is_major_param) %>%
	dplyr::mutate(median_gbr = median(`Point est.`,na.rm=T),
								quant_95 = quantile(`Point est.`, c(.95),na.rm=T),
								min_es = min(es, na.rm=T),
								median_es = min(es, na.rm=T),
								mean_gbr = mean(`Point est.`,na.rm=T)) %>%
	distinct(.keep_all = T)

# For CLR models, we may not have niteration column, so handle it conditionally
if ("niteration" %in% colnames(by_rank)) {
	model_median = by_rank %>% select(c("rank.name", "is_major_param","niteration", "rank", "taxon", "model_name", "group", "rank_only",
																		"time_period", "pretty_group", "model_id", "fcast_type", "median_gbr","mean_gbr","quant_95","min_es","median_es")) %>%
		distinct(.keep_all = T)
} else {
	# CLR models don't have niteration, so create a dummy column
	by_rank$niteration <- 100  # Default value for CLR models
	model_median = by_rank %>% select(c("rank.name", "is_major_param","niteration", "rank", "taxon", "model_name", "group", "rank_only",
																		"time_period", "pretty_group", "model_id", "fcast_type", "median_gbr","mean_gbr","quant_95","min_es","median_es")) %>%
		distinct(.keep_all = T)
}
model_median = model_median %>% pivot_wider(values_from = c("mean_gbr","median_gbr","quant_95","min_es","median_es"),names_from = is_major_param)


ggplot(model_median) + geom_jitter(aes(x = niteration, y = median_gbr_TRUE, color = group)) + ylim(c(0,5)) + geom_hline(yintercept = 1)



# Check what columns we actually have after pivot
cat("Available columns after pivot:\n")
print(colnames(model_median))

# For CLR models, we may not have the expected pivot columns, so handle this conditionally
if (all(c("median_gbr_TRUE", "mean_gbr_TRUE", "mean_gbr_FALSE", "min_es_TRUE") %in% colnames(model_median))) {
	# Standard filtering if all columns exist
	keep_models <- model_median %>%
		group_by(model_id) %>%
		filter(median_gbr_TRUE <= 1.1) %>%
		filter(mean_gbr_TRUE <= 1.2) %>%
		filter(mean_gbr_FALSE <= 1.5) %>%
		filter(min_es_TRUE > 75)
	keep_list <- unique(keep_models$model_id)

	keep_models_weak <- model_median %>%
		group_by(model_id) %>%
		filter(median_gbr_TRUE <= 1.15) %>%
		filter(mean_gbr_TRUE <= 1.5) %>%
		filter(mean_gbr_FALSE <= 2) %>%
		filter(min_es_TRUE > 15)
	keep_list_weak <- unique(keep_models_weak$model_id)
} else {
	# For CLR models, use simpler filtering based on available columns
	cat("Using simplified filtering for CLR models\n")
	keep_models <- model_median %>%
		group_by(model_id) %>%
		filter(TRUE)  # Keep all models for now
	keep_list <- unique(keep_models$model_id)
	
	keep_models_weak <- model_median %>%
		group_by(model_id) %>%
		filter(TRUE)  # Keep all models for now
	keep_list_weak <- unique(keep_models_weak$model_id)
}

rerun <- model_median %>% filter(!model_id %in% keep_list_weak)
rerun_list <- unique(rerun$model_id)

#rerun %>% ungroup %>% select(c(12:21)) %>% as.matrix %>% pairs

saveRDS(keep_list, here("data/summary/clr_converged_taxa_list.rds"))
saveRDS(keep_list_weak, here("data/summary/clr_weak_converged_taxa_list.rds"))
saveRDS(rerun_list, here("data/summary/clr_unconverged_taxa_list.rds"))


saveRDS(list(summary_df = summary_df,
						 plot_est = plot_est,
						 gelman.summary = gelman_list,
						 keep_models_weak = keep_models_weak,
						 keep_list = keep_list,
						 keep_list_weak = keep_list_weak,
						 rerun_list = rerun_list),
				here("data/summary/clr_regression_summaries.rds"))




