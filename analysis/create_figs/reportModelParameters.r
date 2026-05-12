# Summarize all the estimated model parameters for paper supplement.

source("source.R")
library(kableExtra)


sum.in <- readRDS(here("data", paste0("summary/logit_beta_fixed_priors_summaries.rds")))
sum.all <- sum.in$summary_df  %>% mutate(tax_rank = rank)
# Use existing pretty_group column (already has correct Bacteria/Fungi values)
df <- sum.all

# Add prettier data values
df$pretty_name <- recode(df$rank_only, !!!microbialForecast:::pretty_rank_names) %>%
	ordered(levels = c("Genus","Family","Order","Class","Phylum","Functional group","Diversity"))

df$only_rank <- sapply(str_split(df$rank_only, "_",  n = 2), `[`, 1) %>%
	ordered(levels = c("genus","family","order","class","phylum","functional","diversity"))

df = df %>% mutate(parameter_name = ifelse(is.na(beta), rowname, beta))
df_wide = df  %>% filter(time_period == "20130601_20180101" & parameter_name != "site_effect") %>%
	distinct(model_id, taxon, model_name, fcast_type, pretty_group, only_rank, parameter_name, .keep_all = TRUE) %>%
	pivot_wider(id_cols = c("model_id","taxon","model_name","fcast_type","pretty_group","only_rank"),
							names_from = parameter_name,
							values_from = "Mean"
							) %>% as.data.frame() %>%
	mutate(model_name = recode(model_name, !!!model.labs)) %>%
	arrange(model_name, pretty_group, only_rank) %>%
	rename("Covariates included" = "model_name",
				 "Kingdom" = "pretty_group",
				 "Microbial group type" = "fcast_type",
				 "Group" = "taxon",
				 "Rank" = "only_rank") %>%
	{ if ("site_effect_sd" %in% colnames(.)) rename(., "Site effect variation" = "site_effect_sd") else . }


kable(df_wide, "html") %>%
	kable_styling(bootstrap_options = c("striped", "hover")) %>%
	cat(., file = here("figures", "model_params.html"))


# Check number of plots per site - usually 10, some have fewer or more
# hindcast_data not available, skipping this analysis
# aggregate(plotID ~ siteID, hindcast_data, function(x) length(unique(x)))



# df_wide_scores will be defined after scores_list is loaded
# df_wide_scores = merge(df_wide, scores_list$scoring_metrics %>% select(model_id, mean_crps_sample, RSQ, RSQ.1, RMSE.norm), all.y=F)





# Read in forecast scores
scores_list = readRDS(here("data/summary/scoring_metrics_plsr2.rds"))

hindcast_scores_to_merge = scores_list$scoring_metrics %>% ungroup() %>%
	select(model_id, site_prediction, mean_crps_sample, RSQ, RSQ.1,RMSE.norm) %>%
	filter(!grepl("random", site_prediction)) %>%
	rename("Mean hindcast CRPS"="mean_crps_sample",
				 "Overall hindcast R-squared" = "RSQ",
				 "Overall hindcast R-squared (1:1 line)" = "RSQ.1",
				 "Normalized root mean squared error" = "RMSE.norm") %>%
	pivot_wider(id_cols = model_id, names_from = site_prediction,
							values_from = c("Mean hindcast CRPS", "Overall hindcast R-squared",
															"Overall hindcast R-squared (1:1 line)","Normalized root mean squared error"))
	#select(model_id, mean_crps_sample, RSQ)

cal_scores_to_merge = scores_list$calibration_metrics %>% ungroup() %>%
	select(model_id, calibration_RSQ = RSQ)
#select(model_id, mean_crps_sample, RSQ)

df_wide$model_id <- gsub("_beta_regression$", "", df_wide$model_id)
df_wide_scores = merge(df_wide, cal_scores_to_merge, all.x=TRUE)  %>% as.data.frame() %>% merge(hindcast_scores_to_merge, all.x=TRUE)

kable(df_wide_scores, "html") %>%
	kable_styling(bootstrap_options = c("striped", "hover")) %>%
	cat(., file = here("figures", "model_params_scores.html"))
write.csv(df_wide_scores, here("figures", "model_params_scores.csv"))
write.csv(df_wide_scores, here("figures", "table_S2.csv"))

cat("Saved: figures/model_params.html, model_params_scores.html, model_params_scores.csv, table_S2.csv\n")

