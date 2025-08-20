# Summarize all the estimated model parameters for paper supplement.

source("source.R")
library(kableExtra)


sum.in <- readRDS(here("data", paste0("summary/logit_beta_regression_summaries.rds")))
sum.all <- sum.in$summary_df  %>% mutate(tax_rank = rank,
																				 time_period = recode(time_period, !!!microbialForecast:::date_recode))
df <- sum.all %>%
	mutate(pretty_group = ifelse(group %in% c("16S","bac"), "Bacteria", "Fungi"))

# Add prettier data values
df$pretty_name <- recode(df$rank_only, !!!microbialForecast:::pretty_rank_names) %>%
	ordered(levels = c("Genus","Family","Order","Class","Phylum","Functional group","Diversity"))

df$only_rank <- sapply(str_split(df$rank_only, "_",  n = 2), `[`, 1) %>%
	ordered(levels = c("genus","family","order","class","phylum","functional","diversity"))

df = df %>% mutate(parameter_name = ifelse(is.na(beta), rowname, beta))
df_wide = df  %>% filter(time_period ==  "2015-11_2018-01" & parameter_name != "site_effect") %>%
	pivot_wider(id_cols = c("model_id","taxon","model_name","fcast_type","pretty_group","only_rank"),
							names_from = parameter_name,
							values_from = "Mean"#, values_fn = list
							) %>% as.data.frame() %>%
	mutate(model_name = recode(model_name, !!!model.labs)) %>%
	arrange(model_name, pretty_group, only_rank) %>%
	rename("Covariates included" = "model_name",
				 "Site effect variation" = "sig",
				 "Kingdom" = "pretty_group",
				 "Microbial group type" = "fcast_type",
				 "Group" = "taxon",
				 "Rank" = "only_rank")


kable(df_wide, "html") %>%
	kable_styling(bootstrap_options = c("striped", "hover")) %>%
	cat(., file = here("figures", "model_params.html"))


# Check number of plots per site - usually 10, some have fewer or more
# hindcast_data not available, skipping this analysis
# aggregate(plotID ~ siteID, hindcast_data, function(x) length(unique(x)))



# df_wide_scores will be defined after scores_list is loaded
# df_wide_scores = merge(df_wide, scores_list$scoring_metrics %>% select(model_id, mean_crps_sample, RSQ, RSQ.1, RMSE.norm), all.y=F)





# Read in forecast scores
#converged <- readRDS(here("data/summary/weak_converged_taxa_list.rds"))

scores_list = readRDS(here("data/summary/scoring_metrics_plsr2.rds"))
converged = scores_list$converged_list
#converged  = scores_list$converged_strict_list

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

df_wide_scores = merge(df_wide, cal_scores_to_merge, all.y=F)  %>% as.data.frame() %>% merge(hindcast_scores_to_merge)

kable(df_wide_scores, "html") %>%
	kable_styling(bootstrap_options = c("striped", "hover")) %>%
	cat(., file = here("figures", "model_params_scores.html"))
write.csv(df_wide_scores, here("figures", "model_params_scores.csv"))
write.csv(df_wide_scores, here("figures", "table_S2.csv"))


# Below is an unsuccessful attempt to bold only the significant predictors for each row

# this function comes from here: https://stackoverflow.com/questions/28166168/how-to-change-fontface-bold-italics-for-a-cell-in-a-kable-table-in-rmarkdown
format_cells <- function(df, rows ,cols, value = c("italics", "bold", "strikethrough")){

	# select the correct markup
	# one * for italics, two ** for bold
	map <- setNames(c("*", "**", "~~"), c("italics", "bold", "strikethrough"))
	markup <- map[value]

	for (r in rows){
		for(c in cols){

			# Make sure values are not factors
			df[[c]] <- as.character( df[[c]])

			# Update formatting
			df[r, c] <- paste0(markup, df[r, c], markup)
		}
	}

	return(df)
}

# Subset to significant predictors
cells_to_bold <- df  %>% filter(time_period ==  "2015-11_2018-01" & parameter_name != "site_effect" & significant==1) %>%
	pivot_wider(id_cols = c("model_id","taxon","model_name","fcast_type","pretty_group","only_rank"),
							names_from = parameter_name,
							values_from = "Mean"#, values_fn = list
	) %>% as.data.frame() %>%
	mutate(model_name = recode(model_name, !!!model.labs)) %>%
	arrange(model_name, pretty_group, only_rank) %>%
	rename("Covariates included" = "model_name",
				 "Kingdom" = "pretty_group",
				 "Microbial group type" = "fcast_type",
				 "Group" = "taxon",
				 "Rank" = "only_rank")

# Get the cell ids for each model_id x predictor combo
cell_ids_to_bold <- list()
for (i in 1:nrow(df_wide_scores)){
	model_id = df_wide_scores[i,]$model_id
	all_pred_cols = cells_to_bold[cells_to_bold$model_id==model_id,] %>% select("cos", "sin", "LAI", "Ectomycorrhizal trees",
																																							 "pH", "Temperature", "Moisture", "pC")
	col_names_to_bold = names(all_pred_cols)[!is.na(all_pred_cols)]
	cols_to_bold <- which(colnames(df_wide_scores[i,]) %in% col_names_to_bold)
	cell_ids_to_bold[[i]] <- cols_to_bold
}

# HOW DO I VECTORIZE
bolded <- apply(df_wide_scores, 1, format_cells, rows=1:10, cols=cell_ids_to_bold[1:10])
# TO DO: Apply the bolding to every row to highlight significant predictors.

df_wide_scores %>%
	format_cells(1, cell_ids_to_bold[[1]], "bold") %>%
	format_cells(2, cell_ids_to_bold[[2]], "bold") %>%
	format_cells(3, cell_ids_to_bold[[3]], "bold") %>%
	format_cells(4, cell_ids_to_bold[[4]], "bold") %>%
	format_cells(5, cell_ids_to_bold[[5]], "bold") %>%
	format_cells(50, cell_ids_to_bold[[50]], "bold") %>%
	knitr::kable()  %>%
	kable_styling(bootstrap_options = c("striped", "hover")) %>%
	cat(., file = here("figures", "model_params_scores_bolded.html"))
	#save_kable(here("figures", "model_params_scores_bolded.pdf")) %>%
	#kable_styling(bootstrap_options = c("striped", "hover")) %>%
	#cat(., file = here("figures", "model_params_scores_bolded.pdf"))


# Merge the data that has been read in so far
scores_df_env_cycl = scores_list$scoring_metrics %>%
	filter(model_id %in% converged &
				 	site_prediction == "New time (observed site)" &
				 	model_name=="env_cycl")

scores_df_env_cov = scores_list$scoring_metrics %>%
	filter(model_id %in% converged &
				 	#site_prediction != "New time x site (modeled effect)" &
				 	site_prediction == "New time (observed site)" &
				 	model_name=="env_cov")

scores_df_cycl_only = scores_list$scoring_metrics %>%
	filter(model_id %in% converged &
				 	#site_prediction != "New time x site (modeled effect)" &
				 	site_prediction == "New time (observed site)" &
				 	model_name=="cycl_only")

length(unique(scores_df_env_cycl$taxon))
length(unique(scores_df_env_cov$taxon))
length(unique(scores_df_cycl_only$taxon))

reasonable_acc_cycl_only <- scores_df_cycl_only[scores_df_cycl_only$RSQ.1 > .1,]
reasonable_cycl_only <- length(unique(reasonable_acc_cycl_only$model_id))
total_cycl_only <- length(unique(scores_df_cycl_only$model_id))
reasonable_cycl_only/total_cycl_only

reasonable_acc_env_cycl <- scores_df_env_cycl[scores_df_env_cycl$RSQ.1 > .1,]
reasonable_env_cycl <- length(unique(reasonable_acc_env_cycl$model_id))
total_env_cycl <- length(unique(scores_df_env_cycl$model_id))
reasonable_env_cycl/total_env_cycl

table(scores_df_env_cycl$fcast_type)



# Merge the data that has been read in so far
scores_site_env_cycl = scores_list$scoring_metrics_site %>%
	filter(model_id %in% converged &
				 	model_name=="env_cycl" &
				 	!siteID %in% "MLBS")
length(unique(scores_site_env_cycl$siteID))
oldsites <- scores_site_env_cycl[scores_site_env_cycl$site_prediction=="New time (observed site)",]
newsites <- scores_site_env_cycl[scores_site_env_cycl$site_prediction=="New time x site (random effect)",]
newsites_modeled <- scores_site_env_cycl[scores_site_env_cycl$site_prediction=="New time x site (modeled effect)",]


mean(oldsites$RSQ, na.rm=T)
mean(newsites$RSQ, na.rm=T)
mean(newsites_modeled$RSQ, na.rm=T)

mean(newsites$RSQ.1, na.rm=T)
mean(newsites_modeled$RSQ.1, na.rm=T)

mean(newsites$mean_crps_sample, na.rm=T)
mean(newsites_modeled$mean_crps_sample, na.rm=T)


ggplot(scores_site_env_cycl %>%
			 	filter(site_prediction != "New time (observed site)" &
			 				 	!siteID %in% "MLBS" &

			 				 	fcast_type=="Functional")) +
	geom_point(aes(x = site_prediction, y = RSQ.1, color = siteID), position = position_jitter(height=.01, width=.1)) +
	facet_wrap(~taxon)



most_predictable <- scores_list$scoring_metrics_long %>% filter(metric == "RSQ.1" & score > .5)
slightly_predictable <- scores_list$scoring_metrics_long %>% filter(metric == "RSQ.1" & score > .1)


