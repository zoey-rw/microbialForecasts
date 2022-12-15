# This script has yet to be written.
# Should summarize all the model parameters for paper supplement.


sum.all <- readRDS(here("data/summary/all_fcast_effects.rds"))

sum.all_wide <- sum.all  %>% filter(time_period ==  "2015-11_2018-01"  & fcast_type == "Taxonomic" &
																			model_name == "all_covariates" & beta %in% c("sin","cos")) %>%
	pivot_wider(id_cols = c("taxon","model_name","fcast_type",
																																"pretty_name","pretty_group","only_rank","siteID","time_period"),
																										names_from = beta,
																										values_from = "Mean", values_fn = list)



df_all_cov <- sum.all %>% filter(time_period == "2015-11_2020-01" &
																 	model_name == "all_covariates")
all_cov_vals <- df_all_cov %>% pivot_wider(id_cols = c("taxon","model_name","fcast_type",
																																"pretty_name","pretty_group","only_rank"),
																										names_from = beta,
																										values_from = "Mean")


df %>% filter(time_period == "2015-11_2020-01" &
							 	model_name == "all_covariates") %>%
	pivot_wider(id_cols = c("taxon","model_name","fcast_type","siteID",
															 "pretty_name","pretty_group","only_rank"),
									 names_from = rowname,
									 values_from = "Mean")
