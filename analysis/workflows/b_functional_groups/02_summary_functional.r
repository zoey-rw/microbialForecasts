# Summarize output from functional group models
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

# time_period = "20151101_20180101"
# time_period = "20151101_20200101"
#
# # From before the beta-regression fits.
#
# file.list <- list.files(path = here("data", "model_outputs/functional_groups/"),
# 												pattern = time_period,
# 												recursive = T, full.names = T)
# # testing <- summarize_fg_div_model(file.list[[1]], drop_other=T)
#
# cl <- makeCluster(36, outfile="")
# registerDoParallel(cl)
#
# file_summaries = foreach(f=file.list, .errorhandling = "pass") %dopar% {
# 	pacman::p_load(microbialForecast)
# 		out <- summarize_fg_div_model(f, drop_other=T)
# 												 	return(out)
# 												 }
#
# summary_df  <- map_df(file_summaries, 1) %>%
# 	mutate(time_period = recode(as.character(time_period), !!!microbialForecast:::date_recode))
# plot_est    <- map_df(file_summaries, 2) %>%
# 	mutate(time_period = recode(as.character(time_period), !!!microbialForecast:::date_recode))
# gelman_list <- map_df(file_summaries, 3) %>%
# 	mutate(time_period = recode(as.character(time_period), !!!microbialForecast:::date_recode))
#
# out <- list(plot_est = plot_est,
# 						summary_df = summary_df,
# 						gelman_list = gelman_list)
#
# saveRDS(out, here("data", paste0("summary/fg_summaries_", time_period, ".rds")))


cl <- makeCluster(36, outfile="")
registerDoParallel(cl)

time_period="20151101_20200101"
time_period="20151101_20180101"

for (time_period in c("20151101_20180101", "20151101_20200101")) {
file.list <- list.files(path = here("data", "model_outputs/functional_groups/"),
												pattern = time_period,
												recursive = T, full.names = T)
file.list <- file.list[grepl("beta", file.list)]
if (time_period=="20151101_20200101") {
file.list <- file.list[grepl("summary", file.list)]
}
# Summarize the MCMC output into 3 data frames

# Checks if existing summary file is older than the MCMC file; if so, re-summarizes
file_summaries = foreach(f=file.list, .errorhandling = "pass") %dopar% {
	pacman::p_load(microbialForecast)

	# Do we want to keep the "other" category from the beta regression? yes, because needed for site effect analysis
	out <- summarize_fg_beta_model(f, drop_other=F)
	return(out)
}


summary_df  <- map_df(file_summaries, 1) %>%
	mutate(time_period = recode(as.character(time_period), !!!microbialForecast:::date_recode))
plot_est    <- map_df(file_summaries, 2) %>%
	mutate(time_period = recode(as.character(time_period), !!!microbialForecast:::date_recode)) %>% filter(taxon != "other")
gelman_list <- map_df(file_summaries, 3) %>%
	mutate(time_period = recode(as.character(time_period), !!!microbialForecast:::date_recode)) %>% distinct(.keep_all = T)

out <- list(plot_est = plot_est,
						summary_df = summary_df,
						gelman_list = gelman_list)

saveRDS(out, here("data", paste0("summary/bychain_beta_fg_summaries_", time_period, ".rds")))
#saveRDS(out, here("data", "beta_fg_summaries_20151101_20180101.rds")


}

stopCluster(cl)
