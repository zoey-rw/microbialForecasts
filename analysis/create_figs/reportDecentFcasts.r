# Report the number of decent (RSQ > .1) forecasts
source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")
hindcast_data <- readRDS(here("data/summary/all_hindcasts.rds"))


scores_list = readRDS(here("data/summary/scoring_metrics_cv.rds"))
converged = scores_list$converged_list

scores = scores_list$scoring_metrics %>% 
	filter(model_id %in% converged)

# How many groups had decent forecasts using any model? #120
good_rsq = scores %>% filter(RSQ.1 > .1)
unique(good_rsq$taxon)


# How many groups had decent forecasts using env_cycl model? #71
scores_env_cycl = scores %>% filter(model_name == "env_cycl" &
																			grepl("observed", site_prediction)) 
good_rsq = scores_env_cycl %>% filter(RSQ.1 > .1)
unique(good_rsq$taxon)
nrow(good_rsq)/nrow(scores_env_cycl) %>% print


# How many groups had decent forecasts using env_cov model? #83
scores_env_cov = scores %>% filter(model_name == "env_cov" &
																			grepl("observed", site_prediction)) 
good_rsq = scores_env_cov %>% filter(RSQ.1 > .1)
unique(good_rsq$taxon)
nrow(good_rsq)/nrow(scores_env_cov) %>% print
