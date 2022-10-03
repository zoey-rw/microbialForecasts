# Fitting shannon diversity forecasting models for various Nimble models and time periods

pacman::p_load(reshape2, parallel)
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
#source("./functions/MCMC_diversity.r")

# Create parameters to pass
params = data.frame(group = c("ITS", "16S"),
										min.date = c("20130601","20151101", "20130601","20151101"),
										max.date = c("20151101","20180101", "20170101","20200101"),
										scenario =  c("Legacy only","2 year current methods",
																	"Legacy + 1 year current methods",
																	"Full dataset"),
										model_name = c("cycl_only","all_covariates")) %>%
	expand(nesting(min.date, max.date, scenario),
				 nesting(group),
				 model_name) %>%
	mutate(temporalDriverUncertainty =  T,
				 spatialDriverUncertainty =  T)


# Create function that calls run_MCMC for each uncertainty scenario
run_scenarios <- function(j) {
	out <- run_MCMC_diversity(iter = 750000, burnin = 300000, thin = 10,
									#iter = 600000, burnin = 300000, thin = 10,
		#iter = 200000, burnin = 100000, thin = 5,
									test=F, scenario = params$scenario,
									group = params$group[[j]],
									model_name = params$model_name[[j]],
									min.date = params$min.date[[j]],
									max.date = params$max.date[[j]])
	return()
}
#
# print(params[j,])
# run_scenarios(j)

# NOT PARALLEL
# for(j in c(5:8)){
# 												print(params[j,])
# 												run_scenarios(j)
# 											}

to_run1 <- grep("20160101", params$min.date)
#to_run2 <- grep("refit", params$time_period)
#to_run <- intersect(to_run1, to_run2)

# # Create cluster and pass it everything in the workspace
library(doParallel)
cl <- makeCluster(8, type="PSOCK", outfile="")
registerDoParallel(cl)
#
output.list = foreach(j= c(9:16),
#output.list = foreach(j= to_run,
											.errorhandling = 'pass') %dopar% {

	run_scenarios(j)
}
message("Diversity models fit.")
