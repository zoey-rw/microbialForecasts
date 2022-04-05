# Running shannon diversity model for diff uncertainty scenarios and covariate combos

pacman::p_load(reshape2, parallel, nimble, coda, tidyverse) 

source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/MCMC_diversity.r")

# Create parameters to pass	
params = data.frame(index = 1:8,
										scenario = c("no_uncertainty_ITS", "spatial_uncertainty_ITS",
																 "temporal_uncertainty_ITS", "full_uncertainty_ITS",
																 "no_uncertainty_16S", "spatial_uncertainty_16S",
																 "temporal_uncertainty_16S", "full_uncertainty_16S"),
										group = c(rep("ITS", 4),rep("16S", 4)),
										temporalDriverUncertainty = c(F, F, T, T, F, F, T, T),
										spatialDriverUncertainty = c(F, T, F, T, F, T, F, T))

params <- rbind(cbind(params, model_name = "cycl_only"),
								cbind(params, model_name = "all_covariates"))
params <- rbind(cbind(params, time_period = "calibration"),
								cbind(params, time_period = "refit"))


# Create function that calls run_MCMC for each uncertainty scenario
run_scenarios <- function(j) {
	out <- run_MCMC(iter = 1000000, burnin = 500000, thin = 15, 
									#iter = 600000, burnin = 300000, thin = 10, 
									test=F, 
									group = params$group[[j]], 
									temporalDriverUncertainty = params$temporalDriverUncertainty[[j]], 
									spatialDriverUncertainty = params$spatialDriverUncertainty[[j]], 
									scenario = params$scenario[[j]],
									model_name = params$model_name[[j]],
									time_period = params$time_period[[j]])
	return()
}
# 
# j <- 8
# print(params[j,])
# run_scenarios(j)

# NOT PARALLEL
# for(j in c(5:8)){
# 												print(params[j,])
# 												run_scenarios(j)
# 											}
# 

# Only run scenarios with covariate uncertainty
#to_run <- grep("full_uncertainty", params$scenario)


to_run1 <- grep("full_uncertainty", params$scenario) 
to_run2 <- grep("refit", params$time_period) 
#to_run2 <- grep("cal", params$time_period) 
to_run <- intersect(to_run1, to_run2)
#to_run <- to_run1

# # Create cluster and pass it everything in the workspace
library(doParallel)
cl <- makeCluster(4, type="PSOCK", outfile="")
#cl <- makeCluster(8, type="PSOCK", outfile="")
registerDoParallel(cl)
# 
# # output.list = foreach(j= c(3,4,5,6),
output.list = foreach(j= c(28,32),
#output.list = foreach(j= to_run,
#output.list = foreach(j= to_run,
											.errorhandling = 'pass') %dopar% {
	run_scenarios(j)
}
# # saveRDS(output.list, out.path)


