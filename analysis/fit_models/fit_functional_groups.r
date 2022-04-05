#Get arguments from the command line
# argv <- commandArgs(TRUE)
# # Check if the command line is not empty and convert values to numerical values
# if (length(argv) > 0){
# 	k <- as.numeric( argv[1] )
# } else {
# 	k=1
# }

#### Run on all groups ----


# Running dirichlet model 
pacman::p_load(reshape2, parallel, nimble, coda, tidyverse) 
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/MCMC_functional_groups.r")

# Create parameters to pass	
n.groups <- length(keep_fg_names)
params = data.frame(group = rep(1:n.groups),
										scenario = c(rep("no_uncertainty", n.groups), 
																 rep("spatial_uncertainty", n.groups), 
																 rep("temporal_uncertainty", n.groups), 
																 rep("full_uncertainty", n.groups)),
										temporalDriverUncertainty = c(rep(FALSE, n.groups*2),
																									rep(TRUE, n.groups*2)),
										spatialDriverUncertainty = c(rep(FALSE, n.groups),
																								 rep(TRUE, n.groups),
																								 rep(FALSE, n.groups),
																								 rep(TRUE, n.groups)),
										ranks = rep(keep_fg_names, 4))

params <- rbind(cbind(params, model_name = "cycl_only"),
								cbind(params, model_name = "all_covariates"))
params <- rbind(cbind(params, time_period = "calibration"),
								cbind(params, time_period = "refit"))

# Create function that calls run_MCMC for each uncertainty scenario
run_scenarios <- function(j) {
	print(params[j,])
	#run_MCMC(k = params$group[[j]], iter = 60000, burnin = 20000, thin = 3, 
	#run_MCMC(k = params$group[[j]], iter = 1000, burnin = 500, thin = 1, test = F,
  #run_MCMC(k = params$group[[j]], iter = 30000, burnin = 15000, thin = 1, test = F,
	run_MCMC(k = params$group[[j]], iter = 800000, burnin = 400000, thin = 20, test=F,
									temporalDriverUncertainty = params$temporalDriverUncertainty[[j]], 
									spatialDriverUncertainty = params$spatialDriverUncertainty[[j]], 
									scenario = params$scenario[[j]],
					 model_name = params$model_name[[j]],
					 time_period = params$time_period[[j]]
	)
	return()
}


model_outputs <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/fg_summaries.rds")
gelman <- model_outputs$full_uncertainty$gelman_list
names(gelman) <- c(sort(keep_fg_names), sort(keep_fg_names), sort(keep_fg_names))
max_psrf <- lapply(gelman, function(x) max(x[,1])) %>% unlist()
unconverged <- names(gelman[max_psrf > 1.15])


library(doParallel)
cl <- makeCluster(28, type="PSOCK", outfile="")
registerDoParallel(cl)


to_run1 <- grep("full_uncertainty", params$scenario) 

#to_run2 <- grep("refit", params$time_period) 
#to_run <- intersect(to_run1, to_run2)

#to_run3 <- grep("cycl", params$model_name) 
#to_run <- intersect(to_run, to_run3)

to_run4 <- which(params$ranks %in% unconverged) 
to_run <- intersect(to_run1, to_run4)


output.list = foreach(j=to_run,
											.errorhandling = 'pass') %dopar% {
												run_scenarios(j)
	return()
											}
stopCluster(cl)






# # get missing scenarios for re-running
# missing_scenarios <- read.csv("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/analysis/missing_scenarios.csv")
# params$specific <- paste(params$ranks, params$scenario)
# missing_scenarios$specific <- paste(missing_scenarios$scenario, missing_scenarios$rank.name)
# to_rerun <- as.numeric(rownames(params[params$specific %in% missing_scenarios$specific,]))
# 
# fixed_fg_names <- c("cellulolytic", "assim_nitrite_reduction", "dissim_nitrite_reduction",
# 										"assim_nitrate_reduction", "n_fixation", "dissim_nitrate_reduction",
# 										"nitrification", "denitrification", "chitinolytic", "lignolytic",
# 										"methanotroph", "copiotroph", "oligotroph", "glucose_simple", "glycine_simple",
# 										"streptomycin_antibiotic","fsfeso4_anaerobic", "ironcitrate_anaerobic",
# 										"light_stress", "endophyte", "plant_pathogen",
# 										"animal_pathogen", "ectomycorrhizal", "lichenized", "wood_saprotroph",
# 										"soil_saprotroph", "litter_saprotroph", "saprotroph")
# to_rerun <- as.numeric(rownames(params[params$scenario=="full_uncertainty" & !params$ranks %in% fixed_fg_names,]))
# output.list = foreach(j=to_rerun,
# 											.errorhandling = 'pass') %dopar% {
# 	run_scenarios(j)
# 	return()
# }



