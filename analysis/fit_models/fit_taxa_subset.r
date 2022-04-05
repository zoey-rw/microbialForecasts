# Running dirichlet regression on 5 taxonomic ranks for soil fungi and bacteria
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/MCMC_taxa.r")
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

#Get arguments from the command line
argv <- commandArgs(TRUE)
# Check if the command line is not empty and convert values to numerical values
if (length(argv) > 0){
	k <- as.numeric( argv[1] )
}

# Create parameters to pass	
n.groups <- 10 # Number of taxonomic ranks
# params = data.frame(group = rep(1:n.groups, 2),
# 										scenario = c(rep("full_uncertainty", n.groups), 
# 																 rep("no_uncertainty", n.groups)),
# 										temporalDriverUncertainty = c(rep(TRUE, n.groups),
# 																									rep(FALSE, n.groups)),
# 										spatialDriverUncertainty = c(rep(TRUE, n.groups),
# 																								 rep(FALSE, n.groups)))

params = data.frame(group = 1:n.groups)
params <- rbind(cbind(params, model_name = "cycl_only"),
								cbind(params, model_name = "all_covariates"))
params <- rbind(cbind(params, time_period = "calibration"),
								cbind(params, time_period = "refit"))
params$temporalDriverUncertainty <- T
params$spatialDriverUncertainty <- T

to_run1 <- grep("full_uncertainty", params$scenario) 
to_run2 <- grep("calibration", params$time_period) 
#to_run2 <- grep("refit", params$time_period) 
to_run <- intersect(to_run1, to_run2)


# Running dirichlet model 
pacman::p_load(doParallel, reshape2, nimble, coda, tidyverse) 



# Create function that calls run_MCMC for each uncertainty scenario
run_scenarios <- function(j, chain) {
	print(params[j,])
	out <- run_MCMC(k = params$group[[j]], iter_per_chunk = 50000, init_iter = 100000,
									iter = 350000, burnin = 50000, thin = 20, chain_no = chain,
									test=F, 
									temporalDriverUncertainty = params$temporalDriverUncertainty[[j]], 
									spatialDriverUncertainty = params$spatialDriverUncertainty[[j]],
									scenario = params$scenario[[j]],
									model_name = params$model_name[[j]],
									time_period = params$time_period[[j]])
	return(out)
}



# 
# 
logfile <- paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/analysis/fit_models/log/", tax_names[params$group[[k]]], "_", params$model_name[[k]], "_", params$time_period[[k]], ".log")
cl <- makeCluster(4, type="PSOCK", outfile=logfile)
registerDoParallel(cl)
#Run for multiple chains, in parallel (via PSOCK)
output.list = foreach(chain=c(1:4),
											.errorhandling = 'pass') %dopar% {
												set.seed(chain)
												out_chain <- run_scenarios(j = k, chain = chain)
												return(out_chain)
											}

#run_scenarios(k)
stopCluster(cl)
