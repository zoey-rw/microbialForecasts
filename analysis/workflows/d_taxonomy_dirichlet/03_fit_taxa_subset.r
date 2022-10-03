# Running dirichlet regression on 5 taxonomic ranks for soil fungi and bacteria
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

#Get arguments from the command line
argv <- commandArgs(TRUE)
# Check if the command line is not empty and convert values to numerical values
if (length(argv) > 0){
	k <- as.numeric( argv[1] )
} else k = 1

# Create parameters to pass
n.groups <- 10 # Number of taxonomic ranks
# params = data.frame(group = rep(1:n.groups, 2),
# 										scenario = c(rep("full_uncertainty", n.groups),
# 																 rep("no_uncertainty", n.groups)),
# 										temporalDriverUncertainty = c(rep(TRUE, n.groups),
# 																									rep(FALSE, n.groups)),
# 										spatialDriverUncertainty = c(rep(TRUE, n.groups),
# 																								 rep(FALSE, n.groups)))

params = data.frame(group = 1:n.groups, rank.names = microbialForecast:::tax_names)
params <- rbind(cbind(params, model_name = "cycl_only"),
								cbind(params, model_name = "all_covariates"))
# params <- rbind(cbind(params, time_period = "calibration"),
# 								cbind(params, time_period = "refit"))
params$temporalDriverUncertainty <- T
params$spatialDriverUncertainty <- T
params <- rbind(cbind(params, min.date = "20151101", max.date = "20180101", scenario =  "2 year current methods"),
								cbind(params, min.date = "20151101", max.date = "20200101", scenario =  "Current only"))
								#cbind(params, min.date = "20130601", max.date = "20170101", scenario = "Legacy + 1 year current methods"),
								#cbind(params, min.date = "20130601", max.date = "20150101", scenario = "Legacy only"),
								#cbind(params, min.date = "20130601", max.date = "20200101", scenario = "Legacy + current (full dataset)"))


# Running dirichlet model
pacman::p_load(doParallel, reshape2, nimble, coda, tidyverse)




# to_rerun <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/taxa/to_rerun.rds")
#
# params$id <- paste(params$rank.names, params$time_period, params$model_name)
# to_rerun$id <- paste(to_rerun$rank.names, to_rerun$time_period, to_rerun$model)
# params <- params[which(params$id %in% to_rerun$id),]


params <- params[(params$scenario == "2 year current methods"),]

# params <- params[(params$rank.names=="order_bac" & params$model_name == "all_covariates")|
# 								 	(params$rank.names=="family_fun" & params$model_name == "all_covariates")|
# 								 	(params$rank.names=="class_bac" & params$model_name == "all_covariates"),]

# Create function that calls run_MCMC for each uncertainty scenario
run_scenarios <- function(j, chain) {
	print(params[j,])
	out <- run_MCMC_taxa(k = params$group[[j]],
											 iter_per_chunk = 50000,
											 init_iter = 150000,
											 iter = 350000,
											 burnin = 100000, thin = 10, chain_no = chain,
											 #iter = 2000, burnin = 1000, thin = 3, chain_no = chain,
											 test=F,
									temporalDriverUncertainty = params$temporalDriverUncertainty[[j]],
									spatialDriverUncertainty = params$spatialDriverUncertainty[[j]],
									scenario = params$scenario[[j]],
									model_name = params$model_name[[j]],
									time_period = params$time_period[[j]],
									min.date = params$min.date[[j]],
									max.date = params$max.date[[j]])
	return(out)
}






logfile <- paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/analysis/workflows/qsub/", params$rank.names[[k]], "_", params$model_name[[k]], "_", params$min.date[[k]], "_", params$max.date[[k]], ".log")
cl <- makeCluster(4, type="PSOCK", outfile=logfile)
registerDoParallel(cl)
#Run for multiple chains, in parallel (via PSOCK)
output.list = foreach(chain=c(1:4),
											.errorhandling = 'pass') %dopar% {
												source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
												set.seed(chain)
												out_chain <- run_scenarios(j = k, chain = chain)
												return(out_chain)
											}

#run_scenarios(k)
stopCluster(cl)

