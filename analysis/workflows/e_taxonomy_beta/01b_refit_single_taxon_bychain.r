# Fit taxa one-by-one to evaluate convergence.
#install.packages("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/microbialForecast_0.1.0.tar.gz", repos = NULL, type ="source")
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")


#Get arguments from the command line
argv <- commandArgs(TRUE)
# Check if the command line is not empty and convert values to numerical values
if (length(argv) > 0){
	j <- as.numeric( argv[1] )
} else j = 1

# Running single-taxon model for 236 taxa, on full NEON time-series

# Create parameters to pass
grp_key <- cbind.data.frame(group = 1:10, rank.names = microbialForecast:::tax_names)
params <- stack(microbialForecast:::rank_spec_names) %>% select(species = values, rank.names = ind)
params <- merge(params, grp_key)
params <- rbind(cbind(params, model_name = "cycl_only"),
								cbind(params, model_name = "all_covariates"))
params$temporalDriverUncertainty <- T
params$spatialDriverUncertainty <- T
params <- rbind(cbind(params, min.date = "20151101", max.date = "20180101", scenario =  "2 year current methods"),
								cbind(params, min.date = "20151101", max.date = "20200101", scenario =  "Current only"))

#params <- params[(params$scenario == "Current only"  & params$model_name == "all_covariates"),]
params <- params[(params$scenario == "Current only"  & params$model_name == "cycl_only"),]


# Create function that calls run_MCMC for each uncertainty scenario
run_scenarios <- function(j, chain, ...) {
	pacman::p_load(doParallel, reshape2)

	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

													print(params[j,])

													out <- run_MCMC_single_taxon_bychain(k = params$group[[j]],
																															 #iter = 500000, burnin = 300000,
																															 thin = 20,
																															 #iter = 1001,
																															 burnin = 150000, init_iter=200000, iter_per_chunk=10000,
																															 test=F,

																															 #init_iter = 500000, burnin = 300000, iter_per_chunk = 50000,
																															 scenario ="201511_202001",
																															 model_name = params$model_name[[j]],
																															 #time_period = params$time_period[[j]],
																															 min.date = params$min.date[[j]],
																															 max.date = params$max.date[[j]],
																															 chain_no = chain,
																																 s = params$species[[j]])
													return(out)
}


#test_out <- run_scenarios(k)

logfile <- paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/analysis/workflows/qsub/", params$rank.names[[j]], "_", params$model_name[[j]], "_", params$min.date[[j]], "_", params$max.date[[j]], "_", params$species[[j]], ".log")

cl <- makeCluster(4, type="PSOCK", outfile=logfile)
registerDoParallel(cl)
#Run for multiple chains, in parallel (via PSOCK)
output.list = foreach(chain=c(1:4),
											.errorhandling = 'pass') %dopar% {
												source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
												set.seed(chain)
												out_chain <- run_scenarios(j = j, chain = chain)
												return(out_chain)
											}

#run_scenarios(k)
stopCluster(cl)