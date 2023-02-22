

source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/microbialForecast/R/run_MCMC_bychain.r")


#Get arguments from the command line
argv <- commandArgs(TRUE)
# Check if the command line is not empty and convert values to numerical values
if (length(argv) > 0){
	k <- as.numeric( argv[1] )
} else {
	k=1
}

#### Run on all groups ----


# Running dirichlet model
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")


grp_key <- cbind.data.frame(group = 1:10, rank.names = microbialForecast:::tax_names)
params <- stack(microbialForecast:::rank_spec_names) %>% select(species = values, rank.names = ind)
params <- merge(params, grp_key)
params <- rbind(cbind(params, model_name = "cycl_only"),
								cbind(params, model_name = "all_covariates"))
params$temporalDriverUncertainty <- T
params$spatialDriverUncertainty <- T
params <- rbind(cbind(params, min.date = "20151101", max.date = "20180101", scenario =  "2 year current methods"),
								cbind(params, min.date = "20151101", max.date = "20200101", scenario =  "Current only"))


params <- params %>% filter(`max.date` == "20180101" &
															model_name == "all_covariates" &
															#model_name == "cycl_only" &
															rank.names %in% c("phylum_fun","phylum_bac"))

# params <- params %>% filter(`max.date` == "20180101" &
# 															model_name == "all_covariates" &
# 															rank.names == "phylum_bac")

nchains = 4

#params <- crossing(params, chain_no = 1:nchains)

# Create function that calls run_MCMC for each uncertainty scenario
run_scenarios <- function(j, chain_no) {
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/microbialForecast/R/run_MCMC_bychain.r")

	print(params[j,])
	run_MCMC_bychain(k = params$group[[j]],
									# 				 burnin = 75000,
									#  thin = 10,
									#  iter_per_chunk = 10000,
									#  init_iter = 100000,
									 burnin = 200000,
									 thin = 25,
									 iter_per_chunk = 50000,
									 init_iter = 250000,
													 test=F,
													 scenario = params$scenario[[j]],
													 model_name = params$model_name[[j]],
													 chain_no=chain_no,
													 min.date = params$min.date[[j]],
													 max.date = params$max.date[[j]],																													s = params$species[[j]],
									min_eff_size_perchain = 500
)
	return()
}

# Run using command line inputs
#run_scenarios(k, chain_no=1)


# Run in parallel
cl <- makeCluster(4, type="PSOCK", outfile="")
registerDoParallel(cl)
output.list = foreach(chain_no=1:4,
											.errorhandling = 'pass') %dopar% {
												run_scenarios(j = k, chain_no = chain_no)
												return()
											}
stopCluster(cl)

