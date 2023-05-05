

#Get arguments from the command line (run with qsub script & OGE scheduler)
argv <- commandArgs(TRUE)
# Check if the command line is not empty and convert values to numerical values
if (length(argv) > 0){
	k <- as.numeric( argv[1] )
} else {
	k=1
}

# Run with at least 4 cores available (one MCMC chain per core)
nchains = 4

#### Run on all groups ----

source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")


params_in = read.csv(here("data/clean/model_input_df.csv"),
										 colClasses = c(rep("character", 4),
										 							 rep("logical", 2),
										 							 rep("character", 4)))

rerun_list = readRDS(here("data/summary/unconverged_taxa_list.rds"))
#converged_list = readRDS(here("data/summary/converged_taxa_list.rds"))
converged_list  = readRDS(here("data/summary/weak_converged_taxa_list.rds"))

# # Subset to specific params for running


#  Get models that have not converged
params <- params_in %>% filter(!model_id %in% c(converged_list) &
															 	`max.date` == "20180101")# &
															 #	model_name == "env_cycl")



# Double check that none have converged
params <- params %>% filter(!model_id %in% converged_list)



# Check that none are currently running/being written
chain_list <- 
	list.files(path = here("data/model_outputs/logit_beta_regression/"),pattern = "_chain",
												recursive = T, full.names = T)
# Subset to newest output files
info <- file.info(chain_list)
# Files still potentially being written - less than 30 min old
newest <- rownames(info[which(info$mtime > (Sys.time()-18000)), ])
currently_running <- basename(chain_list[chain_list %in% newest])
currently_running_id <- gsub("_chain[1234567].rds|samples_", "", currently_running) %>% unique

params <- params %>% filter(!model_id %in% currently_running_id)

params <- params %>% filter(!model_id %in% c("cycl_only_thelephorales_20151101_20180101","cycl_only_tremellales_20151101_20180101","cycl_only_sordariales_20151101_20180101","cycl_only_thelebolales_20151101_20180101"))



# first run fullest models
#params <- params %>% arrange(desc(model_name))

# first run sparsest models
params <- params %>% arrange(model_name)

# first run fg models
params <- params %>% arrange(fcast_type)

# Create function that calls run_MCMC for each uncertainty scenario
run_scenarios <- function(j, chain_no) {
	source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")
	source("/projectnb/dietzelab/zrwerbin/microbialForecasts/microbialForecast/R/run_MCMC_bychain.r")

	print(params[j,])
	run_MCMC_bychain(k = params$group[[j]],
									 rank.name = params$rank.name[[j]],
									 s = params$species[[j]],
									 model_id = params$model_id[[j]],
									 burnin = 1100000,
									 #thin = 10,
									 iter_per_chunk = 100000,
									 #init_iter = 200000,
									 init_iter = 1300000,
									 #  burnin = 200000,
									 thin = 30,
									 thin2 = 50,
									 #  iter_per_chunk = 50000,
									 #  init_iter = 250000,
									 # 				 test=F,

									 # burnin = 500,
									 # thin = 5,
									 # iter_per_chunk = 500,
									 # init_iter = 1000,
									 # test=T,
									 scenario = params$scenario[[j]],
									 model_name = params$model_name[[j]],
									 chain_no=chain_no,
									 min.date = params$min.date[[j]],
									 max.date = params$max.date[[j]],
									 min_eff_size_perchain = 500
	)
	return()
}

# Run using command line inputs
#run_scenarios(k, chain_no=1)

# Run in parallel - must use PSOCK due to NIMBLE compilation
cl <- makeCluster(nchains, type="PSOCK", outfile="")
registerDoParallel(cl)
output.list = foreach(chain_no=1:nchains,
											.errorhandling = 'pass') %dopar% {
												run_scenarios(j = k, chain_no = chain_no)
												return()
											}

stopCluster(cl)








# j = k = 29
# chain_no=4
#
# run_scenarios <- function(j, chain_no) {
# 	source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")
# 	source("/projectnb/dietzelab/zrwerbin/microbialForecasts/microbialForecast/R/run_MCMC_bychain.r")
#
# 	print(params[j,])
# 	run_MCMC_bychain(k = params$group[[j]],
# 									 rank.name = params$rank.name[[j]],
# 									 s = params$species[[j]],
# 									 burnin = 1000,
# 									 iter_per_chunk = 2000,
# 									 init_iter = 2000,
# 									 thin = 25,
# 									 test=T,
# 									 scenario = params$scenario[[j]],
# 									 model_name = "env_cov",
# 									 chain_no=chain_no,
# 									 min.date = "20151101",
# 									 max.date = "20180101",
# 									 min_eff_size_perchain = 500,
# 									 temporalDriverUncertainty = F,
# 									 spatialDriverUncertainty = F
# 	)
# 	return()
# }
#
# test = readRDS("/projectnb/dietzelab/zrwerbin/microbialForecasts/data/model_outputs/logit_beta_regression/env_cov/samples_env_cov_lignolytic_20151101_20180101_chain4.rds")


