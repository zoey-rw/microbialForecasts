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
pacman::p_load(doParallel) 
#source(here::here("source.R"))
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
source(here("functions", "MCMC_functional_groups.r"))


# Create parameters to pass	
params = data.frame(group = rep(1:length(keep_fg_names)), 
										ranks = keep_fg_names,
										min.date = c("20151101", "20130601"),
										 max.date = c("20180101", "20170101"),
										 scenario =  c("2 year current methods", 
										 							"Legacy + 1 year current methods"),
										model_name = c("cycl_only","all_covariates")) %>% 
	expand(nesting(min.date, max.date, scenario), 
								 nesting(group, ranks),
								 model_name) %>% 
	mutate(temporalDriverUncertainty =  T, 
									spatialDriverUncertainty =  T)


# Create function that calls run_MCMC for each uncertainty scenario
run_scenarios <- function(j) {
	print(params[j,])
	#run_MCMC(k = params$group[[j]], iter = 1000, burnin = 500, thin = 1, test = F,
  #run_MCMC(k = params$group[[j]], iter = 30000, burnin = 15000, thin = 1, test = F,
	run_MCMC(k = params$group[[j]], iter = 300000, burnin = 150000, thin = 10, test=F,
									scenario = params$scenario[[j]],
					 model_name = params$model_name[[j]],
					 time_period = params$time_period[[j]],
  min.date = params$min.date[[j]],
  max.date = params$max.date[[j]])
	return()
}

# Run using command line inputs
run_scenarios(k)


# model_outputs <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/fg_summaries.rds")
# gelman <- model_outputs$full_uncertainty$gelman_list
# names(gelman) <- c(sort(keep_fg_names), sort(keep_fg_names), sort(keep_fg_names))
# max_psrf <- lapply(gelman, function(x) max(x[,1])) %>% unlist()
# unconverged <- names(gelman[max_psrf > 1.15])


# Run in parallel
cl <- makeCluster(8, type="PSOCK", outfile="")
registerDoParallel(cl)
to_run <- grep("2 year", params$scenario) 
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



