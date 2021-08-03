#Get arguments from the command line
# argv <- commandArgs(TRUE)
# # Check if the command line is not empty and convert values to numerical values
# if (length(argv) > 0){
# 	k <- as.numeric( argv[1] )
# } else {
# 	k=1
# }


iter <- 1000
burnin <- 500
thin <- 1
test = T
k = 58
temporalDriverUncertainty <- TRUE
spatialDriverUncertainty <- TRUE

run_MCMC <- function(k = 17, 	
														 iter = 1000,
														 burnin = 500,
														 thin = 1,
														 test = F,
														 temporalDriverUncertainty = TRUE,
														 spatialDriverUncertainty = TRUE,
														scenario = NULL,
														 ...) {
	pacman::p_load(reshape2, parallel, nimble, coda, tidyverse) 
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepFunctionalData.r")
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
	
	# Read in microbial abundances
	d <- c(readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_16S_2021.rds"), 
				 readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_ITS_2021.rds"))
	
	# Only keep functional groups; subset to one rank.
	ranks.keep <- names(d)
	ranks.keep <- ranks.keep[!grepl("bac|fun$", ranks.keep)]
	rank.name <- ranks.keep[k]
	rank.df <- d[[rank.name]] 
	
	# Reduce size for testing
	if (test == T) {
		rank.df <- rank.df %>% arrange(siteID, plotID, dateID)
		rank.df = rank.df[1:500,]
	}
	
	# Prep model inputs/outputs.
	model.dat <- prepFunctionalData(rank.df = rank.df, min.prev = 3)

	# create data object
	constants <- list(N.plot =  length(unique(model.dat$plotID)), 
										N.spp = ncol(model.dat$y), 
										N.core = nrow(model.dat$y), 
										N.date = model.dat$N.date,
										N.site = length(unique(model.dat$siteID)),
										timepoint = model.dat$timepoint,
										mois = model.dat[["mois"]],
										temp = model.dat[["temp"]],
										mois_sd = model.dat[["mois_sd"]],
										temp_sd = model.dat[["temp_sd"]],
										pH = model.dat[["pH"]],
										pC = model.dat[["pC"]],
										pH_sd = model.dat[["pH_sd"]],
										pC_sd = model.dat[["pH_sd"]],
										nspp = model.dat[["nspp"]],
										rc_grass = model.dat[["rc_grass"]],
										plotID = model.dat$plotID,
										plot_site = model.dat$plot_site,
										plot_num = model.dat$plot_num,
										plot_site_num = model.dat$plot_site_num,
										plot_start = model.dat[["plot_start"]],
										plot_index = model.dat[["plot_index"]],
										site_start = model.dat[["site_start"]],
										N.beta = 6)
	truth <- model.dat$truth.plot.long # for output.
	
	initsList <- initsFun(constants)
	initsList$beta <- initsList$beta[1,]
	initsList$rho <- initsList$rho[1]
	initsList$plot_mu <-  matrix(rep(.55, constants$N.plot*constants$N.date), constants$N.plot, constants$N.date)
	initsList$Ex <-  matrix(rep(.55, constants$N.plot*constants$N.date), constants$N.plot, constants$N.date)
	initsList$sigma <- .1
	initsList$core_sd <- .1
	initsList$site_effect <- rep(.1, constants$N.site)
	
	
	## Configure & compile model
	Rmodel <- nimbleModel(code = nimbleModFunctional, 
												constants = constants, 
												data = list(y=model.dat$y), 
												inits = initsList) # Set initial values using function
	
#	Rmodel$checkConjugacy()
	cModel <- compileNimble(Rmodel)
	# Configure & compile MCMC
	myMCMC <- buildMCMC(conf = cModel, monitors = c("beta","sigma","rho","sig",
																									"site_effect"), monitors2 = c("plot_mu"))
	compiled <- compileNimble(myMCMC, project = cModel, resetFunctions = T)

	# Remove large objects due to data usage
	# rm(model.dat, d, prepModelData, rank.df); gc()
	 
	# Sample from MCMC
	samples <- runMCMC(compiled, niter = iter, nburnin = burnin, 
										 nchains = 3, samplesAsCodaMCMC = T, thin = thin)
	samples2 <- rm.NA.mcmc(samples$samples2)
	samples <- rm.NA.mcmc(samples$samples)
	
	# Calculate summary and save output.
	param_summary <- summary(samples)
	plot_summary <- summary(samples2)
	metadata <- c("rank.name" = rank.name,
								"niter" = iter,
								"nburnin" = burnin,
								"thin" = thin,
								"model_data" = truth)
	out <- list(samples = samples, 
							param_summary = param_summary, 
							metadata = metadata, 
							plot_summary = plot_summary)
	out.path <- paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/samples_", rank.name, "_", scenario,".rds")
	saveRDS(out, out.path)
	return(out)
}

#### Run on all groups ----


# Running dirichlet model 
if (!require("pacman")) install.packages("pacman") 
pacman::p_load(reshape2, parallel, nimble, coda, tidyverse) 


# Create parameters to pass	
n.groups <- 66 # Number of functional groups
params = data.frame(group = rep(1:n.groups, 4),
										scenario = c(rep("no_uncertainty", n.groups), 
																 rep("spatial_uncertainty", n.groups), 
																 rep("temporal_uncertainty", n.groups), 
																 rep("full_uncertainty", n.groups)),
										temporalDriverUncertainty = c(rep(FALSE, n.groups*2),
																									rep(TRUE, n.groups*2)),
										spatialDriverUncertainty = c(rep(FALSE, n.groups),
																								 rep(TRUE, n.groups),
																								 rep(FALSE, n.groups),
																								 rep(TRUE, n.groups)))

# Create function that calls run_MCMC for each uncertainty scenario
run_scenarios <- function(j) {
	print(params[j,])
	out <- run_MCMC(k = params$group[[j]], iter = 60000, burnin = 20000, thin = 3, 
									test=F, 
									temporalDriverUncertainty = params$temporalDriverUncertainty[[j]], 
									spatialDriverUncertainty = params$spatialDriverUncertainty[[j]], 
									scenario = params$scenario[[j]])
	return(out)
}


#### First two secnarios
# out.path <- paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/samples_fg_min3_part1.rds")
# # Run for all functional groups, in parallel (via PSOCK)
# output.list <- parLapply(cl, c(1:132), run_scenarios)
# saveRDS(output.list, out.path)
# 
# #### Latter two scenarios
# out.path <- paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/samples_fg_min3_part2.rds")
# # Run for all functional groups, in parallel (via PSOCK)
# output.list <- parLapply(cl, c(133:264), run_scenarios)
# saveRDS(output.list, out.path)


library(doParallel)
cl <- makeCluster(12, type="PSOCK", outfile="")
#cl <- makeCluster(8, type="PSOCK", outfile="")
registerDoParallel(cl)

#output.list <- run_scenarios(200)


d <- c(readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_16S_2021.rds"), 
			 readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_ITS_2021.rds"))
ranks.keep <- names(d)
ranks.keep <- ranks.keep[!grepl("bac|fun$", ranks.keep)]

params <- cbind(params, ranks = rep(ranks.keep, 4))

# get missing scenarios for re-running
missing_scenarios <- read.csv("missing_scenarios.csv")
params$specific <- paste(params$ranks, params$scenario)
missing_scenarios$specific <- paste(missing_scenarios$scenario, missing_scenarios$rank.name)
to_rerun <- params[params$specific %in% missing_scenarios$specific,]$group

# #### Latter two scenarios
output.list = foreach(j=to_rerun,
											.errorhandling = 'pass') %dopar% {
	run_scenarios(j)
	return()										
}
# saveRDS(output.list, out.path)
# 
# 
# 
# out.path <- paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/samples_fg_min3_part2_2.rds")
# output.list = foreach(j=233:264, 
# 											#.export=c("run_scenarios","params","run_MCMC"), 
# 											.errorhandling = 'pass') %dopar% {
# 												run_scenarios(j)
# 											}
# saveRDS(output.list, out.path)
# 
# 





# 
# sample.list <- lapply(output.list, "[[", 1)
# param.summary.list <- lapply(output.list, "[[", 2)
# metadata.list <- lapply(output.list, "[[", 3)
# plot.summary.list <- lapply(output.list, "[[", 4)
# 
# names(sample.list) <- params$scenario
# names(param.summary.list) <- params$scenario
# names(metadata.list) <- params$scenario
# names(plot.summary.list) <- params$scenario
# 
# saveRDS(list(sample.list = sample.list,
# 						 param.summary.list = param.summary.list, 
# 						 metadata.list = metadata.list,
# 						 plot.summary.list = plot.summary.list),
# 				"/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/div_samples_min3.rds")

stopCluster(cl)
