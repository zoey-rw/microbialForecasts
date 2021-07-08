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
														 ...) {
	pacman::p_load(reshape2, parallel, nimble, coda, tidyverse) 
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepDirichletData.r")
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
	
	# Set model run parameters 
	# iter <- 200
	# burnin <- 100
	# thin <- 1

	# Read in microbial abundances
	d <- c(readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_16S_2021.rds"), 
				 readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_ITS_2021.rds"))
	
	# Only keep functional groups
	ranks.keep <- names(d)
	ranks.keep <- ranks.keep[!grepl("bac|fun$", ranks.keep)]
	
	# Subset to one rank.
	rank.name <- ranks.keep[k]
	rank.df <- d[[rank.name]] 
	
	# Reduce size for testing
	if (test == T) {
		rank.df <- rank.df %>% arrange(siteID, plotID, dateID)
		
		rank.df = rank.df[1:500,]
	}
	
	# Prep model inputs/outputs.
	model.dat <- prepModelData(rank.df = rank.df, min.prev = 3)

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
										N.beta = 6,
										alpha0 =  rep(0, ncol(model.dat$y)))
	truth <- model.dat$truth.plot.long # for output.

	## Configure & compile model
	Rmodel <- nimbleModel(code = nimbleModLong, 
												constants = constants, 
												data = list(y=model.dat$y), 
												inits = initsFun(constants)) # Set initial values using function
	cModel <- compileNimble(Rmodel)
	
	# Configure & compile MCMC
	myMCMC <- buildMCMC(conf = cModel, monitors = c("beta","sigma","rho",
																									"site_effect","intercept"),
											monitors2 = c("plot_rel"), useConjugacy = F,
											control = list(scale=.1))
	compiled <- compileNimble(myMCMC, project = cModel, resetFunctions = T)

	
	# Remove large objects due to data usage
	rm(model.dat, d, prepModelData, rank.df)
	gc()

	
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
	
	out.path <- paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/samples_", rank.name, "_min3_10tax.rds")
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
	out <- run_MCMC(k = params$group[[j]], iter = 275000, burnin = 200000, thin = 10, 
									test=F, 
									temporalDriverUncertainty = params$temporalDriverUncertainty[[j]], 
									spatialDriverUncertainty = params$spatialDriverUncertainty[[j]])
	return(out)
}

# Create cluster and pass it everything in the workspace
# cl <- makeCluster(24, outfile="")
# clusterExport(cl, ls())

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
# 

library(doParallel)
cl <- makeCluster(36, type="PSOCK", outfile="")
registerDoParallel(cl)


# #### Latter two scenarios
out.path <- paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/samples_fg_min3_part2_1.rds")
output.list = foreach(j=199:264,
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
