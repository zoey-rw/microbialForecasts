# Running dirichlet regression on 5 taxonomic ranks for soil fungi and bacteria

# For testing
iter <- 5000
burnin <- 2000
thin <- 1
test = T
k = 1
temporalDriverUncertainty <- TRUE
spatialDriverUncertainty <- TRUE
scenario <- "full_uncertainty"

run_MCMC <- function(k = 17, iter = 1000,  burnin = 500, thin = 1,
										 test = F,
										 temporalDriverUncertainty = TRUE, spatialDriverUncertainty = TRUE,
										 scenario=NULL,
										 ...) {
	pacman::p_load(reshape2, parallel, nimble, coda, tidyverse) 
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepTaxonomicData.r")
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
	
	# Read in microbial data
	d <- c(readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_16S_2021.rds"), 
				 readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_ITS_2021.rds"))
	
	# Subset to one rank.
	ranks.keep <- c("phylum_bac", "class_bac", "order_bac", "family_bac", "genus_bac", 
									"phylum_fun", "class_fun", "order_fun", "family_fun", "genus_fun")
	rank.name <- ranks.keep[k]
	rank.df <- d[[rank.name]] 
	
	# Reduce size for testing
	if (test == T) {
		rank.df <- rank.df %>% arrange(siteID, plotID, dateID)
		rank.df = rank.df[1:500,]
	}
	
	# Prep model inputs/outputs.
	print(paste0("Preparing model data for ", rank.name))
	model.dat <- prepTaxonomicData(rank.df = rank.df, min.prev = 3)
	out.path <- paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/samples_", rank.name, "_", scenario, ".rds")
	print(paste0("Completed model data for ", rank.name))

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
											constants = constants, data = list(y=model.dat$y),
											inits = initsFun(constants))
# Compile model
cModel <- compileNimble(Rmodel)
# Configure & compile MCMC
mcmcConf <- configureMCMC(cModel, monitors = c("beta","sigma","site_effect","intercept", "rho"),
													monitors2 = c("plot_rel"), useConjugacy = T#, control = c(scale = .1)
													)
mcmcConf$removeSamplers(c('site_effect'))
mcmcConf$addSampler(target = c('site_effect'), type = 'AF_slice')
myMCMC <- buildMCMC(mcmcConf)
# myMCMC <- buildMCMC(conf = cModel, monitors = c("beta","sigma","site_effect","intercept", "rho"),
# 													monitors2 = c("plot_rel"), useConjugacy = F)
compiled <- compileNimble(myMCMC, project = Rmodel, resetFunctions = TRUE)
# Remove large objects due to data usage
 rm(model.dat, d, rank.df); gc()

# Sample from MCMC
samples <- runMCMC(compiled, niter = iter, nburnin = burnin, 
									 nchains = 3, samplesAsCodaMCMC = T, thin = thin)
samples2 <- rm.NA.mcmc(samples$samples2)
samples <- rm.NA.mcmc(samples$samples)

# Calculate summary and save output.
param_summary <- summary(samples)
plot_summary <- summary(samples2)
metadata <- c("niter" = iter,
						 "nburnin" = burnin,
						 "thin" = thin, "scenario" = scenario,
						 "model_data" = truth)
out <- list(samples = samples, 
						param_summary = param_summary, 
						metadata = metadata, 
						plot_summary = plot_summary)
saveRDS(out, out.path)
print(paste0("Saved ", rank.name, " to ", out.path))
return(out)
}


# Running dirichlet model 
pacman::p_load(doParallel, reshape2, nimble, coda, tidyverse) 


# Create parameters to pass	
n.groups <- 10 # Number of taxonomic ranks
params = data.frame(group = rep(1:n.groups, 2),
										scenario = c(rep("no_uncertainty", n.groups), 
																 rep("full_uncertainty", n.groups)),
										temporalDriverUncertainty = c(rep(FALSE, n.groups),
																									rep(TRUE, n.groups)),
										spatialDriverUncertainty = c(rep(FALSE, n.groups),
																								 rep(TRUE, n.groups)))

# Create function that calls run_MCMC for each uncertainty scenario
run_scenarios <- function(j) {
	print(params[j,])
	out <- run_MCMC(k = params$group[[j]], iter = 3000, burnin = 1000, thin = 5, 
									test=F, 
									temporalDriverUncertainty = params$temporalDriverUncertainty[[j]], 
									spatialDriverUncertainty = params$spatialDriverUncertainty[[j]],
									scenario = params$scenario[[j]])
	return(out)
}


cl <- makeCluster(2, type="PSOCK", outfile="")
registerDoParallel(cl)

# Run for all 10 taxonomic ranks, in parallel (via PSOCK)
output.list = foreach(j=13:20,
											#.export=c("run_scenarios","params","run_MCMC"), 
											.errorhandling = 'pass') %dopar% {
												run_scenarios(j)
											}

stopCluster(cl)

