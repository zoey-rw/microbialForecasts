# Fit taxa one-by-one to evaluate convergence.


# # For testing
iter <- 500
burnin <- 200
thin <- 1
test =F
k = 1
temporalDriverUncertainty <- TRUE
spatialDriverUncertainty <- TRUE
scenario <- "full_uncertainty"
s = "acidobacteriota"


run_MCMC <- function(k = 1, iter = 1000,  burnin = 500, thin = 1,
										 test = F,
										 temporalDriverUncertainty = TRUE, spatialDriverUncertainty = TRUE,
										 scenario=NULL,
										 s = "acidobacteriota",
										 ...) {
	pacman::p_load(reshape2, parallel, nimble, coda, tidyverse) 
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepTaxonomicData.r")
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/cyclical_model_source.r")
	
	# Read in microbial data
	d <- c(readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_16S_2021.rds"), 
				 readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_ITS_2021.rds"))
	
	# Subset to one rank.
	rank.name <- tax_names[k]
	rank.df <- d[[rank.name]] 
	rank.df <- rank.df[which(!rank.df$siteID %in% c("ABBY","LAJA")),]
	
	# Reduce size for testing
	if (test == T) {
		rank.df <- rank.df %>% arrange(siteID, plotID, dateID)
		rank.df = rank.df[1:500,]
	}
	
	# Prep model inputs/outputs.
	print(paste0("Preparing model data for ", rank.name))
	
	spec_names <- colnames(rank.df)[!colnames(rank.df) %in% c("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date","other")] 
	rank.df_spec <- rank.df %>% 
			select("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date", !!s) 
	rank.df_spec$other <- 1-rank.df_spec[[s]]
	
	model.dat <- prepTaxonomicData(rank.df = rank.df_spec, min.prev = 3)
	out.path <- paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/convergence_testing/samples_", rank.name, "_", s, ".rds")
	print(paste("Completed model data for", rank.name, ", group:", s))
	
	mo <- month(as.Date(paste0(colnames(model.dat$mois), "01"), format="%Y%m%d"))
	y_sin = sin((2*pi*mo)/12)
	y_cos = cos((2*pi*mo)/12)
	constants <- list(N.plot =  length(unique(model.dat$plotID)), 
										N.spp = ncol(model.dat$y), 
										N.core = nrow(model.dat$y), 
										N.date = model.dat$N.date,
										N.site = length(unique(model.dat$siteID)),
										timepoint = model.dat$timepoint,
										sin_mo = y_sin,
										cos_mo = y_cos,
										mois = model.dat[["mois"]],
										temp = model.dat[["temp"]],
										mois_sd = model.dat[["mois_sd"]],
										temp_sd = model.dat[["temp_sd"]],
										pH = model.dat[["pH"]],
										pC = model.dat[["pC"]],
										pH_sd = model.dat[["pH_sd"]],
										pC_sd = model.dat[["pC_sd"]],
										nspp = model.dat[["nspp"]],
										relEM = model.dat[["relEM"]],
										LAI = model.dat[["LAI"]],
										plotID = model.dat$plotID,
										plot_site = model.dat$plot_site,
										plot_num = model.dat$plot_num,
										plot_site_num = model.dat$plot_site_num,									
										plot_start = model.dat[["plot_start"]],
										plot_index = model.dat[["plot_index"]],
										site_start = model.dat[["site_start"]],
										N.beta = 2
	)
	
	# for output.
	metadata <- list("rank.name" = rank.name,
									 "niter" = iter,
									 "nburnin" = burnin,
									 "thin" = thin,
									 "model_data" = model.dat$truth.plot.long)
	
	inits <- initsFun(constants, type = "tax")
	## Configure & compile model
	Rmodel <- nimbleModel(code = nimbleModTaxa_cycl_only,
												constants = constants, data = list(y=model.dat$y),
												inits = inits)
	# Compile model
	cModel <- compileNimble(Rmodel)
	nimbleOptions(multivariateNodesAsScalars = TRUE)
	# Configure & compile MCMC
	mcmcConf <- configureMCMC(cModel, monitors = c("beta","sigma","site_effect",
																								 "sig","intercept", 
																								 "rho"),
														monitors2 = c("plot_rel"), thin2 = 25,
														useConjugacy = T)

	myMCMC <- buildMCMC(mcmcConf)
	compiled <- compileNimble(myMCMC, project = Rmodel, resetFunctions = TRUE)

	
	# Sample from MCMC
	samples <- runMCMC(compiled, niter = iter, nburnin = burnin, 
										 nchains = 3, samplesAsCodaMCMC = T, thin = thin)
	
	cat(paste("Finished sampling for ", rank.name, ", group: ", s)); gc()
	
	# Calculate summary and save output.
	param_summary <- fast.summary.mcmc(samples$samples)
	plot_summary <- fast.summary.mcmc(samples$samples2)
	
	out <- list(samples = samples$samples,
							param_summary = param_summary,
							metadata = metadata,
							plot_summary = plot_summary,
							gelman = gelman.diag(samples$samples))
	saveRDS(out, out.path, compress = F)
	cat(paste("Saved summary output for ", rank.name, ", group: ", s, "to: \n", out.path))
	return("ok")
}


# Running dirichlet model 
pacman::p_load(doParallel, reshape2, nimble, coda, tidyverse) 

	
# Create function that calls run_MCMC for each uncertainty scenario
run_scenarios <- function(k) {
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/cyclical_model_source.r")
	
	# Read in microbial data
	d <- c(readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_16S_2021.rds"), 
				 readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_ITS_2021.rds"))
	# Subset to one rank.
	rank.name <- tax_names[k]
	rank.df <- d[[rank.name]] 
	spec_names <- colnames(rank.df)[!colnames(rank.df) %in% c("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date","other")] 
	
	start_no <- ifelse(k == 5, 23, 29)
	
	for (s in spec_names[start_no:length(spec_names)]) {
		print(s)
	out <- run_MCMC(k = k, 
									#iter = 600000, burnin = 300000, thin = 10, 
									iter = 50001, burnin = 20002, thin = 3, 
									test=F, 
									temporalDriverUncertainty = TRUE, 
									spatialDriverUncertainty = TRUE,
									scenario ="full_uncertainty",
									s = s)
	}
	return(out)

}


# 
# 
# cl <- makeCluster(2, type="PSOCK", outfile="fit_taxa_convergence.log")
# registerDoParallel(cl)
# #Run for all 10 taxonomic ranks, in parallel (via PSOCK)
# output.list = foreach(k=c(1:2),
# 											.errorhandling = 'pass') %dopar% {
# 												run_scenarios(k)
# 											}



cl <- makeCluster(5, type="PSOCK", outfile="")
registerDoParallel(cl)
#Run for all 10 taxonomic ranks, in parallel (via PSOCK)
output.list = foreach(k=c(1:5),
											.errorhandling = 'pass') %dopar% {
												run_scenarios(k)
											}
stopCluster(cl)
