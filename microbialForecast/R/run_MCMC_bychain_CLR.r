#' @title run_spatial_model
#' @description run_spatial_model
#' @export
#
# For testing
# iter <- 500
# burnin <- 1000
# thin <- 1
# test =F
# k = 1
# scenario <- "full_uncertainty"
# s = "acidobacteriota"
# k = 6
# s = "ascomycota"
# model_name = "all_covariates"
# chain_no = 1
# min.date = "20151101"
# max.date = "20200101"
# max.date = "20180101"
# iter_per_chunk = 500
# init_iter = 5000
# model_name = "cycl_only"
# min_eff_size_perchain = 10
# max_loops=100
# thin2 = 20
# rank.name = "phylum_bac"
# s = "armatimonadota"
# #rank.name = s = species = "nitrification"
# model_id = NULL

run_MCMC_bychain_CLR <- function(rank.name = "phylum_bac",
														 s = "acidobacteriota",
														 temporalDriverUncertainty = TRUE,
														 spatialDriverUncertainty = TRUE,
														 scenario=NULL,
														 min.date = "20151101",
														 max.date = "20180101",
														 model_name = "spatial",
														 model_id = NULL,
														 test = F,
														 
														 chain_no = 1,
														 burnin = 500,
														 thin = 1,
														 thin2 = 20,
														 init_iter = 10000,
														 iter_per_chunk = 5000,
														 min_eff_size_perchain = 5,
														 max_loops=100,
														 max_save_size = 40000,
														 ...) {
	if (is.null(model_id)){
		model_id <- paste(model_name, s, min.date, max.date, sep = "_")
	}
	
	# Create savepaths
	out.path <- here("data",paste0("model_outputs/clr/",
																 model_name,"/samples_", model_id, ".rds"))
	
	# Create directory structure
	dir.create(dirname(out.path), recursive = TRUE, showWarnings = FALSE)
	
	if (test==T){
		out.path2 <- gsub(".rds", paste0("_test_chain", chain_no, ".rds"), out.path)
	} else {
		out.path2 <- gsub(".rds", paste0("_chain", chain_no, ".rds"), out.path)
	}
	
	
	CLR_modelCode_env_cycl <-  nimble::nimbleCode({
		# Loop through core observations ----
		
		# Observation model (cores ~ plot means)
		for (i in 1:N.core) {
			y[i, 1] ~ dnorm(plot_mu[plot_num[i], timepoint[i]],  core_sd)
		}
		
		# Plot-level process model ----
		for (p in 1:N.plot) {
			
			for (t in plot_start[p]) {
				Ex[p, t] ~ dnorm(0,sd=2)
				# Add process error (sigma)
				plot_mu[p, t] ~ dnorm(Ex[p, t], sd = sigma)
			}
			
			for (t in plot_index[p]:N.date) {
				# Dynamic linear model
				Ex[p, t] <- rho * plot_mu[p, t - 1] +
					site_effect[plot_site_num[p]] +
					beta[1] * temp_est[plot_site_num[p], t] +
					beta[2] * mois_est[plot_site_num[p], t] +
					beta[3] * pH_est[p, plot_start[p]] +
					beta[4] * pC_est[p, plot_start[p]] +
					beta[5] * relEM[p, t] +
					beta[6] * LAI[plot_site_num[p], t] +
					beta[7] * sin_mo[t] +
					beta[8] * cos_mo[t] +
					intercept
				plot_mu[p, t] ~ dnorm(Ex[p, t], sd = sigma)
			}
		}
		# Priors for site effects----
		for (k in 1:N.site) {
			site_effect[k] ~ dnorm(0,  sd = sig)
		}
		# Priors for everything else ----
		
		
		sig ~ dgamma(.5, 1)
		sigma ~ dgamma(.5, 1)
		core_sd ~ dgamma(.5, 1)
		
		intercept ~ dnorm(0, sd = 1)
		rho ~ dnorm(0, sd = 1) # Autocorrelation
		
		beta[1:8] ~ dmnorm(zeros[1:8], omega[1:8, 1:8])
		
		
		
		# Add driver uncertainty if desired ----
		if (temporalDriverUncertainty) {
			for (k in 1:N.site) {
				for (t in site_start[k]:N.date) {
					mois_est[k, t] ~ dnorm(mois[k, t], sd = mois_sd[k, t])
					temp_est[k, t] ~ dnorm(temp[k, t], sd = temp_sd[k, t])
				}
			}
		} else {
			for (k in 1:N.site) {
				for (t in site_start[k]:N.date) {
					mois_est[k, t] <- mois[k, t]
					temp_est[k, t] <- temp[k, t]
				}
			}
		}
		
		# Using 40th time point (values are constant over time)
		if (spatialDriverUncertainty) {
			for (p in 1:N.plot) {
				pH_est[p, plot_start[p]] ~ dnorm(pH[p, plot_start[p]], sd = pH_sd[p, plot_start[p]])
				pC_est[p, plot_start[p]] ~ dnorm(pC[p, plot_start[p]], sd = pC_sd[p, plot_start[p]])
			}
		} else {
			for (p in 1:N.plot) {
				pH_est[p, plot_start[p]] <- pH[p, plot_start[p]]
				pC_est[p, plot_start[p]] <- pC[p, plot_start[p]]
			}
		}
		
	}) #end NIMBLE model.
	
	CLR_modelCode_cycl_only <- nimble::nimbleCode({
		
		
		# Observation model (cores ~ plot means)
		for (i in 1:N.core) {
			y[i, 1] ~ dnorm(plot_mu[plot_num[i], timepoint[i]],  core_sd)
		}
		
		# Plot-level process model ----
		for (p in 1:N.plot) {
			
			for (t in plot_start[p]) {
				Ex[p, t] ~ dnorm(0,sd=1)
				# Add process error (sigma)
				plot_mu[p, t] ~ dnorm(Ex[p, t], sd = sigma)
			}
			
			for (t in plot_index[p]:N.date) {

				# Dynamic linear model
				Ex[p, t] <- rho * plot_mu[p, t - 1] +
					beta[1] * sin_mo[t] +
					beta[2] * cos_mo[t] +
					site_effect[plot_site_num[p]] +
					intercept
				
				# Add process error (sigma)
				plot_mu[p, t] ~ dnorm(Ex[p, t], sd = sigma)
			}
		}
		
		
		# Priors for site effects----
		for (k in 1:N.site) {
			site_effect[k] ~ dnorm(0,  sd = sig)
		}
		
		
		sig ~ dgamma(.5, 1)
		sigma ~ dgamma(.5, 1)
		core_sd ~ dgamma(.5, 1)
		
		
		intercept ~ dnorm(0, sd = 1)
		rho ~ dnorm(0, sd = 1) # Autocorrelation
		
		beta[1:2] ~ dmnorm(zeros[1:2], omega[1:2, 1:2])
		
	}) #end NIMBLE model.
	
	
	nimbleMod_clr_cycl_only <- nimble::nimbleCode({
		# Observation model (cores ~ plot means)
		for (i in 1:N.core) {
			y[i, 1] ~ dnorm(plot_mu[plot_num[i], timepoint[i]],  core_sd)
		}
		
		# Process model
		for (p in 1:N.plot) {
			for (t in plot_start[p]) {
				plot_mu[p, t] ~ dgamma(2, 1) # Plot means for first date
			}
			
			for (t in plot_index[p]:N.date) {
				# Starts from second date
				# Previous value * rho + covariates
				Ex[p, t] <- rho * plot_mu[p, t - 1] +
					beta[1] * sin_mo[t] + beta[2] * cos_mo[t] +
					site_effect[plot_site_num[p]] +
					intercept
				# Add process error, sigma
				plot_mu[p, t] ~ dnorm(Ex[p, t], sigma)
			}
		}
		
		rho ~ dnorm(0, sd = 1)
		core_sd ~ dgamma(.1, 1)
		sigma ~ dgamma(.5, .1)
		intercept ~ dnorm(0, sd = 1)
		
		# Priors for covariates:
		for (n in 1:2) {
			beta[n] ~ dnorm(0, sd = 1)
		}
		# Priors for site effects----
		for (k in 1:N.site) {
			site_effect[k] ~ dnorm(0,  sig)
		}
		# Priors for site effect variance ----
		sig ~ dgamma(.5, 1)
		
	}) #end NIMBLE model.
	
	
	# Read in microbial abundances
	bacteria <- readRDS(here("data", "clean/groupAbundances_16S_2023.rds"))
	fungi <- readRDS(here("data", "clean/groupAbundances_ITS_2023.rds"))
	all_ranks = c(bacteria, fungi)
	
	# Subset to one rank.
	rank.df <- all_ranks[[rank.name]]
	message("Preparing model data for ", rank.name)
	
	
	# Reduce size if testing
	if (test == T) {
		rank.df <- rank.df %>% arrange(siteID, plotID, dateID)
		rank.df = rank.df[1:500,]
	}
	
	# Subset to one group/species, only keep essential columns, create "other" column
	rank.df_spec <- rank.df %>%
		select("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date", !!s) %>%
		mutate(other = 1-rank.df[[s]])
	
	
	# Prep model inputs/outputs.
	model.dat <- prepCLRData(rank.df = rank.df,
															 min.prev = 3,
															  min.date = min.date,
															  max.date = max.date,
													 s=s)
	
	
	
	message("Completed model data for ", rank.name, ", group: ", s)
	
	constants <- model.dat[c("plotID",  "timepoint","plot_site",
													 "site_start", "plot_start", "plot_index",
													 "plot_num", "plot_site_num",
													 "N.plot", "N.spp", "N.core", "N.site", "N.date",
													 "mois", "mois_sd", "temp", "temp_sd", "pH", "pH_sd",
													 "pC", "pC_sd", "LAI", "relEM", "sin_mo", "cos_mo")]
	
	if (model_name == "cycl_only"){
		constants$N.beta = 2
		Nimble_model = CLR_modelCode_cycl_only
	} else if(model_name == "env_cov"){
		constants$N.beta = 6
		Nimble_model = CLR_modelCode_env_covariates
	} else if(model_name == "env_cycl"){
		constants$N.beta = 8
		constants$temporalDriverUncertainty=TRUE
		constants$spatialDriverUncertainty=TRUE
		Nimble_model = CLR_modelCode_env_cycl
	} else message("Missing specification of Nimble model.")
	
	# Some model hyperparameters
	constants$omega <- 0.0001 * diag(constants$N.beta)
	constants$zeros = rep(0, constants$N.beta)
	
	# This date scaling should be done within the prepBetaRegData function instead
	constants$sin_mo = scale(constants$sin_mo, center = F) %>% as.numeric()
	constants$cos_mo = scale(constants$cos_mo, center = F) %>% as.numeric()
	
	
	# for output.
	metadata <- list("rank.name" = rank.name,
									 "niter" = init_iter,
									 "nburnin" = burnin,
									 "thin" = thin,
									 "model_data" = model.dat$truth.plot.long,
									 "sample_values" = model.dat$sample_values)
	
	inits <- createInits(constants)
	inits$y = inits$y[,1,drop=F]
	
	# mcmc.out <- nimbleMCMC(code = Nimble_model,
	# 											 constants = constants,
	# 											 nchains = 3, niter = 100000,
	# 											 nburnin = 50000,
	# 											 #nburnin = 25000,
	# 											 thin = 10,
	# 											 summary = TRUE,
	# 											 monitors = c("beta","sigma","site_effect",
	# 											 						 "sig","intercept",
	# 											 						 "rho","core_sd","plot_mu"),
	# 											 data = list(y=model.dat$y),
	# 											 inits = inits,
	# 											 samplesAsCodaMCMC = T)
	
	
	#plot(mcmc.out$samples)
	
	## Configure & compile model
	Rmodel <- nimbleModel(code = Nimble_model,
												constants = constants,
												data = list(y=model.dat$y),
												inits = inits)
	#Rmodel$Ex <- inits$Ex
	#Rmodel$site_effect <- rep(0, constants$N.site)
	#Rmodel$y <- inits$y
	# Compile model
	cModel <- compileNimble(Rmodel)
	nimbleOptions(multivariateNodesAsScalars = TRUE)
	# Configure & compile MCMC
	mcmcConf <- configureMCMC(cModel, monitors = c("beta","sigma","site_effect",
																								 "sig","intercept",
																								 "rho","core_sd"),
														monitors2 = c("plot_mu"), thin2 = thin2,
														useConjugacy = F)
	
	myMCMC <- buildMCMC(mcmcConf)
	compiled <- compileNimble(myMCMC, project = Rmodel, resetFunctions = TRUE)
	
	
	compiled$run(niter=init_iter, thin=thin, thin2 = thin2, nburnin = burnin)
	cat(paste0("\nInitial run finished for chain", chain_no))
	
	# Sample from MCMC
	
	out.run<-as.mcmc(as.matrix(compiled$mvSamples))
	out.run2<-as.mcmc(as.matrix(compiled$mvSamples2))
	out.run2 <- na.omit(out.run2)
	# Remove timepoints that had no sampling at all (only zeros)
	out.run2 <- mcmc(out.run2[, which(colSums(out.run2) != 0)])
	
	saveRDS(list(samples = out.run, samples2 = out.run2, metadata = metadata),
					out.path2)
	
	continue <- check_continue(out.run, min_eff_size = min_eff_size_perchain)
	loop_counter = 0
	while (continue){
		message("\nEffective sample size too low; running for another ", iter_per_chunk, " iterations\n")
		if (loop_counter < max_loops) {
			
			compiled$run(niter=iter_per_chunk, thin=thin, reset=F)
			
			# Shorten if more than 10k samples have accumulated
			out.run = window_chain(compiled$mvSamples, max_size = max_save_size)
			out.run2 <- window_chain(compiled$mvSamples2, max_size = max_save_size)
			out.run2 <- na.omit(out.run2)
			# Remove timepoints that had no sampling at all (only zeros)
			out.run2 <- mcmc(out.run2[, which(colSums(out.run2) != 0)])
			
			metadata$niteration = iter_per_chunk * loop_counter
			
			saveRDS(list(samples = out.run, samples2 = out.run2, metadata = metadata),
							out.path2)
			continue <- check_continue(out.run, min_eff_size = min_eff_size_perchain)
			loop_counter = loop_counter + 1
		} else {
			message("Exceeded the max number of loops (max_loops). No more sampling!")
			continue <- FALSE
		}
	}
	
	return("ok")
}

