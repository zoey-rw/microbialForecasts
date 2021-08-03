# New source script (started Mar 2021 for Dirichlet regression in Nimble)
# Contains miscellaneous functions and objects to make model scripts clearer.
library(nimble)


### CONVERT PRECISION TO SD
prec_to_sd <- function(x) 1/sqrt(x)


### CONVERT PRECISION TO SD for all samples 
tau_to_sd <- function(samples, var.list = c("tau_proc", "tau_obs","sigma"), mean = TRUE){
	if (!is.mcmc(samples) && !is.mcmc.list(samples)) {
		cat("Samples must be MCMC or MCMC.list objects"); stop()
	}
	out.list <- list()
	for (i in 1:length(var.list)){
		varname <- var.list[[i]]
		tau_samples <- samples[,varname]
		samples_sd <- unlist(lapply(tau_samples, function(y) lapply(y, function(x) 1/sqrt(x))))
		samples_mean_sd <- mean(samples_sd)
		out.list[[i]] <- samples_mean_sd
	}
	names(out.list) <- var.list
	# 
	# tau_samples_mean <-   summary(tau_samples)[[1]][1]
	# samples_mean_sd <- 1/sqrt(tau_samples_mean)
	return(out.list)
}

rbind.named.dfs <- function(df.list){ 
	# solution from https://stackoverflow.com/questions/15162197/combine-rbind-data-frames-and-create-column-with-name-of-original-data-frames
	dfs <- df.list[sapply(df.list, function(x) !is.null(dim(x)))]
	all.out <- cbind.data.frame(do.call(rbind,dfs), 
															name = rep(names(dfs), vapply(dfs, nrow, numeric(1))))
	return(all.out)
}


interval_transform <- function(x,C = ncol(x), N = nrow(x)){
	        out <- (x * (N - 1) + (1/C)) / N
     return(out)
}


lmIndicatorCode <- nimbleCode({ 
	
	for(i in 1:n.core.by.date){
		for (t in site_start[core_site[i]]:N.date){
			# Core observations (3 per plot):
			y[i,1:N.spp,t] ~ ddirch(plot_mu[plotID[i],1:N.spp,t])
		}
	}
	#	ypred[1:6,1:6] <- y[,,1]
	
	
	for(s in 1:N.spp){
		for(p in 1:N.plot){
			#for (t in 1:1) {
			Ex[p,s,site_start[plot_site[p]]] ~ dgamma(0.01, 0.01) # Plot means for first date
			# Add process error, sigma
			plot_mu[p,s,site_start[plot_site[p]]] ~ dnorm(Ex[p,s,site_start[plot_site[p]]], 
																										sd = sigma[s])
			# Convert back to relative abundance
			plot_rel[p,s,site_start[plot_site[p]]] <- plot_mu[p,s,site_start[plot_site[p]]] / 
				sum(plot_mu[p,1:N.spp,site_start[plot_site[p]]])
			#}
			for (t in site_index[plot_site[p]]:N.date) {
				# Previous value * rho
				log(Ex[p,s,t]) <- rho[s] * log(plot_mu[p,s,t-1]) + 
					zbeta[s,1]*temp[p,t] + 
					zbeta[s,2]*mois[p,t] + 
					zbeta[s,3]*pH[p,t] + 
					zbeta[s,4]*pC[p,t] +
					zbeta[s,5]*nspp[p,t] + 
					zbeta[s,6]*rc_grass[p,t] +
					#zbeta[s,7]*beta[s,7]*rc_exotic[p,t] +
					site_effect[plot_site[p],s] +
					intercept[s]
				# Add process error, sigma
				plot_mu[p,s,t] ~ dnorm(Ex[p,s,t], sd = sigma[s])
				# Convert back to relative abundance
				plot_rel[p,s,t] <- plot_mu[p,s,t] / sum(plot_mu[p,1:N.spp,t])
			}
		}
	}
	
	psi ~ dunif(0,1)    ## prior on inclusion probability
	for (s in 1:N.spp){
		rho[s] ~ dnorm(0, sd = 3)
		sigma[s] ~ dgamma(.1, .1)
		intercept[s] ~ dgamma(.1, .1)
		for (n in 1:N.beta){
			beta[s,n] ~ dnorm(0, sd = 3)
			z[s,n] ~ dbern(psi) ## indicator variable for each coefficient
			zbeta[s,n] <- z[s,n] * beta[s,n] 
		}
	}
	
	# Priors for site random effects:
	for(s in 1:N.spp){
		for(k in 1:N.site){
			site_effect[k,s] ~ dnorm(0, 3)
		}
	}
	
}) #end NIMBLE model.


initsFun <- function(constants){
	#	core_per_plot <- 3
	y_init <- matrix(rep(rep(1/constants$N.spp, constants$N.spp), constants$N.core),
				 ncol = constants$N.spp, nrow = constants$N.core)
	
	plot_mu_init <- array(rnorm(constants$N.plot  * constants$N.spp * constants$N.date,
															1,.1), 
												dim = c(constants$N.plot, constants$N.spp, constants$N.date))
	plot_rel_init <- array(rep(rep(rep(1/constants$N.spp, constants$N.spp), constants$N.plot), constants$N.date),
												dim = c(constants$N.plot, constants$N.spp, constants$N.date))
	beta_init <- matrix(rep(rep(0, constants$N.beta), constants$N.spp), constants$N.spp, constants$N.beta)
	rho_init <- rep(0, constants$N.spp)
	sigma_init <- rep(1, constants$N.spp)
	intercept_init <- rep(.3, constants$N.spp)
	Ex_init <- plot_mu_init
	mois_est <- constants$mois
	temp_est <- constants$temp
	pH_est <- constants$pH
	pC_est <- constants$pC
	sig_init <- 1
	
	if (constants$N.spp == 1){ # For diversity models 
		plot_mu_init <- matrix(rnorm(constants$N.plot * constants$N.date,
																1,.1), 
													nrow = constants$N.plot, ncol = constants$N.date)
		beta_init <- rep(0, constants$N.beta)
		site_effect_init <- rep(0, constants$N.site)
		Ex_init <- plot_mu_init
		out <- list(
			y = y_init,
			plot_mu = plot_mu_init,
			intercept = intercept_init,
			sig = sig_init,
			beta = beta_init,
			rho = rho_init,
			core_sd = 1,
			sigma = sigma_init,
			site_effect = site_effect_init,
			Ex = plot_mu_init,
			mois_est = mois_est,
			temp_est = temp_est,
			pH_est = pH_est,
			pC_est = pC_est 
		)
	} else { # For taxa/functional groups
		SIGMA <- diag(rep(.1, constants$N.spp))		

		out <- list(
			y = y_init,
			plot_mu = plot_mu_init,
			intercept = intercept_init,
			sig = sig_init,
			beta = beta_init,
			rho = rho_init,
			sigma = sigma_init,
			plot_rel = plot_rel_init,
			site_effect <- diag(rep(sig_init, constants$N.spp)),
			Ex = plot_mu_init,
			SIGMA = SIGMA,
			mois_est = mois_est,
			temp_est = temp_est,
			pH_est = pH_est,
			pC_est = pC_est
		)
	}
	return(out)
}

# Remove NAs from MCMC objects
rm.NA.mcmc <- function(samples){
	sample.list  <- mcmc.list(mcmc(na.omit(samples[[1]])), 
														mcmc(na.omit(samples[[2]])), 
														mcmc(na.omit(samples[[3]])))
	return(sample.list)
}


# 
# #Prep abundance data (just run once, for this testing dataset)
# legacy_ps <- readRDS("/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/data/NEON_16S_phyloseq_legacy.rds")
# tax_rank <- "Phylum"
# ps.filt <- prune_samples(rowSums(otu_table(legacy_ps)) > 5000, legacy_ps)
# ps.rel <- transform_sample_counts(ps.filt, function(x) x/sum(x))
# ps.phy <- tax_glom(ps.rel, tax_rank)
# 
# glom_melt <- speedyseq::psmelt(ps.phy)
# form <- as.formula(paste0("sampleID ~ ", tax_rank))
# glom_wide <- reshape2::dcast(glom_melt, form, value.var = "Abundance", fun.aggregate = sum)
# out_abun <- transform(glom_wide, row.names=sampleID, sampleID=NULL)
# 
# most_abundant_taxa <- names(sort(colSums(out_abun), decreasing = T)[1:5])
# seqDepth <- rowSums(out_abun)
# out_top10 <- out_abun[,colnames(out_abun) %in% most_abundant_taxa]
# out_top10$Other <- 1-rowSums(out_top10)
# rank.df <- cbind(sample_data(ps.phy)[,c("siteID","plotID","dateID","sampleID","dates","plot_date")], out_top10)
# saveRDS(rank.df,"/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/phylum_dat_testing.rds")


nimbleMod_shannon <- nimbleCode({ 
	
	# Observation model (cores ~ plot means)
	for(i in 1:N.core){
		y[i,1] ~ dnorm(plot_mu[plot_num[i],timepoint[i]],  core_sd)
	}
	
	# Process model
	for(p in 1:N.plot){
		for (t in plot_start[p]) {
			plot_mu[p,t] ~ dgamma(2, 1) # Plot means for first date
		}
		
		for (t in plot_index[p]:N.date) { # Starts from second date
			# Previous value * rho + covariates
			log(Ex[p,t]) <- rho * log(plot_mu[p,t-1]) + 
				beta[1]*temp_est[plot_site_num[p],t] + 
				beta[2]*mois_est[plot_site_num[p],t] + 
				beta[3]*pH_est[p,1] + 
				beta[4]*pC_est[p,1] +
				beta[5]*nspp[p,t] +
				beta[6]*rc_grass[p,t] +
				#	beta[7]*sin_mo[t] + beta[8]*cos_mo[t] +
				site_effect[plot_site_num[p]] #+
#				intercept
			# Add process error, sigma
			plot_mu[p,t] ~ dnorm(Ex[p,t], sigma)
		}
	}
	
	# Add driver uncertainty if desired ----
	if(temporalDriverUncertainty) {
		for(k in 1:N.site){
			for (t in site_start[k]:N.date) {
				mois_est[k,t] ~ dnorm(mois[k,t], sd = mois_sd[k,t])
				temp_est[k,t] ~ dnorm(temp[k,t], sd = temp_sd[k,t])
			}
		}
	} else {
		for(k in 1:N.site){
			for (t in site_start[k]:N.date) {
				mois_est[k,t] <- mois[k,t]
				temp_est[k,t] <- temp[k,t]
			}
		} 
	}
	
	# Add spatial uncertainty (values are constant over time)
	if(spatialDriverUncertainty) {
		for(p in 1:N.plot){
			pH_est[p,1] ~ dnorm(pH[p,1], sd = pH_sd[p,1])
			pC_est[p,1] ~ dnorm(pC[p,1], sd = pC_sd[p,1])
		}
	} else {
		for(p in 1:N.plot){
			pH_est[p,1] <- pH[p,1]
			pC_est[p,1] <- pC[p,1]
		}
	}
	
	rho ~ dnorm(0, sd = 1)
	#core_sd ~ dinvgamma(3, .5)
	core_sd ~ dgamma(.1, 1)
	sigma ~ dgamma(.5, .1)
	#	intercept ~ dgamma(.1, .1)
#	intercept ~ dnorm(0, sd = 3)
	
	# Priors for covariates:
	for (n in 1:N.beta){
		beta[n] ~ dnorm(0, sd = 1)
	}
	# Priors for site effects----
	for(k in 1:N.site){
		site_effect[k] ~ dnorm(0,  sig)
	}
	# Priors for site effect variance ----
	sig ~ dgamma(.5,1)
	
}) #end NIMBLE model.


nimbleMod_shannon_nolog <- nimbleCode({ 
	
	# Observation model (cores ~ plot means)
	for(i in 1:N.core){
		y[i,1] ~ dnorm(plot_mu[plot_num[i],timepoint[i]],  core_sd)
	}
	
	# Process model
	for(p in 1:N.plot){
		for (t in plot_start[p]) {
			plot_mu[p,t] ~ dgamma(2, 1) # Plot means for first date
		}
		
		for (t in plot_index[p]:N.date) { # Starts from second date
			# Previous value * rho + covariates
			Ex[p,t] <- rho * plot_mu[p,t-1] + 
				beta[1]*temp_est[plot_site_num[p],t] + 
				beta[2]*mois_est[plot_site_num[p],t] + 
				beta[3]*pH_est[p,1] + 
				beta[4]*pC_est[p,1] +
				beta[5]*nspp[p,t] +
				beta[6]*rc_grass[p,t] +
				#	beta[7]*sin_mo[t] + beta[8]*cos_mo[t] +
				site_effect[plot_site_num[p]] #+
			#				intercept
			# Add process error, sigma
			plot_mu[p,t] ~ dnorm(Ex[p,t], sigma)
		}
	}
	
	# Add driver uncertainty if desired ----
	if(temporalDriverUncertainty) {
		for(k in 1:N.site){
			for (t in site_start[k]:N.date) {
				mois_est[k,t] ~ dnorm(mois[k,t], sd = mois_sd[k,t])
				temp_est[k,t] ~ dnorm(temp[k,t], sd = temp_sd[k,t])
			}
		}
	} else {
		for(k in 1:N.site){
			for (t in site_start[k]:N.date) {
				mois_est[k,t] <- mois[k,t]
				temp_est[k,t] <- temp[k,t]
			}
		} 
	}
	
	# Add spatial uncertainty (values are constant over time)
	if(spatialDriverUncertainty) {
		for(p in 1:N.plot){
			pH_est[p,1] ~ dnorm(pH[p,1], sd = pH_sd[p,1])
			pC_est[p,1] ~ dnorm(pC[p,1], sd = pC_sd[p,1])
		}
	} else {
		for(p in 1:N.plot){
			pH_est[p,1] <- pH[p,1]
			pC_est[p,1] <- pC[p,1]
		}
	}
	
	rho ~ dnorm(0, sd = 1)
	#core_sd ~ dinvgamma(3, .5)
	core_sd ~ dgamma(.1, 1)
	sigma ~ dgamma(.5, .1)
	#	intercept ~ dgamma(.1, .1)
	#	intercept ~ dnorm(0, sd = 3)
	
	# Priors for covariates:
	for (n in 1:N.beta){
		beta[n] ~ dnorm(0, sd = 1)
	}
	# Priors for site effects----
	for(k in 1:N.site){
		site_effect[k] ~ dnorm(0,  sig)
	}
	# Priors for site effect variance ----
	sig ~ dgamma(.5,1)
	
}) #end NIMBLE model.







nimbleModLong <- nimbleCode({ 
	
	# Loop through core observations ----
	for(i in 1:N.core){
		y[i,1:N.spp] ~ ddirch(plot_mu[plot_num[i], 1:N.spp, timepoint[i]])
	}
	
	# Plot-level process model ----
	for(s in 1:N.spp){
		for(p in 1:N.plot){
			for (t in plot_start[p]) {
				plot_mu[p,s,t] ~ dgamma(0.01, 0.01) # Plot means for first date
				# Convert back to relative abundance
				plot_rel[p,s,t] <- plot_mu[p,s,t] / sum(plot_mu[p,1:N.spp,t])
			}
			
			for (t in plot_index[p]:N.date) {
				# Previous value * rho
				log(Ex[p,s,t]) <- rho[s] * log(plot_mu[p,s,t-1]) + 
					beta[s,1]*temp_est[plot_site_num[p],t] +
					beta[s,2]*mois_est[plot_site_num[p],t] +
					beta[s,3]*pH_est[p,1] +
					beta[s,4]*pC_est[p,1] +
					beta[s,5]*nspp[p,t] +
					beta[s,6]*rc_grass[p,t] +
					site_effect[plot_site_num[p],s] +
					intercept[s]
				# Add process error (sigma)
				plot_mu[p,s,t] ~ dnorm(Ex[p,s,t], sd = sigma[s])
				# Convert back to relative abundance
				plot_rel[p,s,t] <- plot_mu[p,s,t] / sum(plot_mu[p,1:N.spp,t])
			}
		}
	}
	
	# Add driver uncertainty if desired ----
	if(temporalDriverUncertainty) {
		for(k in 1:N.site){
			for (t in site_start[k]:N.date) {
				mois_est[k,t] ~ dnorm(mois[k,t], sd = mois_sd[k,t])
				temp_est[k,t] ~ dnorm(temp[k,t], sd = temp_sd[k,t])
			}
		}
	} else {
		for(k in 1:N.site){
			for (t in site_start[k]:N.date) {
				mois_est[k,t] <- mois[k,t]
				temp_est[k,t] <- temp[k,t]
			}
		} 
	}
	
	# Using 40th time point (values are constant over time)
	if(spatialDriverUncertainty) {
		for(p in 1:N.plot){
				pH_est[p,1] ~ dnorm(pH[p,1], sd = pH_sd[p,1])
				pC_est[p,1] ~ dnorm(pC[p,1], sd = pC_sd[p,1])
		}
	} else {
		for(p in 1:N.plot){
				pH_est[p,1] <- pH[p,1]
				pC_est[p,1] <- pC[p,1]
		} 
	}
	
	# Priors for site effect covariance matrix ----
	#sig ~ dgamma(3,1)
	# SIGMA[1:N.spp,1:N.spp] <- diag(rep(sig^2, N.spp))		
	# 
	# # Priors for site random effects:
	# for(k in 1:N.site){
	# 	site_effect[k,1:N.spp] ~ dmnorm(alpha0[1:N.spp], # vector of zeros
	# 																	cov = SIGMA[1:N.spp,1:N.spp])
	# }
	
	# Priors for site random effects:
	for(s in 1:N.spp){
		for(k in 1:N.site){
			site_effect[k,s] ~ dnorm(0, sig)
		}
	}
	
	
	# Priors for everything else ----
	for (s in 1:N.spp){
		rho[s] ~ dnorm(0, sd = 1)
		sigma[s] ~ dgamma(.1, .1)
		intercept[s] ~ dgamma(.1, .1)
		for (n in 1:N.beta){
			beta[s,n] ~ dnorm(0, sd = 1)
		}
	}
	
	
}) #end NIMBLE model.





nimbleModFunctional <- nimbleCode({ 
	
	# Loop through core observations ----
	for(i in 1:N.core){
		y[i,1] ~ dbeta(mean = plot_mu[plot_num[i],timepoint[i]], 
									 sd = core_sd)
	}
	
	# Plot-level process model ----
	for(p in 1:N.plot){
		for (t in plot_start[p]) {
			plot_mu[p,t] ~ dbeta(mean=.55, sd=.1) # Plot means for first date
		}
		
		for (t in plot_index[p]:N.date) {
			# Previous value * rho
			logit(Ex[p,t]) <- rho * logit(plot_mu[p,t-1]) + 
				
				beta[1]*temp_est[plot_site_num[p],t] +
				beta[2]*mois_est[plot_site_num[p],t] +
				beta[3]*pH_est[p,1] +
				beta[4]*pC_est[p,1] +
				beta[5]*nspp[p,t] +
				beta[6]*rc_grass[p,t] +
				site_effect[plot_site_num[p]] #+
			#intercept
			# Add process error (sigma)
			#	plot_mu[p,t] ~ dnorm(Ex[p,t], sigma)
			plot_mu[p,t] ~ dnorm(mean = Ex[p,t], sigma)
		}
	}
	

	# Add driver uncertainty if desired ----
	if(temporalDriverUncertainty) {
		for(k in 1:N.site){
			for (t in site_start[k]:N.date) {
				mois_est[k,t] ~ dnorm(mois[k,t], sd = mois_sd[k,t])
				temp_est[k,t] ~ dnorm(temp[k,t], sd = temp_sd[k,t])
			}
		}
	} else {
		for(k in 1:N.site){
			for (t in site_start[k]:N.date) {
				mois_est[k,t] <- mois[k,t]
				temp_est[k,t] <- temp[k,t]
			}
		} 
	}
	
	# Using 40th time point (values are constant over time)
	if(spatialDriverUncertainty) {
		for(p in 1:N.plot){
			pH_est[p,1] ~ dnorm(pH[p,1], sd = pH_sd[p,1])
			pC_est[p,1] ~ dnorm(pC[p,1], sd = pC_sd[p,1])
		}
	} else {
		for(p in 1:N.plot){
			pH_est[p,1] <- pH[p,1]
			pC_est[p,1] <- pC[p,1]
		} 
	}
	
	
	# Priors for site effect covariance matrix ----
	sig ~ dgamma(.5,1)
	
	# Priors for site effects----
	for(k in 1:N.site){
		site_effect[k] ~ dnorm(0,  sig)
	}
	
	# Priors for everything else ----
	rho ~ dnorm(0, sd = 1)
	sigma ~ dgamma(.5, .1)
	#intercept ~ dgamma(.1, .1)
	for (n in 1:N.beta){
		beta[n] ~ dnorm(0, sd = 1)
	}
}) #end NIMBLE model.
