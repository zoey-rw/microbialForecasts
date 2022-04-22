# New source script (started Mar 2021 for Dirichlet regression in Nimble)
# Contains miscellaneous functions and objects to make model scripts clearer.
if (!require("pacman")) install.packages("pacman") 
library(nimble)
library(lubridate)

#### global variables ####
keep_fg_names <- c("cellulolytic", "assim_nitrite_reduction", "dissim_nitrite_reduction", 
									 "assim_nitrate_reduction", "n_fixation", "dissim_nitrate_reduction", 
									 "nitrification", "denitrification", "chitinolytic", "lignolytic", 
									 "copiotroph", "oligotroph", "benomyl_antibiotic", "glucose_simple",  "pyruvate_simple", 
									 "streptomycin_antibiotic", "sucrose_complex", "acetogen_anaerobic", 
									 "chloramphenicol_antibiotic", "erythromycin_antibiotic", 
									 "gentamycin_antibiotic", "glycerol_simple", 
									 "acetate_simple",
									 "acidic_stress", "cellobiose_complex", 
									 "cellulose_complex", "chitin_complex", "galactose_simple", 
									 "xylose_simple", "salt_stress", "herbicide_stress", "osmotic_stress", 
									 "heat_stress", "light_stress", "endophyte", "plant_pathogen", 
									 "animal_pathogen", "ectomycorrhizal", "lichenized", "saprotroph")

fg_names <- c("cellulolytic", "assim_nitrite_reduction", "dissim_nitrite_reduction", 
							"assim_nitrate_reduction", "n_fixation", "dissim_nitrate_reduction", 
							"nitrification", "denitrification", "chitinolytic", "lignolytic", 
							"methanotroph", "copiotroph", "oligotroph", "benomyl_antibiotic", 
							"citrate_simple", "glucose_simple", "glycine_simple", "pyruvate_simple", 
							"streptomycin_antibiotic", "sucrose_complex", "acetogen_anaerobic", 
							"fsfeso4_anaerobic", "ironcitrate_anaerobic", "nonfsfeso4_anaerobic", 
							"potassiumnitrate_anaerobic", "chloramphenicol_antibiotic", "erythromycin_antibiotic", 
							"gentamycin_antibiotic", "nystatin_antibiotic", "glycerol_simple", 
							"acetate_simple", "glutamate_simple", "late_stress", "propionate_simple", 
							"acidic_stress", "alkaline_stress", "d_galacturonicacid_simple", 
							"d_glucuronicacid_simple", "arabinose_simple", "cellobiose_complex", 
							"cellulose_complex", "chitin_complex", "galactose_simple", "glucosamine_simple", 
							"mannose_simple", "n_acetylglucosamine_simple", "pectin_complex", 
							"rhamnose_simple", "trehalose_complex", "xylan_complex", "xylose_simple", 
							"salt_stress", "herbicide_stress", "lowcarbon_stress", "osmotic_stress", 
							"heat_stress", "light_stress", "endophyte", "plant_pathogen", 
							"animal_pathogen", "ectomycorrhizal", "lichenized", "wood_saprotroph", 
							"soil_saprotroph", "litter_saprotroph", "saprotroph")

tax_names <- c("phylum_bac", "class_bac", "order_bac", "family_bac", "genus_bac", 
								"phylum_fun", "class_fun", "order_fun", "family_fun", "genus_fun")

div_scenarios <- c("no_uncertainty_ITS", "spatial_uncertainty_ITS", "temporal_uncertainty_ITS", 
									 "full_uncertainty_ITS", "no_uncertainty_16S", "spatial_uncertainty_16S", 
									 "temporal_uncertainty_16S", "full_uncertainty_16S")

#####


#### misc functions ####
assign_fg_categories <- function(vector) {
	out <- rep(NA, length(vector))
	out[which(grepl("simple", vector))] <- "Simple substrates"
	out[which(grepl("complex|lign|cellu|chitin", vector))] <- "Complex substrates"
	out[which(grepl("stress", vector))] <- "Stresses"
	out[which(grepl("antibiotic", vector))] <- "Antibiotic resistance"
	out[which(grepl("anaerobic", vector))] <- "Anaerobic"
	out[which(grepl("nitr|fixa", vector))] <- "N-cycling"
	out[which(grepl("sapr|path|arbusc|ecto|endo|lichen", vector))] <- "Trophic guild"
	out[which(grepl("copio|oligo", vector))] <- "Life-history"
	return(out)
}

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

fixDate <- function(datesToFix){
	as.Date(paste0(datesToFix, "01"), format = "%Y%m%d")
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


initsFun <- function(constants, type = NULL){
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
	mois_est[is.na(mois_est)] <- 0
	temp_est <- constants$temp
	temp_est[is.na(temp_est)] <- 0
	
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
	} else if (type == "fg") { # for functional groups
		out <- list(
			y = y_init,
			plot_mu = matrix(rep(.1, constants$N.plot*constants$N.date), 
											 constants$N.plot, constants$N.date),
			intercept = 0,
			core_sd = .1,
			sig = sig_init,
			beta = beta_init[1,],
			rho = rho_init[1],
			sigma = .1,
			plot_rel = plot_rel_init,
			site_effect = rep(.1, constants$N.site),
			Ex = matrix(rep(.1, constants$N.plot*constants$N.date), 
									constants$N.plot, constants$N.date),
			#SIGMA = SIGMA,
			mois_est = mois_est,
			temp_est = temp_est,
			pH_est = pH_est,
			pC_est = pC_est
		)
	} else { # For taxa 
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

#####

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
			Ex[p,t] <- rho * plot_mu[p,t-1] + 
				beta[1]*temp_est[plot_site_num[p],t] + 
				beta[2]*mois_est[plot_site_num[p],t] + 
				beta[3]*pH_est[p,1] + 
				beta[4]*pC_est[p,1] +
				beta[5]*relEM[p,t] +
				beta[6]*LAI[plot_site_num[p],t] +
				beta[7]*sin_mo[t] + beta[8]*cos_mo[t] +
				site_effect[plot_site_num[p]] +
				intercept
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
 intercept ~ dnorm(0, sd = 1)
	
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




nimbleModTaxa <- nimbleCode({ 
	
	# Loop through core observations ----
	for(i in 1:N.core){
		y[i,1:N.spp] ~ ddirch(plot_mu[plot_num[i], 1:N.spp, timepoint[i]])
	}
	
	# Plot-level process model ----
	for(s in 1:N.spp){
		for(p in 1:N.plot){
			for (t in plot_start[p]) {
				plot_mu[p,s,t] ~ dgamma(0.5, 1) # Plot means for first date
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
					beta[s,5]*relEM[p,t] +
					beta[s,6]*LAI[plot_site_num[p],t] +
					beta[s,7]*sin_mo[t] + beta[s,8]*cos_mo[t] +
					site_effect[plot_site_num[p],s] +
					intercept[s]
				# Add process error (sigma)
				plot_mu[p,s,t] ~ dnorm(Ex[p,s,t], sigma[s])
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
	sig ~ dgamma(3,1)
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
	
	
	# # Priors for everything else ----
	# for (s in 1:N.spp){
	# 	rho[s] ~ dnorm(0, sd = 1)
	# 	sigma[s] ~ dgamma(.1, .1)
	# 	intercept[s] ~ dnorm(0, sd = 1)
	# 	for (n in 1:N.beta){
	# 		beta[s,n] ~ dnorm(0, sd = 1)
	# 	}
	# }
	
	# Priors for everything else ----
	for (s in 1:N.spp){
		sigma[s] ~ dgamma(.1, .1)
		beta[s,1:N.beta] ~ dmnorm(zeros[1:N.beta], omega[1:N.beta, 1:N.beta])
	}
	rho[1:N.spp] ~ dmnorm(zeros[1:N.spp], omega[1:N.spp, 1:N.spp])
	intercept[1:N.spp] ~ dmnorm(zeros[1:N.spp], omega[1:N.spp, 1:N.spp])
	
	
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
			plot_mu[p,t] ~ dbeta(mean=.3, sd=.3) # Plot means for first date
		}
		
		for (t in plot_index[p]:N.date) {
			# Previous value * rho
			logit(Ex[p,t]) <- rho * logit(plot_mu[p,t-1]) + 
				beta[1]*temp_est[plot_site_num[p],t] +
				beta[2]*mois_est[plot_site_num[p],t] +
				beta[3]*pH_est[p,1] +
				beta[4]*pC_est[p,1] +
				#beta[5]*nspp[p,t] +
				beta[5]*relEM[p,t] +
				beta[6]*LAI[plot_site_num[p],t] +
				beta[7]*sin_mo[t] + beta[8]*cos_mo[t] +
				site_effect[plot_site_num[p]] +
			intercept
			# Add process error (sigma)
			#	plot_mu[p,t] ~ dnorm(Ex[p,t], sigma)
			plot_mu[p,t] ~ dnorm(mean = Ex[p,t], sigma)
			#plot_mu[p,t] <- Ex[p,t]
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
	core_sd ~ dgamma(.5,1)
	rho ~ dnorm(0, sd = 1)
	sigma ~ dgamma(.5, .1)
	#intercept ~ dgamma(.1, .1)
	intercept ~ dnorm(0, sd = 1)
	
	for (n in 1:N.beta){
		beta[n] ~ dnorm(0, sd = 1)
	}
}) #end NIMBLE model.




nimbleModFunctional_trunc <- nimbleCode({ 
	
	# Loop through core observations ----
	for(i in 1:N.core){
		y[i,1] ~ dbeta(mean = plot_mu[plot_num[i],timepoint[i]], 
									 sd = core_sd)
	}
	
	# Plot-level process model ----
	for(p in 1:N.plot){
		for (t in plot_start[p]) {
			plot_mu[p,t] ~ dbeta(mean=.3, sd=.3) # Plot means for first date
		}
		
		for (t in plot_index[p]:N.date) {
			# Previous value * rho
			logit(Ex[p,t]) <- rho * logit(plot_mu[p,t-1]) + 
				beta[1]*temp_est[plot_site_num[p],t] +
				beta[2]*mois_est[plot_site_num[p],t] +
				beta[3]*pH_est[p,1] +
				beta[4]*pC_est[p,1] +
				beta[5]*relEM[p,t] +
				beta[6]*LAI[plot_site_num[p],t] +
				beta[7]*sin_mo[t] + beta[8]*cos_mo[t] +
				site_effect[plot_site_num[p]] +
				intercept
			# Add process error (sigma)
			#	plot_mu[p,t] ~ dnorm(Ex[p,t], sigma)
			plot_mu[p,t] ~ T(dnorm(mean = Ex[p,t], sigma), 0, Inf)
			#plot_mu[p,t] <- Ex[p,t]
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
	core_sd ~ dgamma(.5,1)
	rho ~ dnorm(0, sd = 1)
	sigma ~ dgamma(.5, .1)
	intercept ~ dnorm(0, sd = 1)
	
	for (n in 1:N.beta){
		beta[n] ~ dnorm(0, sd = 1)
	}
}) #end NIMBLE model.





# object <- read_in$samples$samples2
# quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975)

fast.summary.mcmc <- function (object, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), 
															 ...) {
	require(matrixStats)
	require(data.table)
	setDTthreads(threads = 8)
	
	x <- mcmc.list(object)
	statnames <- c("Mean", "SD", "Naive SE", "Time-series SE")
	varstats <- matrix(nrow = nvar(x), ncol = length(statnames), 
										 dimnames = list(varnames(x), statnames))
	xtsvar <- matrix(nrow = nchain(x), ncol = nvar(x))
	
	if (is.matrix(x[[1]])) {
		for (i in 1:nchain(x)) {
			print(paste0("Summarizing chain ", i))
			pb <- txtProgressBar(min = 0, max = nvar(x), style = 3)
			for (j in 1:nvar(x)) {
				setTxtProgressBar(pb, j)
				if (all(na.omit(x[[i]][, j])==0)) xtsvar[i,j] <- 0; next()
				xtsvar[i,j] <- fast.spectrum0.ar(x[[i]][, j])$spec
			}}
		cat("\nCombining MCMC chains...\n")
		#xlong <- as.matrix(data.table::rbindlist(lapply(x,as.data.frame)))
		xlong <- as.matrix(x)
	} else {
		for (i in 1:nchain(x)) {
			xtsvar[i, ] <- fast.spectrum0.ar(x[[i]])$spec
		}
		xlong <- as.matrix(x)
	}
	rm(object)
	cat("\nWrapping up output...\n")
	xmean <- matrixStats::colMeans2(xlong)
	xvar <- matrixStats::colVars(xlong)
	xtsvar <- matrixStats::colMeans2(xtsvar)
	varquant <- matrixStats::colQuantiles(xlong, probs = quantiles)
	
	varstats[, 1] <- xmean
	varstats[, 2] <- sqrt(xvar)
	varstats[, 3] <- sqrt(xvar/(niter(x) * nchain(x)))
	varstats[, 4] <- sqrt(xtsvar/(niter(x) * nchain(x)))
	varquant <- drop(varquant)
	varstats <- drop(varstats)
	out <- list(statistics = varstats, quantiles = varquant, 
							start = start(x), end = end(x), thin = thin(x), nchain = nchain(x))
	class(out) <- "summary.mcmc"
	return(out)
}
	
fast.spectrum0.ar <- function (x) {
	x <- as.matrix(x)
	v0 <- order <- numeric(ncol(x))
	names(v0) <- names(order) <- colnames(x)
	z <- 1:nrow(x)
	for (i in 1:ncol(x)) {
			ar.out <- ar(na.omit(x[, i]))
			v0[i] <- ar.out$var.pred/(1 - sum(ar.out$ar))^2
			order[i] <- ar.out$order
	}
	return(list(spec = v0, order = order))
}



# Calculate summary and save output.
# system.time(v1 <- coda:::spectrum0.ar(y))
# # system.time(v2 <- fast.spectrum0.ar(y))
# 
# parameter_names <- varnames(mcmc_list)
# saved_steps <- as.integer(row.names(mcmc_list[[1]]))
# out <- data.frame("chain" = factor(rep(1 : length(mcmc_list), each = length(saved_steps))),
# 									"step" = rep(saved_steps, length(mcmc_list)) )
# out <- cbind(out, as.data.frame(as.matrix(chain_samples)))

create_div_constants <- function(model.dat){
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
									relEM = model.dat[["relEM"]],
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
return(constants)
}






# nimbleModTaxa <- nimbleCode({ 
# 	
# 	# Loop through core observations ----
# 	for(i in 1:N.core){
# 		y[i,1:N.spp] ~ ddirch(plot_mu[plot_num[i], 1:N.spp, timepoint[i]])
# 	}
# 	
# 	# Plot-level process model ----
# 	for(s in 1:N.spp){
# 		for(p in 1:N.plot){
# 			for (t in plot_start[p]) {
# 				plot_mu[p,s,t] ~ dgamma(0.5, 1) # Plot means for first date
# 				# Convert back to relative abundance
# 				plot_rel[p,s,t] <- plot_mu[p,s,t] / sum(plot_mu[p,1:N.spp,t])
# 			}
# 			
# 			for (t in plot_index[p]:N.date) {
# 				# Previous value * rho
# 				log(Ex[p,s,t]) <- rho[s] * log(plot_mu[p,s,t-1]) + 
# 					beta[s,1]*temp_est[plot_site_num[p],t] +
# 					beta[s,2]*mois_est[plot_site_num[p],t] +
# 					beta[s,3]*pH_est[p,1] +
# 					beta[s,4]*pC_est[p,1] +
# 					beta[s,5]*relEM[p,t] +
# 					beta[s,6]*LAI[plot_site_num[p],t] +
# 					beta[s,7]*sin_mo[t] + beta[s,8]*cos_mo[t] +
# 					site_effect[plot_site_num[p],s] +
# 					intercept[s]
# 				# Add process error (sigma)
# 				plot_mu[p,s,t] ~ dnorm(Ex[p,s,t], sigma[s])
# 				# Convert back to relative abundance
# 				plot_rel[p,s,t] <- plot_mu[p,s,t] / sum(plot_mu[p,1:N.spp,t])
# 			}
# 		}
# 	}
# 	
# 	# Add driver uncertainty if desired ----
# 	if(temporalDriverUncertainty) {
# 		for(k in 1:N.site){
# 			for (t in site_start[k]:N.date) {
# 				mois_est[k,t] ~ dnorm(mois[k,t], sd = mois_sd[k,t])
# 				temp_est[k,t] ~ dnorm(temp[k,t], sd = temp_sd[k,t])
# 			}
# 		}
# 	} else {
# 		for(k in 1:N.site){
# 			for (t in site_start[k]:N.date) {
# 				mois_est[k,t] <- mois[k,t]
# 				temp_est[k,t] <- temp[k,t]
# 			}
# 		} 
# 	}
# 	
# 	# Using 40th time point (values are constant over time)
# 	if(spatialDriverUncertainty) {
# 		for(p in 1:N.plot){
# 			pH_est[p,1] ~ dnorm(pH[p,1], sd = pH_sd[p,1])
# 			pC_est[p,1] ~ dnorm(pC[p,1], sd = pC_sd[p,1])
# 		}
# 	} else {
# 		for(p in 1:N.plot){
# 			pH_est[p,1] <- pH[p,1]
# 			pC_est[p,1] <- pC[p,1]
# 		} 
# 	}
# 	
# 	# Priors for site effect covariance matrix ----
# 	sig ~ dgamma(3,1)
# 	
# 	# Priors for site random effects:
# 	for(s in 1:N.spp){
# 		for(k in 1:N.site){
# 			site_effect[k,s] ~ dnorm(0, sig)
# 		}
# 	}
# 	
# 	
# 	# Priors for everything else ----
# 	for (s in 1:N.spp){
# 		rho[s] ~ dnorm(0, sd = 1)
# 		sigma[s] ~ dgamma(.1, .1)
# 		intercept[s] ~ dnorm(0, sd = 1)
# 		for (n in 1:N.beta){
# 			beta[s,n] ~ dnorm(0, sd = 1)
# 		}
# 	}
# 	
# 	
# }) #end NIMBLE model.



 
check_continue <- function(run1, min_eff_size = 50) {
	require(coda)
	
	cat(paste0("\n Current size: ", nrow(run1)))
	
	effsize <- effectiveSize(run1)
	
	# Get lowest non-zero effective sample size
	lowest_eff_size <- min(effsize[effsize != 0])
	
	# If lower than our preset, continue sampling
	if(lowest_eff_size < min_eff_size){
		cat("\n Effective samples sizes too low:", min(lowest_eff_size))
		return(TRUE)
		#	}
	} else {
		cat("\n Effective samples sizes is sufficient:", min(lowest_eff_size))
		return(FALSE)
	}
}

# 
# combine_chains <- function(chain_paths, save = FALSE){
# 	
# 	# initialize
# 	samples <- samples2 <- metadata <- list()
# 	for(i in 1:length(chain_paths)){
# 		print(i)
# 		# paste model file path to chain number
# 		chain <- readRDS(chain_paths[[i]])
# 		samples[[i]] <- chain[[1]]
# 		samples2[[i]] <- chain[[2]]
# 		metadata[[i]] <- chain[[3]]
# 	}
# 	
# 	nrows <- lapply(samples, nrow) %>% unlist()
# 	min_nrow <- min(nrows)
# 	
# 	for(i in 1:length(chain_paths)){
# 		current_nrow <- nrow(samples[[i]])
# 		if (min_nrow < current_nrow){
# 			print(i)
# 			samples[[i]] <- as.mcmc(samples[[i]][(current_nrow-min_nrow+1):current_nrow,])
# 		}
# 	}
# 	nrows <- lapply(samples2, nrow) %>% unlist()
# 	min_nrow <- min(nrows)
# 	for(i in 1:length(chain_paths)){
# 		current_nrow <- nrow(samples2[[i]])
# 		if (min_nrow < current_nrow){
# 			print(i)
# 			samples2[[i]] <- as.mcmc(samples2[[i]][(current_nrow-min_nrow+1):current_nrow,])
# 		}
# 	}
# 	
# 	out <- list(samples = as.mcmc.list(samples),
# 							samples2 = as.mcmc.list(samples2),
# 							metadata = metadata[[1]])
# 	
# 	if(!isFALSE(save)){
# 		saveRDS(out, file = save)
# 	}
# 	return(out)
# }





combine_chains <- function(chain_paths, 
													 save = FALSE, 
													 cut_size1 = NULL, 
													 cut_size2 = NULL){
	require(coda)
	require(tidyverse)
	
	
	if (is.null(cut_size1)) cut_size1 <- 19999 
	if (is.null(cut_size2)) cut_size2 <- 9999 
	
	
	# initialize
	samples <- samples2 <- metadata <- list()
	for(i in 1:length(chain_paths)){
		
		print(i)
		# paste model file path to chain number
		chain <- readRDS(chain_paths[[i]])
		nrow_samples <- nrow(chain[[1]])
		nrow_samples2 <- nrow(chain[[2]])
		
		if (nrow_samples < cut_size1) {
			cut_size1 <- nrow_samples
		}
		
		samples[[i]] <- as.mcmc(as.matrix(window(chain[[1]], nrow_samples-cut_size1, nrow_samples, 1)))
		
		
		if (nrow_samples2 < cut_size2) {
			cut_size2 <- nrow_samples2
		}
		samples2[[i]] <- as.mcmc(as.matrix(window(chain[[2]], nrow_samples2-cut_size1, nrow_samples2, 1)))
		
		metadata[[i]] <- chain[[3]]
	}
	
	
	nrows <- lapply(samples, nrow) %>% unlist()
	min_nrow <- min(nrows)
	for(i in 1:length(chain_paths)){
		current_nrow <- nrow(samples[[i]])
		if (min_nrow < current_nrow){
			print(i)
			samples[[i]] <- as.mcmc(as.matrix(window(samples[[i]], (current_nrow-min_nrow+1), current_nrow, 1)))
		}
	}
	
	
	nrows <- lapply(samples2, nrow) %>% unlist()
	min_nrow <- min(nrows)
	for(i in 1:length(chain_paths)){
		current_nrow <- nrow(samples2[[i]])
		if (min_nrow < current_nrow){
			print(i)
			samples2[[i]] <- as.mcmc(as.matrix(window(samples2[[i]], (current_nrow-min_nrow+1), current_nrow, 1)))
		}
	}
	
	
	out <- list(samples = as.mcmc.list(samples),
							samples2 = as.mcmc.list(samples2),
							metadata = metadata[[1]])
	
	if(!isFALSE(save)){
		saveRDS(out, file = save)
	}
	return(out)
}


norm_sample <- function(x, y) Rfast::Rnorm(1, x, y)
