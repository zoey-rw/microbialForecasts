# Helper functions and global variables for soil microbial forecasts

message("Loading helper functions")


pacman::p_load(nimble, coda, lubridate, tidyverse, ggpubr, reshape2) 



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


all_covariates_key <- c("1" = "Temperature",
												"2" = "Moisture",
												"3" = "pH",
												"4" = "pC",
												"5" = "Ectomycorrhizal trees",
												"6" = "LAI",
												"7" = "sin",
												"8" = "cos",
												"NA" = "NA")

cycl_only_key <- list("1" = "sin",
											"2" = "cos")

# Combine MCMC chains using paths, shortening each chain due to RAM constraints
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


# Determine whether MCMC should continue running, based on the minimum effective sample size
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


# Version of coda::summary.mcmc(), except using datatable for speed
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

# Needed for faster summary function above
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


# Create input constants for diversity models
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


# Create initial values for MCMC runs
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


# Take a vector of dates in dateID format ("201401") and return as date with first of month
fixDate <- function(datesToFix){
  as.Date(paste0(datesToFix, "01"), format = "%Y%m%d")
}


# Take a vector of values and put them in the (0,1) interval
interval_transform <- function(x,C = ncol(x), N = nrow(x)){
	out <- (x * (N - 1) + (1/C)) / N
	return(out)
}




gel <- function(samples){
  if ("WAIC" %in% names(samples)) {
    samples <- samples$samples
  }
  allvars <- c("beta[1]","beta[2]","beta[3]","beta[4]","beta[5]","alpha","plot_var","time_var","site_var","tau_proc","tau_obs")
  #time_effs <- grep("time_effect", colnames(samples[[1]]), fixed=T)
  vars <- allvars[allvars %in% colnames(samples[[1]])]
  coda::gelman.diag(samples[,vars])
}

var_list <- c("tau_obs","tau_proc","plot_var","site_var", "plot_effect", "time_var", "beta", "site_effect","plot_mean_hat","alpha")
var_listSimple <- c("tau_obs","tau_proc","plot_var", "plot_effect",#"glob_mean",
                    #"time_var", 
                    "plot_mean",
                    "beta","plot_mean_hat","alpha")


# Remove NAs from MCMC objects
rm.NA.mcmc <- function(samples){
	sample.list  <- mcmc.list(mcmc(na.omit(samples[[1]])), 
														mcmc(na.omit(samples[[2]])), 
														mcmc(na.omit(samples[[3]])))
	return(sample.list)
}


### CONVERT PRECISION TO SD
prec_to_sd <- function(x) 1/sqrt(x)


### CONVERT PRECISION TO SD for all samples 
var_to_sd <- function(samples, var.list = c("tau_proc", "tau_obs"), mean = TRUE){
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


plot_betas <- function(samples_var = samples){
  samples <- samples_var$samples
  beta_names <- colnames(samples[[1]])[grep("beta", colnames(samples[[1]]))]
  par(ask = TRUE)
  plot(samples[,beta_names])
}

traceplots <- function(samples = samples, var = "betas"){
  if ("WAIC" %in% names(samples)) {
    samples <- samples$samples
  }
  if (var == "betas"){
    to_plot <- colnames(samples[[1]])[grep("beta", colnames(samples[[1]]))]
  } else if (var == "plot_effect"){
    to_plot <- colnames(samples[[1]])[grep("plot_effect", colnames(samples[[1]]))]
  } else if (var == "alpha"){
    to_plot <- colnames(samples[[1]])[grep("alpha", colnames(samples[[1]]))]
  } else if (var == "tau"){
    to_plot <- colnames(samples[[1]])[grep("tau", colnames(samples[[1]]))]
  } else if (var == "site_effect"){
    to_plot <- colnames(samples[[1]])[grep("site_effect", colnames(samples[[1]]))]
  } else {
    to_plot <- colnames(samples[[1]])[grep(var, colnames(samples[[1]]), fixed=T)]
  }
  par(ask = TRUE)
  plot(samples[,to_plot])
}


# tic/toc functions from Colin Averill:
#' Two clock functions.
#' Place tic() at the line in the code where you want to start timing.
#' Place toc() at the position in the code where you want to stop timing and report.
tic <- function() {assign("timer", Sys.time(), envir=.GlobalEnv)}
toc <- function() print(Sys.time()-timer)



# Assign functional group categories based on group name
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


# Assign functional group categories based on group name
assign_fg_kingdoms <- function(vector) {
	out <- rep(NA, length(vector))
	out[which(grepl("Simple|Complex", vector))] <- "16S"
	out[which(grepl("Stress|Antibiotic|Anaerobic|cycling", vector))] <- "16S"
	out[which(grepl("Troph", vector))] <- "ITS"
	out[which(grepl("Life", vector))] <- "16S"
	return(out)
}


# create sample information data.frame from NEON sample names
parseNEONsampleIDs <- function(sampleID){
  df <- data.frame(siteID = substr(sampleID, 1, 4), sampleID = sampleID, stringsAsFactors = F) %>% 
  	mutate(sample = sapply(strsplit(sampleID, "-GEN|-gen"),  "[[" , 1)) %>% 
  	mutate(geneticSampleID = sapply(strsplit(sampleID, "-DNA"),  "[[" , 1)) %>% 
  	mutate(sampleID = sapply(strsplit(sampleID, "-gen.fastq"),  "[[" , 1)) %>% 
    mutate(dates = sapply(strsplit(sample, "-"), function(x) x[grep("[2]\\d\\d\\d\\d\\d\\d\\d", x)])) %>% 
    mutate(dates = ifelse(dates == "21040514", "20140514", dates)) %>% 
    mutate(asDate = as.Date(as.character(dates), "%Y%m%d")) %>% 
    mutate(dateID = substr(as.character(dates), 1, 6)) %>% 
    mutate(plotID = substr(sample, 1, 8)) %>% 
    mutate(site_date = paste0(siteID, "-", dateID)) %>% 
    mutate(horizon = ifelse(grepl("-M-", sample), "M", "O")) %>% 
    mutate(without_horizon = gsub("-[M|O]-", "-", sample)) %>% 
    mutate(plot_date = paste0(plotID, "-", dateID)) %>% 
    as.data.frame()
  rownames(df) <- make.unique(sampleID)
  return(df)
}





rbind.named.dfs <- function(df.list){ 
  # solution from https://stackoverflow.com/questions/15162197/combine-rbind-data-frames-and-create-column-with-name-of-original-data-frames
  dfs <- df.list[sapply(df.list, function(x) !is.null(dim(x)))]
  all.out <- cbind.data.frame(do.call(rbind,dfs), 
                              name = rep(names(dfs), vapply(dfs, nrow, numeric(1))))
  return(all.out)
}




base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}




pretty_breaks <- function(n = 5, ...) {
  function(x) {
    breaks <- pretty(x, n, ...)
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
}

rbind.df.list <- function(pl.out){
  lapply(pl.out, function(x){
    do.call(rbind, x)
  })
}


pretty_names_old <- list("cellulolytic" = "Cellulolytic bacteria",
                     "assim_nitrite_reduction" = "Assimilatory nitrite reducing bacteria",
                     "dissim_nitrite_reduction" = "Dissimilatory nitrite reducing bacteria",
                     "assim_nitrate_reduction" = "Assimilatory nitrate reducing bacteria",
                     "n_fixation" = "Nitrogen fixing bacteria",
                     "dissim_nitrate_reduction" = "Dissimilatory nitrate reducing bacteria",
                     "nitrification" = "Nitrifying bacteria",
                     "denitrification" = "Denitrifying bacteria",
                     "chitinolytic" = "Chitinolytic bacteria",
                     "lignolytic" = "Ligninolytic bacteria",
                     "methanotroph" = "Methanotrophic bacteria",
                     "copiotroph" = "Copiotrophic bacteria",
                     "oligotroph" = "Oligotrophic bacteria",
                     "Arbuscular" = "Arbuscular mycorrhizal fungi",
                     "Animal_Pathogen" = "Animal pathogen fungi",
                     "Plant_Pathogen" = "Plant pathogen fungi",
                     "Saprotroph" = "Saprotrophs fungi",
                     "Wood_Saprotroph" = "Wood saprotrophs (fungi)",
                     "Ectomycorrhizal" = "Ectomycorrhizal fungi",
                     "Bdellovibrio" = "Genus Bdellovibrio",
                     "Pseudogymnoascus" = "Genus Pseudogymnoascus",
                     "Cyanobacteria" = "Phylum Cyanobacteria"
)

pretty_rank_names <- list("genus_bac" = "Genus",
													 "family_bac" = "Family",
													 "order_bac" = "Order", 
													 "class_bac" = "Class", 
													 "phylum_bac" = "Phylum",
													 "genus_fun" = "Genus",
													 "family_fun" = "Family",
													 "order_fun" = "Order", 
													 "class_fun" = "Class", 
													 "phylum_fun" = "Phylum",
													 "functional_group" = "Functional group",
													"diversity_16S" = "Diversity",
													"diversity_ITS" = "Diversity"
													)

pretty_names <- list("cellulolytic" = "Cellulose degraders",
                     "assim_nitrite_reduction" = "Assimilatory nitrite reducers",
                     "dissim_nitrite_reduction" = "Dissimilatory nitrite reducers",
                     "assim_nitrate_reduction" = "Assimilatory nitrate reducers",
                     "n_fixation" = "Nitrogen fixers",
                     "dissim_nitrate_reduction" = "Dissimilatory nitrate reducers",
                     "nitrification" = "Nitrifiers",
                     "denitrification" = "Denitrifiers",
                     "chitinolytic" = "Chitin degraders",
                     "lignolytic" = "Lignin degraders",
                     "methanotroph" = "Methanotrophs",
                     "copiotroph" = "Copiotrophs",
                     "oligotroph" = "Oligotrophs",
                     "Arbuscular" = "Arbuscular mycorrhizae",
                     "Animal_Pathogen" = "Animal pathogens",
                     "Plant_Pathogen" = "Plant pathogens",
                     "Saprotroph" = "Saprotrophs",
                     "Wood_Saprotroph" = "Wood saprotrophs",
                     "Ectomycorrhizal" = "Ectomycorrhizae",
                     "Bdellovibrio" = "Genus Bdellovibrio",
                     "Pseudogymnoascus" = "Genus Pseudogymnoascus",
                     "Cyanobacteria" = "Phylum Cyanobacteria"
)


FG_kingdoms <- list("cellulolytic" = "Bacterial_functional_group",
                     "assim_nitrite_reduction" = "Bacterial_functional_group",
                     "dissim_nitrite_reduction" = "Bacterial_functional_group",
                     "assim_nitrate_reduction" = "Bacterial_functional_group",
                     "n_fixation" = "Bacterial_functional_group",
                     "dissim_nitrate_reduction" = "Bacterial_functional_group",
                     "nitrification" = "Bacterial_functional_group",
                     "denitrification" = "Bacterial_functional_group",
                     "chitinolytic" = "Bacterial_functional_group",
                     "lignolytic" = "Bacterial_functional_group",
                     "methanotroph" = "Bacterial_functional_group",
                     "copiotroph" = "Bacterial_functional_group",
                     "oligotroph" = "Bacterial_functional_group",
                     "Arbuscular" = "Fungal_functional_group",
                     "Animal_Pathogen" = "Fungal_functional_group",
                     "Plant_Pathogen" = "Fungal_functional_group",
                     "Saprotroph" = "Fungal_functional_group",
                     "Wood_Saprotroph" = "Fungal_functional_group",
                     "Ectomycorrhizal" = "Fungal_functional_group")


N_cyclers <- c("assim_nitrite_reduction", "dissim_nitrite_reduction","assim_nitrate_reduction","n_fixation","dissim_nitrate_reduction","nitrification","denitrification")

tukey <- function(x, y, extra_info = NULL, y.offset = .3){
	new.df <- cbind.data.frame("x" = x, "y" = y)
	abs_max <- max(new.df[,2])
  maxs <- new.df %>%
    group_by(x) %>%
    summarise(tot=max(y)+ y.offset * abs_max)
  Tukey_test <- aov(y ~ x, data=new.df) %>%
    agricolae::HSD.test("x", group=TRUE) %>%
    .$groups %>%
    as_tibble(rownames="x") %>%
    rename("Letters_Tukey"="groups") %>% 
    dplyr::select(-y) %>%
    left_join(maxs, by="x") 
  if (!is.null(extra_info)){
    Tukey_test <- cbind.data.frame(Tukey_test)
  }
  return(Tukey_test)
}





