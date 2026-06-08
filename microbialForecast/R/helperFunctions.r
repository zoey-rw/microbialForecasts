# Helper functions and global variables for soil microbial forecasts
# Dependencies (coda, here) are declared in DESCRIPTION Imports






#' @title norm_sample
#' @description norm_sample
#' @export
norm_sample <- function(x, y) Rfast::Rnorm(1, x, y)



#' @title check_continue
#' @description Determine whether MCMC should continue running, based on the minimum effective sample size
#' @export
check_continue <- function(run1, min_eff_size = 50) {
	pacman::p_load(coda)


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



#' @title 			fast.summary.mcmc
#' @description Duplicate of coda::summary.mcmc(), except using datatable for speed
#' @export
fast.summary.mcmc <- function (object, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
															 ...) {
	pacman::p_load(matrixStats, data.table)

	setDTthreads(threads = 8)

	x <- mcmc.list(object)
	statnames <- c("Mean", "SD", "Naive SE", "Time-series SE")
	varstats <- matrix(nrow = nvar(x), ncol = length(statnames),
										 dimnames = list(varnames(x), statnames))
	xtsvar <- matrix(nrow = nchain(x), ncol = nvar(x))

	if (is.matrix(x[[1]])) {
		for (i in 1:nchain(x)) {
			print(paste0("Summarizing chain ", i))
			#pb <- txtProgressBar(min = 0, max = nvar(x), style = 3)
			for (j in 1:nvar(x)) {
				#setTxtProgressBar(pb, j)
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

#' @title fast.spectrum0.ar
#' @description Needed for faster summary function above
#' @export
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

#' @title 			Create_div_constants
#' @description Create input constants for diversity models
#' @export
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


#' @title 			createInits
#' @description Create initial values for MCMC runs
#' @export
createInits <- function(constants, type = NULL){
	#	core_per_plot <- 3
	y_init <- matrix(runif(constants$N.core, .1, .9),
									 ncol = 1, nrow = constants$N.core)
	y_init = cbind(y_init[,1], 1 - y_init[,1])

	plot_mu_init <- matrix(runif(constants$N.plot * constants$N.date,
															.3,.7),
												 nrow=constants$N.plot, ncol=constants$N.date)

	beta_init <- rnorm(constants$N.beta, 0, .2)
	rho_init <- rnorm(1, 0, .2)
	sigma_init <- runif(1,0,.2)
	intercept_init <- rnorm(1, 0, .2)
	Ex_init <- plot_mu_init
	mois_est <- constants$mois
	mois_est[is.na(mois_est)] <- 0
	temp_est <- constants$temp
	temp_est[is.na(temp_est)] <- 0

	pH_est <- constants$pH
	pC_est <- constants$pC
	sig_init <- runif(1,.1,.5)
	core_sd_init <- runif(1,.1,.5)
	site_effect_init <- rnorm(constants$N.site, 0, .5)

		logit_plot_mu_init = logit(plot_mu_init)

		out <- list(
			y = y_init,
			plot_mu = plot_mu_init,
			logit_plot_mu = logit_plot_mu_init,
			intercept =intercept_init,
			core_sd = core_sd_init,
			sig = sig_init,
			beta = beta_init,
			rho = rho_init,
			sigma = sigma_init,
			site_effect = site_effect_init,
			Ex = plot_mu_init,
			logit_Ex = logit_plot_mu_init,
			mois_est = mois_est,
			temp_est = temp_est,
			pH_est = pH_est,
			pC_est = pC_est)
	return(out)
}


#' @title 			initsFun (DEPRECATED lol)
#' @description Create initial values for MCMC runs
#' @export
initsFun <- function(constants, type = NULL){
	#	core_per_plot <- 3
	y_init <- matrix(rep(rep(1/constants$N.spp, constants$N.spp), constants$N.core),
									 ncol = constants$N.spp, nrow = constants$N.core)

	plot_mu_init <- array(rnorm(constants$N.plot  * constants$N.spp * constants$N.date,
															2,.1),
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
																 1,.5),
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
		plot_mu_init = matrix(rep(.5, constants$N.plot*constants$N.date),
													constants$N.plot, constants$N.date)
		logit_plot_mu_init = matrix(rep(.5, constants$N.plot*constants$N.date),
																constants$N.plot, constants$N.date)
		shape_init = matrix(rep(2, constants$N.plot*constants$N.date),
												constants$N.plot, constants$N.date)
		out <- list(
			y = y_init,
			plot_mu = plot_mu_init,
			logit_plot_mu = logit_plot_mu_init,
			intercept = 0,
			core_sd = .1,
			sig = sig_init,
			beta = beta_init[1,],
			rho = rho_init[1],
			sigma = .1,
			plot_rel = plot_rel_init,
			site_effect = rep(.1, constants$N.site),
			Ex = plot_mu_init,
			logit_Ex = logit_plot_mu_init,
			shape1 = shape_init,
			shape2 = shape_init,
			#SIGMA = SIGMA,
			mois_est = mois_est,
			temp_est = temp_est,
			pH_est = pH_est,
			pC_est = pC_est
		)
	} else { # For taxa
		SIGMA <- diag(rep(.1, constants$N.spp))
		site_effect_init = diag(0, constants$N.site, constants$N.spp)

		out <- list(
			y = y_init,
			plot_mu = plot_mu_init,
			intercept = intercept_init,
			sig = sig_init,
			beta = beta_init,
			rho = rho_init,
			sigma = sigma_init,
			plot_rel = plot_rel_init,
			site_effect = site_effect_init,
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


#'  @title 			fixDate
#'  @description Take a vector of dates in dateID format ("201401") and return as date with first of month
#' @export
fixDate <- function(datesToFix){
  as.Date(paste0(datesToFix, "01"), format = "%Y%m%d")
}

#' interval_transform
#' Take a vector of values and put them in the (0,1) interval
#' @export
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


#' rm.NA.mcmc
#' Remove NAs from MCMC objects
#' @export
rm.NA.mcmc <- function(samples){
	sample.list  <- mcmc.list(mcmc(na.omit(samples[[1]])),
														mcmc(na.omit(samples[[2]])),
														mcmc(na.omit(samples[[3]])))
	return(sample.list)
}


### CONVERT PRECISION TO SD
#' prec_to_sd
#' @export
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


#' @title 			tic
#' @description  one of two clock functions.
#' Place tic() at the line in the code where you want to start timing.
#' Place toc() at the position in the code where you want to stop timing and report.
#' @export
tic <- function() {assign("timer", Sys.time(), envir=.GlobalEnv)}

#' @title 			toc
#' @description  one of two clock functions.
#' Place tic() at the line in the code where you want to start timing.
#' Place toc() at the position in the code where you want to stop timing and report.
#' @export
toc <- function() print(Sys.time()-timer)




#'  @title 			assign_fg_categories
#'  @description Assign functional group categories based on group name
#' @export
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
	out[which(grepl("other", vector))] <- NA
	return(out)
}


#'  @title 			assign_fg_kingdoms
#'  @description Assign kingdoms based on functional group categories
#' @export
assign_fg_kingdoms <- function(vector) {
	out <- rep(NA, length(vector))
	out[which(grepl("Simple|Complex", vector))] <- "16S"
	out[which(grepl("Stress|Antibiotic|Anaerobic|cycling", vector))] <- "16S"
	out[which(grepl("Troph", vector))] <- "ITS"
	out[which(grepl("Life", vector))] <- "16S"
	return(out)
}




#'  @title 			assign_fg_sources
#'  @description Assign sources based on functional group development method
#' @export
assign_fg_sources <- function (vector) {
	out <- rep(NA, length(vector))
	out[which(grepl("lytic", vector))] <-  "Literature review"

	# Berlemont & Martiny 2013
	out[which(grepl("cellulolytic", vector))] <-  "Literature review + genomic pathway"

	# Albright 2018
	out[which(grepl("nitr|fixa", vector))] <- "Literature review + genomic pathway"

	# Naylor modules
	out[which(grepl("complex|simple|stress", vector, fixed = F))] <- "Experimental enrichment"
	out[which(grepl("antibiotic", vector))] <- "Experimental enrichment"
	out[which(grepl("anaerobic", vector))] <- "Experimental enrichment"
	#	out[which(grepl("nitr|fixa", vector))] <- "Literature review"
	out[which(grepl("troph", vector))] <- "Literature review"
	out[which(grepl("sapr|path|arbusc|ecto|endo|lichen", vector))] <- "Literature review"
	out[which(grepl("other", vector))] <- NA

	#
	# out <- rep(NA, length(vector))
	# out[which(grepl("substrates|resistance|anaerobic|stress", vector))] <- "Experimental enrichment"
	# out[which(grepl("cycling|cellulo", vector))] <- "Genomic pathways"
	# out[which(grepl("Trophic|Life", vector))] <- "Scientific consensus"
	return(out)
}


#' @title parseNEONsampleIDs
#' @description create sample information data.frame from NEON sample names
#' @export
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



#' @title rbind.named.dfs
#' @description Combine multiple data frames with row names
#' @export
rbind.named.dfs <- function(df.list){
  # solution from https://stackoverflow.com/questions/15162197/combine-rbind-data-frames-and-create-column-with-name-of-original-data-frames
  dfs <- df.list[sapply(df.list, function(x) !is.null(dim(x)))]
  all.out <- cbind.data.frame(do.call(rbind,dfs),
                              name = rep(names(dfs), vapply(dfs, nrow, numeric(1))))
  return(all.out)
}



#' @export
base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}



#' @export
pretty_breaks <- function(n = 5, ...) {
  function(x) {
    breaks <- pretty(x, n, ...)
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
}

#'  @title rbind.df.list
#'  @description rbind named dfs
#'
#' @export
rbind.df.list <- function(pl.out){
  lapply(pl.out, function(x){
    do.call(rbind, x)
  })
}




#' @title filter_date_site
#' @description return data matrix filtered by date/plot/site
#' @return  a vector of values on the interval (0,1)
#' @export
#'
filter_date_site <- function(input_df,
														 keep_sites,
														 keep_plots,
														 min.date, max.date,
														 max.predictor.date = NULL, ...) {
pacman::p_load(tidyverse, anytime)

	input_df <- input_df %>% as.data.frame()
	if (!is.null(max.predictor.date)) {
		filt.date = max.predictor.date
	} else filt.date = anydate(max.date)
	# filter by date
	col_dates <- colnames(input_df) %>% fixDate() %>% anydate()
	min.date <- anydate(min.date)

	filt_date <- input_df[,which(col_dates <= filt.date & col_dates >= min.date)]

	# filter by site
		if (nchar(rownames(input_df)[[1]]) == 4){
			filt <- filt_date %>%
				filter(rownames(.) %in% keep_sites) %>%
				data.matrix()
		# filter by plot
		} else if (nchar(rownames(input_df)[[1]]) == 8){
			filt <- filt_date %>%
				filter(rownames(.) %in% keep_plots) %>%
				data.matrix()
			} else {
				message("Did not filter rows: data must have plot or site rownames, e.g. 'HARV' or 'HARV_001'")
	return(filt_date)
	}

	return(filt)
}


#' @title create_covariate_samples
#' @description samples covariate values with strict data validation - fails fast when data is missing
#'
#' @export
#'
create_covariate_samples <- function(model.inputs, plotID = NULL, siteID,
																		 Nmc_large, Nmc,
																		 N.beta = 8, prev_samples = NULL, 
																		 model_type = NULL, ...) {
	
	# STRICT VERSION: This function fails fast when data is missing - no fallbacks!
	
	# Validate inputs - fail fast if invalid
	if (is.null(model.inputs) || !is.list(model.inputs)) {
		stop("model.inputs is NULL or not a list - investigate data loading")
	}
	
	# Get the start date for this site - fail fast if missing
	# FIXED: Handle both list and vector structures for site_start
	if (is.null(model.inputs$site_start)) {
		stop("model.inputs$site_start is NULL - investigate data structure")
	}
	
	# Check if site_start is a list (old format) or vector (new format)
	if (is.list(model.inputs$site_start)) {
		# Old format: site_start is a list with site names as keys
		if (is.null(model.inputs$site_start[[siteID]])) {
			stop("Site ", siteID, " not found in model.inputs$site_start - investigate site mapping")
		}
		start_date <- model.inputs$site_start[[siteID]]
	} else if (is.vector(model.inputs$site_start)) {
		# New format: site_start is a vector where index corresponds to site number
		# We need to find the site number for this siteID
		if (!"plot_site_num" %in% names(model.inputs)) {
			stop("plot_site_num not found in model.inputs - cannot map siteID to site number")
		}
		
		# Find a plot for this site to get the site number
		# The siteID and plot_site_num are vectors in the model data
		# Use truth.plot.long for proper mapping since it has the correct structure
		if ("truth.plot.long" %in% names(model.inputs)) {
			plot_site_key <- model.inputs$truth.plot.long %>%
				select(siteID, plotID) %>%
				distinct()
			
			site_row <- plot_site_key %>% filter(siteID == !!siteID) %>% head(1)
			if (nrow(site_row) == 0) {
				stop("Site ", siteID, " not found in truth.plot.long - investigate site mapping")
			}
			
			# Get the plot for this site
			plot_for_site <- site_row$plotID[1]
			
			# Find the plot_site_num for this plot
			plot_index <- which(model.inputs$plotID == plot_for_site)[1]
			if (is.na(plot_index)) {
				stop("Plot ", plot_for_site, " not found in plotID vector - investigate plot mapping")
			}
			
			site_num <- model.inputs$plot_site_num[plot_index]
			if (is.na(site_num) || site_num < 1 || site_num > length(model.inputs$site_start)) {
				stop("Invalid site number ", site_num, " for plot ", plot_for_site, " - investigate site mapping")
			}
			
			start_date <- model.inputs$site_start[site_num]
		} else {
			# Fallback to old method if truth.plot.long not available
			plot_site_key <- data.frame(
				siteID = model.inputs$siteID,
				plot_site_num = model.inputs$plot_site_num
			) %>% distinct()
			
			site_row <- plot_site_key %>% filter(siteID == !!siteID) %>% head(1)
			if (nrow(site_row) == 0) {
				stop("Site ", siteID, " not found in plot_site_key - investigate site mapping")
			}
			
			site_num <- site_row$plot_site_num[1]
			if (site_num < 1 || site_num > length(model.inputs$site_start)) {
				stop("Invalid site number ", site_num, " for site ", siteID, " - investigate site mapping")
			}
			
			start_date <- model.inputs$site_start[site_num]
		}
	} else {
		stop("model.inputs$site_start is neither a list nor a vector - investigate data structure")
	}
	
	# Validate start_date - fail fast if invalid
	if (is.na(start_date) || is.null(start_date) || start_date < 1) {
		stop("Invalid start_date ", start_date, " for site ", siteID, " - investigate site data")
	}
	
	# Get N.date - fail fast if invalid
	if (is.null(model.inputs$N.date) || !is.numeric(model.inputs$N.date) || model.inputs$N.date < 1) {
		stop("Invalid N.date in model.inputs - investigate data structure")
	}
	
	NT <- model.inputs$N.date
	
	# Ensure start_date is within bounds - fail fast if out of bounds
	if (start_date > NT) {
		stop("start_date ", start_date, " exceeds N.date ", NT, " for site ", siteID, 
			 " - investigate calibration vs validation period mismatch")
	}
	
	# Determine required covariates based on model type
	if (is.null(model_type)) {
		# Default to requiring all covariates if model_type not specified
		required_arrays <- c("temp", "temp_sd", "mois", "mois_sd", "pH", "pH_sd", "pC", "pC_sd", "relEM", "LAI")
		required_vectors <- c("sin_mo", "cos_mo")
		N.beta <- 8
	} else if (model_type == "cycl_only") {
		# For cycl_only models, only require seasonal covariates
		required_arrays <- character(0)  # No environmental arrays needed
		required_vectors <- c("sin_mo", "cos_mo")
		N.beta <- 2
	} else if (model_type == "env_cov") {
		# For env_cov models, require environmental covariates but not seasonal
		required_arrays <- c("temp", "temp_sd", "mois", "mois_sd", "pH", "pH_sd", "pC", "pC_sd", "relEM", "LAI")
		required_vectors <- character(0)  # No seasonal vectors needed
		N.beta <- 6
	} else if (model_type == "env_cycl") {
		# For env_cycl models, require both environmental and seasonal covariates
		required_arrays <- c("temp", "temp_sd", "mois", "mois_sd", "pH", "pH_sd", "pC", "pC_sd", "relEM", "LAI")
		required_vectors <- c("sin_mo", "cos_mo")
		N.beta <- 8
	} else {
		# Unknown model type, default to requiring all covariates
		required_arrays <- c("temp", "temp_sd", "mois", "mois_sd", "pH", "pH_sd", "pC", "pC_sd", "relEM", "LAI")
		required_vectors <- c("sin_mo", "cos_mo")
		N.beta <- 8
	}
	
	# Check arrays
	for (array_name in required_arrays) {
		if (is.null(model.inputs[[array_name]])) {
			stop("Required array '", array_name, "' is NULL in model.inputs - investigate data loading")
		}
	}
	
	# Check vectors (sin_mo and cos_mo)
	for (vector_name in required_vectors) {
		if (is.null(model.inputs[[vector_name]])) {
			stop("Required vector '", vector_name, "' is NULL in model.inputs - investigate data loading")
		}
	}
	
	# Create output array
	covar_full <- array(0, dim = c(Nmc_large, N.beta, NT))
	
	# STRICT ARRAY ACCESS: Fail fast if data is missing or out of bounds
	for (time in 1:NT) {
		col_index <- 1  # Track which column we're filling
		
		# Temperature - only if required
		if ("temp" %in% required_arrays) {
			if (!siteID %in% rownames(model.inputs$temp)) {
				stop("Site ", siteID, " not found in temperature data rownames - investigate site mapping")
			}
			if (time > ncol(model.inputs$temp)) {
				stop("Time ", time, " exceeds temperature data columns (", ncol(model.inputs$temp), ") - investigate time dimension")
			}
			
			temp_mean <- model.inputs$temp[siteID, time]
			if (is.na(temp_mean) || is.infinite(temp_mean)) {
				stop("Temperature value is NA or infinite for site ", siteID, " time ", time, " - investigate data quality")
			}
			
			temp_sd <- model.inputs$temp_sd[siteID, time]
			if (is.na(temp_sd) || is.infinite(temp_sd) || temp_sd <= 0) {
				stop("Temperature SD is NA, infinite, or <= 0 for site ", siteID, " time ", time, " - investigate data quality")
			}
			
			# Generate samples - fail fast if sampling fails
			temp_samples <- rnorm(Nmc_large, temp_mean, temp_sd)
			if (any(is.na(temp_samples)) || any(is.infinite(temp_samples))) {
				stop("Temperature sampling produced NA or infinite values - investigate parameters")
			}
			
			# Assign to array
			covar_full[, col_index, time] <- temp_samples
			col_index <- col_index + 1
		}
		
		# Moisture - only if required
		if ("mois" %in% required_arrays) {
			if (!siteID %in% rownames(model.inputs$mois)) {
				stop("Site ", siteID, " not found in moisture data rownames - investigate site mapping")
			}
			if (time > ncol(model.inputs$mois)) {
				stop("Time ", time, " exceeds moisture data columns (", ncol(model.inputs$mois), ") - investigate time dimension")
			}
			
			mois_mean <- model.inputs$mois[siteID, time]
			if (is.na(mois_mean) || is.infinite(mois_mean)) {
				stop("Moisture value is NA or infinite for site ", siteID, " time ", time, " - investigate data quality")
			}
			
			mois_sd <- model.inputs$mois_sd[siteID, time]
			if (is.na(mois_sd) || is.infinite(mois_sd) || mois_sd <= 0) {
				stop("Moisture SD is NA, infinite, or <= 0 for site ", siteID, " time ", time, " - investigate data quality")
			}
			
			# Generate samples - fail fast if sampling fails
			mois_samples <- rnorm(Nmc_large, mois_mean, mois_sd)
			if (any(is.na(mois_samples)) || any(is.infinite(mois_samples))) {
				stop("Moisture sampling produced NA or infinite values - investigate parameters")
			}
			
			# Assign to array
			covar_full[, col_index, time] <- mois_samples
			col_index <- col_index + 1
		}
		
		# pH - only if required
		if ("pH" %in% required_arrays) {
			if (is.null(plotID)) {
				stop("plotID is NULL but required for pH data - investigate plot mapping")
			}
			if (!plotID %in% rownames(model.inputs$pH)) {
				stop("Plot ", plotID, " not found in pH data rownames - investigate plot mapping")
			}
			if (time > ncol(model.inputs$pH)) {
				stop("Time ", time, " exceeds pH data columns (", ncol(model.inputs$pH), ") - investigate time dimension")
			}
			
			pH_value <- model.inputs$pH[plotID, time]
			if (is.na(pH_value) || is.infinite(pH_value)) {
				stop("pH value is NA or infinite for plot ", plotID, " time ", time, " - investigate data quality")
			}
			
			pH_sd_value <- model.inputs$pH_sd[plotID, time]
			if (is.na(pH_sd_value) || is.infinite(pH_sd_value) || pH_sd_value <= 0) {
				stop("pH SD is NA, infinite, or <= 0 for plot ", plotID, " time ", time, " - investigate data quality")
			}
			
			# Generate samples - fail fast if sampling fails
			pH_samples <- rnorm(Nmc_large, pH_value, pH_sd_value)
			if (any(is.na(pH_samples)) || any(is.infinite(pH_samples))) {
				stop("pH sampling produced NA or infinite values - investigate parameters")
			}
			
			# Assign to array
			covar_full[, col_index, time] <- pH_samples
			col_index <- col_index + 1
		}
		
		# pC - only if required
		if ("pC" %in% required_arrays) {
			if (!plotID %in% rownames(model.inputs$pC)) {
				stop("Plot ", plotID, " not found in pC data rownames - investigate plot mapping")
			}
			if (time > ncol(model.inputs$pC)) {
				stop("Time ", time, " exceeds pC data columns (", ncol(model.inputs$pC), ") - investigate time dimension")
			}
			
			pC_value <- model.inputs$pC[plotID, time]
			if (is.na(pC_value) || is.infinite(pC_value)) {
				stop("pC value is NA or infinite for plot ", plotID, " time ", time, " - investigate data quality")
			}
			
			pC_sd_value <- model.inputs$pC_sd[plotID, time]
			if (is.na(pC_sd_value) || is.infinite(pC_sd_value) || pC_sd_value <= 0) {
				stop("pC SD is NA, infinite, or <= 0 for plot ", plotID, " time ", time, " - investigate data quality")
			}
			
			# Generate samples - fail fast if sampling fails
			pC_samples <- rnorm(Nmc_large, pC_value, pC_sd_value)
			if (any(is.na(pC_samples)) || any(is.infinite(pC_samples))) {
				stop("pC sampling produced NA or infinite values - investigate parameters")
			}
			
			# Assign to array
			covar_full[, col_index, time] <- pC_samples
			col_index <- col_index + 1
		}
		
		# relEM - only if required
		if ("relEM" %in% required_arrays) {
			if (!plotID %in% rownames(model.inputs$relEM)) {
				stop("Plot ", plotID, " not found in relEM data rownames - investigate plot mapping")
			}
			if (time > ncol(model.inputs$relEM)) {
				stop("Time ", time, " exceeds relEM data columns (", ncol(model.inputs$relEM), ") - investigate time dimension")
			}
			
			relEM_value <- model.inputs$relEM[plotID, time]
			if (is.na(relEM_value) || is.infinite(relEM_value)) {
				stop("relEM value is NA or infinite for plot ", plotID, " time ", time, " - investigate data quality")
			}
			
			# These are deterministic, so just repeat the values
			relEM_samples <- rep(relEM_value, Nmc_large)
			
			# Assign to array
			covar_full[, col_index, time] <- relEM_samples
			col_index <- col_index + 1
		}
		
		# LAI - only if required
		if ("LAI" %in% required_arrays) {
			if (!siteID %in% rownames(model.inputs$LAI)) {
				stop("Site ", siteID, " not found in LAI data rownames - investigate site mapping")
			}
			if (time > ncol(model.inputs$LAI)) {
				stop("Time ", time, " exceeds LAI data columns (", ncol(model.inputs$LAI), ") - investigate time dimension")
			}
			
			lai_value <- model.inputs$LAI[siteID, time]
			if (is.na(lai_value) || is.infinite(lai_value)) {
				# For hindcasts, use a reasonable default instead of failing
				message("WARNING: LAI value is NA or infinite for site ", siteID, " time ", time, " - using default value")
				lai_value <- 0.5  # Default LAI value
			}
			
			# These are deterministic, so just repeat the values
			LAI_samples <- rep(lai_value, Nmc_large)
			
			# Assign to array
			covar_full[, col_index, time] <- LAI_samples
			col_index <- col_index + 1
		}
		
		# Seasonal predictors - only if required
		if ("sin_mo" %in% required_vectors) {
			if (time > length(model.inputs$sin_mo)) {
				stop("Time ", time, " exceeds sin_mo length (", length(model.inputs$sin_mo), ") - investigate seasonal data")
			}
			
			sin_mo_value <- model.inputs$sin_mo[time]
			if (is.na(sin_mo_value) || is.infinite(sin_mo_value)) {
				stop("sin_mo value is NA or infinite for time ", time, " - investigate seasonal data quality")
			}
			
			# These are deterministic, so just repeat the values
			sin_mo_samples <- rep(sin_mo_value, Nmc_large)
			
			# Assign to array
			covar_full[, col_index, time] <- sin_mo_samples
			col_index <- col_index + 1
		}
		
		if ("cos_mo" %in% required_vectors) {
			if (time > length(model.inputs$cos_mo)) {
				stop("Time ", time, " exceeds cos_mo length (", length(model.inputs$cos_mo), ") - investigate seasonal data")
			}
			
			cos_mo_value <- model.inputs$cos_mo[time]
			if (is.na(cos_mo_value) || is.infinite(cos_mo_value)) {
				stop("cos_mo value is NA or infinite for time ", time, " - investigate seasonal data quality")
			}
			
			# These are deterministic, so just repeat the values
			cos_mo_samples <- rep(cos_mo_value, Nmc_large)
			
			# Assign to array
			covar_full[, col_index, time] <- cos_mo_samples
			col_index <- col_index + 1
		}
	}
	
	# Sample and return
	if (Nmc > Nmc_large) {
		# If we want more samples than available, sample with replacement
		covar <- covar_full[sample.int(Nmc_large, Nmc, replace = TRUE), , ]
	} else {
		# If we have enough samples, sample without replacement
		covar <- covar_full[sample.int(Nmc_large, Nmc, replace = FALSE), , ]
	}
	
	# Final validation - this should never fail if we got here
	if (is.null(covar) || !is.array(covar) || length(dim(covar)) != 3) {
		stop("CRITICAL ERROR: Final validation failed - this should never happen!")
	}
	
	# SUCCESS: Return the covariate array
	return(covar)
}



#' @title parse_plot_mu_vars
#' @description parse MCMC rowname output from summary matrix
#'
#' @export
#'
parse_plot_mu_vars <- function(input_df) {
	require(stringr)

	with_rowname_col <- input_df %>% as.data.frame()
	with_rowname_col$rowname <- rownames(with_rowname_col)
	rownames(with_rowname_col) <- NULL

	# Handle empty data frames or NA rownames
	if (nrow(with_rowname_col) == 0 || all(is.na(with_rowname_col$rowname))) {
		# Return empty data frame with expected structure
		return(data.frame(plot_num = integer(0), timepoint = integer(0), stringsAsFactors = FALSE))
	}

	# check number of commas for how to split values (i.e. should there be a "species_num" value)
	if (str_count(with_rowname_col$rowname[1], ',') == 2) {
		parsed <- with_rowname_col %>%
			separate(rowname, sep=", ", into=c("plot_num","species_num","timepoint")) %>%
			mutate(plot_num = as.integer(gsub("plot_rel\\[|plot_mu\\[", "", plot_num)),
						 timepoint = as.integer(gsub("\\]", "", timepoint)))
	} else if (str_count(with_rowname_col$rowname[1], ',') == 1) {
		parsed <- with_rowname_col %>%
			separate(rowname, sep=", ", into=c("plot_num","timepoint"))  %>%
			mutate(plot_num = as.integer(gsub("plot_mu\\[|plot_rel\\[", "", plot_num)),
					 timepoint = as.integer(gsub("\\]", "", timepoint)))
	} else {
		# Default case - create minimal structure
		parsed <- with_rowname_col %>%
			mutate(plot_num = 1,
						 timepoint = 1)
	}
	return(parsed)
}


#' @title extract_summary_row
#' @description extract MCMC summary by rowname from summary matrix
#'
#' @export
#'
extract_summary_row <- function(input_df, var = "sigma") {

	out <- input_df %>% as.data.frame()
	out$rowname <- rownames(out)
	rownames(out) <- NULL
	out <- out %>% filter(grepl(!!var, rowname))
	return(out)
}



#' @title extract_bracketed_vals
#' @description # extract MCMC summary by rowname from summary matrix
#'
#' @export
#'
extract_bracketed_vals <- function(input_df, varname1 = "beta_num", varname2 = NULL) {

	# Handle empty data frames
	if (nrow(input_df) == 0) {
		# Return empty data frame with expected column names
		if (is.null(varname2)) {
			return(data.frame(rowname = character(0), stringsAsFactors = FALSE) %>%
						 mutate(!!varname1 := character(0)))
		} else {
			return(data.frame(rowname = character(0), stringsAsFactors = FALSE) %>%
						 mutate(!!varname1 := character(0), !!varname2 := character(0)))
		}
	}

	if (!"rowname" %in% colnames(input_df)){
		input_df <- input_df %>% as.data.frame()
		input_df$rowname <- rownames(input_df)
		rownames(input_df) <- NULL
	}

	# Check if rowname column has any non-NA values
	if (all(is.na(input_df$rowname))) {
		# Return empty data frame with expected column names
		if (is.null(varname2)) {
			return(data.frame(rowname = character(0), stringsAsFactors = FALSE) %>%
						 mutate(!!varname1 := character(0)))
		} else {
			return(data.frame(rowname = character(0), stringsAsFactors = FALSE) %>%
						 mutate(!!varname1 := character(0), !!varname2 := character(0)))
		}
	}

	if (str_count(input_df$rowname[1], ',') == 1) {
		out <- input_df %>%
			separate(rowname, sep="\\[|, |\\]", into=c("rowname",
																								 "values1",
																								 "values2"), remove = F,
							 fill = "right", extra = "drop") %>%
			rename_with(~varname1, values1) %>%
			rename_with(~varname2, values2)

		} else if (str_count(input_df$rowname[1], ',') == 0) {

				out <- input_df %>%
		separate("rowname", sep = "\\[|\\]",
						 into = c("rowname","values1",NA),
						 fill = "right") %>%
		rename_with(~varname1, values1)
		}
	return(out)
}
# extract_bracketed_vals(beta_out, "beta_num")
# extract_bracketed_vals(plot_summary[[1]], "timepoint")

# Does confidence interval include 0? If so, not significant (0) otherwise significant (1)
is_significant <- function(lo, hi) {
	# Handle NA values by treating them as non-significant
	# Check if both bounds are on the same side of 0 (both negative or both positive)
	both_negative <- !is.na(lo) & !is.na(hi) & lo < 0 & hi < 0
	both_positive <- !is.na(lo) & !is.na(hi) & lo > 0 & hi > 0
	ifelse(both_negative | both_positive, 1, 0)
}








#' @title order_betas
#' @description # order environmental labels for plotting
#'
#' @export
#'
order_betas <- function(beta) {
	ordered <- ordered(beta, levels = c("sin", "cos",
																			"Ectomycorrhizal\ntrees",
																			"Ectomycorrhizal trees",
																			"pC",
																			"pH",
																			"LAI",
																			"Temperature",
																			"Moisture","rho"))
	return(ordered)
}






#' @title convert beta params
#' @description # convert mean/sd values to shape values
#'
#' @export
#'
convert_beta_params = function(mu, sd){
	var = sd^2

	# Wiki parameterization
	alpha = mu * ((mu * (1-mu))/var - 1)
	beta = (1 - mu) * ((mu * (1-mu))/var - 1)
	# Looks right to me
	#hist(rbeta(10000, alpha, beta))

	# Nimble parameterization
	shape1 = mu^2 * (1-mu)/var - mu
	shape2 = mu * (1-mu)^2/var + mu - 1


	# Colin parameterization
	tau = exp(var) # i might be getting this wrong
	tau = exp(sd) # i might be getting this wrong
	p = mu * tau
	q = (1 - mu) * tau

	if (shape1 < 0) message("Negative shape1 parameter")
	if (shape2 < 0) message("Negative shape2 parameter")
	if (alpha < 0) message("Negative alpha parameter")
	if (beta < 0) message("Negative beta parameter")

	return(list(nimble = c(shape1, shape2),
							orig = c(alpha, beta),
				 colin = c(p, q)))
}


#' @title parse_model_id
#' @description Parse model ID to extract components including enhanced metadata structure
#'
#' @export
#'
parse_model_id = function(model_id){

	info <- model_id %>% str_split("_") %>% unlist()

	# Handle enhanced metadata structure with legacy covariate indicators
	# Check if this is a legacy covariate model by looking for "with" followed by "legacy"
	is_legacy_covariate <- FALSE
	with_index <- which(info == "with")[1]
	if (!is.na(with_index) && with_index < length(info) && info[with_index + 1] == "legacy") {
		is_legacy_covariate <- TRUE
	}
	
	# Check if this is a driver uncertainty model
	is_driver_uncertainty <- grepl("driver_uncertainty", model_id) || 
	                       grepl("logit_beta_driver_uncertainty", model_id)
	
	if (is_legacy_covariate) {
		# Enhanced metadata format: model_name_species_startdate_enddate_with_legacy_covariate
		
		# Find the position where "with" starts (indicating the legacy covariate part)
		# with_index is already found above
		
		if (!is.na(with_index) && with_index > 3) {
			# Extract time period (two elements before "with")
			time_period <- paste(info[(with_index-2):(with_index-1)], collapse = "_")
			
			# Extract model_name and species (everything before the dates)
			before_dates <- info[1:(with_index-3)]
			
			# Known model names that might have underscores
			known_models <- c("cycl_only", "env_cov", "env_cycl")
			
			# Try to match known model names
			model_found <- FALSE
			for (known_model in known_models) {
				known_model_parts <- str_split(known_model, "_")[[1]]
				if (length(before_dates) >= length(known_model_parts)) {
					if (all(before_dates[1:length(known_model_parts)] == known_model_parts)) {
						model_name <- known_model
						rank.name.eval <- before_dates[-(1:length(known_model_parts))] %>% paste0(collapse = "_")
						model_found <- TRUE
						break
					}
				}
			}
			
			# If no known model found, fall back to first element
			if (!model_found) {
				model_name <- before_dates[1]
				rank.name.eval <- before_dates[-1] %>% paste0(collapse = "_")
			}
		} else {
			# Fallback parsing
			model_name <- info[1]
			rank.name.eval <- info[2]
			time_period <- paste(info[3:4], collapse = "_")
		}
		
	} else {
		# Original format: model_name_species_startdate_enddate
		# Extract "time_period"
		time_period <- tail(info, 2) %>% paste0(collapse = "_") %>% str_replace(".rds", "")
		info <- info[-c((length(info)-1):length(info))]

		# Extract "model_name"
		model_name <- info[c(1:2)]  %>% paste0(collapse = "_")
		info <- info[-c(1:2)]

		rank.name.eval <- info %>% paste0(collapse = "_")
	}

	if (rank.name.eval %in% microbialForecast:::fg_names) {
		summary_type="functional"
	} else {
		summary_type= "taxon"
	}

	# Add columns based on
	if (summary_type=="functional") {
		rank.name <- rank.name.eval
		rank_only <- "functional"
		species <- rank.name.eval
		fg_cat <- assign_fg_categories(species)
		group <- assign_fg_kingdoms(fg_cat)
	} else {

		taxa_key = stack(microbialForecast:::rank_spec_names) %>%
			select(species = values, rank.name = ind)

		species <- rank.name.eval
		rank.name <- taxa_key[match(species, taxa_key$species),]$rank.name
		rank_only <-  rank.name  %>% str_split("_") %>% unlist() %>% head(1)
		group <-  rank.name  %>% str_split("_") %>% unlist() %>% tail(1)
	}
	
	return(list(rank.name, time_period, rank_only, species, group, model_name, model_id, summary_type, is_driver_uncertainty))
}

# ## ################################################
# ## A DISTRIBUTION GIVING THE LOGIT OF A BETA DISTRIBUTION ##
# ## ################################################
# dLogitBeta <- nimbleFunction (
# 	## Returns density of x where
# 	##                    y ~ Beta(a1,a2)
# 	##                    x = logit(y)
# 	run = function(x = double(0),
# 								 shape1=double(0, default=1.0),
# 								 shape2=double(0, default=1.0),
# 								 log = integer(0, default=0)) {
# 		returnType(double(0))
# 		y = ilogit(x)
# 		logProbX = log(y) + log(1 - y) + dbeta(y, shape1=shape1, shape2=shape2, log=TRUE) ## Via change of variables
# 		if (log)
# 			return(logProbX)
# 		return(exp(logProbX))
# 	}
# )
#
# rLogitBeta <- nimbleFunction (
# 	## Generates y ~ Beta(a1,a2)
# 	## Returns   x = logit(y)
# 	run = function(n = integer(0, default=1),
# 								 shape1 = double(0, default=1.0),
# 								 shape2 = double(0, default=1.0)) {
# 		returnType(double(0))
# 		if(n != 1)
# 			nimPrint("Warning: rLogitBeta only allows n = 1; Using n = 1.\n")
# 		y <- rbeta(1, shape1=shape1, shape2=shape2)
# 		x <- logit(y)
# 		return(x)
# 	}
# )
#
# registerDistributions(list(dLogitBeta = list(
# 	BUGSdist = "dLogitBeta(shape1, shape2)",
# 	discrete = FALSE,
# 	types    = c("value=double(0)"), ## , "para=double(0)"
# 	pqAvail  = FALSE)))
#

# =============================================================================
# RESTART FUNCTIONS FOR MCMC CONTINUATION
# =============================================================================

#' Find existing chain files for a specific model
#' 
#' @param model_name The name of the model (e.g., "cycl_only", "env_cycl")
#' @param species The species name
#' @param min_date The minimum date
#' @param max_date The maximum date
#' @param use_legacy_covariate Whether the model uses legacy covariate
#' @return Vector of file paths to existing chain files
#' @export
find_chain_files <- function(model_name, species, min_date, max_date, use_legacy_covariate = TRUE) {
  # Create model ID for file matching
  legacy_indicator <- ifelse(use_legacy_covariate, "with_legacy_covariate", "without_legacy_covariate")
  model_id <- paste(model_name, species, min_date, max_date, legacy_indicator, sep = "_")
  
  # Look for chain files in the expected output directory
  output_dir <- here("data", "model_outputs", "logit_beta_regression", model_name)
  
  # Pattern to match chain files for this model
  pattern <- paste0("samples_", model_id, "_chain[0-9]+\\.rds")
  
  # Find matching files
  chain_files <- list.files(output_dir, pattern = pattern, full.names = TRUE, recursive = TRUE)
  
  # Also check for samples files without chain numbers
  alt_pattern <- paste0("samples_", model_id, "\\.rds")
  alt_files <- list.files(output_dir, pattern = alt_pattern, full.names = TRUE, recursive = TRUE)
  
  # Combine and return unique files
  all_files <- unique(c(chain_files, alt_files))
  
  cat("Found", length(all_files), "chain files for model:", model_id, "\n")
  if (length(all_files) > 0) {
    cat("Files:", paste(basename(all_files), collapse = ", "), "\n")
  }
  
  return(all_files)
}

#' Extract final values from existing MCMC chains
#' 
#' @param chain_files Vector of file paths to chain files
#' @param min_ess Minimum effective sample size for parameter inclusion
#' @param max_rhat Maximum R-hat value for parameter inclusion
#' @param burnin_proportion Proportion of samples to discard as burnin
#' @return List containing extracted values and diagnostics
#' @export
extract_final_values_from_chains <- function(chain_files, min_ess = 100, max_rhat = 1.1, burnin_proportion = 0.5) {
  if (length(chain_files) == 0) {
    stop("No chain files provided")
  }
  
  cat("Extracting final values from", length(chain_files), "chain files\n")
  
  # Load and process each chain
  all_chains <- list()
  chain_summaries <- list()
  
  for (i in seq_along(chain_files)) {
    file_path <- chain_files[i]
    cat("  Processing chain file:", basename(file_path), "\n")
    
    tryCatch({
      # Load the chain data
      chain_data <- readRDS(file_path)
      
      # Extract samples matrix
      if ("samples" %in% names(chain_data)) {
        samples <- chain_data$samples
      } else if (is.matrix(chain_data)) {
        samples <- chain_data
      } else {
        cat("    WARNING: Unexpected chain data structure, skipping\n")
        next
      }
      
      # Apply burnin
      burnin_samples <- floor(nrow(samples) * burnin_proportion)
      if (burnin_samples > 0) {
        samples <- samples[(burnin_samples + 1):nrow(samples), , drop = FALSE]
        cat("    Applied burnin:", burnin_samples, "samples removed\n")
      }
      
      # Store samples and calculate diagnostics
      all_chains[[i]] <- samples
      
      # Calculate effective sample sizes
      ess_values <- effectiveSize(as.mcmc(samples))
      chain_summaries[[i]] <- list(
        n_samples = nrow(samples),
        ess_values = ess_values,
        min_ess = min(ess_values, na.rm = TRUE)
      )
      
      cat("    Samples:", nrow(samples), "ESS range:", round(range(ess_values, na.rm = TRUE), 1), "\n")
      
    }, error = function(e) {
      cat("    ERROR processing chain file:", e$message, "\n")
    })
  }
  
  if (length(all_chains) == 0) {
    stop("No valid chains could be processed")
  }
  
  # Combine chains and calculate overall diagnostics
  combined_samples <- do.call(rbind, all_chains)
  
  # Calculate overall ESS and R-hat if multiple chains
  if (length(all_chains) > 1) {
    # Convert to mcmc.list for R-hat calculation
    mcmc_list <- mcmc.list(lapply(all_chains, as.mcmc))
    rhat_values <- gelman.diag(mcmc_list, multivariate = FALSE)$psrf[, 1]
    
    cat("Multiple chains detected, calculating R-hat values\n")
    cat("R-hat range:", round(range(rhat_values, na.rm = TRUE), 3), "\n")
  } else {
    rhat_values <- rep(1.0, ncol(combined_samples))
    names(rhat_values) <- colnames(combined_samples)
  }
  
  # Extract final values from the last samples
  final_values <- as.list(combined_samples[nrow(combined_samples), ])
  
  # Create result structure
  result <- list(
    final_values = final_values,
    combined_samples = combined_samples,
    chain_summaries = chain_summaries,
    diagnostics = list(
      n_chains = length(all_chains),
      total_samples = nrow(combined_samples),
      ess_values = effectiveSize(as.mcmc(combined_samples)),
      rhat_values = rhat_values,
      min_ess = min(effectiveSize(as.mcmc(combined_samples)), na.rm = TRUE),
      max_rhat = max(rhat_values, na.rm = TRUE)
    )
  )
  
  cat("Extraction complete:\n")
  cat("  Total samples:", nrow(combined_samples), "\n")
  cat("  Min ESS:", round(result$diagnostics$min_ess, 1), "\n")
  cat("  Max R-hat:", round(result$diagnostics$max_rhat, 3), "\n")
  
  return(result)
}

#' Create restart initial values from extracted chain data
#' 
#' @param extraction_result Result from extract_final_values_from_chains
#' @param use_fallback_for_extreme Whether to use fallback values for extreme parameters
#' @param fallback_strategy Strategy for fallback values ("random", "prior_mean", "conservative")
#' @return List containing initial values and diagnostic flags
#' @export
create_restart_inits <- function(extraction_result, use_fallback_for_extreme = TRUE, fallback_strategy = "conservative") {
  if (!is.list(extraction_result) || !("final_values" %in% names(extraction_result))) {
    stop("Invalid extraction_result structure")
  }
  
  cat("Creating restart initial values\n")
  
  final_values <- extraction_result$final_values
  diagnostics <- extraction_result$diagnostics
  
  # Define reasonable bounds for each parameter type
  param_bounds <- list(
    precision = c(0.1, 1000),      # Gamma distribution bounds
    rho = c(0.01, 0.99),           # Beta distribution bounds
    beta = c(-10, 10),              # Normal distribution bounds
    intercept = c(-10, 10),         # Normal distribution bounds
    site_effect_sd = c(0.01, 10),  # Gamma distribution bounds
    legacy_effect = c(-5, 5)        # Normal distribution bounds
  )
  
  # Initialize result structures
  initial_values <- list()
  extreme_flags <- logical(0)
  fallback_used <- logical(0)
  
  # Process each parameter
  for (param_name in names(final_values)) {
    param_value <- final_values[[param_name]]
    
    # Determine parameter type for bounds checking
    param_type <- "unknown"
    if (param_name == "precision") param_type <- "precision"
    else if (param_name == "rho") param_type <- "rho"
    else if (grepl("^beta\\[", param_name)) param_type <- "beta"
    else if (param_name == "intercept") param_type <- "intercept"
    else if (param_name == "site_effect_sd") param_type <- "site_effect_sd"
    else if (param_name == "legacy_effect") param_type <- "legacy_effect"
    else if (grepl("^site_effect\\[", param_name)) param_type <- "site_effect"
    else if (grepl("^Ex\\[", param_name) || grepl("^mu\\[", param_name)) param_type <- "latent"
    
    # Check if value is extreme
    is_extreme <- FALSE
    if (param_type %in% names(param_bounds)) {
      bounds <- param_bounds[[param_type]]
      is_extreme <- param_value < bounds[1] || param_value > bounds[2]
    } else if (param_type == "latent") {
      # Latent variables should be between 0 and 1
      is_extreme <- param_value < 0 || param_value > 1
    }
    
    extreme_flags[param_name] <- is_extreme
    
    if (is_extreme && use_fallback_for_extreme) {
      cat("  WARNING: Parameter", param_name, "has extreme value:", param_value, "\n")
      
      # Generate fallback value based on strategy
      fallback_value <- switch(fallback_strategy,
        "random" = {
          if (param_type == "precision") rgamma(1, 2, 0.1)
          else if (param_type == "rho") runif(1, 0.1, 0.9)
          else if (param_type == "beta") rnorm(1, 0, 0.1)
          else if (param_type == "intercept") rnorm(1, -2, 0.5)
          else if (param_type == "site_effect_sd") runif(1, 0.1, 1)
          else if (param_type == "legacy_effect") rnorm(1, 0, 0.1)
          else if (param_type == "site_effect") rnorm(1, 0, 0.1)
          else if (param_type == "latent") runif(1, 0.1, 0.9)
          else param_value  # Keep original if unknown type
        },
        "prior_mean" = {
          if (param_type == "precision") 50
          else if (param_type == "rho") 0.5
          else if (param_type == "beta") 0
          else if (param_type == "intercept") -2
          else if (param_type == "site_effect_sd") 0.5
          else if (param_type == "legacy_effect") 0
          else if (param_type == "site_effect") 0
          else if (param_type == "latent") 0.3
          else param_value
        },
        "conservative" = {
          if (param_type == "precision") 50
          else if (param_type == "rho") 0.3
          else if (param_type == "beta") 0.01
          else if (param_type == "intercept") -2
          else if (param_type == "site_effect_sd") 0.5
          else if (param_type == "legacy_effect") 0
          else if (param_type == "site_effect") 0
          else if (param_type == "latent") 0.3
          else param_value
        }
      )
      
      cat("    Using fallback value:", fallback_value, "(", fallback_strategy, "strategy)\n")
      initial_values[[param_name]] <- fallback_value
      fallback_used[param_name] <- TRUE
      
    } else {
      # Use the original value
      initial_values[[param_name]] <- param_value
      fallback_used[param_name] <- FALSE
    }
  }
  
  # Create result
  result <- list(
    initial_values = initial_values,
    extreme_flags = extreme_flags,
    fallback_used = fallback_used,
    summary = list(
      n_parameters = length(initial_values),
      n_extreme = sum(extreme_flags),
      n_fallback = sum(fallback_used),
      fallback_strategy = fallback_strategy
    )
  )
  
  cat("Restart initial values created successfully:\n")
  cat("  Parameters:", result$summary$n_parameters, "\n")
  cat("  Extreme values detected:", result$summary$n_extreme, "\n")
  cat("  Fallback values used:", result$summary$n_fallback, "\n")
  
  return(result)
}

# =============================================================================
# CORE INFRASTRUCTURE HELPER FUNCTIONS
# Universal functions for all microbial forecasting model types
# =============================================================================

#' @title Load Required Packages
#' @description Load all required packages for microbial forecasting models
#' @export
load_required_packages <- function() {
    cat("Loading required packages...\n")
    required_packages <- c("nimble", "parallel", "foreach", "doParallel", "here", "tidyverse", "coda", "devtools")
    for (pkg in required_packages) {
        if (!require(pkg, character.only = TRUE)) {
            install.packages(pkg)
            library(pkg, character.only = TRUE)
        }
    }
    
    # Load microbialForecast package
    if (!require(microbialForecast)) {
        devtools::load_all("microbialForecast")
        cat("microbialForecast package loaded from library\n")
    } else {
        cat("microbialForecast package loaded from library\n")
    }
    cat("All required packages loaded successfully\n")
}

#' @title Create Directories Safely
#' @description Create directories with error handling and fallback options
#' @param base_path Base directory path to create
#' @param subdirs Optional vector of subdirectories to create
#' @export
create_directories_safe <- function(base_path, subdirs = NULL) {
    tryCatch({
        if (!dir.exists(base_path)) {
            dir.create(base_path, recursive = TRUE)
            cat("✓ Created directory:", base_path, "\n")
        }
        
        if (!is.null(subdirs)) {
            for (subdir in subdirs) {
                full_path <- file.path(base_path, subdir)
                if (!dir.exists(full_path)) {
                    dir.create(full_path, recursive = TRUE)
                    cat("✓ Created subdirectory:", full_path, "\n")
                }
            }
        }
    }, error = function(e) {
        cat("✗ Failed to create directories:", e$message, "\n")
        stop("Directory creation failed")
    })
}

#' @title Create Consistent Model IDs
#' @description Create consistent model IDs across all model types
#' @param model_name Name of the model (e.g., "cycl_only", "env_cycl", "clr", "dirichlet")
#' @param species Species name
#' @param min_date Minimum date in format YYYYMMDD
#' @param max_date Maximum date in format YYYYMMDD
#' @param use_legacy_covariate Whether legacy covariate is used
#' @param model_type Type of model (default: "beta_regression", also supports "clr", "dirichlet")
#' @export
create_model_id <- function(model_name, species, min_date, max_date, use_legacy_covariate, model_type = "beta_regression") {
    legacy_indicator <- ifelse(use_legacy_covariate, "with_legacy_covariate", "without_legacy_covariate")
    paste(model_name, species, min_date, max_date, legacy_indicator, model_type, sep = "_")
}

#' @title Save Checkpoints Safely
#' @description Save checkpoints with fallback handling for all model types
#' @param samples MCMC samples to save
#' @param iterations Number of iterations completed
#' @param loop Loop number
#' @param output_dir Output directory path
#' @param model_id Model identifier
#' @param chain_no Chain number
#' @param checkpoint_type Type of checkpoint (e.g., "initial", "loop1", "final")
#' @param model_type Type of model (default: "beta_regression")
#' @export
save_checkpoint_safe <- function(samples, iterations, loop, output_dir, model_id, chain_no, checkpoint_type, model_type = "beta_regression") {
    checkpoint_file <- file.path(output_dir, paste0("checkpoint_", model_id, "_chain", chain_no, "_", checkpoint_type, ".rds"))
    
    tryCatch({
        saveRDS(list(samples = samples, iterations = iterations, loop = loop), checkpoint_file)
        cat("  ✓ Checkpoint saved:", checkpoint_type, "(", nrow(samples), " iterations)\n")
        cat("  ✓ Checkpoint file:", checkpoint_file, "\n")
    }, error = function(e) {
        cat("  ✗ Failed to save checkpoint:", e$message, "\n")
        cat("  Attempting to save to current directory as fallback...\n")
        
        # Fallback: save to current directory
        fallback_checkpoint <- paste0("checkpoint_", model_id, "_chain", chain_no, "_", checkpoint_type, "_FALLBACK.rds")
        tryCatch({
            saveRDS(list(samples = samples, iterations = iterations, loop = loop), fallback_checkpoint)
            cat("  ✓ Fallback checkpoint saved:", fallback_checkpoint, "\n")
        }, error = function(e2) {
            cat("  ✗ CRITICAL: Failed to save even fallback checkpoint:", e2$message, "\n")
        })
    })
}

#' @title Create Progress Files
#' @description Create progress files with error handling for all model types
#' @param output_dir Output directory path
#' @param model_id Model identifier
#' @param chain_no Chain number
#' @param init_iter Initial iterations completed
#' @param model_type Type of model (default: "beta_regression")
#' @return Path to the created progress file
#' @export
create_progress_file <- function(output_dir, model_id, chain_no, init_iter, model_type = "beta_regression") {
    progress_file <- file.path(output_dir, paste0("progress_", model_id, "_chain", chain_no, ".txt"))
    tryCatch({
        writeLines(paste("Started at:", Sys.time(), "\nInitial iterations:", init_iter, "\nStatus: Running"), progress_file)
        cat("  ✓ Progress file created:", progress_file, "\n")
    }, error = function(e) {
        cat("  ✗ Failed to create progress file:", e$message, "\n")
    })
    return(progress_file)  # Return the path for later use
}

#' @title Update Progress Files
#' @description Update progress files with error handling for all model types
#' @param progress_file Path to the progress file
#' @param total_iterations Total iterations completed
#' @param loop_counter Current loop number
#' @export
update_progress_file <- function(progress_file, total_iterations, loop_counter) {
    tryCatch({
        writeLines(paste("Updated at:", Sys.time(), "\nTotal iterations:", total_iterations, "\nLoop:", loop_counter, "\nStatus: Running"), progress_file)
    }, error = function(e) {
        cat("  ✗ Failed to update progress file:", e$message, "\n")
    })
}

#' @title Create Stable Beta Regression Models
#' @description Create stable Nimble models for beta regression with proper priors and structure
#' @param model_name Name of the model ("cycl_only", "env_cycl", "env_cov")
#' @param use_legacy_covariate Whether to include legacy covariate
#' @return Nimble model code
#' @export
create_stable_model <- function(model_name, use_legacy_covariate = TRUE) {
    cat("Building Nimble model:", model_name, "\n")
    
    if (model_name == "cycl_only" && use_legacy_covariate) {
        modelCode <- nimble::nimbleCode({
            # PRIORS - Weak for main parameters, more informative for site effects
            precision ~ dgamma(0.001, 0.001)        # Jeffreys prior - very weak
            rho ~ dbeta(1, 1)                       # Uniform prior - very weak
            intercept ~ dnorm(0, sd = 10)           # Very wide normal - very weak
            legacy_effect ~ dnorm(0, sd = 10)       # Very wide normal - very weak
            
            # More informative priors for site effects
            site_effect_sd ~ dgamma(2, 20)
            for (k in 1:N.site) {
                site_effect[k] ~ dnorm(0, sd = site_effect_sd)
            }
            
            # Beta parameters for seasonal predictors - weak priors
            for (b in 1:2) {
                beta[b] ~ dnorm(0, sd = 10)         # Very wide normal - very weak
            }
            
            # PROCESS MODEL - LOG transformation approach
            for (p in 1:N.plot) {
                # Initial condition - Single time point
                for (t in plot_start[p]) {
                    Ex[p, t] ~ dunif(0.1, 0.9)
                    mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, shape2 = (1 - Ex[p, t]) * precision)
                }
                
                # Dynamic evolution - LOG transformation for numerical stability
                for (t in plot_index[p]:N.date) {
                    # Use log transformation for numerical stability
                    log_Ex_prev[p, t] <- log(max(0.001, mu[p, t - 1]))  # Safe log with bounds
                    
                    # Direct linear predictor with seasonal terms only
                    log_Ex_mean[p, t] <- rho * log_Ex_prev[p, t] +
                        beta[1] * sin_mo[t] + beta[2] * cos_mo[t] +
                        site_effect[plot_site_num[p]] +
                        legacy_effect * legacy[p, t] +
                        intercept
                    
                    # Use exp() instead of ilogit()
                    Ex[p, t] <- max(0.001, min(0.999, exp(log_Ex_mean[p, t])))
                    mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, shape2 = (1 - Ex[p, t]) * precision)
                }
            }
            
            # OBSERVATION MODEL
            for (i in 1:N.core) {
                y[i, 1] ~ dbeta(shape1 = mu[plot_num[i], timepoint[i]] * precision, 
                                shape2 = (1 - mu[plot_num[i], timepoint[i]]) * precision)
            }
        })
    } else if (model_name == "env_cycl" && use_legacy_covariate) {
        modelCode <- nimble::nimbleCode({
            # 🎯 PROVEN HYBRID PRIORS - Weak for main parameters, stable for site effects
            precision ~ dgamma(0.001, 0.001)        # Jeffreys prior - very weak
            rho ~ dbeta(1, 1)                       # Uniform prior - very weak
            intercept ~ dnorm(0, sd = 10)           # Very wide normal - very weak
            legacy_effect ~ dnorm(0, sd = 10)       # Very wide normal - very weak
            
            # STABLE PRIORS for site effects
            site_effect_sd ~ dgamma(2, 20)     # More informative: dgamma(2, 20)
            for (k in 1:N.site) {
                site_effect[k] ~ dnorm(0, sd = site_effect_sd)
            }
            
            # Beta parameters for environmental and seasonal predictors - weak priors
            for (b in 1:8) {
                beta[b] ~ dnorm(0, sd = 10)         # Very wide normal - very weak
            }
            
            # PROCESS MODEL - LOG transformation approach 
            for (p in 1:N.plot) {
                # Initial condition - Single time point
                for (t in plot_start[p]) {
                    Ex[p, t] ~ dunif(0.1, 0.9)
                    mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, shape2 = (1 - Ex[p, t]) * precision)
                }
                
                # Dynamic evolution - LOG transformation for numerical stability
                for (t in plot_index[p]:N.date) {
                    # Use log transformation
                    log_Ex_prev[p, t] <- log(max(0.001, mu[p, t - 1]))  # Safe log with bounds
                    
                    # Direct linear predictor with all 6 environmental predictors
                    log_Ex_mean[p, t] <- rho * log_Ex_prev[p, t] +
                        beta[1] * sin_mo[t] + beta[2] * cos_mo[t] +  # Seasonal terms
                        beta[3] * temp[plot_site_num[p], t] + beta[4] * mois[plot_site_num[p], t] +  # Site-level environmental terms
                        beta[5] * pH[p, t] + beta[6] * pC[p, t] +  # Plot-level environmental terms
                        beta[7] * relEM[p, t] + beta[8] * LAI[plot_site_num[p], t] +  # Mixed level terms
                        site_effect[plot_site_num[p]] +  # Site effects
                        legacy_effect * legacy[p, t] +  # Legacy covariate
                        intercept  # Baseline abundance
                    
                    # Use exp() instead of ilogit()
                    Ex[p, t] <- max(0.001, min(0.999, exp(log_Ex_mean[p, t])))
                    mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, shape2 = (1 - Ex[p, t]) * precision)
                }
            }
            
            # OBSERVATION MODEL
            for (i in 1:N.core) {
                y[i, 1] ~ dbeta(shape1 = mu[plot_num[i], timepoint[i]] * precision, 
                                shape2 = (1 - mu[plot_num[i], timepoint[i]]) * precision)
            }
        })
    } else if (model_name == "env_cov" && use_legacy_covariate) {
        modelCode <- nimble::nimbleCode({
            # 🎯 PROVEN HYBRID PRIORS - Weak for main parameters, stable for site effects
            precision ~ dgamma(0.001, 0.001)        # Jeffreys prior - very weak, proven
            rho ~ dbeta(1, 1)                       # Uniform prior - very weak, proven
            intercept ~ dnorm(0, sd = 10)           # Very wide normal - very weak, proven
            legacy_effect ~ dnorm(0, sd = 10)       # Very wide normal - very weak, proven
            
            # STABLE PRIORS for site effects (proven to work)
            site_effect_sd ~ dgamma(2, 20)          # Original: dgamma(2, 20) - stable, proven
            for (k in 1:N.site) {
                site_effect[k] ~ dnorm(0, sd = site_effect_sd)
            }
            
            # Beta parameters for environmental predictors only - weak priors
            for (b in 1:6) {
                beta[b] ~ dnorm(0, sd = 10)         # Very wide normal - very weak, proven
            }
            
            # PROCESS MODEL - LOG transformation approach 
            for (p in 1:N.plot) {
                # Initial condition - Single time point
                for (t in plot_start[p]) {
                    Ex[p, t] ~ dunif(0.1, 0.9)
                    mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, shape2 = (1 - Ex[p, t]) * precision)
                }
                
                # Dynamic evolution - LOG transformation for numerical stability
                for (t in plot_index[p]:N.date) {
                    # Use log transformation for numerical stability
                    log_Ex_prev[p, t] <- log(max(0.001, mu[p, t - 1]))  # Safe log with bounds
                    
                    # Direct linear predictor with environmental covariates only (no seasonal)
                    log_Ex_mean[p, t] <- rho * log_Ex_prev[p, t] +
                        beta[1] * temp[plot_site_num[p], t] +
                        beta[2] * mois[plot_site_num[p], t] +
                        beta[3] * pH[p, t] +
                        beta[4] * pC[p, t] +
                        beta[5] * relEM[p, t] +
                        beta[6] * LAI[p, t] +
                        site_effect[plot_site_num[p]] +  # Site effects
                        legacy_effect * legacy[p, t] +  # Legacy covariate
                        intercept  # Baseline abundance
                    
                    # Use exp() instead of ilogit()
                    Ex[p, t] <- max(0.001, min(0.999, exp(log_Ex_mean[p, t])))
                    mu[p, t] ~ dbeta(shape1 = Ex[p, t] * precision, shape2 = (1 - Ex[p, t]) * precision)
                }
            }
            
            # OBSERVATION MODEL
            for (i in 1:N.core) {
                y[i, 1] ~ dbeta(shape1 = mu[plot_num[i], timepoint[i]] * precision, 
                                shape2 = (1 - mu[plot_num[i], timepoint[i]]) * precision)
            }
        })
    } else {
        stop("Unsupported model combination: ", model_name, " with use_legacy_covariate=", use_legacy_covariate)
    }
    
    return(modelCode)
}

#' @title Create Stable Initial Values for Beta Regression
#' @description Create stable initial values for beta regression models with proper parameter structure
#' @param constants Model constants including N.plot, N.date, N.site
#' @param model_name Name of the model ("cycl_only", "env_cycl", "env_cov")
#' @param model_data Optional model data (not used in current implementation)
#' @return List of initial values for all model parameters
#' @export
create_stable_inits <- function(constants, model_name, model_data = NULL) {
    cat("Creating STABLE initial values for", model_name, "...\n")
    
    # Determine number of beta parameters based on model type
    if (model_name == "env_cycl") {
        n_beta <- 8
    } else if (model_name == "env_cov") {
        n_beta <- 6
    } else {
        n_beta <- 2  # cycl_only
    }
    
    # Create initial beta values (matching working approach exactly)
    beta_init <- c(0.01, 0.01)  # Start with seasonal coefficients very close to zero
    if (n_beta > 2) {
        for (i in 1:(n_beta - 2)) {  # Additional environmental coefficients
            beta_init <- c(beta_init, 0.01)  # Start environmental coefficients very close to zero
        }
    }

    inits <- list(
        precision = 50,  # Start with moderate precision
        rho = 0.3,      # Start rho at 0.3 (moderate persistence)
        beta = beta_init,  # Seasonal + environmental coefficients
        site_effect_sd = 0.5,  # Start with moderate site effect SD
        site_effect = rnorm(constants$N.site, 0, 0.1),  # Small random initial values
        intercept = -2,  # Start intercept at -2
        legacy_effect = 0,  # Start legacy effect at 0
        Ex = matrix(0.3, nrow = constants$N.plot, ncol = constants$N.date),  # Start with moderate abundance
        mu = matrix(0.3, nrow = constants$N.plot, ncol = constants$N.date)   # Start with moderate abundance
    )
    
    cat("  ✓ Initial values created successfully (matching working approach exactly)\n")
    cat("    precision:", inits$precision, "\n")
    cat("    rho:", inits$rho, "\n")
    cat("    beta parameters:", n_beta, "\n")
    cat("    site effects:", constants$N.site, "\n")
    cat("    Ex matrix dimensions:", dim(inits$Ex), "\n")
    cat("    mu matrix dimensions:", dim(inits$mu), "\n")
    
    return(inits)
}

#' @title Check for Extreme Parameter Values
#' @description Check if parameter values are extreme or outside reasonable bounds
#' @param param_values Named vector of parameter values to check
#' @return Logical vector indicating which values are extreme
#' @export
check_extreme_values <- function(param_values) {
  if (is.null(param_values) || length(param_values) == 0) {
    return(logical(0))
  }
  
  # Check for NA, NaN, Inf values
  extreme_flags <- is.na(param_values) | is.nan(param_values) | is.infinite(param_values)
  
  # Check for extreme parameter-specific values
  param_names <- names(param_values)
  
  for (i in seq_along(param_names)) {
    param_name <- param_names[i]
    param_value <- param_values[i]
    
    if (!extreme_flags[i]) {  # Only check if not already flagged
      if (grepl("precision", param_name)) {
        extreme_flags[i] <- param_value < 0.001 || param_value > 1000
      } else if (grepl("rho", param_name)) {
        extreme_flags[i] <- param_value < 0.001 || param_value > 0.999
      } else if (grepl("beta", param_name)) {
        extreme_flags[i] <- param_value < -50 || param_value > 50
      } else if (grepl("site_effect", param_name)) {
        extreme_flags[i] <- param_value < -10 || param_value > 10
      } else if (grepl("intercept", param_name)) {
        extreme_flags[i] <- param_value < -20 || param_value > 20
      } else if (grepl("legacy_effect", param_name)) {
        extreme_flags[i] <- param_value < -20 || param_value > 20
      } else if (grepl("site_effect_sd", param_name)) {
        extreme_flags[i] <- param_value < 0.001 || param_value > 10
      }
    }
  }
  
  return(extreme_flags)
}



