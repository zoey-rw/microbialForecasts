# Helper functions and global variables for soil microbial forecasts






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



#'  @title 			parseNEONsampleIDs
#'  @description create sample information data.frame from NEON sample names
#'  @export
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
#' @description
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
#' @description samples covariate values
#'
#' @export
#'
create_covariate_samples <- function(model.inputs, plotID = NULL, siteID,
																		 Nmc_large, Nmc,
																		 N.beta = 8, prev_samples = NULL, ...) {
if (!is.null(prev_samples)){


}
	else {
	start_date <- model.inputs$site_start[siteID]
	NT = model.inputs$N.date
	covar_full <- array(NA, dim = c(Nmc_large, N.beta, NT))

	set.seed(1)

	for (time in start_date:NT) {
	covar_full[,,time] <- c(Rfast::Rnorm(Nmc_large, model.inputs$temp[siteID, time],
																			 model.inputs$temp_sd[siteID, time]),
													Rfast::Rnorm(Nmc_large, model.inputs$mois[siteID, time],
																			 model.inputs$mois_sd[siteID, time]),
													Rfast::Rnorm(Nmc_large, model.inputs$pH[plotID,start_date],
																			 model.inputs$pH_sd[plotID,start_date]),
													Rfast::Rnorm(Nmc_large, model.inputs$pC[plotID,start_date],
																			 model.inputs$pC_sd[plotID,start_date]),
													rep(model.inputs$relEM[plotID, time], Nmc_large),
													rep(model.inputs$LAI[siteID, time], Nmc_large),
													rep(model.inputs$sin_mo[time], Nmc_large),
													rep(model.inputs$cos_mo[time], Nmc_large))
}
covar <- covar_full[sample.int(Nmc_large,Nmc),,]
return(covar)
}
}



#' @title parse_plot_mu_vars
#' @description parse MCMC rowname output from summary matrix
#'
#' @export
#'
parse_plot_mu_vars <- function(input_df) {
	require(stringr)

	with_rowname_col <-input_df %>% as.data.frame() %>%
		rownames_to_column()

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
	}
	return(parsed)
}


#' @title extract_summary_row
#' @description extract MCMC summary by rowname from summary matrix
#'
#' @export
#'
extract_summary_row <- function(input_df, var = "sigma") {

	out <- input_df %>% as.data.frame() %>%
		rownames_to_column("rowname") %>%
		filter(grepl(!!var, rowname))
	return(out)
}



#' @title extract_bracketed_vals
#' @description # extract MCMC summary by rowname from summary matrix
#'
#' @export
#'
extract_bracketed_vals <- function(input_df, varname1 = "beta_num", varname2 = NULL) {


	if (!"rowname" %in% colnames(input_df)){
		input_df <-input_df %>% as.data.frame() %>%
			rownames_to_column()
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
	ifelse(lo < 0 & hi < 0 |
				 	lo > -0 & lo > -0, 1, 0)
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
#' @description  parse_model_id
#'
#' @export
#'
parse_model_id = function(model_id){

	info <- model_id %>% str_split("_") %>% unlist()

	# Extract "time_period"
		time_period <- tail(info, 2) %>% paste0(collapse = "_") %>% str_replace(".rds", "")
		info <- info[-c((length(info)-1):length(info))]

	# Extract "model_name"
		model_name <- info[c(1:2)]  %>% paste0(collapse = "_")
		info <- info[-c(1:2)]

	rank.name.eval <- info %>% paste0(collapse = "_")

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
return(list(rank.name, time_period, rank_only, species, group, model_name, model_id, summary_type))
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



