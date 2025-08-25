# Function to forecast functional  groups at all NEON sites, using parameters estimated from model samples
#
##### Use Nmc samples to make predictions, returns a dataframe with CIs and observed truth values (plot means)
# Nmc = 1000
# drop_other = T
# predict_site_effects = pred_effects
# model.inputs = full.ts.model.inputs
# model_id="cycl_only_chitinolytic_20151101_20180101"
#
# Nmc = 1000
# predict_site_effects = NULL
# model.inputs = full.ts.model.inputs
#' @title fg_fcast_beta
#' @description Forecast functional groups at NEON plots and sites, using beta regression output
#' @export
fcast_logit_beta <- function(plotID,
													model.inputs,
													param_samples,
													truth.plot.long,
													plot_summary,
													Nmc = 1000,
													drop_other = T,
													predict_site_effects = NULL,
													rank.name=NULL,
													model_id,
													...) {



	siteID <- substr(plotID, 1, 4)

	# Debug: Check if site exists in model inputs
	if (!siteID %in% names(model.inputs$site_start)) {
		message("Warning: Site ", siteID, " not found in model.inputs$site_start")
		message("Available sites: ", paste(names(model.inputs$site_start), collapse = ", "))
		return(NULL)
	}

	# Prep MCMC sampling IDs
	# Initial condition uncertainty
	#	ic <- truncnorm::rtruncnorm(Nmc, mean = .5, sd = .2, a = 0, b = 1)

	# Check if we need to set initial condition
	if (!exists("ic")) {
		# This happens for sites not in original model data
		ic <- truncnorm::rtruncnorm(Nmc, mean(plot_summary$truth, na.rm=T), sd(plot_summary$truth, na.rm=T), a = 0, b = 1)
	}
	
	# Validate initial condition
	if (any(is.na(ic)) || any(is.infinite(ic))) {
		message("Warning: Invalid initial condition for plot ", plotID, " - skipping")
		return(NULL)
	}

	NT = model.inputs$N.date
	Nmc_large <- max(nrow(param_samples)) #20000 # Larger sample number for covariate/IC set of values
	#Nmc_large <- 20000
	# Remove verbose print statement
	# print(Nmc_large)
	# Handle case where we want more samples than available
	if (Nmc > Nmc_large) {
		# If we want more samples than available, sample with replacement
		row_samples <- sample.int(Nmc_large, Nmc, replace = TRUE)
	} else {
		# If we have enough samples, sample without replacement
		row_samples <- sample.int(Nmc_large, Nmc, replace = FALSE)
	}

	date_key <- model.inputs$truth.plot.long %>%
		select(dateID, date_num) %>% distinct()

	#Sample covariate data - FAIL FAST: This will stop execution if data is missing
	tryCatch({
		covar <- create_covariate_samples(model.inputs, plotID, siteID,
																		Nmc_large, Nmc)
	}, error = function(e) {
		message("ERROR in create_covariate_samples for plot ", plotID, " site ", siteID, ":")
		message("  ", e$message)
		message("  This indicates missing or invalid data that needs investigation.")
		message("  Check:")
		message("    1. Site mapping between calibration and validation periods")
		message("    2. Array dimensions and rownames")
		message("    3. Data quality (NA, infinite values)")
		message("    4. Time period compatibility")
		stop("Hindcasting stopped due to data issues - investigate and fix root cause")
	})
	
	# Validation - covar should always be valid now since function fails fast
	if (!is.array(covar) || length(dim(covar)) != 3) {
		stop("CRITICAL ERROR: create_covariate_samples returned invalid array - this should never happen!")
	}
	
	# Validate we have enough covariates for the model type
	expected_covars <- if (model_name == "cycl_only") 2 else if (model_name == "env_cov") 6 else 8
	if (dim(covar)[2] < expected_covars) {
		stop("CRITICAL ERROR: covar has ", dim(covar)[2], " covariates but model expects ", expected_covars, " for plot ", plotID)
	}

	plot_obs <- model.inputs$truth.plot.long %>%
		filter(plotID==!!plotID) %>%
		select(-c(plot_num,site_num))

	taxon_names <- model.inputs$truth.plot.long %>% select(species) %>% distinct() %>% unlist()
	taxon_name <- taxon_names[1]  # Use first taxon name for output

	# Check whether there's already an estimated site effect. If not, we'll sample!
	is_new_site <- ifelse(siteID %in% truth.plot.long$siteID, FALSE, TRUE)

	if (is_new_site) {
		# For new sites, we need to predict site effects
		if (!is.null(predict_site_effects)) {
			# Use predicted site effects if available
			site_effect <- predict_site_effects$fit
			# Ensure we have the right number of samples
			if (length(site_effect) == 1) {
				site_effect <- rep(site_effect, Nmc)
			} else if (length(site_effect) != Nmc) {
				# Sample from the available predictions
				site_effect <- sample(site_effect, Nmc, replace = TRUE)
			}
		} else {
			# Sample from site effect variance
			site_effect_sd <- param_samples[row_samples,] %>% select(grep("site_effect_sd", colnames(.))) %>% unlist()
			site_effect_sd <- mean(site_effect_sd, na.rm = TRUE)
			site_effect <- rnorm(Nmc, 0, site_effect_sd)
		}

		# Take initial condition & start forecast from mean observed value if possible
		plot_start_date <- model.inputs$plot_start[plotID]
		if (is.na(plot_start_date)) {
			plot_start_date <- 1  # Default to first timepoint
		}

	} else {
		site_num <- unique(truth.plot.long[truth.plot.long$siteID==siteID,]$site_num)
		site_param <- paste0("site_effect[", site_num, "]")
		site_effect <- 	param_samples[row_samples,] %>% select(!!site_param) %>% unlist()

		# add model estimates if possible
		plot_est <- plot_summary %>%
			filter(plotID == !!plotID) %>%
			select(-c(plot_num, site_num, dateID, dates, truth, rank))
		plot_obs <- left_join(plot_obs, plot_est, by = intersect(colnames(plot_obs), colnames(plot_est)))

		# Take initial condition & start forecast from last observed value if possible
		last_obs <- plot_est %>% filter(timepoint==max(timepoint))
		plot_start_date <- last_obs$timepoint
		ic <- last_obs$`50%`
	}

	# CRITICAL FIX: Validate plot_start_date before using it
	if (is.na(plot_start_date) || is.null(plot_start_date) || plot_start_date < 1) {
		message("Warning: Invalid plot_start_date ", plot_start_date, " for plot ", plotID, " - using default")
		plot_start_date <- 1
	}

	# Get NT from model.inputs for bounds checking
	NT <- tryCatch({
		if (is.null(model.inputs$N.date) || !is.numeric(model.inputs$N.date) || model.inputs$N.date < 1) {
			warning("Invalid N.date in model.inputs, using default")
			12
		} else {
			model.inputs$N.date
		}
	}, error = function(e) {
		warning("Error accessing N.date: ", e$message, " - using default")
		12
	})

	# Ensure plot_start_date is within bounds
	if (plot_start_date > NT) {
		message("Warning: plot_start_date ", plot_start_date, " exceeds N.date ", NT, " for plot ", plotID, " - adjusting")
		plot_start_date <- 1
	}

	if (length(site_effect)==0) {
		message("No site effect found!")
		return(NULL)  # Skip this plot entirely
	}
	
	# Check for valid site effect values
	if (any(is.na(site_effect))) {
		message("Warning: NA values in site effects for plot ", plotID, " - skipping")
		return(NULL)
	}

	### Get other parameter estimates
	### Rho
	rho <- param_samples[row_samples,] %>% select(grep("rho", colnames(.))) %>% unlist()
	### Betas
	betas <- param_samples[row_samples,] %>% select(grep("beta", colnames(.)))
	### Intercept
	intercept <- param_samples[row_samples,] %>% select(grep("intercept", colnames(.))) %>% unlist()
	### Process error (precision parameter)
	precision <- param_samples[row_samples,] %>% select(grep("precision", colnames(.))) %>% unlist()
	### Legacy effect (if using legacy covariate)
	legacy_effect <- param_samples[row_samples,] %>% select(grep("legacy_effect", colnames(.))) %>% unlist()
	# If no legacy effect parameter found, set to 0
	if (length(legacy_effect) == 0) {
		legacy_effect <- rep(0, Nmc)
	}
#	precision <- lapply(precision_samp, function(x) lapply(x, prec_to_sd)) %>% unlist()

	# Validate that we have all required parameters
	if (length(rho) == 0 || length(betas) == 0 || length(intercept) == 0 || length(precision) == 0) {
		message("Warning: Missing required parameters for plot ", plotID, " - skipping")
		return(NULL)
	}
	
	# Check for valid parameter values
	if (any(is.na(rho)) || any(is.na(betas)) || any(is.na(intercept)) || any(is.na(precision))) {
		message("Warning: NA values in parameters for plot ", plotID, " - skipping")
		return(NULL)
	}

	# Parse model_id to determine model type and legacy covariate usage
	model_info <- tryCatch({
		parse_model_id(model_id)
	}, error = function(e) {
		message("Warning: Error parsing model_id '", model_id, "': ", e$message)
		# Fallback parsing
		info <- strsplit(model_id, "_")[[1]]
		if (length(info) >= 2) {
			if (info[1] == "cycl" && info[2] == "only") {
				model_name <- "cycl_only"
			} else if (info[1] == "env" && info[2] == "cov") {
				model_name <- "env_cov"
			} else if (info[1] == "env" && info[2] == "cycl") {
				model_name <- "env_cycl"
			} else {
				model_name <- info[1]
			}
		} else {
			model_name <- "unknown"
		}
		list(model_name = model_name)
	})
	
	# Extract model name safely
	if (is.list(model_info) && length(model_info) >= 6) {
		model_name <- model_info[[6]]
	} else if (is.list(model_info) && "model_name" %in% names(model_info)) {
		model_name <- model_info$model_name
	} else {
		model_name <- "unknown"
	}
	
	use_legacy_covariate <- grepl("with_legacy_covariate", model_id)
	
	# Handle covariate selection based on model type
	if (model_name == "cycl_only") {
		# Cyclical only: just sin/cosine terms (indices 7,8)
		covar <- covar[,c(7,8),]
	} else if (model_name == "env_cov") {
		# Environmental covariates only: temp, mois, pH, pC, relEM, LAI (indices 1:6)
		covar <- covar[,c(1:6),]
	} else if (model_name == "env_cycl") {
		# Environmental + cyclical: all 8 covariates (indices 1:8)
		covar <- covar[,c(1:8),]
	} else {
		# Default fallback: assume all covariates
		message("Warning: Unknown model type '", model_name, "', using all covariates")
	}

	all_tax_abs <- matrix(NA, Nmc, NT)
	# Initialize first timepoint with initial condition
	all_tax_abs[,plot_start_date] <- ic
	## simulate
	x <- ic

	# Only forecast for valid timepoints (where we have data or can reasonably predict)
	valid_timepoints <- plot_start_date:NT
	
	# CRITICAL FIX: Ensure we don't try to forecast beyond available data
	if (max(valid_timepoints) > NT) {
		valid_timepoints <- valid_timepoints[valid_timepoints <= NT]
	}
	
	# plot_start_date bounds checking already done earlier in the function

	for (time in valid_timepoints) {
		if (time == plot_start_date) {
			# Skip the initial condition timepoint
			next
		}
		
		# CRITICAL FIX: Bounds checking for covariate access
		if (time > dim(covar)[3]) {
			message("Warning: Time ", time, " exceeds covariate array bounds for plot ", plotID, " - skipping")
			next
		}
		
		Z  <- covar[, ,time]
		
		# CRITICAL FIX: Validate Z dimensions
		if (ncol(Z) != ncol(betas)) {
			message("Warning: Covariate dimension mismatch at time ", time, " for plot ", plotID, " - skipping")
			next
		}
		
		Ex <- rho * logit(x) + apply(Z * betas, 1, sum) + site_effect + intercept
		
		# Add legacy effect if using legacy covariate
		if (use_legacy_covariate) {
			# Create legacy covariate using the same logic as model fitting
			# Get the date for this timepoint
			date_key_row <- date_key %>% filter(date_num == time)
			if (nrow(date_key_row) > 0) {
				date_val <- date_key_row$dateID
				# Convert to date and check if it's in legacy period (2013-06-27 to 2015-11-30)
				date_obj <- as.Date(date_val, format = "%Y%m%d")
				legacy_start <- as.Date("2013-06-27")
				legacy_end <- as.Date("2015-11-30")
				legacy_indicator <- as.numeric(date_obj >= legacy_start & date_obj <= legacy_end)
				Ex <- Ex + legacy_effect * legacy_indicator
			}
		}
		
		# norm_sample <- function(x, y) truncnorm::rtruncnorm((1, x, y)
		# x <- unlist(Map(norm_sample, mu, precision))

		mu = invlogit(Ex)

		#prevent NAs using interval transform:
		# transformed_mu = interval_transform(cbind(mu, 1-mu))

		to_transform = cbind(mu, 1-mu)
		N = 200 # Making the interval a little squeezier since values get too close to 1
		N = 2000 # Making the interval a little squeezier since values get too close to 1
		N = Nmc
		transformed_mu <- (to_transform * (N - 1) + .5)/N

		mu = transformed_mu[,1]
		# Add process error using logit beta
		alpha = mu * ((mu * (1 - mu))/precision^2 - 1)
		beta = (1 - mu) * ((mu * (1 - mu))/precision^2 - 1)

		# logit_mu = rLogitBeta(1, alpha, beta)
		# x = invlogit(logit_mu)

		alpha_beta = cbind(alpha, beta)
		logit_mu_vec = lapply(1:Nmc, function(x) rLogitBeta(1, alpha_beta[x,1], alpha_beta[x,2])) %>% unlist()
		x = invlogit(logit_mu_vec)

		# Save to array
		all_tax_abs[,time] <- x
	}


	# Create output dataframe
	predict <- all_tax_abs
	ci <- as.data.frame(t(apply(predict, 2, quantile, c(0.025,.25,0.5,.75,0.975), na.rm=T))) %>%
		mutate(mean = apply(predict, 2, mean, na.rm=T),
												sd = apply(predict, 2, sd, na.rm=T),
												date_num = as.numeric(1:NT),
												plotID = plotID,
												siteID = siteID,
												taxon_name = taxon_name,
												#rank = rank.name,
												species = taxon_name,
												taxon = taxon_name,
												taxon_name = taxon_name,
												new_site = ifelse(is_new_site, T, F))
	colnames(ci)[1:5] <- c("lo","lo_25","med","hi_75","hi")

	ci <- left_join(ci, date_key, by=c("date_num"))
	ci$dates <- fixDate(ci$dateID)
	ci <- ci %>% filter(!is.na(med))

	coalesce_by_column <- function(df) {
		return(coalesce(df[1], df[2]))
	}

	# Add on truth values and/or model estimates
	if (!is_new_site){
		plot_obs_simple = plot_obs %>%
			select(siteID,plotID, dateID,species,
							 truth,lo,lo_25,med,hi_75,hi,
							 mean = Mean, sd = SD) %>% mutate(truth=as.numeric(truth))
		ci <- full_join(ci, plot_obs_simple, by = intersect(colnames(ci), colnames(plot_obs_simple)))
		ci = ci %>%
			group_by(dateID) %>%
			summarise_all(coalesce_by_column)
	} else {
		# For new sites, we don't have truth values, so just return the CI
	}

	# Add model metadata
	ci <- ci %>% mutate(
		model_id = model_id,
		model_name = model_name,
		use_legacy_covariate = use_legacy_covariate,
		legacy_effect_mean = mean(legacy_effect, na.rm = TRUE)
	)

	return(ci)
}



#
# ggplot(ci) +
# 	facet_grid(rows=vars(species), drop=T, scales="free") +
# 	geom_line(aes(x = dates, y = med), show.legend = F) +
# 	geom_ribbon(aes(x = dates, ymin = lo, ymax = hi), alpha=0.6, fill="blue") +
# 	geom_ribbon(aes(x = dates, ymin = lo_25, ymax = hi_75),fill="red", alpha=0.6) +
# 	theme_bw()+
# 	scale_fill_brewer(palette = "Paired") +
# 	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
# 				legend.position = "bottom",legend.title = element_text(NULL),
# 				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) #+
# 	#geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='')
