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

	# Prep MCMC sampling IDs
	# Initial condition uncertainty
	ic <- truncnorm::rtruncnorm(Nmc, mean = .5, sd = .2, a = 0, b = 1)
	NT = model.inputs$N.date
	Nmc_large <- max(nrow(param_samples)) #20000 # Larger sample number for covariate/IC set of values
	#Nmc_large <- 20000
	print(Nmc_large)
	row_samples <- sample.int(Nmc_large,Nmc)

	date_key <- model.inputs$truth.plot.long %>%
		select(dateID, date_num) %>% distinct()

	#Sample covariate data
	covar <- create_covariate_samples(model.inputs, plotID, siteID,
																		Nmc_large, Nmc)
	plot_obs <- model.inputs$truth.plot.long %>%
		filter(plotID==!!plotID) %>%
		select(-c(plot_num,site_num))# %>% rename(species = name)


	taxon_names <- model.inputs$truth.plot.long %>% select(species) %>% distinct() %>% unlist()


	all_tax_abs <- matrix(nrow = Nmc, ncol = NT)


		taxon_name <- taxon_names[1]
		print(paste0("Forecasting for taxon: ", taxon_name))

		# Check whether there's already an estimated site effect. If not, we'll sample!
		is_new_site <- ifelse(siteID %in% truth.plot.long$siteID, FALSE, TRUE)

		if (is_new_site) {


			if (!is.null(predict_site_effects)){
				message("Using predicted site effect")

				pred_site_taxon = predict_site_effects %>%
					filter(siteID==!!siteID & model_id ==!!model_id)
				new_site_effect <- data.frame(rnorm(Nmc, pred_site_taxon$fit, pred_site_taxon$se.fit))
				#new_site_effect <- data.frame(runif(Nmc, pred_site_taxon$ci_lo, pred_site_taxon$ci_hi))
				site_effect <- unlist(new_site_effect)

			} else {
				message("Sampling random site effect")
				# Sample from site effect variance
				site_effect_sd <- param_samples[row_samples,] %>%
					select(grep("sig$", colnames(.))) %>% unlist() %>% mean()
				new_site_effect <- data.frame(rnorm(Nmc, 0, site_effect_sd))
				site_effect <- unlist(new_site_effect)

			}

			# # Take initial condition & start forecast from mean observed value if possible
			# plot_start_date <- model.inputs$plot_index[plotID]

			# No! For new sites, begin forecast the timepoint before the first observation, still with random ic
			plot_start_date <- model.inputs$plot_start[plotID] - 1

		} else {

			message("Using known site effect")

			site_num <- unique(truth.plot.long[truth.plot.long$siteID==siteID,]$site_num)
			site_param <- paste0("site_effect[", site_num, "]")
			site_effect <- 	param_samples[row_samples,] %>% select(!!site_param) %>% unlist()

				# add model estimates if possible
				plot_est <- plot_summary %>%
					filter(plotID == !!plotID) %>%
					select(-c(plot_num, site_num, dateID, dates, truth, rank)) %>%
					dplyr::rename(lo = `2.5%`,
								 hi = `97.5%`,
								 lo_25 = `25%`,
								 hi_75 = `75%`,
								 med = `50%`)
				plot_obs <- left_join(plot_obs, plot_est,
															by = intersect(colnames(plot_obs), colnames(plot_est))) %>% select(-taxon)

				# Take initial condition & start forecast from last observed value if possible
				last_obs <- plot_est %>% filter(timepoint==max(timepoint))
				plot_start_date <- last_obs$timepoint
				print(last_obs)
				ic <- rep(last_obs$med, Nmc)
				plot_obs <- plot_obs  %>% select(-c(timepoint,date_num))
				if (length(ic)==0) {
					message("No initial condition found!")
				} else message("Using initial condition: ", ic[[1]])

		}

		if (length(site_effect)==0) {
			message("No site effect found!")
		}

		### Get other parameter estimates
		### Rho
		rho <- param_samples[row_samples,] %>% select(grep("rho", colnames(.))) %>% unlist()
		### Betas
		betas <- param_samples[row_samples,] %>% select(grep("beta", colnames(.)))
		### Intercept
		intercept <- param_samples[row_samples,] %>% select(grep("intercept", colnames(.))) %>% unlist()
		### Process error
		sigma <- param_samples[row_samples,] %>% select(grep("sigma", colnames(.))) %>% unlist()
	#	sigma <- lapply(sigma_samp, function(x) lapply(x, prec_to_sd)) %>% unlist()



		# If the model only had sin/cosine, remove the other covariate data
		if (model_name=="cycl_only"){
			covar <- covar[,c(7,8),]
		} else if (model_name=="env_cov"){
			covar <- covar[,c(1:6),]
		}

		predict <- matrix(NA, Nmc, NT)
		## simulate
		x <- ic


		for (time in (plot_start_date+1):NT) {
			Z  <- covar[, ,time]
			Ex <- rho * logit(x) + apply(Z * betas, 1, sum) + site_effect + intercept
			# norm_sample <- function(x, y) truncnorm::rtruncnorm((1, x, y)
			# x <- unlist(Map(norm_sample, mu, sigma))

			mu = invlogit(Ex)

			#prevent NAs using interval transform:
			# transformed_mu = interval_transform(cbind(mu, 1-mu))

			to_transform = cbind(mu, 1-mu)
			N = 200 # Making the interval a little squeezier since values get too close to 1
			N = 2000 # Making the interval a little squeezier since values get too close to 1
			transformed_mu <- (to_transform * (N - 1) + .5)/N

			mu = transformed_mu[,1]
			# Add process error using logit beta
			alpha = mu * ((mu * (1 - mu))/sigma^2 - 1)
			beta = (1 - mu) * ((mu * (1 - mu))/sigma^2 - 1)

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
			plot_obs_simple = plot_obs %>%
				select(siteID,plotID, dateID,species,
							 truth) %>% mutate(truth=as.numeric(truth))
			ci <- full_join(ci, plot_obs_simple, by = intersect(colnames(ci), colnames(plot_obs_simple)))
			#ci <- merge(ci, plot_obs, by = c("plotID", "siteID", "species", "dateID" ))
			ci = ci %>%
				group_by(dateID) %>%
				summarise_all(coalesce_by_column)
		}

		ci$date_num <- coalesce(ci$date_num, date_key$date_num)

		obs_date_num = ci[which(!is.na(ci$truth) & !is.na(ci$date_num)),]$date_num
		obs_date_num = obs_date_num[which(obs_date_num > plot_start_date)]
		if (length(obs_date_num) > 0){
		crps_df = data.frame()
		for (i in obs_date_num){
			obs = ci[which(ci$date_num==i),]$truth
			# mu = ci[which(ci$date_num==i),]$mean
			# sd = ci[which(ci$date_num==i),]$sd
			X = predict[,i]
			X = X[complete.cases(X)]
			crps_val = cbind(crps = crps_sample(obs, X), date_num = i)
			crps_df = rbind(crps_df, crps_val)
		}
		ci = left_join(ci, crps_df)
		}
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