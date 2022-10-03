# Function to forecast taxonomic groups at all NEON sites, using parameters estimated from model samples
#
##### Use Nmc samples to make predictions, returns a dataframe with CIs and observed truth values (plot means)

#' @title taxa_fcast
#' @description Forecast taxa at NEON plots and sites
#' @export
taxa_fcast <- function(plotID,
											 model.inputs,
											 param_samples,
											 truth.plot.long,
											 plot_summary,
											 Nmc = 1000, ...) {


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

	if (length(taxon_names) == 1) taxon_names <- c(taxon_names[[1]], "other")

	n.taxa <- length(taxon_names)
	all_tax_abs <- array(dim = c(Nmc, n.taxa, NT))

	plot_est_list <- list()
	i = 1
	for (i in 1:n.taxa){
		taxon_name <- taxon_names[i]
		print(paste0("Forecasting for taxon: ", taxon_name))

		# Check whether there's already an estimated site effect. If not, we'll sample!
		is_new_site <- ifelse(siteID %in% truth.plot.long$siteID, FALSE, TRUE)

		if (is_new_site) {

			# Sample from site effect variance
			site_effect_tau <- param_samples[row_samples,] %>% select(grep("sig$", colnames(.))) %>% unlist()
			# Convert precision to SD
			site_effect_tau <- unlist(lapply(site_effect_tau,
																			 function(x) lapply(x, prec_to_sd)))
			site_tau <- mean(site_effect_tau)
			new_site_effect <- data.frame(rnorm(Nmc, 0, site_tau))
			site_effect <- unlist(new_site_effect)

			# Take initial condition & start forecast from mean observed value if possible
			plot_start_date <- model.inputs$plot_index[plotID]

		} else {


			site_num <- unique(truth.plot.long[truth.plot.long$siteID==siteID,]$site_num)
			site_param <- paste0("site_effect[", site_num, ", ", i, "]")
			site_effect <- 	param_samples[row_samples,] %>% select(!!site_param) %>% unlist()

			if (taxon_name != "other"){
			# add model estimates if possible
			plot_est <- plot_summary %>%
				filter(plotID == !!plotID & rank_name == !!rank.name & species == taxon_name) %>%
				select(-c(plot_num, site_num, dateID, dates, truth, rank)) #%>% rename(species = name)
			plot_obs <- left_join(plot_obs, plot_est,
														by = intersect(colnames(plot_obs), colnames(plot_est))) %>% select(-taxon)
			#plot_est_list[[i]] <- plot_obs %>% select(-taxon)

			# Take initial condition & start forecast from last observed value if possible
			last_obs <- plot_est %>% filter(timepoint==max(timepoint))
			plot_start_date <- last_obs$timepoint
			print(last_obs)
			ic <- last_obs$`50%`
			} else {
				ic = 1 - last_obs$`50%`
			}
		}

		### Get other parameter estimates
		### Rho
		rho <- param_samples[row_samples,] %>% select(grep("rho", colnames(.))) %>%
			select(grep(paste0("[",i,"]"), colnames(.), fixed = T)) %>% unlist()
		### Betas
		betas <- param_samples[row_samples,] %>% select(grep("beta", colnames(.))) %>%
			select(grep(paste0("[",i,","), colnames(.), fixed = T))
		### Intercept
		intercept <- param_samples[row_samples,] %>% select(grep("intercept", colnames(.))) %>%
			select(grep(paste0("[",i,"]"), colnames(.), fixed = T)) %>% unlist()
		### Process error
		sigma_samp <- param_samples[row_samples,] %>% select(grep("sigma", colnames(.))) %>%
			select(grep(paste0("[",i,"]"), colnames(.), fixed = T)) %>% unlist()
		sigma <- lapply(sigma_samp, function(x) lapply(x, prec_to_sd)) %>% unlist()



		# If the model only had sin/cosine, remove the other covariate data
		if (ncol(betas)==2) {
			if(ncol(covar)==8) {
				covar <- covar[,c(7,8),]
			}
		}

		predict <- matrix(NA, Nmc, NT)
		## simulate
		x <- exp(ic)
		#x <- ic

		for (time in (plot_start_date+1):NT) {
			Z  <- covar[, ,time]
			mu <- rho * log(x) + apply(Z * betas, 1, sum) + site_effect + intercept
			# norm_sample <- function(x, y) truncnorm::rtruncnorm((1, x, y)
			# x <- unlist(Map(norm_sample, mu, sigma))

			# Truncated to prevent negative values
			x <- truncnorm::rtruncnorm(Nmc, mean = exp(mu), sd = sigma, a = 0, b = Inf)
			# Save to array
			all_tax_abs[,i,time] <- x
		}
	}

	# Make abundances relative within each MCMC sample
	all_tax_rel <- array(dim = c(Nmc, n.taxa, NT))
	for (time in (plot_start_date+1):NT) {
		all_tax_rel[,,time] <- all_tax_abs[,,time]/rowSums(all_tax_abs[,,time])
	}
	colnames(all_tax_rel) <- taxon_names

	# Create output dataframe
	all_tax_single_plot <- list()
	for (i in 1:n.taxa){
		taxon_name <- taxon_names[i]
		predict <- all_tax_rel[,i,]
		ci <- as.data.frame(t(apply(predict, 2, quantile, c(0.025,0.5,0.975), na.rm=T)))
		ci <- ci %>% mutate(mean = apply(predict, 2, mean, na.rm=T),
												sd = apply(predict, 2, sd, na.rm=T),
												date_num = as.numeric(1:NT),
												plotID = plotID,
												siteID = siteID,
												taxon_name = taxon_name,
												rank = rank.name,
												species = taxon_name,
												taxon = taxon_name,
												taxon_name = taxon_name,
												new_site = ifelse(is_new_site, T, F))
		colnames(ci)[1:3] <- c("lo","med","hi")


		ci <- left_join(ci, date_key, by=c("date_num"))
		ci$dates <- fixDate(ci$dateID)

		if (taxon_name != "other"){
		ci <- left_join(ci, plot_obs, by = intersect(colnames(ci), colnames(plot_obs)))
		}
		all_tax_single_plot[[i]] <- ci
	}
	ci_single_plot <-  plyr::rbind.fill(all_tax_single_plot)
	return(ci_single_plot)
}




# ggplot(ci_single_plot) +
# 	facet_grid(rows=vars(species), drop=T, scales="free") +
# 	geom_line(aes(x = dates, y = mean), show.legend = F, linetype=2) +
# 	geom_line(aes(x = dates, y = `50%`), show.legend = F) +
# 	geom_ribbon(aes(x = dates, ymin = lo, ymax = hi), alpha=0.6, fill="blue") +
# 	geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill="red", alpha=0.6) +
# 	theme_bw()+
# 	scale_fill_brewer(palette = "Paired") +
# 	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
# 				legend.position = "bottom",legend.title = element_text(NULL),
# 				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
# 	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='')