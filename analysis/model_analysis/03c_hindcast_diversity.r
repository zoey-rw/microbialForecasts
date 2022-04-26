
# Create forecasts for functional groups, using structure from SOBOL code
#source("./source.R")
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
source("./functions/prepDiversityData.r")
source("./functions/forecastDiversity.r")

pacman::p_load(readxl, rjags, Rfast, moments, scales, data.table, doParallel)

# Set some global parameters
Nmc_AB <- 5000 # Number of samples for subset
Nmc <- 5000 
N.beta = 8 

summaries <- readRDS("./data/summary/div_summaries.rds")

# for testing
model_name = "cycl_only"
#model_name = "all_covariates"
scenario = "full_uncertainty_ITS"
scenario = "full_uncertainty_16S"

#Run for multiple groups at once, in parallel (via PSOCK)
cl <- makeCluster(2, type="FORK", outfile="")
registerDoParallel(cl)

#for (scenario in c("full_uncertainty_ITS","full_uncertainty_16S")) {
# Loop through each group
output.list = foreach(scenario = c("full_uncertainty_ITS","full_uncertainty_16S"), 
											.errorhandling = 'pass') %dopar% {
												
																								
	message("Beginning forecast loop for: ", scenario)
	
	# Read in observations and covariate data
	
	if (grepl("ITS", scenario)){
		group <- "ITS"
		div_in <- readRDS("./data/clean/alpha_div_ITS.rds")
	} else {
		group <- "16S"
		div_in <- readRDS("./data/clean/alpha_div_16S.rds")
	}
	# Prep validation data
	rank.df = div_in$full %>% filter(!siteID %in% c("ABBY","LAJA"))
	model.inputs <- prepDivData(rank.df = rank.df, min.prev = 3, max.date = "20200101", full_timeseries = T)

	# Loop through linear model versions
	model_output_list <- list()
	for (model_name in c("all_covariates", "cycl_only")){
		message("Forecasting with model: ", model_name)
		
		# Filter model estimates for each plot abundance
		plot_summary <- summaries$plot_est %>% filter(time_period == "calibration" & 
																										model_name == !!model_name &
																										scenario == !!scenario)
		
		# Get actual model MCMC samples & inputs
		f <- file.path("./data/model_outputs/diversity/", model_name, paste0("/calibration_samples_div_", scenario, ".rds"))
		read_in <- readRDS(f)
		param_samples_orig <- read_in$samples
		truth.plot.long <- model.dat <- read_in$metadata$model_data
		plot_site_key <- model.dat %>% select(siteID, plotID, dateID, date_num, plot_num, site_num) %>% distinct()
		site_list <- unique(plot_site_key$siteID)
		
		# Use new model inputs for full date, site, and plot keys
		date_key <- model.inputs$truth.plot.long %>% select(dateID, date_num) %>% distinct()
		new_plot_site_key <- model.inputs$truth.plot.long %>% select(siteID, plotID, dateID, date_num, plot_num, site_num) %>% 
			distinct() %>% filter(!siteID %in% plot_site_key$siteID)
		new_site_list <- unique(new_plot_site_key$siteID)
		
		# Prep MCMC sampling IDs
		param_samples <- as.data.frame(as.matrix(param_samples_orig))
		Nmc_large <- max(nrow(param_samples)) #20000 # Larger sample number for covariate/IC set of values
		row_samples <- sample.int(max(nrow(param_samples)),Nmc_AB)
		ic = Rfast::Rnorm(Nmc_AB,0,1) # Initial condition uncertainty
		
		full_site_list <- c(site_list, new_site_list)
		site_output_list <- list()
		
		siteID <- site_list[[1]] #testing
		siteID <- new_site_list[[1]] #testing
		for (siteID in full_site_list){
			message("SiteID: ", siteID)
			
			# Change based on each site
			start_date <- model.inputs$site_start[siteID]
			NT = model.inputs$N.date
			
			if (siteID %in% new_plot_site_key$siteID) {
				plot_key <- new_plot_site_key %>% filter(siteID == !!siteID)
				plot_list <- unique(plot_key$plotID)
			} else {
				plot_key <- plot_site_key %>% filter(siteID == !!siteID)
				plot_list <- unique(plot_key$plotID)
			}
			plot_output_list <- list()
			plotID <- plot_list[[1]] #testing
			for (plotID in plot_list){
				message("PlotID: ", plotID)
				#Sample covariate data
				covar_full <- array(NA, dim = c(Nmc_large, N.beta, NT))
				set.seed(1)
				for (time in start_date:NT) {
					covar_full[,,time] <- c(Rfast::Rnorm(Nmc_large, model.inputs$temp[siteID, time],
																							 model.inputs$temp_sd[siteID, time]),
																	Rfast::Rnorm(Nmc_large, model.inputs$mois[siteID, time],
																							 model.inputs$mois_sd[siteID, time]),
																	Rfast::Rnorm(Nmc_large, model.inputs$pH[plotID,],
																							 model.inputs$pH_sd[plotID,]),
																	Rfast::Rnorm(Nmc_large, model.inputs$pC[plotID,],
																							 model.inputs$pC_sd[plotID,]),
																	rep(model.inputs$relEM[plotID, time], Nmc_large),
																	rep(model.inputs$LAI[siteID, time], Nmc_large),
																	rep(model.inputs$y_sin[time], Nmc_large),
																	rep(model.inputs$y_cos[time], Nmc_large))}
				
				covar <- covar_full[row_samples,,]
				#go for it!!!
				hindcast.plot <- diversity_fcast(plotID, covar, param_samples, scenario,
																	ic, truth.plot.long, Nmc = 5000,  model.inputs,
																	plot_summary, plot_start_date, date_key)
				hindcast.plot <- hindcast.plot %>% mutate(dates = fixDate(dateID),
																									model_name = !!model_name,
																									scenario = !!scenario,
																									group = !!group,
																									time_period = "calibration")
				plot_output_list[[plotID]] <- hindcast.plot
			}
			site_output_list[[siteID]] <- rbindlist(plot_output_list)	
		}
		model_output_list[[model_name]] <- rbindlist(site_output_list, fill = T)	
	}
	group_output <- rbindlist(model_output_list)	
	#rank_output_list[[rank.name]] = rank_output
	#output.list[[2]] <- group_output
	return(group_output)
}				
out <- rbindlist(output.list)	
out$fcast_period <- ifelse(out$dates < "2017-01-01", "calibration", "hindcast")
saveRDS(out, "./data/summary/hindcast_div.rds")

out <- readRDS("./data/summary/hindcast_div.rds")

# Merging in some truth values that weren't added correctly. Not ideal but it works.
group <- "16S"
div_in <- readRDS("./data/clean/alpha_div_16S.rds")
rank.df = div_in$full %>% filter(!siteID %in% c("ABBY","LAJA"))
model.inputs_16S <- prepDivData(rank.df = rank.df, min.prev = 3, max.date = "20200101", full_timeseries = T)
truth_16S <- model.inputs_16S$truth.plot.long  %>% select(-c(plot_num, site_num, timepoint))
truth_16S$scenario = "full_uncertainty_16S"

group <- "ITS"
div_in <- readRDS("./data/clean/alpha_div_ITS.rds")
rank.df = div_in$full %>% filter(!siteID %in% c("ABBY","LAJA"))
model.inputs_ITS <- prepDivData(rank.df = rank.df, min.prev = 3, max.date = "20200101", full_timeseries = T)
truth_ITS <- model.inputs_ITS$truth.plot.long %>% select(-c(plot_num, site_num, timepoint))
truth_ITS$scenario = "full_uncertainty_ITS"

out <- out %>% select(-c(Shannon_orig, Shannon_scale_site, truth))
truth <- rbind(truth_16S, truth_ITS)
out_fixed <- left_join(out, truth)
saveRDS(out_fixed, "./data/summary/hindcast_div.rds")

# View example output
ggplot(out %>% filter(plotID=="BART_002" & scenario == "full_uncertainty_ITS")) + 
	facet_grid(#rows=vars(taxon), 
		cols = vars(model_name), drop=T, scales="free") +
	geom_line(aes(x = dates, y = med), show.legend = F, linetype=2) +
	geom_line(aes(x = dates, y = `50%`), show.legend = F) +
	geom_ribbon(aes(x = dates, ymin = lo, ymax = hi), alpha=0.6, fill="blue") +
	geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill="red", alpha=0.6) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='')

