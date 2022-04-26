
# Create forecasts for taxonomic groups, using structure from SOBOL code
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
source("./functions/prepTaxonomicData.r")
source("./functions/forecastTaxonomic.r")

pacman::p_load(readxl, rjags, Rfast, moments, scales, data.table, doParallel)

# Set some global parameters
Nmc_AB <- 5000 # Number of samples for subset
Nmc <- 5000 
N.beta = 8 
k = 10

k = 1

# Read in microbial abundances
cal <- c(readRDS("./data/clean/cal_groupAbundances_16S_2021.rds"), 
				 readRDS("./data/clean/cal_groupAbundances_ITS_2021.rds"))
val <- c(readRDS("./data/clean/val_groupAbundances_16S_2021.rds"), 
				 readRDS("./data/clean/val_groupAbundances_ITS_2021.rds"))

summaries <- readRDS("./data/summary/taxa_summaries.rds")
# Loop through each model

max.date = "20200101"

model_name = "cycl_only"
#model_name = "all_covariates"
scenario = "full_uncertainty"


cl <- makeCluster(10, type="FORK", outfile="")
registerDoParallel(cl)

#rank_output_list <- list()

#Run for multiple chains, in parallel (via PSOCK)
#for (k in 2:10){
rank_output_list = foreach(k=1:10, .errorhandling = 'pass') %dopar% {
	rank.name <- tax_names[[k]]
	message("Beginning forecast loop for: ", rank.name)
	
	cal.rank.df <- cal[[rank.name]] 
	val.rank.df <- val[[rank.name]] 
	rank.df <- rbind(cal.rank.df, val.rank.df)	
	rank.df <- rank.df[!rank.df$siteID %in% c("ABBY","LAJA"),]
	
	# Grab names of taxa to keep - this code is ugly but oh well
	keep_list <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/converged_taxa_list.rds")
	keep_names <- keep_list[[rank.name]]$taxon.name
	keep_vec <- c(keep_names, "siteID", "plotID", "dateID", "sampleID", "dates", "plot_date")
	colnames(rank.df) <- lapply(strsplit(colnames(rank.df), "\\."), "[[", 1)
	rank.df_spec <- rank.df[,colnames(rank.df) %in% keep_vec] 
	tots <- rowSums(rank.df_spec[,7:ncol(rank.df_spec)])
	rank.df_spec$other <- 1-tots
	# Remove samples where more than 99% of reads are "Other"
	rank.df_spec <- rank.df_spec %>% filter(tots > .01)
	
	# Prep validation data
	model.inputs <- prepTaxonomicData(rank.df = rank.df_spec, min.prev = 3, max.date = max.date,	full_timeseries = T)
	
	
	model_output_list <- list()
	for (model_name in c("all_covariates", "cycl_only")){
		message("Forecasting with model: ", model_name)
		
		# Filter model estimates for each plot abundance
		plot_summary <- summaries$plot_est %>% filter(rank == rank.name &
																										model_name == !!model_name &
																										scenario == !!scenario & 
																										time_period == "calibration")
		
		# Get model outputs
		f <- file.path("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/taxa/", model_name,  paste0("/calibration_samples_", rank.name, "_full_uncertainty_summary.rds"))
		read_in <- readRDS(f)
		param_samples_orig <- read_in$samples
		model.dat <- read_in$metadata$model_data
		truth.plot.long <- model.dat
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
		#siteID = new_site_list[[1]]
		#full_site_list <- c(site_list[[1]], new_site_list[[1]])
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
			plotID <- plot_list[[1]]
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
																	rep(model.inputs$y_cos[time], Nmc_large))
				}
				
				covar <- covar_full[row_samples,,]
				#go for it!!!
				hindcast.plot <- taxa_fcast(plotID,
																		covar,
																		param_samples,
																		ic,
																		truth.plot.long,
																		Nmc = 5000,  plot_summary,plot_start_date)
				
				hindcast.plot$dateID <- date_key[match(hindcast.plot$date_num, date_key$date_num),]$dateID
				hindcast.plot <- hindcast.plot %>% mutate(dates = fixDate(dateID),
																									model_name = !!model_name,
																									time_period = "calibration")
				plot_output_list[[plotID]] <- hindcast.plot
				
			}
			site_output_list[[siteID]] <- rbindlist(plot_output_list)	
		}
		model_output_list[[model_name]] <- rbindlist(site_output_list, fill = T)	
	}
	#rank_output_list[[rank.name]] <- rbindlist(model_output_list)	
	rank_output <- rbindlist(model_output_list)	
	#rank_output_list[[rank.name]] = rank_output
	return(rank_output)
}				


keep <- list()									
for (i in 1:10){
	if (is.data.frame(rank_output_list[[i]])){
		keep[[i]] <- T
	} else keep[[i]] <- F
}			
rank_output_save <- rank_output_list[unlist(keep)] 
all_out <- rbindlist(rank_output_save, fill = T)	
all_out$group <- ifelse(grepl("_bac", all_out$rank, fixed = T), "16S", "ITS")
saveRDS(all_out, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/hindcast_tax_test.rds")

all_out <- readRDS("./data/summary/hindcast_tax_test.rds")


ggplot(all_out %>% filter(plotID=="BART_002" & rank == "phylum_fun")) + 
	facet_grid(rows=vars(taxon), cols = vars(model_name), drop=T, scales="free") +
	geom_line(aes(x = dates, y = mean), show.legend = F, linetype=2) +
	geom_line(aes(x = dates, y = `50%`), show.legend = F) +
	geom_ribbon(aes(x = dates, ymin = lo, ymax = hi), alpha=0.6, fill="blue") +
	geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill="red", alpha=0.6) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='')
