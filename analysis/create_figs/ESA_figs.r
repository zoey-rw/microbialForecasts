
library(hrbrthemes)
library(ggplot2)
library(tidyr)
library(coda)
library(dplyr)
library(tibble)

# # Read in samples for visualization
##### Load libraries and input/output data ----------------------

pacman::p_load(scoringRules, tidyr, dplyr, reshape2, parallel, 
							 lubridate, nimble, coda, tidyverse, runjags) 

read_in <- list()
read_in[[1]] <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/nolog_div_no_uncertainty_ITS.rds")
read_in[[2]] <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/nolog_div_spatial_uncertainty_ITS.rds")
read_in[[3]] <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/nolog_div_temporal_uncertainty_ITS.rds")
read_in[[4]] <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/nolog_div_full_uncertainty_ITS.rds")
read_in[[5]] <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/nolog_div_no_uncertainty_16S.rds")
read_in[[6]] <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/nolog_div_spatial_uncertainty_16S.rds")
read_in[[7]] <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/nolog_div_temporal_uncertainty_16S.rds")
read_in[[8]] <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/nolog_div_full_uncertainty_16S.rds")
names(read_in) <- c("no_uncertainty_ITS", "spatial_uncertainty_ITS", "temporal_uncertainty_ITS", 
										"full_uncertainty_ITS", "no_uncertainty_16S", "spatial_uncertainty_16S", 
										"temporal_uncertainty_16S", "full_uncertainty_16S")
read_in$sample.list <- lapply(read_in, "[[", 1)
read_in$param.summary.list <- lapply(read_in, "[[", 2)
read_in$metadata.list <- lapply(read_in, "[[", 3)
read_in$plot.summary.list <- lapply(read_in, "[[", 4)


div_in = switch(group,
								"ITS" = readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_ITS.rds"),
								"16S" = readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_16S.rds"))


rank.df = div_in$cal
rank.df$Shannon_scaled <- scale(rank.df$Shannon, scale = F)

# Custom function for organizing model data.
model.dat <- prepDivData(rank.df = rank.df, min.prev = 3)
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



data_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/div_summaries.rds")


#### PARAMETER ESTIMATES ####
# ## GGPLOT
sum.all <- data_in$summary_df

beta_out <- sum.all[which(!is.na(sum.all$beta)),]
# By rank with every taxon - cluttered
ggplot(data=beta_out,
			 aes(x = reorder(beta, effSize),y = effSize)) +
	facet_grid(rows = vars(group), cols = vars(scenario), drop = T) +
	geom_point(aes(shape = as.factor(significant), color = beta), size = 4) +
	labs(col = "Parameter", title = "Absolute effect size") + 
	xlab("Parameter")+ 
	ylab(NULL) +
	theme_bw() + theme(#axis.ticks.x=element_blank(),
		text = element_text(size = 16),
		axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
			angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=22,face="bold")
		#strip.text.y = element_text(size=24,hjust=0,vjust = 1,angle=180,face="bold")
	) + scale_shape_manual(values = c(21, 16), name = NULL, 
												 labels = c("Not significant","Significant")) 

# Only full-uncertainty scenario
ggplot(data=beta_out[beta_out$scenario=="full uncertainty",],
			 aes(x = reorder(beta, effSize),y = effSize)) +
	facet_grid(rows = vars(group), drop = T) +
	geom_point(aes(shape = as.factor(significant), color = beta), size = 4) +
	labs(col = "Parameter", title = "Absolute effect size") + 
	xlab("Parameter")+ 
	ylab(NULL) +
	theme_bw() + theme(#axis.ticks.x=element_blank(),
		text = element_text(size = 16),
		axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
			angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=22,face="bold")
		#strip.text.y = element_text(size=24,hjust=0,vjust = 1,angle=180,face="bold")
	) + scale_shape_manual(values = c(21, 16), name = NULL, 
												 labels = c("Not significant","Significant")) 


### FUNGI
## Violin plots of parameter estimates
samps <- read_in$sample.list$full_uncertainty_ITS
samps <- do.call(rbind.data.frame, samps[,grep("beta|rho", colnames(samps[[1]]))])
samps_long_fungi <- samps %>% 
	rownames_to_column("rowname") %>% pivot_longer(2:8) %>% 
	mutate(beta_num = as.numeric(gsub("beta\\[|\\]", "", name))) %>% 
	mutate(beta = recode(beta_num,
											 "1" = "Temperature",
											 "2" = "Moisture",
											 "3" = "pH",
											 "4" = "pC",
											 "5" = "Plant species richness",
											 "6" = "% grasses",
											 .missing = "Autocorrelation"),
				 pretty_group = "Fungi")
samps <- read_in$sample.list$full_uncertainty_16S
samps <- do.call(rbind.data.frame, samps[,grep("beta|rho", colnames(samps[[1]]))])
samps_long_bacteria <- samps %>% 
	rownames_to_column("rowname") %>% pivot_longer(2:8) %>% 
	mutate(beta_num = as.numeric(gsub("beta\\[|\\]", "", name))) %>% 
	mutate(beta = recode(beta_num,
											 "1" = "Temperature",
											 "2" = "Moisture",
											 "3" = "pH",
											 "4" = "pC",
											 "5" = "Plant species richness",
											 "6" = "% grasses",
											 .missing = "Autocorrelation"),
				 pretty_group = "Bacteria")
samps_long <- rbind(samps_long_fungi, samps_long_bacteria)



params_plot <- ggplot(data=samps_long,
			 aes(x = reorder(beta, value),y = value)) +
	geom_violin(aes(fill = beta), trim=FALSE, show.legend = F) + 
	facet_grid(rows=vars(pretty_group)) +
	ylab("Effect size") + 
	xlab("Parameter")+ ggtitle("Drivers of Shannon evenness") + theme_minimal(base_size = 16) + geom_hline(aes(yintercept=0), linetype=2) + 
	theme(
		axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
			angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=18,face="bold")
	) 
params_plot

# Without rho
params_no_rho <- ggplot(data=samps_long[samps_long$beta != "Autocorrelation",],
			 aes(x = reorder(beta, value),y = value)) +
	geom_violin(aes(fill = beta), trim=FALSE, show.legend = F) + 
	facet_grid(cols=vars(pretty_group)) +
	ylab("Effect size") + 
	xlab("Parameter")+ ggtitle("Drivers of Shannon evenness") + theme_minimal(base_size = 20) + geom_hline(aes(yintercept=0), linetype=2) + 
	theme(
		axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
			angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=18,face="bold")
	) 
 params_no_rho




### BACTERIA
## Violin plots of parameter estimates

bacteria_params_plot <- ggplot(data=samps_long,
			 aes(x = reorder(beta, value),y = value)) +
	geom_violin(aes(fill = beta), trim=FALSE) + 
	ylab("Effect size") + 
	xlab("Parameter")+ ggtitle("Drivers of Shannon evenness of soil bacteria") + theme_minimal(base_size = 16) + geom_hline(aes(yintercept=0), linetype=2) + 
	theme(
		axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
			angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=18,face="bold")
	) 
bacteria_params_plot

# Without rho
bacteria_params_no_rho <- ggplot(data=samps_long[samps_long$beta != "Autocorrelation",],
			 aes(x = reorder(beta, value),y = value)) +
	geom_violin(aes(fill = beta), trim=FALSE) + 
	ylab("Effect size") + 
	xlab("Parameter")+ ggtitle("Drivers of Shannon evenness of soil bacteria") + theme_minimal(base_size = 16) + geom_hline(aes(yintercept=0), linetype=2) + 
	theme(
		axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
			angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=18,face="bold")
	) 
bacteria_params_no_rho


fungi_params_plot
fungi_params_no_rho
bacteria_params_plot
bacteria_params_no_rho


#### end parameter estimates ------------------



#### Example time-series -----



rank.df 
model.dat 
constants 
predictor_data <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/all_predictor_data.rds")


mois <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/monthly_soil_moisture.rds")
# Soil moisture
mois_df <- mois %>% 
	select(siteID, date, low, hi, dateID=month, source, mois=moisture_out, mois_sd=moisture_sd_out) %>% 
	filter(siteID=="CPER")
mois_plot <- ggplot(mois_df, aes(x = date, y = mois, color = source)) + 
	geom_point(show.legend = F) + #facet_wrap(~plotID) +#, nrow = 8) + 
	geom_errorbar(aes(ymin=low,ymax=hi),alpha=0.8) + theme_minimal() + 
	theme(legend.position = c(0.55, 0.9)) +
	guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 5, alpha = 1))) + ylab("Soil moisture")
mois_plot




temp <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/monthly_soil_temperature.rds")
# Soil moisture
temp_df <- temp %>% 
	select(siteID, date, low, hi, dateID=month, source, temp=temperature_out, temp_sd=temperature_sd_out) %>% 
	filter(siteID=="CPER")
temp_plot <- ggplot(temp_df, aes(x = date, y = temp, color = source)) + 
	geom_point(show.legend = F) + #facet_wrap(~plotID) +#, nrow = 8) + 
	geom_errorbar(aes(ymin=low,ymax=hi),alpha=0.8) + theme_minimal() + 
	theme(legend.position = c(0.55, 0.9)) +
	guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 5, alpha = 1))) + ylab("Soil temperature") + ylim(c(-10,50)) + xlab(NULL)
temp_plot



dat_soil <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/soilSample_data_allsites.rds")
chem_df <- dat_soil %>% select(siteID, plotID, soilInCaClpH, organicCPercent, year)
pH_plot <- ggplot(chem_df[chem_df$plotID %in% c("CPER_006"),]) + 
	geom_jitter(aes(x = year, y = soilInCaClpH, color = year), height=0, size=5, width=.15, show.legend = F) + 
	#facet_grid(~siteID, scales = "free_x", drop = T) + 
	theme(axis.text.x = element_text(angle = 200)) + theme_minimal(base_size=24) + ylab("pH") + xlab(NULL) + ylim(c(5, 7))
pH_plot
####


#### Old site/new site plots




#### SHannon diversity
library(phyloseq)
recent_ps <- readRDS("/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/data/NEON_16S_phyloseq_subset.rds")
# shannon <- estimate_richness(recent_ps, measures = "Shannon")

low_shannon <- shannon[shannon$Shannon < 4,,drop=F]
hi_shannon <- shannon[shannon$Shannon > 6.5,,drop=F]

low_ps <- prune_samples(sample_names(recent_ps) %in% rownames(low_shannon),recent_ps)
hi_ps <- prune_samples(sample_names(recent_ps) %in% rownames(hi_shannon),recent_ps)


plot_bar(low_ps)

p = plot_bar(low_ps, "Family", fill="Family", facet_grid=~sampleID) + guides(color = FALSE, fill = FALSE)
p + geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack", show.legend = F) + theme(legend.position="bottom")
