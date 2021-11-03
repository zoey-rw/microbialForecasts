# Summarize outputs from diversity models
library(tidyverse)
library(coda)

# # Read in samples for visualization
# read_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/div_samples_min3.rds")


read_list <- list()
read_in <- list()
read_list[[1]] <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/refit_nolog_div_no_uncertainty_ITS.rds")
read_list[[2]] <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/refit_nolog_div_full_uncertainty_ITS.rds")
read_list[[3]] <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/refit_nolog_div_no_uncertainty_16S.rds")
read_list[[4]] <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/refit_nolog_div_full_uncertainty_16S.rds")
names(read_list) <- c("no_uncertainty_ITS", 
										"full_uncertainty_ITS", "no_uncertainty_16S", "full_uncertainty_16S")
read_in$sample.list <- lapply(read_list, "[[", 1)
read_in$param.summary.list <- lapply(read_list, "[[", 2)
read_in$metadata.list <- lapply(read_list, "[[", 3)
read_in$plot.summary.list <- lapply(read_list, "[[", 4)


# Create parameters to pass	
params = data.frame(index = 1:8,
										scenario = c("no_uncertainty_ITS", "spatial_uncertainty_ITS",
																 "temporal_uncertainty_ITS", "full_uncertainty_ITS",
																 "no_uncertainty_16S", "spatial_uncertainty_16S",
																 "temporal_uncertainty_16S", "full_uncertainty_16S"),
										group = c(rep("ITS", 4),rep("16S", 4)),
										temporalDriverUncertainty = c(F, F, T, T, F, F, T, T),
										spatialDriverUncertainty = c(F, T, F, T, F, T, F, T))

plot_est_df_all <- list()
summary_df_all <- list()
gelman_list <- list()
allplots.scores.list <- list()

for (scenario in names(read_in$sample.list)){
	
	cat(paste0("\nSummarizing ", scenario, "..."))
	samples <- read_in$sample.list[[scenario]]
	param_summary <- read_in$param.summary.list[[scenario]]
	plot_summary <- read_in$plot.summary.list[[scenario]]
	truth.plot.long <- read_in$metadata.list[[scenario]]$model_data 
	
	
	# Calculate plot means and actual values per rank
	plot_est <- plot_summary[[2]]
	pred.plot <- plot_est %>% as.data.frame() %>%  rownames_to_column() %>%
		separate(rowname, sep=", ", into=c("plot_num","timepoint")) %>%
		mutate(plot_num = as.integer(gsub("plot_mu\\[", "", plot_num)),
					 timepoint = as.integer(gsub("\\]", "", timepoint)),
					 rank = "bacterial_richness")
	truth.plot.long <- truth.plot.long %>% 
		mutate(dates = as.Date(paste0(dateID, "01"), "%Y%m%d"),
					 timepoint = as.integer(timepoint))
	allplots <- merge(truth.plot.long, pred.plot, by = c("plot_num","timepoint"), all=T)
	allplots$scenario <- scenario
	plot_est_df_all[[scenario]] <- allplots
	
	# For scoring the predictions (need mean and SD)
	pred.plot.scores <- plot_summary[[1]] %>% as.data.frame() %>%  rownames_to_column() %>%
		separate(rowname, sep=", ", into=c("plot_num","timepoint")) %>%
		mutate(plot_num = as.integer(gsub("plot_mu\\[", "", plot_num)),
					 timepoint = as.integer(gsub("\\]", "", timepoint)))
	allplots.scores <- merge(truth.plot.long, pred.plot.scores, by = c("plot_num","timepoint"), all=T)
	allplots.scores$scenario <- scenario
	allplots.scores.list[[scenario]] <- allplots.scores
	
	# Get mean values for parameters
	means <- param_summary[[1]]
	beta_out <- means[grep("beta", rownames(means)),]
	rho_out <- means[grep("rho", rownames(means)),,drop=F] %>% as.data.frame() %>% 
		rownames_to_column("rowname") %>% mutate(beta = "rho")
	intercept_out <- means[grep("intercept", rownames(means)),,drop=F] %>% as.data.frame() %>% 
		rownames_to_column("rowname")
	sigma_out <- means[grep("sigma", rownames(means)),,drop=F] %>% as.data.frame() %>% 
		rownames_to_column("rowname")
	sig_out <- means[grep("sig$", rownames(means)),,drop=F] %>% as.data.frame() %>% 
		rownames_to_column("rowname")
	core_sd_out <- means[grep("core", rownames(means)),,drop=F] %>% as.data.frame() %>% 
		rownames_to_column("rowname")
	
	
	# Get beta sizes per rank
	beta_out <-  beta_out %>% as.data.frame() %>% 
		rownames_to_column("rowname") %>% 
		mutate(beta_num = as.numeric(gsub("beta\\[|\\]", "", rowname))) %>% 
		mutate(beta = recode(beta_num,
												 "1" = "Temperature",
												 "2" = "Moisture",
												 "3" = "pH",
												 "4" = "pC",
												 "5" = "Plant species richness",
												 "6" = "Ectomycorrhizal trees"))
	
	# Use quantiles to assign significance to beta parameters.
	beta_ci <-  param_summary[[2]] %>% as.data.frame() %>% 
		rownames_to_column("rowname") %>% filter(grepl("beta", rowname)) %>% 
		mutate(beta_num = as.numeric(gsub("beta\\[|\\]", "", rowname)))
	
	beta_out$significant <- ifelse(beta_ci$`2.5%` < 0 & beta_ci$`97.5%` < 0 |
																 	beta_ci$`2.5%` > -0 & beta_ci$`97.5%` > -0,
																 1, 0)
	beta_out$effSize <- abs(beta_out$Mean)
	
	# Get site effect sizes per rank
	site_eff_out <- means %>% as.data.frame() %>% 
		rownames_to_column("rowname") %>% filter(grepl("site", rowname)) %>% 
		mutate(site_num = as.character(gsub("site_effect\\[|\\]", "", rowname))) %>% 
		mutate(siteID = truth.plot.long[match(site_num, truth.plot.long$site_num),]$siteID)
	
	summary_df <- plyr::rbind.fill(site_eff_out, beta_out, intercept_out, sigma_out, sig_out, rho_out, core_sd_out)
	summary_df$scenario <- scenario
	summary_df_all[[scenario]] <- summary_df
	
	
	## Calculate gelman diagnostics to assess convergence
	gd <- gelman.diag(samples, multivariate = FALSE)
	gelman_list[[scenario]] <- cbind(gd[[1]], effSize = effectiveSize(samples))
	
}

plot_est_df_all <- plyr::rbind.fill(plot_est_df_all) %>% 
	tidyr::separate(scenario, sep = "_", into = c("uncert1", "uncert2", "group"), remove = F) %>% 
	mutate(uncert = paste(uncert1, uncert2, sep = "_")) %>% select(-c(uncert1, uncert2)) %>% 
	mutate(fcast_type = "Diversity",  
				 fcast_period = "refit")

scores.list <- plyr::rbind.fill(allplots.scores.list) %>% 
	tidyr::separate(scenario, sep = "_", into = c("uncert1", "uncert2", "group"), remove = F) %>% 
	mutate(uncert = paste(uncert1, uncert2, sep = "_")) %>% select(-c(uncert1, uncert2)) %>% 
	mutate(fcast_type = "Diversity",  
				 fcast_period = "refit")

summary_df_all <- do.call(rbind, summary_df_all) %>% 
	tidyr::separate(scenario, sep = "_", into = c("uncert1", "uncert2", "group"), remove = F) %>% 
	mutate(uncert = paste(uncert1, uncert2, sep = "_")) %>% select(-c(uncert1, uncert2)) %>% 
	mutate(fcast_type = "Diversity",  
				 fcast_period = "refit")

saveRDS(list(plot_est = plot_est_df_all,
						 summary_df = summary_df_all,
						 samples = read_in$sample.list,
						 gelman_list = gelman_list,
						 scores.list = scores.list),
				"/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/refit_div_summaries.rds")


# Copy over outputs
sum.all <- summary_df_all

# Or read in outputs instead
data_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/refit_div_summaries.rds")
sum.all <- data_in$summary_df

# ## GGPLOT

beta_out <- sum.all[which(!is.na(sum.all$beta)),]
# By rank with every taxon - cluttered
ggplot(data=beta_out,
			 aes(x = reorder(beta, effSize),y = effSize)) +
	facet_grid(rows = vars(group), cols = vars(uncert), drop = T) +
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


## Violin plots of parameter estimates
samps <- read_in$sample.list$full_uncertainty_ITS
samps <- do.call(rbind.data.frame, samps[,grep("beta|rho", colnames(samps[[1]]))])
samps_long <- samps %>% rownames_to_column("rowname") %>% pivot_longer(2:8)

samps_long <- samps_long %>% 
	mutate(beta_num = gsub("beta\\[|\\]", "", name)) %>% 
	mutate(beta = recode(beta_num,
											 "1" = "Temperature",
											 "2" = "Moisture",
											 "3" = "pH",
											 "4" = "pC",
											 "5" = "Plant species richness",
											 "6" = "% ectomycorrhizal trees",
											 .missing = "Autocorrelation"))
ggplot(data=samps_long,
			 aes(x = reorder(beta, value),y = value)) +
	geom_violin(aes(fill = beta), trim=FALSE) + 
	ylab("Effect size") + 
	xlab("Parameter")+ ggtitle("Drivers of Shannon evenness of soil fungi") + theme_minimal(base_size = 16) + geom_hline(aes(yintercept=0), linetype=2) + 
	theme(
		axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
			angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=18,face="bold")
	) 


# Takes a long time to render...
ggplot(data=samps_long,
			 aes(x = reorder(beta, value),y = value)) +
	geom_jitter(aes(color = beta), width = .2, height = 0, alpha = .3) + 
	ylab("Effect size") + 
	xlab("Parameter")+ ggtitle("Drivers of Shannon evenness of soil fungi") + theme_minimal(base_size = 16) + geom_hline(aes(yintercept=0), linetype=2) + 
	theme(
		axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
			angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=18,face="bold")
	) 
