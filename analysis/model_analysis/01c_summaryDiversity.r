# Summarize outputs from diversity models
pacman::p_load(coda, tidyverse) 
source("./source.R")


# # Read in samples for visualization
# read_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/div_samples_min3.rds")



plot_est_df_all <- list()
summary_df_all <- list()
gelman_list <- list()
allplots.scores.list <- list()

file.list <- list.files(path = "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/diversity/", 
												recursive = T,
												pattern = "samples_div_[I1]",
												full.names = T)
#f <- file.list[[8]]
for (f in file.list){
#for (scenario in names(read_in$sample.list)){
	
	# Extract run information
	info <- basename(f) %>% str_split("_") %>% unlist()
	model_name <- basename(dirname(f))
	group <- info[3]
	time_period <- paste0(info[4:5], collapse = "_") %>% str_replace(".rds", "")

	message("\nSummarizing ", group, ", ", time_period, ", ", model_name)

	read_in <- readRDS(f)
	samples <- read_in$samples
	param_summary <- read_in$param_summary
	plot_summary <- read_in$plot_summary
	truth.plot.long <- read_in$metadata$model_data
	
	# Add some info to observational data for merging
	truth.plot.long <- truth.plot.long %>% 
		mutate(dates = fixDate(dateID),
					 group = !!group,
					 model_name = !!model_name,
					 time_period = !!time_period,
					 fcast_type = "Diversity")
	
	# Calculate plot means and actual values per rank
	plot_est <- plot_summary[[2]]
	pred.plot <- plot_est %>% as.data.frame() %>%  rownames_to_column() %>%
		separate(rowname, sep=", ", into=c("plot_num","timepoint")) %>%
		mutate(plot_num = as.integer(gsub("plot_mu\\[", "", plot_num)),
					 timepoint = as.integer(gsub("\\]", "", timepoint)))
	allplots <- merge(truth.plot.long, pred.plot, by = c("plot_num","timepoint"), all=T)
	plot_est_df_all[[f]] <- allplots
	
	# For scoring the predictions (need mean and SD)
	pred.plot.scores <- plot_summary[[1]] %>% as.data.frame() %>%  rownames_to_column() %>%
		separate(rowname, sep=", ", into=c("plot_num","timepoint")) %>%
		mutate(plot_num = as.integer(gsub("plot_mu\\[", "", plot_num)),
					 timepoint = as.integer(gsub("\\]", "", timepoint)))
	allplots.scores <- merge(truth.plot.long, pred.plot.scores, by = c("plot_num","timepoint"), all=T)
	allplots.scores.list[[f]] <- allplots.scores
	
	# Get mean values for parameters
	means <- param_summary[[1]]
	beta_out <- means[grep("beta|rho", rownames(means)),]
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
	
	if(model_name == "all_covariates") cov_key <- all_covariates_key
	if(model_name == "cycl_only") cov_key <- cycl_only_key
	
	# Get beta sizes per rank
	beta_out <-  beta_out %>% as.data.frame() %>% 
		rownames_to_column("rowname") %>% 
		mutate(beta_num = as.numeric(gsub("beta\\[|\\]", "", rowname))) %>% 
		mutate(beta = recode(as.character(beta_num),!!!cov_key))
	beta_out[grep("rho", beta_out$rowname),]$beta = "rho"
	beta_out[grep("rho", beta_out$rowname),]$beta_num = "0"
	
	
	# Use quantiles to assign significance to beta parameters.
	beta_ci <-  param_summary[[2]] %>% as.data.frame() %>% 
		rownames_to_column("rowname") %>% filter(grepl("beta|rho", rowname)) %>% 
		mutate(beta_num = as.numeric(gsub("beta\\[|\\]", "", rowname)))
	beta_ci[grep("rho", beta_ci$rowname),]$beta_num = "0"
	beta_out$significant <- ifelse(beta_ci$`2.5%` < 0 & beta_ci$`97.5%` < 0 |
																 	beta_ci$`2.5%` > -0 & beta_ci$`97.5%` > -0,
																 1, 0)

	# Get site effect sizes per rank
	site_eff_out <- means %>% as.data.frame() %>% 
		rownames_to_column("rowname") %>% filter(grepl("site", rowname)) %>% 
		mutate(site_num = as.character(gsub("site_effect\\[|\\]", "", rowname))) %>% 
		mutate(siteID = truth.plot.long[match(site_num, truth.plot.long$site_num),]$siteID)
	
	summary_df <- plyr::rbind.fill(site_eff_out, beta_out, intercept_out, sigma_out, sig_out, core_sd_out)
	summary_df <- summary_df %>% mutate(group = !!group,
																			model_name = !!model_name,
																			time_period = !!time_period,
																			effSize = abs(Mean),
																			fcast_type = "Diversity")
	summary_df_all[[f]] <- summary_df
	
	
	## Calculate gelman diagnostics to assess convergence
	gd <- gelman.diag(samples, multivariate = FALSE)
	gelman_list[[f]] <- cbind(gd[[1]], effSize = effectiveSize(samples))
	
}

plot_est_df_all <- plyr::rbind.fill(plot_est_df_all) 
scores.list <- plyr::rbind.fill(allplots.scores.list)
summary_df_all <- do.call(rbind, summary_df_all)
plot_est_df_all$pretty_group <- ifelse(plot_est_df_all$group=="16S", "Bacteria", "Fungi")
scores.list$pretty_group <- ifelse(scores.list$group=="16S", "Bacteria", "Fungi")
summary_df_all$pretty_group <- ifelse(summary_df_all$group=="16S", "Bacteria", "Fungi")
rownames(summary_df_all) <- NULL

scores.list[scores.list$time_period=="20160101_20180101",]$time_period <- "2016-01_2020-01"
plot_est_df_all[plot_est_df_all$time_period=="20160101_20180101",]$time_period <- "2016-01_2020-01"
summary_df_all[summary_df_all$time_period=="20160101_20180101",]$time_period <- "2016-01_2020-01"

scores.list[scores.list$time_period=="20151101_20180101",]$time_period <- "2015-11_2018-01"
plot_est_df_all[plot_est_df_all$time_period=="20151101_20180101",]$time_period <- "2015-11_2018-01"
summary_df_all[summary_df_all$time_period=="20151101_20180101",]$time_period <- "2015-11_2018-01"


saveRDS(list(plot_est = plot_est_df_all,
						 summary_df = summary_df_all,
						 gelman_list = gelman_list,
						 scores.list = scores.list),
				"/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/div_summaries.rds")


# read in outputs instead
data_in_legacy <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/div_summaries_legacy_incl.rds")
sum.all_legacy <- data_in_legacy$summary_df
beta_legacy <- sum.all_legacy[which(!is.na(sum.all_legacy$beta)),]
rownames(beta_legacy) <- NULL

data_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/div_summaries.rds")
sum.all <- data_in$summary_df
beta_recent <- sum.all[which(!is.na(sum.all$beta)),]

combined <- plyr::rbind.fill(beta_legacy, beta_recent)
combined[combined$time_period=="20151101_20180101",]$time_period <- "2015-11_2018-01"
combined[combined$time_period=="refit",]$time_period				     <- "2013-06_2017-01"
combined[combined$time_period=="calibration",]$time_period       <- "2013-06_2016-01"

# ## GGPLOT

cycl_combined <- combined %>% filter(model_name == "cycl_only")

##### Effect sizes & directions for refit vs calibration ####
ggplot(data=beta_recent,
			 aes(x = reorder(beta, Mean),y = Mean)) +
	facet_grid(cols = vars(model_name), 
						 rows = vars(group), drop = T) +
	geom_point(aes(shape = as.factor(time_period), 
								 color = beta), size = 4) +
	labs(col = "Parameter", title = "Absolute effect size") + 
	xlab("Parameter")+ 
	ylab(NULL) +
	theme_bw() + theme(#axis.ticks.x=element_blank(),
		text = element_text(size = 16),
		axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
			angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=22,face="bold")
		#strip.text.y = element_text(size=24,hjust=0,vjust = 1,angle=180,face="bold")
	) + #scale_shape_manual(values = c(21, 16), name = NULL, 
			#									 labels = c("Calibration","Full dataset")) + 
	geom_hline(yintercept = 0)

ggplot(data=combined,
			 aes(x = reorder(beta, Mean),y = Mean)) +
	facet_grid(cols = vars(model_name), 
						 rows = vars(group), drop = T) +
	geom_point(aes(shape = as.factor(time_period), 
								 color = as.factor(time_period)), size = 4) +
	labs(col = "Parameter", title = "Calibration period effect size") + 
	xlab("Parameter")+ 
	ylab(NULL) +
	theme_bw() + theme(#axis.ticks.x=element_blank(),
		text = element_text(size = 16),
		axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
			angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=22,face="bold")
		#strip.text.y = element_text(size=24,hjust=0,vjust = 1,angle=180,face="bold")
	) + geom_hline(yintercept = 0) + 
	scale_colour_manual(name = "Time period",
										labels = c("2015-11_2018-01",
															 "2016-01_2018-01",
															 "2013-06_2017-01",
															 "2013-06_2016-01"), values = c(1,2,3,4)) +   
	scale_shape_manual(name = "Time period",
										 labels = c("2015-11_2018-01",
										 					 "2016-01_2018-01",
										 										"2013-06_2017-01",
										 										"2013-06_2016-01"), values = c(19,22,21,23))


# Supplemental fig 1?
# Abs. effect sizes for all covariates vs cyc alone
ggplot(data=refit,
			 aes(x = reorder(beta, effSize),y = effSize)) +
	facet_grid(#cols = vars(model_name), 
						 rows = vars(group), drop = T) +
	geom_point(aes(shape = as.factor(model_name), 
								 color = beta), size = 4) +
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
												 labels = c("All covariates","No environmental covariates")) 


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
