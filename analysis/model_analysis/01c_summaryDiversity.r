# Summarize outputs from diversity models
pacman::p_load(coda, tidyverse) 


# # Read in samples for visualization
# read_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/div_samples_min3.rds")



plot_est_df_all <- list()
summary_df_all <- list()
gelman_list <- list()
allplots.scores.list <- list()



file.list <- list.files(path = "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/diversity/", recursive = T,
												pattern = "samples",
												full.names = T)
f <- file.list[[8]]
for (f in file.list){

	
#for (scenario in names(read_in$sample.list)){
	
	info <- basename(f) %>% str_split("_")
	model_name <- basename(dirname(f))
	time_period <- info[[1]][1]
	scenario <- paste0(info[[1]][4:6], collapse = "_") %>% str_replace(".rds", "")
	
	cat(paste0("\nSummarizing ", scenario, ", ", time_period, ", ", model_name))
	
	read_in <- readRDS(f)
	
	
	samples <- read_in$samples
	param_summary <- read_in$param_summary
	plot_summary <- read_in$plot_summary
	truth.plot.long <- read_in$metadata$model_data
	
	
	# Calculate plot means and actual values per rank
	plot_est <- plot_summary[[2]]
	pred.plot <- plot_est %>% as.data.frame() %>%  rownames_to_column() %>%
		separate(rowname, sep=", ", into=c("plot_num","timepoint")) %>%
		mutate(plot_num = as.integer(gsub("plot_mu\\[", "", plot_num)),
					 timepoint = as.integer(gsub("\\]", "", timepoint)))
	truth.plot.long <- truth.plot.long %>% 
		mutate(dates = as.Date(paste0(dateID, "01"), "%Y%m%d"),
					 timepoint = as.integer(timepoint))
	allplots <- merge(truth.plot.long, pred.plot, by = c("plot_num","timepoint"), all=T)
	allplots$scenario <- scenario
	allplots$time_period <- time_period
	allplots$model_name <- model_name
	plot_est_df_all[[f]] <- allplots
	
	# For scoring the predictions (need mean and SD)
	pred.plot.scores <- plot_summary[[1]] %>% as.data.frame() %>%  rownames_to_column() %>%
		separate(rowname, sep=", ", into=c("plot_num","timepoint")) %>%
		mutate(plot_num = as.integer(gsub("plot_mu\\[", "", plot_num)),
					 timepoint = as.integer(gsub("\\]", "", timepoint)))
	allplots.scores <- merge(truth.plot.long, pred.plot.scores, by = c("plot_num","timepoint"), all=T)
	allplots.scores$scenario <- scenario
	allplots.scores$time_period <- time_period
	allplots.scores$model_name <- model_name
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
	
	all_covariates_key <- c("1" = "Temperature",
												"2" = "Moisture",
												"3" = "pH",
												"4" = "pC",
												"5" = "Ectomycorrhizal trees",
												"6" = "LAI",
												"7" = "sin",
												"8" = "cos",
												"NA" = "NA")
	
	cycl_only_key <- list("1" = "sin",
												"2" = "cos")
	
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
	
	summary_df <- plyr::rbind.fill(site_eff_out, beta_out, intercept_out, sigma_out, sig_out, 
																 #rho_out, 
																 core_sd_out)
	summary_df$scenario <- scenario
	summary_df$time_period <- time_period
	summary_df$model_name <- model_name
	summary_df$effSize <- abs(summary_df$Mean)
	
	summary_df_all[[f]] <- summary_df
	
	
	## Calculate gelman diagnostics to assess convergence
	gd <- gelman.diag(samples, multivariate = FALSE)
	gelman_list[[f]] <- cbind(gd[[1]], effSize = effectiveSize(samples))
	
}

plot_est_df_all <- plyr::rbind.fill(plot_est_df_all) %>% 
	tidyr::separate(scenario, sep = "_", into = c("uncert1", "uncert2", "group"), remove = F) %>% 
	mutate(uncert = paste(uncert1, uncert2, sep = "_")) %>% select(-c(uncert1, uncert2)) %>% 
	mutate(fcast_type = "Diversity")

scores.list <- plyr::rbind.fill(allplots.scores.list) %>% 
	tidyr::separate(scenario, sep = "_", into = c("uncert1", "uncert2", "group"), remove = F) %>% 
	mutate(uncert = paste(uncert1, uncert2, sep = "_")) %>% select(-c(uncert1, uncert2)) %>% 
	mutate(fcast_type = "Diversity")

summary_df_all <- do.call(rbind, summary_df_all) %>% 
	tidyr::separate(scenario, sep = "_", into = c("uncert1", "uncert2", "group"), remove = F) %>% 
	mutate(uncert = paste(uncert1, uncert2, sep = "_")) %>% select(-c(uncert1, uncert2)) %>% 
	mutate(fcast_type = "Diversity")

summary_df_all$pretty_group <- ifelse(summary_df_all$group=="16S", "Bacteria", "Fungi")


saveRDS(list(plot_est = plot_est_df_all,
						 summary_df = summary_df_all,
						 samples = read_in$sample.list,
						 gelman_list = gelman_list,
						 scores.list = scores.list),
				"/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/div_summaries.rds")


# Copy over outputs
sum.all <- summary_df_all

# Or read in outputs instead
data_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/div_summaries.rds")
sum.all <- data_in$summary_df

# ## GGPLOT

beta_out <- sum.all[which(!is.na(sum.all$beta)),]
rownames(beta_out) <- NULL

refit <- beta_out[which(beta_out$time_period=="refit"),]


# Effect sizes & directions for refit vs calibration
ggplot(data=beta_out,
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
	) + scale_shape_manual(values = c(21, 16), name = NULL, 
												 labels = c("Calibration","Full dataset")) + geom_hline(yintercept = 0)


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
