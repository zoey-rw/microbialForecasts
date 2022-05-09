# Summarize outputs from diversity models
pacman::p_load(coda, tidyverse) 
source("./source.R")

plot_est_df_all <- list()
summary_df_all <- list()
gelman_list <- list()
allplots.scores.list <- list()

file.list <- list.files(path = "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/diversity/", 
												recursive = T,
												pattern = "samples_div_[I1]",
												full.names = T)

file_path = file.list[[1]]
file_summaries <- lapply(file.list, summarize_fg_div_model)

summary_df <- map(file_summaries, 1) %>% plyr::rbind.fill() %>% 
	mutate(time_period = recode(as.character(time_period), !!!date_recode))
plot_est <- map(file_summaries, 3) %>% plyr::rbind.fill() %>% 
	mutate(time_period = recode(as.character(time_period), !!!date_recode))
gelman_list <- map(file_summaries, 4) 
names(gelman_list) <- paste(basename(dirname(file.list)),basename(file.list), sep = "_")
scores.list <- map(file_summaries, 2) %>% plyr::rbind.fill() %>% 
	mutate(time_period = recode(as.character(time_period), !!!date_recode))



out <- list(plot_est = plot_est,
						summary_df = summary_df,
						gelman_list = gelman_list,
						scores.list = scores.list)

saveRDS(out, here("data", "/summary/div_summaries.rds"))




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
combined[combined$time_period=="2015-11_2018-01",]$time_period <- "2015-11_2020-01"
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
															 "2016-01_2020-01",
															 "2013-06_2017-01",
															 "2013-06_2016-01"), values = c(1,2,3,4)) +   
	scale_shape_manual(name = "Time period",
										 labels = c("2015-11_2018-01",
										 					 "2016-01_2020-01",
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
