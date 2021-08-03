# Summarize outputs from diversity models
library(hrbrthemes)
library(ggplot2)
library(tidyr)
library(coda)
library(dplyr)
library(tibble)

# # Read in samples for visualization
# read_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/div_samples_min3.rds")

read_in <- list()
read_in[[1]] <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/refit_nolog_div_no_uncertainty_ITS.rds")
read_in[[2]] <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/refit_nolog_div_full_uncertainty_ITS.rds")
read_in[[3]] <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/refit_nolog_div_no_uncertainty_16S.rds")
read_in[[4]] <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/refit_nolog_div_full_uncertainty_16S.rds")
names(read_in) <- c("no_uncertainty_ITS", 
										"full_uncertainty_ITS", "no_uncertainty_16S", 
										 "full_uncertainty_16S")
read_in$sample.list <- lapply(read_in, "[[", 1)
read_in$param.summary.list <- lapply(read_in, "[[", 2)
read_in$metadata.list <- lapply(read_in, "[[", 3)
read_in$plot.summary.list <- lapply(read_in, "[[", 4)


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
					 rank = "diversity")
	truth.plot.long <- truth.plot.long %>% 
		mutate(dates = as.Date(paste0(dateID, "01"), "%Y%m%d"),
					 timepoint = as.integer(timepoint))
	allplots <- merge(truth.plot.long, pred.plot, by = c("plot_num","timepoint"), all=T)
	allplots$scenario <- scenario
	plot_est_df_all[[scenario]] <- allplots
	
	
	pred.plot.scores <- plot_summary[[1]] %>% as.data.frame() %>%  rownames_to_column() %>%
		separate(rowname, sep=", ", into=c("plot_num","timepoint")) %>%
		mutate(plot_num = as.integer(gsub("plot_mu\\[", "", plot_num)),
					 timepoint = as.integer(gsub("\\]", "", timepoint)))
	allplots.scores <- merge(truth.plot.long, pred.plot.scores, by = c("plot_num","timepoint"), all=T)
	allplots.scores$scenario <- scenario
	allplots.scores.list[[scenario]] <- allplots.scores
	
	means <- param_summary[[2]]
	beta_out <- means[grep("beta", rownames(means)),]
	rho_out <- means[grep("rho", rownames(means)),,drop=F] %>% as.data.frame() %>% 
		rownames_to_column("rowname") %>% mutate(beta = "rho")
	intercept_out <- means[grep("intercept", rownames(means)),,drop=F] %>% as.data.frame() %>% 
		rownames_to_column("rowname")
	sigma_out <- means[grep("sigma", rownames(means)),,drop=F] %>% as.data.frame() %>% 
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
												 "6" = "% grasses",
												 "7" = "sin(month)",
												 "8" = "cosin(month)"))
	beta_out <- plyr::rbind.fill(beta_out, rho_out)
	beta_out$significant <- ifelse(beta_out$`2.5%` < 0 & beta_out$`97.5%` < 0 |
																 	beta_out$`2.5%` > -0 & beta_out$`97.5%` > -0,
																 1, 0)
	beta_out$effSize <- abs(beta_out$`50%`)
	
	
	# Get site effect sizes per rank
	site_eff_out <- means %>% as.data.frame() %>% 
		rownames_to_column("rowname") %>% filter(grepl("site", rowname)) %>% 
		mutate(site_num = as.character(gsub("site_effect\\[|\\]", "", rowname))) %>% 
		mutate(siteID = truth.plot.long[match(site_num, truth.plot.long$site_num),]$siteID)
	
	summary_df <- plyr::rbind.fill(site_eff_out, beta_out, intercept_out, sigma_out)
	summary_df$scenario <- scenario
	summary_df_all[[scenario]] <- summary_df
	
	
	## Calculate gelman diagnostics to assess convergence
	gd <- gelman.diag(samples, multivariate = FALSE)
	gelman_list[[scenario]] <- gd
}

plot_est_df_all <- plyr::rbind.fill(plot_est_df_all)
plot_est_df_all <- plot_est_df_all %>% tidyr::separate(scenario, sep = "_", into = c("uncert1", "uncert2", "group"))
plot_est_df_all$scenario <- paste(plot_est_df_all$uncert1, plot_est_df_all$uncert2)

scores.list <- plyr::rbind.fill(allplots.scores.list)
scores.list <- scores.list %>% tidyr::separate(scenario, sep = "_", into = c("uncert1", "uncert2", "group"))
scores.list$scenario <- paste(scores.list$uncert1, scores.list$uncert2)

summary_df_all <- do.call(rbind, summary_df_all)
summary_df_all <- summary_df_all %>% tidyr::separate(scenario, sep = "_", into = c("uncert1", "uncert2", "group"))
summary_df_all$scenario <- paste(summary_df_all$uncert1, summary_df_all$uncert2)

saveRDS(list(plot_est = plot_est_df_all,
						 summary_df = summary_df_all,
						 samples = read_in$sample.list,
						 gelman_list = gelman_list,
						 scores.list = scores.list),
				"/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/refit_div_summaries.rds")




data_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/refit_div_summaries.rds")
# ## GGPLOT
sum.all <- data_in$summary_df

beta_out <- sum.all[which(!is.na(sum.all$beta)),]


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

samps_long <- samps_long %>% 
	mutate(beta_num = as.numeric(gsub("beta\\[|\\]", "", name))) %>% 
	mutate(beta = recode(beta_num,
											 "1" = "Temperature",
											 "2" = "Moisture",
											 "3" = "pH",
											 "4" = "pC",
											 "5" = "Plant species richness",
											 "6" = "% grasses",
											 .missing = "Autocorrelation"))




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





out.data <- data_in$plot_est[data_in$plot_est$scenario == "full uncertainty",]
site_df <-  out.data[out.data$siteID %in% c("HARV","OSBS","CPER"),]
ggplot(site_df) +
	#	ggplot(out.data[out.data$group=="16S",]) +
	
	#facet_grid(rows=vars(siteID), cols=vars(group), drop=T, scales="free", space="free") +
	#geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`, fill = plotID), 
	#					 alpha=0.4, show.legend = F) +
	#theme_ipsum(base_size = 18, strip_text_size = 22) + 
	ggtitle(paste0("Shannon diversity hindcasts for ")) + #scale_x_date() +	
	theme_minimal(base_size=18) +
	theme(panel.spacing = unit(.2, "cm"),
				plot.margin = unit(c(.2, .2, .2, .2), "cm")) +#+ ylim(c(1.5,7)) +
	ylab("Shannon diversity") + geom_point(aes(x = dates, y = truth, color=siteID))
