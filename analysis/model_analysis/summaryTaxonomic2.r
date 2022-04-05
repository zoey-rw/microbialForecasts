library(scales) 
library(coda)
library(tidyverse)

source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepTaxonomicData.r")
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

summary_list <- list()
taxon_key_list <- list()

# Create lists for outputs.
plot_est_df_all <- list()
summary_df_all <- list()
gelman_list <- list()
allplots.scores.list <- list()

# For testing.
scenario <- "full_uncertainty"
rank.name <- "phylum_bac"

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


file.list <- list.files(path = "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/taxa/", recursive = T,
												#pattern = "^refit_samples",
												pattern = "chain",
												full.names = T)
f <- file.list[[1]]

for (f in file.list) {
	# Subset to one rank.
	#rank.df <- d[[rank.name]]  
	group.model.path <- f

	
	info <- basename(f) %>% str_split("_")
	model_name <- basename(dirname(f))
	time_period <- info[[1]][1]
	last <- length(info[[1]])
	scenario <- paste0(info[[1]][(last-1):last], collapse = "_") %>% str_replace(".rds", "")
	rank.name <- paste0(info[[1]][3:(last-2)], collapse = "_") 
	
	
	cat(paste0("\nSummarizing ", model_name, ", ", time_period, ", ", rank.name, "\n"))
	

		# Read in samples
	read_in <- readRDS(group.model.path)
		samples <- read_in$samples
		param_summary <- read_in$param_summary
		plot_summary <- read_in$plot_summary
		# TO DO: remove once output is fixed
		truth.plot.long <- read_in$metadata$model_data
		print(read_in$metadata[1:3])
		taxon_key <- unique(truth.plot.long$species)
		names(taxon_key) <- seq(1, length(taxon_key))
		
		
		# truth.plot.long <- cbind.data.frame(plotID = read_in$metadata$model_data.plotID,
		# 																	plot_num = read_in$metadata$model_data.plot_num,
		# 																	timepoint = read_in$metadata$model_data.timepoint,
		# 																	truth = read_in$metadata$model_data.truth,
		# 																	dateID = read_in$metadata$model_data.dateID,
		# 																	taxon = read_in$metadata$model_data.species,
		# 																	site_num = read_in$metadata$model_data.site_num,
		# 																	siteID = read_in$metadata$model_data.siteID)
		truth.plot.long <- truth.plot.long %>% 
			mutate(dates = as.Date(paste0(dateID, "01"), "%Y%m%d")) %>% 
						 rename(taxon = species)
		
		
		# For scoring the predictions (need mean and SD)
		pred.plot.scores <- plot_summary[[1]] %>% as.data.frame() %>%  rownames_to_column() %>%
			separate(rowname, sep=", ", into=c("plot_num","species_num","timepoint")) %>%
			mutate(plot_num = as.integer(gsub("plot_rel\\[", "", plot_num)),
						 timepoint = as.integer(gsub("\\]", "", timepoint)),
						 rank = !!rank.name,
						 scenario = !!scenario,
						 time_period = !!time_period,
						 model_name = !!model_name,
						 taxon = recode(species_num,
						 							 !!!taxon_key))
		allplots.scores <- merge(truth.plot.long, pred.plot.scores, by = c("plot_num","timepoint","taxon"), all=T)
		allplots.scores.list[[f]] <- allplots.scores
		
		
		# Calculate plot means and actual values per rank
		plot_est <- plot_summary[[2]]
		pred.plot <- plot_est %>% as.data.frame() %>%  rownames_to_column() %>%
			separate(rowname, sep=", ", into=c("plot_num","species_num","timepoint")) %>%
			
			mutate(plot_num = as.integer(gsub("plot_rel\\[", "", plot_num)),
						 timepoint = as.integer(gsub("\\]", "", timepoint)),
						 rank = !!rank.name,
						 scenario = !!scenario,
						 time_period = !!time_period,
						 model_name = !!model_name,
						 taxon = recode(species_num,
						 							 !!!taxon_key))
		allplots <- merge(truth.plot.long, pred.plot, by = c("plot_num","timepoint","taxon"), all=T) %>% arrange(plot_num, date_num)
		plot_est_df_all[[f]] <- allplots
		
		
		
		if(model_name == "all_covariates") cov_key <- all_covariates_key
		if(model_name == "cycl_only") cov_key <- cycl_only_key
		
		
		# Get median/quantile values for parameters
		means <- param_summary[[1]]
		rho_out <- means[grep("rho", rownames(means)),,drop=F] %>% as.data.frame() %>% 
			rownames_to_column("rowname") %>% mutate(beta = "rho") %>% 
			mutate(taxon_num = as.numeric(gsub("rho\\[|\\]", "", rowname)),
						 taxon = recode(taxon_num, !!!taxon_key))
		sigma_out <- means[grep("sigma", rownames(means)),,drop=F] %>% as.data.frame() %>% 
			rownames_to_column("rowname") %>% 
			mutate(taxon_num = as.numeric(gsub("sigma\\[|\\]", "", rowname)),
		taxon = recode(taxon_num, !!!taxon_key))
		
		sig_out <- means[grep("sig$", rownames(means)),,drop=F] %>% as.data.frame() %>% 
			rownames_to_column("rowname")
		int_out <- means[grep("intercept", rownames(means)),,drop=F] %>% as.data.frame() %>% 
			rownames_to_column("rowname") %>% 
			mutate(taxon_num = as.numeric(gsub("intercept\\[|\\]", "", rowname)),
																		taxon = recode(taxon_num, !!!taxon_key),
																		rank = rank.name)
		
		# Get site effect sizes per rank
		site_eff_out <- means %>% as.data.frame() %>% 
			rownames_to_column("rowname") %>% filter(grepl("site", rowname)) %>% 
			separate(rowname, sep=", ", into=c("site_num","taxon_num"), remove=F) %>% 
			mutate(site_num = as.character(gsub("site_effect\\[|\\]", "", site_num)),
						taxon_num = as.character(gsub("\\]", "", taxon_num))) %>% 
			mutate(siteID = truth.plot.long[match(site_num, truth.plot.long$site_num),]$siteID,
						 taxon = recode(taxon_num, !!!taxon_key), 
						 rank = !!rank.name)
		
		# Get beta sizes per rank
		beta_out <-  means %>% as.data.frame() %>% 
			rownames_to_column("rowname") %>% filter(grepl("beta", rowname)) %>% 
			separate(rowname, sep=", ", into=c("taxon_num","beta_num"), remove = F) %>% 
			mutate(taxon_num = as.numeric(gsub("beta\\[", "", taxon_num)),
						 beta_num = as.numeric(gsub("\\]", "", beta_num)),
						 param = "beta") %>% 
			mutate(beta = recode(beta_num, !!!cov_key),
						 taxon = recode(taxon_num, !!!taxon_key),
						 rank = rank.name)
		
		
		# Use quantiles to assign significance to beta parameters.
		beta_ci <-  param_summary[[2]] %>% as.data.frame() %>% 
			rownames_to_column("rowname") %>% filter(grepl("beta", rowname)) %>% 
			separate(rowname, sep=", ", into=c("taxon_num","beta_num"), remove = F) %>% 
			mutate(taxon_num = as.numeric(gsub("beta\\[", "", taxon_num)),
						 beta_num = as.numeric(gsub("\\]", "", beta_num)))
		
		beta_out$significant <- ifelse(beta_ci$`2.5%` < 0 & beta_ci$`97.5%` < 0 |
																	 	beta_ci$`2.5%` > -0 & beta_ci$`97.5%` > -0,
																	 1, 0)
		beta_out$effSize <- abs(beta_out$Mean)
		
		# Combine parameter estimates into summary
		summary_df <- plyr::rbind.fill(int_out, site_eff_out, beta_out, sigma_out, rho_out, sig_out) %>%
			mutate(scenario = scenario,
						 rank = !!rank.name,
						 scenario = !!scenario,
						 time_period = !!time_period,
						 model_name = !!model_name,)
		summary_df_all[[f]] <- summary_df
		
		
		## Calculate gelman diagnostics to assess convergence
		gd <- gelman.diag(samples, multivariate = FALSE)
		gelman_list[[f]] <- cbind(gd[[1]], effSize = effectiveSize(samples))
		
		
}
	
	plot_est_df_all_ranks <- plyr::rbind.fill(plot_est_df_all)
	scores.list_all_ranks <- plyr::rbind.fill(allplots.scores.list)
	summary_df_all_ranks <- do.call(rbind, summary_df_all)
	
	
	summary_df_all_ranks <- summary_df_all_ranks %>% mutate(fcast_type = "Taxonomic",
																	fg_cat = "Bacterial taxa",
																	group = "16S") 
	summary_df_all_ranks[grepl("fun$", summary_df_all_ranks$rank),]$fg_cat = "Fungal taxa"
	summary_df_all_ranks[grepl("fun$", summary_df_all_ranks$rank),]$group = "ITS"

saveRDS(list(plot_est = plot_est_df_all_ranks,
						 summary_df = summary_df_all_ranks,
						 gelman_list = gelman_list,
						 scores.list = scores.list_all_ranks), "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/taxa_summaries.rds")
