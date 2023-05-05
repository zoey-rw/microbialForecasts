# Source script for Nimble models (w/ environmental covariates, & cyclical covariates only)
setwd("/projectnb/dietzelab/zrwerbin/microbialForecasts/")
here::i_am("source.R")

# Pacman package loader is used throughout scripts
if (!require("pacman")) install.packages("pacman")
if (!require("microbialForecast")) install.packages("/projectnb/dietzelab/zrwerbin/microbialForecasts/microbialForecast_0.1.0.tar.gz", repos = NULL, type="source")

library("microbialForecast", lib.loc="/usr3/graduate/zrwerbin/R/x86_64-pc-linux-gnu-library/4.2")

	suppressPackageStartupMessages(library(tidyverse, warn.conflicts = F))
	pacman::p_load(pacman,
							 nimble, coda, lubridate, here,
							 doParallel, data.table, Rfast, moments,
							 scoringRules, Metrics, ggpubr)


	# Create output directory for MCMC runs
	model_output_dir = here("data", "model_outputs", "logit_beta_regression")
	dir.create(model_output_dir, showWarnings = FALSE)
	model_output_dir = here("data", "model_outputs", "logit_beta_regression", "env_cycl")
	dir.create(model_output_dir, showWarnings = FALSE)
	model_output_dir = here("data", "model_outputs", "logit_beta_regression", "cycl_only")
	dir.create(model_output_dir, showWarnings = FALSE)

	model_output_dir = here("data", "model_outputs", "logit_beta_regression", "env_cov")
	dir.create(model_output_dir, showWarnings = FALSE)

	model_summary_dir = here("data", "summary")
	dir.create(model_summary_dir, showWarnings = FALSE)


	# New facet label names for supp variable
	model.labs <- c("Environmental\npredictors", "Seasonality", "Environmental predictors\n& seasonality")
	names(model.labs) <- c("env_cov", "cycl_only","env_cycl")


	
	metric.labs <- c("Relative forecast error (nRMSE)", "Absolute forecast error (CRPS)")
	names(metric.labs) <- c("RMSE.norm", "mean_crps")
	
	# Convert sin and cos effect sizes to a seasonal amplitude parameter
	sin_cos_to_seasonality <- function(sin, cos){
		if (sin==0 & cos==0|is.na(sin)|is.na(cos)) {return(cbind.data.frame(max=NA,
																									amplitude_orig=NA,
																									amplitude = NA))}
		min_max <- getMaxMin(sin, cos, max_only = F)
		amplitude <- sqrt(sin^2 + cos^2)

		t=seq(0,12,0.1)
		monthly_vals = sin*sin(2*pi*t/12)+cos*cos(2*pi*t/12)
		max_val = max(monthly_vals)
		# Average of minimum and maximum wave values
		avg_val <- mean(min_max[[1]], min_max[[2]])
		out <- cbind.data.frame(max=min_max[[1]],
														amplitude_orig=amplitude,
														amplitude = max_val)
		return(out)
	}

	source("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/NEFI_functions/rsq_1.1.r")


	predictive_loss = function(observed, predicted, predicted_sd){
		npred = length(predicted)
		predictive_variance = predicted_sd^2
		residual_variance = (predicted - observed)^2
		P = sum(predictive_variance, na.rm=T)/npred
		G = sum(residual_variance, na.rm=T)/npred
		total_PL = P+G
		data.frame(total_PL = total_PL, predictive_variance=P, residual_variance=G)

	}

	# observed = cal_test$truth
	# mean_predicted = cal_test$Mean
	# sd_predicted = cal_test$SD
	# type=c("RMSE","BIAS","MAE",
	# 			 "CRPS", "RSQ", "RSQ.1",
	# 			 "RMSE.norm",  "residual_variance", "predictive_variance", "total_PL")
	#
	# add_scoring_metrics(observed = cal_test$truth,
	# 										mean_predicted = cal_test$Mean,
	# 										sd_predicted = cal_test$SD)



	first = function(x) x %>% nest %>% ungroup %>% slice(1) %>% unnest(data)



	quick_get_rank_df = function(k = 1,
															 min.date = "20151101",
															 max.date = "20200101"){

		pacman::p_load(reshape2, parallel, nimble, coda, tidyverse)

		# Subset to one rank.
		rank.name <- microbialForecast:::tax_names[k]
	# Read in microbial abundances
	cal <- c(readRDS(here("data", "clean/cal_groupAbundances_16S_2021.rds")),
					 readRDS(here("data", "clean/cal_groupAbundances_ITS_2021.rds")))
	val <- c(readRDS(here("data", "clean/val_groupAbundances_16S_2021.rds")),
					 readRDS(here("data", "clean/val_groupAbundances_ITS_2021.rds")))

	cal.rank.df <- cal[[rank.name]]
	val.rank.df <- val[[rank.name]]
	rank.df <- rbind(cal.rank.df, val.rank.df)

	# Prep model inputs/outputs.
	print(paste0("Preparing model data for ", rank.name))

	# spec_names <- colnames(rank.df)[!colnames(rank.df) %in% c("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date","other")]
	# rank.df_spec <- rank.df %>%
	# 	select("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date", !!spec_names)
	# rank.df_spec$other <- 1-rank.df_spec[[s]]

	model.dat <- prepTaxonomicData(rank.df = rank.df,
																 min.prev = 3,
																 min.date = min.date,
																 max.date = max.date)
	return(model.dat)
	}


	calc_cv <- function(x) sd(x, na.rm = T) / mean(x, na.rm = T) * 100



	pivot_metrics = function(df) {
		df %>% pivot_longer(cols = c(RMSE, BIAS, MAE, CRPS,CRPS_truncated, RSQ, RSQ.1,
																 RMSE.norm, residual_variance, predictive_variance, total_PL),
												names_to = "metric", values_to = "score")
	}







	combine_chains_existing = function(input_list,
																		 save = FALSE,
																		 cut_size1 = NULL,
																		 cut_size2 = NULL){
		require(coda)
		require(tidyverse)

		if (is.null(cut_size1)) cut_size1 <- 19999
		if (is.null(cut_size2)) cut_size2 <- 9999

		readInputRdsFile = function(input_rds){
			input = tryCatch(readRDS(input_rds),
											 error = function(c) {
											 	message("The input *rds is invalid")
											 	return(NA)
											 }
			)
		}

		# initialize
		samples <- metadata <- list()
		first_iter <- last_iter <- list()
		for(i in 1:length(input_list)){
			print(i)

			if (class(input_list[[i]])=="character") {
				# paste model file path to chain number
				chain <- readInputRdsFile(chain_paths[[i]])
				if (any(is.na(chain))) next()
				samples[[i]] <- chain[[1]]
				#samples2[[i]] <- chain[[2]]
			} else {
				samples[[i]]  = input_list[[i]]
				#samples2[[i]]  = input_list[[i]][[2]]
			}
		}

		samples<-samples[!sapply(samples,is.null)]
#		samples2<-samples2[!sapply(samples2,is.null)]

		# Now make them all the same size
		nrows <- lapply(samples, nrow) %>% unlist()
		min_nrow <- min(nrows)
		for(i in 1:length(samples)){
			current_nrow <- nrow(samples[[i]])
			if (min_nrow < current_nrow){
				samples[[i]] <- window_chain(samples[[i]], max_size = (min_nrow-1))
			}
		}


		# Now make them all the same size, v2
		# nrows <- lapply(samples2, nrow) %>% unlist()
		# min_nrow <- min(nrows)
		# for(i in 1:length(samples2)){
		# 	current_nrow <- nrow(samples2[[i]])
		# 	if (min_nrow < current_nrow){
		# 		samples2[[i]] <- window_chain(samples2[[i]], max_size = (min_nrow-1))
		# 	}
		# }

		# Make the attributes match up (sort of arbitrary)
		for (i in 1:length(samples)) {
			attr(samples[[i]], "mcpar") = attr(samples[[1]], "mcpar")
			#attr(samples2[[i]], "mcpar") = attr(samples2[[1]], "mcpar")
		}

		out <- as.mcmc.list(samples)
		#out2 <- as.mcmc.list(samples2)

		return(out)
	}


assign_fg_sources <-function (vector)
	{
	out <- rep(NA, length(vector))
	out[which(grepl("lytic", vector))] <-  "Literature review"
	out[which(grepl("cellulolytic", vector))] <-  "Literature review + genomic pathway"
	out[which(grepl("nitr|fixa", vector))] <- "Literature review + genomic pathway"
	out[which(grepl("complex|simple|stress", vector, fixed = F))] <- "Experimental enrichment"
	out[which(grepl("antibiotic", vector))] <- "Experimental enrichment"
	out[which(grepl("anaerobic", vector))] <- "Experimental enrichment"
#	out[which(grepl("nitr|fixa", vector))] <- "Literature review"
	out[which(grepl("troph", vector))] <- "Literature review"
	out[which(grepl("sapr|path|arbusc|ecto|endo|lichen", vector))] <- "Literature review"
	out[which(grepl("other", vector))] <- NA

		#
		# out <- rep(NA, length(vector))
		# out[which(grepl("substrates|resistance|anaerobic|stress", vector))] <- "Experimental enrichment"
		# out[which(grepl("cycling|cellulo", vector))] <- "Genomic pathways"
		# out[which(grepl("Trophic|Life", vector))] <- "Scientific consensus"
		return(out)
}



rLogitBeta <- nimbleFunction (
	## Generates y ~ Beta(a1,a2)
	## Returns   x = logit(y)
	run = function(n = integer(0, default=1),
								 shape1 = double(0, default=1.0),
								 shape2 = double(0, default=1.0)) {
		returnType(double(0))
		if(n != 1)
			nimPrint("Warning: rLogitBeta only allows n = 1; Using n = 1.\n")
		y <- rbeta(1, shape1=shape1, shape2=shape2)
		x <- logit(y)
		return(x)
	}
)




summarize_nan_iterations = function(problem_param_name, samples){
	problem_param = samples[,problem_param_name]
	problem_iter1 = which(is.nan(problem_param[[1]]))
	problem_iter2 = which(is.nan(problem_param[[2]]))
	problem_iter3 = which(is.nan(problem_param[[3]]))
	fixed_param = problem_param[-c(problem_iter1,problem_iter2,problem_iter3),] %>%
		lapply(as.mcmc) %>%
		as.mcmc.list()
	fixed_param_summary = fast.summary.mcmc(fixed_param)
	mean_val = fixed_param_summary[[1]]["Mean"]
	mean_sd = fixed_param_summary[[1]]["SD"]
	return(list(Mean = mean_val, SD = mean_sd))
}



tukey2 = function (x, y, extra_info = NULL, y.offset = 0.3)
{
	new.df <- cbind.data.frame(x = x, y = y)
	abs_max <- max(new.df[, 2], na.rm=T)
	maxs <- new.df %>% group_by(x) %>% summarise(tot = max(y, na.rm=T) +
																							 	y.offset * abs_max)
	Tukey_test <- aov(y ~ x, data = new.df) %>% agricolae::HSD.test("x",
																																	group = TRUE) %>% .$groups %>% as_tibble(rownames = "x") %>%
		rename(Letters_Tukey = "groups") %>% dplyr::select(-y) %>%
		left_join(maxs, by = "x")
	if (!is.null(extra_info)) {
		Tukey_test <- cbind.data.frame(Tukey_test)
	}
	return(Tukey_test)
}


tag_facet <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf,
											hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {

	gb <- ggplot_build(p)
	lay <- gb$layout$layout
	tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
	p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust,
								vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
}
