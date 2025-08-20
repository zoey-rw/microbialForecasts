# ====================
# Microbial Forecasts Environment Setup
# ====================
# Source script for NIMBLE models and analysis environment
# This script sets up the R environment for microbial forecasting analysis

# Automatically detect project root directory
if (!exists("project_root")) {
  # Try multiple methods to find project root
  if (rstudioapi::isAvailable() && !is.null(rstudioapi::getActiveProject())) {
    # If in RStudio with a project
    project_root <- rstudioapi::getActiveProject()
  } else if (file.exists("source.R")) {
    # If source.R exists in current directory
    project_root <- getwd()
  } else if (file.exists("microbialForecasts/source.R")) {
    # If in parent directory
    project_root <- file.path(getwd(), "microbialForecasts")
  } else {
    # Default to current working directory
    project_root <- getwd()
    warning("Could not automatically detect project root. Using current working directory.")
  }
}

# Set working directory and establish here root
setwd(project_root)
here::i_am("source.R")

cat("Project root set to:", project_root, "\n")

# ====================
# Package Installation and Loading
# ====================
# Install pacman if not available
if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman")
  library("pacman")
}

# Try to load microbialForecast package
microbial_pkg_loaded <- FALSE

# First, try loading from standard library locations
if (require("microbialForecast", quietly = TRUE)) {
  microbial_pkg_loaded <- TRUE
  cat("microbialForecast package loaded from library\n")
} else {
  # Try to install from local package directory
  local_pkg_path <- here::here("microbialForecast")
  if (dir.exists(local_pkg_path)) {
    tryCatch({
      devtools::install_local(local_pkg_path, quiet = TRUE)
      library("microbialForecast")
      microbial_pkg_loaded <- TRUE
      cat("microbialForecast package installed and loaded from local source\n")
    }, error = function(e) {
      cat("Failed to install from local source:", e$message, "\n")
    })
  }
  
  # If still not loaded, try to find tar.gz file
  if (!microbial_pkg_loaded) {
    pkg_files <- list.files(pattern = "microbialForecast.*\\.tar\\.gz$", recursive = TRUE)
    if (length(pkg_files) > 0) {
      tryCatch({
        install.packages(pkg_files[1], repos = NULL, type = "source")
        library("microbialForecast")
        microbial_pkg_loaded <- TRUE
        cat("microbialForecast package installed from:", pkg_files[1], "\n")
      }, error = function(e) {
        cat("Failed to install from tar.gz:", e$message, "\n")
      })
    }
  }
}

if (!microbial_pkg_loaded) {
  warning("microbialForecast package could not be loaded. Some functions may not be available.")
}

# ====================
# Load Required Packages
# ====================
cat("Loading required packages...\n")

# Core packages
suppressPackageStartupMessages({
  library(tidyverse, warn.conflicts = FALSE)
})

# Load additional packages with error handling
required_packages <- c(
  "nimble", "coda", "lubridate", "here",
  "doParallel", "data.table", "Rfast", "moments",
  "scoringRules", "Metrics", "ggpubr", "ggallin"
)

# Check if devtools is available for local package installation
if ("devtools" %in% required_packages && !requireNamespace("devtools", quietly = TRUE)) {
  required_packages <- c(required_packages, "devtools")
}

tryCatch({
  pacman::p_load(char = required_packages, character.only = TRUE)
  cat("All required packages loaded successfully\n")
}, error = function(e) {
  cat("Warning: Some packages could not be loaded:", e$message, "\n")
  cat("Continuing with available packages...\n")
})

# ====================
# Create Directory Structure
# ====================
cat("Setting up directory structure...\n")

# Define directory structure
dirs_to_create <- c(
  here("data", "model_outputs"),
  here("data", "model_outputs", "logit_beta_regression"),
  here("data", "model_outputs", "logit_beta_regression", "env_cycl"),
  here("data", "model_outputs", "logit_beta_regression", "cycl_only"),
  here("data", "model_outputs", "logit_beta_regression", "env_cov"),
  here("data", "model_outputs", "functional_groups"),
  here("data", "model_outputs", "diversity"),
  here("data", "summary"),
  here("figures")
)

# Create directories
for (dir_path in dirs_to_create) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    cat("Created directory:", dir_path, "\n")
  }
}

# ====================
# Define Analysis Parameters and Labels
# ====================

# Model type labels for plotting
model.labs <- c(
  "Environmental\npredictors", 
  "Seasonality", 
  "Environmental predictors\n& seasonality"
)
names(model.labs) <- c("env_cov", "cycl_only", "env_cycl")

# Metric labels for plotting
metric.labs <- c(
  "Relative forecast error (nRMSE)", 
  "Absolute forecast error (CRPS)"
)
names(metric.labs) <- c("RMSE.norm", "mean_crps")

# ====================
# Define Utility Functions
# ====================
	
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

# ====================
# Environment Setup Complete
# ====================
cat(paste0("\n", paste(rep("=", 50), collapse=""), "\n"))
cat("Microbial Forecasts Environment Setup Complete!\n")
cat("Project root:", project_root, "\n")
cat("Package status:", ifelse(microbial_pkg_loaded, "microbialForecast loaded", "microbialForecast not available"), "\n")
cat("Ready for analysis.\n")
cat(paste0(paste(rep("=", 50), collapse=""), "\n"))
