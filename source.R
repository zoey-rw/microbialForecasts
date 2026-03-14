# ====================
# Microbial Forecasts Environment Setup
# ====================
# Source script for NIMBLE models and analysis environment
# This script sets up the R environment for microbial forecasting analysis

# Automatically detect project root directory
if (!exists("project_root")) {
  # Try multiple methods to find project root
  rstudio_available <- requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable() &&
    !is.null(rstudioapi::getActiveProject())
  if (rstudio_available) {
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

tryCatch({
  pacman::p_load(char = required_packages, character.only = TRUE)
  cat("All required packages loaded successfully\n")
}, error = function(e) {
  cat("Warning: Some packages could not be loaded:", e$message, "\n")
  cat("Continuing with available packages...\n")
})

# Arrow package is optional
tryCatch({
  if (!require(arrow, quietly = TRUE)) {
    cat("Note: arrow package not available (optional)\n")
  }
}, error = function(e) {
  # Silently ignore arrow loading errors
})

# ====================
# Create Directory Structure
# ====================
dirs_to_create <- c(
  here("data", "model_outputs"),
  here("data", "model_outputs", "cloglog_beta_driver_uncertainty"),
  here("data", "model_outputs", "CLR_regression"),
  here("data", "summary"),
  here("figures")
)

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
# Utility Functions (not suitable for package)
# ====================

# Load rsq function from external source
source("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/NEFI_functions/rsq_1.1.r")

first = function(x) x %>% nest %>% ungroup %>% slice(1) %>% unnest(data)

quick_get_rank_df = function(k = 1,
                             min.date = "20151101",
                             max.date = "20200101"){
	pacman::p_load(reshape2, parallel, nimble, coda, tidyverse)
	rank.name <- microbialForecast:::tax_names[k]
	cal <- c(readRDS(here("data", "clean/cal_groupAbundances_16S_2021.rds")),
	         readRDS(here("data", "clean/cal_groupAbundances_ITS_2021.rds")))
	val <- c(readRDS(here("data", "clean/val_groupAbundances_16S_2021.rds")),
	         readRDS(here("data", "clean/val_groupAbundances_ITS_2021.rds")))
	cal.rank.df <- cal[[rank.name]]
	val.rank.df <- val[[rank.name]]
	rank.df <- rbind(cal.rank.df, val.rank.df)
	print(paste0("Preparing model data for ", rank.name))
	model.dat <- prepTaxonomicData(rank.df = rank.df,
	                               min.prev = 3,
	                               min.date = min.date,
	                               max.date = max.date)
	return(model.dat)
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
		                 })
	}

	samples <- metadata <- list()
	first_iter <- last_iter <- list()
	for(i in 1:length(input_list)){
		print(i)
		if (class(input_list[[i]])=="character") {
			chain <- readInputRdsFile(chain_paths[[i]])
			if (any(is.na(chain))) next()
			samples[[i]] <- chain[[1]]
		} else {
			samples[[i]]  = input_list[[i]]
		}
	}

	samples<-samples[!sapply(samples,is.null)]

	nrows <- lapply(samples, nrow) %>% unlist()
	min_nrow <- min(nrows)
	for(i in 1:length(samples)){
		current_nrow <- nrow(samples[[i]])
		if (min_nrow < current_nrow){
			samples[[i]] <- window_chain(samples[[i]], max_size = (min_nrow-1))
		}
	}

	for (i in 1:length(samples)) {
		attr(samples[[i]], "mcpar") = attr(samples[[1]], "mcpar")
	}

	out <- as.mcmc.list(samples)
	return(out)
}

# Only define nimbleFunction if nimble is available
if (requireNamespace("nimble", quietly = TRUE)) {
	rLogitBeta <- nimbleFunction (
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
} else {
	rLogitBeta <- function(n = 1, shape1 = 1.0, shape2 = 1.0) {
		warning("nimble not available; rLogitBeta using basic R implementation")
		y <- rbeta(1, shape1=shape1, shape2=shape2)
		x <- qlogis(y)
		return(x)
	}
}

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

# ====================
# Environment Setup Complete
# ====================
cat(paste0("\n", paste(rep("=", 50), collapse=""), "\n"))
cat("Microbial Forecasts Environment Setup Complete!\n")
cat("Project root:", project_root, "\n")
cat("Package status:", ifelse(microbial_pkg_loaded, "microbialForecast loaded", "microbialForecast not available"), "\n")
cat("Ready for analysis.\n")
cat(paste0(paste(rep("=", 50), collapse=""), "\n"))
