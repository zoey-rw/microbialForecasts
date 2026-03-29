# Microbial Forecasts Environment Setup
# Sourced by all analysis and figure scripts

# Establish project root — works whether sourced from project root or subdirectories
here::i_am("source.R")

# Prefer project R_library (e.g. when getwd() is analysis/model_analysis, .Rprofile is not read)
proj_lib <- here::here("R_library")
if (dir.exists(proj_lib)) {
  .libPaths(c(proj_lib, .libPaths()))
}

# Load microbialForecast package (provides helper functions and global variables)
library(microbialForecast)

# Core packages
suppressPackageStartupMessages({
  library(tidyverse, warn.conflicts = FALSE)
  library(here)
  library(data.table)
})

# Analysis packages
pacman::p_load(nimble, coda, lubridate,
               doParallel, Rfast, moments,
               scoringRules, Metrics, ggpubr, ggallin)

# Plotting labels (used by figure and analysis scripts)
model.labs <- c("Environmental\npredictors", "Seasonality",
                "Environmental predictors\n& seasonality")
names(model.labs) <- c("env_cov", "cycl_only", "env_cycl")

metric.labs <- c("Relative forecast error (nRMSE)", "Absolute forecast error (CRPS)")
names(metric.labs) <- c("RMSE.norm", "mean_crps")

# Ensure key output directories exist
for (d in c(here("data", "model_outputs"), here("data", "summary"), here("figures"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
