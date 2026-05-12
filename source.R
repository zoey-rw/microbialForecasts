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

# Canonical color palette (Wong colorblind-friendly).
# Bacteria/Fungi use orange/blue EXCLUSIVELY across the manuscript.
# Taxonomic/Functional use a distinct green/pink pair.
# Model types use sky-blue / green / vermillion.
# Any figure that introduces a new color scale should reuse one of these or
# pick from the unused Wong colors (#F0E442 yellow, #000000 black).
kingdom_colors  <- c("Bacteria"   = "#E69F00", "Fungi"      = "#0072B2")
fcast_type_colors <- c("Taxonomic" = "#009E73", "Functional" = "#CC79A7")
model_colors    <- c("env_cov"    = "#56B4E9",
                     "cycl_only"  = "#009E73",
                     "env_cycl"   = "#D55E00")

# Ensure key output directories exist
for (d in c(here("data", "model_outputs"), here("data", "summary"), here("figures"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
