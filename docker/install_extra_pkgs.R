#!/usr/bin/env Rscript
# Extra CRAN + GitHub packages for analysis scripts (05–07, source.R) not declared in
# microbialForecast DESCRIPTION alone.

options(Ncpus = max(1L, parallel::detectCores() - 1L))

cran_pkgs <- c(
  "pacman",
  "remotes",
  "nanoparquet",
  "duckdb",     # required by microbialForecast::load_hindcasts()
  "MuMIn",
  "pls",
  "forestplot",
  "gridExtra",
  "ggforce",
  "ggallin",
  "plotrix",
  "foreach",
  "anytime",
  "padr",
  "matrixStats",
  "Rfast",
  "RcppZiggurat",
  "logger",
  "truncnorm"
)

# Use the image's pinned repository snapshot (rocker/tidyverse sets a dated Posit
# Package Manager CRAN snapshot) for reproducible versions; do NOT override it with
# the always-latest cloud.r-project.org.
install.packages(cran_pkgs, dependencies = NA)

# PLSR helpers used by microbialForecast/R/spectra_site_eff_permutation_fixed.r (not on CRAN).
# Pinned to spectratrait v1.2.6 for reproducible rebuilds.
remotes::install_github("plantphys/spectratrait@d159c5dda9d6739105d95b357e48f7754daf4978",
                        upgrade = "never", dependencies = NA)

# speedyseq: faster phyloseq accessors (e.g. psmelt) used by the phylogeny scripts.
# The scripts fall back to phyloseq when it is absent, but it is installed here for
# the faster path. Pinned to v0.2.0 for reproducible rebuilds.
remotes::install_github("mikemc/speedyseq@66c527e6e90f70de4a02e4ae5fbb08996ad2febf",
                        upgrade = "never", dependencies = NA)
