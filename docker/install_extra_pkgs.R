#!/usr/bin/env Rscript
# Extra CRAN + GitHub packages for analysis scripts (05–07, source.R) not declared in
# microbialForecast DESCRIPTION alone.

options(Ncpus = max(1L, parallel::detectCores() - 1L))

cran_pkgs <- c(
  "pacman",
  "remotes",
  "nanoparquet",
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

install.packages(cran_pkgs, repos = "https://cloud.r-project.org", dependencies = NA)

# PLSR helpers used by microbialForecast/R/spectra_site_eff_permutation_fixed.r (not on CRAN)
remotes::install_github("plantphys/spectratrait", upgrade = "never", dependencies = NA)
