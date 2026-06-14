#!/usr/bin/env Rscript
# Quick check that the image can load the project stack (no data required).

options(warn = 1)

root <- Sys.getenv("MICROBIAL_FORECASTS_ROOT", "/opt/microbialForecasts")
setwd(root)

source(file.path(root, "source.R"))

req <- c(
  "microbialForecast", "here", "tidyverse", "nimble", "data.table",
  "doParallel", "foreach", "nanoparquet", "duckdb", "MuMIn", "pls",
  "spectratrait", "pacman", "ggallin"
)
missing <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing: ", paste(missing, collapse = ", "))
}

cat("smoke_test: OK\n")
