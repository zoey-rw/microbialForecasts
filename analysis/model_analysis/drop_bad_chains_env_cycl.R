#!/usr/bin/env Rscript
# Drop one bad chain from 24 env_cycl taxa where it fixes convergence (max Rhat < 1.1)
# Chain-drop mapping was computed by leave-one-out Gelman diagnostic analysis.
# After dropping, recompute Gelman diagnostics and param_summary, then re-summarize.

source("../../source.R")
library(coda)
library(dplyr)

drop_map <- readRDS(here("data/model_outputs/cloglog_beta_driver_uncertainty/env_cycl_chain_drop_map.rds"))
cat("Processing", length(drop_map), "taxa with chain drops\n\n")

for (taxon in names(drop_map)) {
  info <- drop_map[[taxon]]
  f <- here(info$file)
  drop_chain <- info$drop

  cat(sprintf("%-35s: dropping chain %d (Rhat %.2f -> %.2f) ... ",
              taxon, drop_chain, info$before, info$after))

  # Read combined samples
  s <- readRDS(f)
  samples <- s$samples

  if (!inherits(samples, "mcmc.list") || length(samples) < 3) {
    cat("SKIP (not mcmc.list or < 3 chains)\n")
    next
  }

  # Drop the bad chain
  s$samples <- samples[-drop_chain]

  # Also drop from samples2 if present
  if ("samples2" %in% names(s) && inherits(s$samples2, "mcmc.list")) {
    s$samples2 <- s$samples2[-drop_chain]
  }

  # Recompute param_summary from remaining chains
  remaining <- s$samples
  param_means <- do.call(rbind, lapply(remaining, function(ch) colMeans(as.matrix(ch))))
  param_sds <- do.call(rbind, lapply(remaining, function(ch) apply(as.matrix(ch), 2, sd)))
  overall_means <- colMeans(param_means)
  overall_sds <- colMeans(param_sds)

  s$param_summary <- list(
    means = data.frame(Mean = overall_means, SD = overall_sds, param = names(overall_means)),
    quantiles = as.data.frame(t(apply(
      do.call(rbind, lapply(remaining, as.matrix)),
      2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)
    )))
  )
  s$param_summary$quantiles$param <- rownames(s$param_summary$quantiles)

  # Recompute Gelman diagnostics
  gd <- tryCatch({
    gelman.diag(remaining, multivariate = FALSE)$psrf
  }, error = function(e) {
    matrix(NA, nrow = ncol(remaining[[1]]), ncol = 2,
           dimnames = list(colnames(remaining[[1]]), c("Point est.", "Upper C.I.")))
  })
  s$gelman <- gd

  # Update metadata
  s$metadata$chains_dropped <- drop_chain
  s$metadata$original_nchain <- length(samples)
  s$metadata$nchain <- length(remaining)

  # Save updated combined file
  saveRDS(s, f)

  # Verify
  new_max <- max(gd[grep("^(beta|intercept|rho|precision|site_effect_sd)", rownames(gd)), "Point est."], na.rm = TRUE)
  cat(sprintf("saved (%d chains, max Rhat %.3f)\n", length(remaining), new_max))
}

cat("\nDone. Now re-run summaries with:\n")
cat("  Rscript 03_summarize_env_cycl_relabel.r  # or equivalent\n")
cat("  Rscript 03b_combine_summaries.r\n")
