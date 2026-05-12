#!/usr/bin/env Rscript
# Combine chain files for a single taxon, summarize, and optionally clean up.
# Usage: Rscript combine_single_taxon.R <chain_dir> <output_file> [--cleanup]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript combine_single_taxon.R <chain_dir> <output_file> [--cleanup]")

chain_dir <- args[1]
output_file <- args[2]
cleanup <- "--cleanup" %in% args

source("../../source.R")
library(coda)

all_chain_files <- sort(list.files(chain_dir, pattern = "samples_.*chain[0-9]", full.names = TRUE))
cat("All chain files:", length(all_chain_files), "\n")

if (length(all_chain_files) < 2) stop("Need at least 2 chain files")

# If multiple batches exist (different dates/sizes), prefer the newer/larger files.
# Group by chain number and keep the largest file for each chain.
chain_nums <- as.integer(gsub(".*chain(\\d+).*", "\\1", basename(all_chain_files)))
chain_sizes <- file.info(all_chain_files)$size
chain_df <- data.frame(file = all_chain_files, chain = chain_nums, size = chain_sizes,
                       stringsAsFactors = FALSE)
# For each chain number, keep the largest file (= most iterations)
chain_files <- character()
for (cn in sort(unique(chain_df$chain))) {
  sub <- chain_df[chain_df$chain == cn, ]
  best <- sub$file[which.max(sub$size)]
  chain_files <- c(chain_files, best)
  if (nrow(sub) > 1) {
    cat(sprintf("  Chain %d: %d versions, keeping %s (%s MB)\n",
                cn, nrow(sub), basename(best), round(max(sub$size)/1e6)))
  }
}
cat("Selected", length(chain_files), "chain files (1 per chain number)\n")

if (length(chain_files) < 2) stop("Need at least 2 chains after deduplication")

# Read all chains as matrices
chains_raw <- lapply(chain_files, function(f) {
  cat("  Reading", basename(f), "\n")
  s <- readRDS(f)
  if ("samples" %in% names(s)) as.matrix(s$samples)
  else if (is.matrix(s)) s
  else NULL
})
chains_raw <- Filter(Negate(is.null), chains_raw)
iters <- sapply(chains_raw, nrow)
cat("Iterations per chain:", iters, "\n")

# Strategy: prefer longer (newer) chains. If >=3 chains share the
# max iteration count, use only those. Otherwise keep all and truncate.
max_iter <- max(iters)
min_iter <- min(iters)
long_idx <- which(iters == max_iter)
short_idx <- which(iters < max_iter)

if (length(long_idx) >= 3) {
  # Enough long chains — drop the shorter (older) ones
  cat("Using", length(long_idx), "chains with", max_iter,
      "iterations (dropping", length(short_idx), "shorter chains)\n")
  chains_raw <- chains_raw[long_idx]
  chain_files <- chain_files[long_idx]
  iters <- iters[long_idx]
  min_iter <- max_iter
} else if (length(short_idx) > 0) {
  cat("Keeping all", length(chains_raw), "chains,",
      "truncating to last", min_iter, "iterations\n")
}
chains <- lapply(chains_raw, function(ch) {
  n <- nrow(ch)
  as.mcmc(ch[(n - min_iter + 1):n, ])
})
combined <- as.mcmc.list(chains)

# Gelman diagnostics on major params
param_cols <- grep("^(beta|intercept|rho|precision|site_effect_sd|sigma_proc|legacy)",
                   colnames(chains[[1]]), value = TRUE)
gd <- tryCatch(
  gelman.diag(combined[, param_cols], multivariate = FALSE)$psrf,
  error = function(e) { cat("Gelman failed:", e$message, "\n"); NULL }
)
if (!is.null(gd)) {
  max_rhat <- max(gd[, "Point est."], na.rm = TRUE)
  cat(sprintf("Max Rhat (major params): %.3f\n", max_rhat))
  cat("Converged (<1.1):", max_rhat < 1.1, "\n")
}

# Build output from first chain file as template
template <- readRDS(chain_files[1])
template$samples <- combined

# Recompute param_summary
all_s <- do.call(rbind, lapply(combined, as.matrix))
template$param_summary <- list(
  means = data.frame(Mean = colMeans(all_s), SD = apply(all_s, 2, sd), param = colnames(all_s)),
  quantiles = as.data.frame(t(apply(all_s, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))))
)
template$param_summary$quantiles$param <- colnames(all_s)
template$gelman <- gd
template$metadata$nchain <- length(chains)
template$metadata$niter <- min_iter

saveRDS(template, output_file)
cat("Saved:", round(file.info(output_file)$size / 1e6), "MB\n")

# Re-summarize using the package function
cat("Summarizing...\n")
tryCatch(
  summarize_beta_model(output_file, save_summary = TRUE, overwrite = TRUE),
  error = function(e) cat("Summary failed:", e$message, "\n")
)

# Clean up chain files if requested
if (cleanup) {
  file.remove(chain_files)
  cat("Cleaned up", length(chain_files), "chain files\n")
}

cat("Done\n")
