#!/usr/bin/env Rscript
# Create warm-start initial values for 86 poorly converged env_cycl taxa
# that can't be fixed by dropping a single chain.
#
# For each taxon, extract posterior means from existing chains and save as
# initial values for continued MCMC sampling on SCC.
#
# Usage: Rscript create_env_cycl_warmstart.R
# Output: data/model_outputs/cloglog_beta_driver_uncertainty/env_cycl/{taxon}/warmstart_inits.rds

source("../../source.R")
library(coda)

# Identify the 86 still-poor taxa (max Rhat >= 1.2 even after best chain drop)
drop_map <- readRDS(here("data/model_outputs/cloglog_beta_driver_uncertainty/env_cycl_chain_drop_map.rds"))
fixable_taxa <- names(drop_map)

sample_files <- list.files(here("data/model_outputs/cloglog_beta_driver_uncertainty/env_cycl"),
                           pattern = "^samples_env_cycl_", full.names = TRUE)
summary_files <- list.files(here("data/model_outputs/cloglog_beta_driver_uncertainty/env_cycl"),
                            pattern = "^summary_", full.names = TRUE)

# Get all poorly converged taxa
poor_taxa <- c()
for (f in summary_files) {
  taxon <- sub("^summary_env_cycl_", "", sub("_20130601.*$", "", basename(f)))
  s <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(s)) next
  for (i in 3:min(4, length(s))) {
    if (is.data.frame(s[[i]]) && "Point est." %in% names(s[[i]])) {
      gd <- s[[i]]
      major <- gd[grep("^(beta|intercept|rho|precision|site_effect_sd|sigma_proc|legacy)", gd$parameter), ]
      if (nrow(major) > 0 && max(major[["Point est."]], na.rm = TRUE) >= 1.2) {
        poor_taxa <- c(poor_taxa, taxon)
      }
      break
    }
  }
}

# Exclude the 24 fixable taxa
rerun_taxa <- setdiff(poor_taxa, fixable_taxa)
cat("Taxa needing warm-start re-runs:", length(rerun_taxa), "\n\n")

created <- 0
for (taxon in rerun_taxa) {
  f <- grep(paste0("samples_env_cycl_", taxon, "_"), sample_files, value = TRUE, fixed = FALSE)
  if (length(f) != 1) next

  s <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(s)) next

  samples <- s$samples
  if (!inherits(samples, "mcmc.list")) next

  # Extract posterior means from the LAST values of each chain (better starting point)
  # Use last 1000 iterations averaged across chains
  n_iter <- nrow(samples[[1]])
  tail_start <- max(1, n_iter - 999)

  all_chain_means <- lapply(samples, function(ch) {
    colMeans(as.matrix(ch)[tail_start:n_iter, , drop = FALSE])
  })
  # Average across chains
  param_names <- names(all_chain_means[[1]])
  param_means <- setNames(
    sapply(param_names, function(p) mean(sapply(all_chain_means, function(m) m[p]))),
    param_names
  )

  # Build inits list matching create_inits_with_uncertainty structure
  inits <- list()

  # Scalar parameters
  for (p in c("precision", "rho", "intercept", "legacy_effect", "sigma_proc",
              "sigma_init", "site_effect_sd")) {
    if (p %in% param_names) {
      inits[[p]] <- param_means[p]
    }
  }

  # Beta vector
  beta_idx <- grep("^beta\\[", param_names)
  if (length(beta_idx) > 0) {
    inits$beta <- param_means[beta_idx]
    names(inits$beta) <- NULL  # NIMBLE wants unnamed vectors
  }

  # Site effects vector
  site_idx <- grep("^site_effect\\[", param_names)
  if (length(site_idx) > 0) {
    inits$site_effect <- param_means[site_idx]
    names(inits$site_effect) <- NULL
  }

  # Eta matrix (plot x time)
  eta_idx <- grep("^eta\\[", param_names)
  if (length(eta_idx) > 0) {
    # Parse dimensions from parameter names like eta[1, 1]
    eta_coords <- regmatches(param_names[eta_idx],
                             regexpr("\\d+,\\s*\\d+", param_names[eta_idx]))
    if (length(eta_coords) > 0) {
      rows <- as.integer(sub(",.*", "", eta_coords))
      cols <- as.integer(sub(".*,\\s*", "", eta_coords))
      eta_mat <- matrix(NA, nrow = max(rows), ncol = max(cols))
      for (k in seq_along(eta_idx)) {
        eta_mat[rows[k], cols[k]] <- param_means[eta_idx[k]]
      }
      inits$eta <- eta_mat
    }
  }

  # Driver uncertainty variables (temp_est, mois_est, pH_est_p, pC_est_p)
  for (driver in c("temp_est", "mois_est")) {
    d_idx <- grep(paste0("^", driver, "\\["), param_names)
    if (length(d_idx) > 0) {
      d_coords <- regmatches(param_names[d_idx],
                             regexpr("\\d+,\\s*\\d+", param_names[d_idx]))
      if (length(d_coords) > 0) {
        rows <- as.integer(sub(",.*", "", d_coords))
        cols <- as.integer(sub(".*,\\s*", "", d_coords))
        d_mat <- matrix(NA, nrow = max(rows), ncol = max(cols))
        for (k in seq_along(d_idx)) {
          d_mat[rows[k], cols[k]] <- param_means[d_idx[k]]
        }
        inits[[driver]] <- d_mat
      }
    }
  }
  for (driver in c("pH_est_p", "pC_est_p")) {
    d_idx <- grep(paste0("^", driver, "\\["), param_names)
    if (length(d_idx) > 0) {
      inits[[driver]] <- param_means[d_idx]
      names(inits[[driver]]) <- NULL
    }
  }

  # Save warm-start
  ws <- list(
    inits = inits,
    source_iterations = s$metadata$niter,
    source_chains = length(samples),
    taxon = taxon,
    model_name = "env_cycl",
    created = Sys.time()
  )

  out_dir <- here("data/model_outputs/cloglog_beta_driver_uncertainty/env_cycl", taxon)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_file <- file.path(out_dir, "warmstart_inits.rds")
  saveRDS(ws, out_file)
  created <- created + 1
}

cat("\nCreated", created, "warm-start files\n")
cat("Re-fit: 01_fitModels.R auto-loads warmstart_inits.rds from each taxon directory\n")
cat("  (env_cycl/<taxon>/warmstart_inits.rds), or set WARMSTART_FILE to an explicit path.\n")
cat("Optional: SKIP_WARMSTART=true cold start; WARMSTART_JITTER_SD (default 0.05).\n")
cat("Suggested: ~50k additional iterations per chain (or use existing ESS continuation).\n")
