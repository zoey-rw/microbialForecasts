#!/usr/bin/env Rscript
# Create warm-start initialization file for the truncated normal ascomycota env_cycl model
# by pooling the two existing chain files and computing parameter means.
#
# Usage (from project root):
#   Rscript analysis/model_analysis/create_truncnorm_warmstart.R

here::i_am("source.R")
library(here)

chain_dir <- here("data", "model_outputs", "truncated_normal", "env_cycl", "ascomycota")

chain1_path <- file.path(chain_dir, "samples_env_cycl_ascomycota_20130601_20180101_with_legacy_covariate_chain1.rds")
chain2_path <- file.path(chain_dir, "samples_env_cycl_ascomycota_20130601_20180101_with_legacy_covariate_chain2.rds")

for (p in c(chain1_path, chain2_path)) {
  if (!file.exists(p)) stop("Chain file not found: ", p)
}

chain1 <- readRDS(chain1_path)
chain2 <- readRDS(chain2_path)

cat("Chain1 samples dims:", dim(chain1$samples), "\n")
cat("Chain2 samples dims:", dim(chain2$samples), "\n")
cat("Chain1 iterations:", chain1$metadata$niter, "\n")
cat("Chain2 iterations:", chain2$metadata$niter, "\n")

# Verify both chains have identical column names
stopifnot(identical(colnames(chain1$samples), colnames(chain2$samples)))

# Pool the two chains row-wise and compute means
pooled <- rbind(chain1$samples, chain2$samples)
cat("Pooled dims:", dim(pooled), "\n")

param_means <- colMeans(pooled, na.rm = TRUE)
cat("Computed param_means for", length(param_means), "parameters\n")
cat("Key parameter means:\n")
key_params <- c("core_sd", "sigma", "rho", "intercept", "legacy_effect", "site_effect_sd",
                "beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", "beta[6]", "beta[7]", "beta[8]")
print(round(param_means[key_params[key_params %in% names(param_means)]], 4))

warmstart <- list(
  param_means       = param_means,
  source_iterations = chain1$metadata$niter + chain2$metadata$niter,
  source_chains     = 2L,
  colnames          = colnames(pooled)
)

out_path <- file.path(chain_dir, "warmstart_inits.rds")
saveRDS(warmstart, out_path)

cat("\nSaved warm-start file to:", out_path, "\n")
cat("source_iterations:", warmstart$source_iterations, "\n")
cat("source_chains:", warmstart$source_chains, "\n")
cat("To use: WARMSTART_FILE=", out_path, " Rscript analysis/model_analysis/01_fitModels_truncNorm.R --species ascomycota --model-name env_cycl\n", sep = "")
