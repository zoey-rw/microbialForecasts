#!/usr/bin/env Rscript

# Dirichlet model convergence and consistency

source("source.R")

# Load the effects data
effects <- readRDS("data/model_outputs/dirichlet_regression/dirichlet_effects.rds")
convergence <- readRDS("data/model_outputs/dirichlet_regression/dirichlet_convergence.rds")

cat("=== DIRICHLET MODEL CONVERGENCE & CONSISTENCY ===\n\n")

# Overall convergence summary
cat("CONVERGENCE SUMMARY:\n")
cat("Total parameters analyzed:", nrow(effects), "\n")
cat("Significant effects:", sum(effects$significant, na.rm=TRUE), 
    "(", round(100*sum(effects$significant, na.rm=TRUE)/nrow(effects), 1), "%)\n")
cat("Parameters with Gelman-Rubin < 1.1:", sum(effects$gelman_rubin < 1.1, na.rm=TRUE), 
    "(", round(100*sum(effects$gelman_rubin < 1.1, na.rm=TRUE)/nrow(effects), 1), "%)\n")
cat("Parameters with ESS > 100:", sum(effects$effective_sample_size > 100, na.rm=TRUE), 
    "(", round(100*sum(effects$effective_sample_size > 100, na.rm=TRUE)/nrow(effects), 1), "%)\n")
cat("Parameters with ESS > 50:", sum(effects$effective_sample_size > 50, na.rm=TRUE), 
    "(", round(100*sum(effects$effective_sample_size > 50, na.rm=TRUE)/nrow(effects), 1), "%)\n\n")

# Model-by-model comparison
cat("MODEL-BY-MODEL COMPARISON:\n")
model_summary <- effects %>% 
  group_by(model_id) %>% 
  summarise(
    n_params = n(),
    n_sig = sum(significant, na.rm=TRUE),
    pct_sig = round(100*n_sig/n_params, 1),
    mean_effect = mean(abs(mean), na.rm=TRUE),
    mean_gelman = mean(gelman_rubin, na.rm=TRUE),
    max_gelman = max(gelman_rubin, na.rm=TRUE),
    mean_ess = mean(effective_sample_size, na.rm=TRUE),
    min_ess = min(effective_sample_size, na.rm=TRUE),
    n_converged = sum(gelman_rubin < 1.1, na.rm=TRUE),
    pct_converged = round(100*n_converged/n_params, 1)
  ) %>%
  arrange(desc(mean_effect))

print(model_summary)

# Parameter consistency analysis
cat("\nPARAMETER CONSISTENCY ANALYSIS:\n")
beta_params <- effects[grepl("beta", effects$parameter), ]

# Extract covariate and taxon numbers from parameter names
beta_params$cov_num <- as.numeric(gsub("beta\\[(\\d+),\\d+\\]", "\\1", beta_params$parameter))
beta_params$taxon_num <- as.numeric(gsub("beta\\[\\d+,(\\d+)\\]", "\\1", beta_params$parameter))

# Analyze by covariate
cat("\nEFFECT SIZES BY COVARIATE:\n")
covariate_summary <- beta_params %>% 
  group_by(cov_num) %>% 
  summarise(
    n_params = n(),
    n_sig = sum(significant, na.rm=TRUE),
    pct_sig = round(100*n_sig/n_params, 1),
    mean_effect = mean(abs(mean), na.rm=TRUE),
    max_effect = max(abs(mean), na.rm=TRUE),
    mean_gelman = mean(gelman_rubin, na.rm=TRUE),
    mean_ess = mean(effective_sample_size, na.rm=TRUE)
  ) %>%
  arrange(desc(mean_effect))

print(covariate_summary)

# Analyze by taxon
cat("\nEFFECT SIZES BY TAXON:\n")
taxon_summary <- beta_params %>% 
  group_by(taxon_num) %>% 
  summarise(
    n_params = n(),
    n_sig = sum(significant, na.rm=TRUE),
    pct_sig = round(100*n_sig/n_params, 1),
    mean_effect = mean(abs(mean), na.rm=TRUE),
    max_effect = max(abs(mean), na.rm=TRUE),
    mean_gelman = mean(gelman_rubin, na.rm=TRUE),
    mean_ess = mean(effective_sample_size, na.rm=TRUE)
  ) %>%
  arrange(desc(mean_effect))

print(taxon_summary)

# Cross-model consistency
cat("\nCROSS-MODEL CONSISTENCY:\n")
# Look at the same parameter across different models
cat("Largest effects across models:\n")
top_effects <- effects %>% 
  arrange(desc(effect_size)) %>% 
  head(15)

print(top_effects[, c("parameter", "model_id", "mean", "significant", "gelman_rubin", "effective_sample_size")])

# Convergence issues
cat("\nCONVERGENCE ISSUES:\n")
cat("Worst converged parameters (highest Gelman-Rubin):\n")
worst_convergence <- effects %>% 
  arrange(desc(gelman_rubin)) %>% 
  head(10)

print(worst_convergence[, c("parameter", "model_id", "gelman_rubin", "effective_sample_size", "significant")])

cat("\nLowest effective sample sizes:\n")
lowest_ess <- effects %>% 
  arrange(effective_sample_size) %>% 
  head(10)

print(lowest_ess[, c("parameter", "model_id", "gelman_rubin", "effective_sample_size", "significant")])

# Model type comparison
cat("\nMODEL TYPE COMPARISON:\n")
model_type_summary <- effects %>% 
  group_by(model_name) %>% 
  summarise(
    n_params = n(),
    n_sig = sum(significant, na.rm=TRUE),
    pct_sig = round(100*n_sig/n_params, 1),
    mean_effect = mean(abs(mean), na.rm=TRUE),
    mean_gelman = mean(gelman_rubin, na.rm=TRUE),
    mean_ess = mean(effective_sample_size, na.rm=TRUE),
    n_converged = sum(gelman_rubin < 1.1, na.rm=TRUE),
    pct_converged = round(100*n_converged/n_params, 1)
  )

print(model_type_summary)

# Time period comparison
cat("\nTIME PERIOD COMPARISON:\n")
time_summary <- effects %>% 
  group_by(time_period) %>% 
  summarise(
    n_params = n(),
    n_sig = sum(significant, na.rm=TRUE),
    pct_sig = round(100*n_sig/n_params, 1),
    mean_effect = mean(abs(mean), na.rm=TRUE),
    mean_gelman = mean(gelman_rubin, na.rm=TRUE),
    mean_ess = mean(effective_sample_size, na.rm=TRUE),
    n_converged = sum(gelman_rubin < 1.1, na.rm=TRUE),
    pct_converged = round(100*n_converged/n_params, 1)
  )

print(time_summary)

cat("\n=== SUMMARY ===\n")
cat("1. Convergence: ", round(100*sum(effects$gelman_rubin < 1.1, na.rm=TRUE)/nrow(effects), 1), "% of parameters converged\n")
cat("Effective sample sizes: ", round(100*sum(effects$effective_sample_size > 100, na.rm=TRUE)/nrow(effects), 1), "% have ESS > 100\n")
cat(round(100*sum(effects$significant, na.rm=TRUE)/nrow(effects), 1), "% of effects are significant\n")