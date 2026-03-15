# Summarize MCMC output from all single-taxon truncated normal models
# Assumes input files have already had MCMC chains combined
# Run from project root: Rscript analysis/model_analysis/03_summarizeModelOutputs_truncNorm.r

source("../../source.R")
library(dplyr)
library(purrr)
library(data.table)
library(parallel)

# =============================================================================
# FIND SAMPLE FILES
# =============================================================================

base_path <- here("data/model_outputs/truncated_normal")

# Find all combined sample files (exclude chain files and deduplicate subdirs)
file.list <- list.files(base_path, pattern = "^samples_.*_with_legacy_covariate\\.rds$",
                        recursive = TRUE, full.names = TRUE)
file.list <- file.list[!grepl("_chain[0-9]", file.list)]

# Keep only the top-level copy when both top-level and subdir copies exist
# (prefer the shorter path, i.e. directly in env_cycl/, cycl_only/, env_cov/)
file.list <- file.list[order(nchar(file.list))]
file.list <- file.list[!duplicated(basename(file.list))]

cat("Truncated normal sample files found:", length(file.list), "\n")
for (f in file.list) cat(" ", basename(f), "\n")

# =============================================================================
# SUMMARIZE EACH MODEL FILE
# =============================================================================

max_cores <- min(4, max(1, detectCores() - 1))
use_sequential <- length(file.list) <= 6

if (use_sequential) {
  cat("Using sequential processing (", length(file.list), " files)\n", sep = "")
  results <- lapply(file.list, function(f) {
    cat("Processing:", basename(f), "\n")
    tryCatch(
      microbialForecast::summarize_beta_model(f, save_summary = TRUE, overwrite = TRUE, drop_other = TRUE),
      error = function(e) { cat("  ERROR:", e$message, "\n"); NULL }
    )
  })
} else {
  cat("Using parallel processing with", max_cores, "cores\n")
  cl <- makeCluster(max_cores)
  clusterEvalQ(cl, { library(microbialForecast); library(here) })
  results <- tryCatch(
    parLapply(cl, file.list, function(f) {
      tryCatch(
        microbialForecast::summarize_beta_model(f, save_summary = TRUE, overwrite = TRUE, drop_other = TRUE),
        error = function(e) NULL
      )
    }),
    error = function(e) {
      cat("Parallel failed, falling back to sequential:", e$message, "\n")
      lapply(file.list, function(f)
        tryCatch(microbialForecast::summarize_beta_model(f, save_summary = TRUE, overwrite = TRUE, drop_other = TRUE),
                 error = function(e) NULL))
    }
  )
  stopCluster(cl)
}

# =============================================================================
# COLLECT SUMMARY FILES
# =============================================================================

summary_file_list <- list.files(base_path, pattern = "^summary_.*_with_legacy_covariate\\.rds$",
                                recursive = TRUE, full.names = TRUE)
summary_file_list <- summary_file_list[!grepl("_chain[0-9]", summary_file_list)]
summary_file_list <- summary_file_list[order(nchar(summary_file_list))]
summary_file_list <- summary_file_list[!duplicated(basename(summary_file_list))]

cat("Summary files found:", length(summary_file_list), "\n")

file_summaries <- lapply(summary_file_list, function(f) {
  tryCatch(readRDS(f), error = function(e) { cat("Error reading:", basename(f), "\n"); NULL })
})
file_summaries <- file_summaries[!sapply(file_summaries, is.null)]

if (length(file_summaries) == 0) {
  stop("No valid summaries produced. Check model files and summarize_beta_model() errors above.")
}

summary_df <- map_df(file_summaries, ~ .x[[1]])
plot_est   <- map_df(file_summaries, ~ if (is.data.frame(.x[[2]]) && nrow(.x[[2]]) > 0) .x[[2]] else data.frame())
gelman_list <- map_df(file_summaries, ~ if (length(.x) >= 3 && is.data.frame(.x[[3]])) .x[[3]] else data.frame())

cat("Summary rows:", nrow(summary_df), "\n")
cat("Plot estimate rows:", nrow(plot_est), "\n")
cat("Parameters found:", paste(unique(gsub("\\[.*", "", summary_df$rowname)), collapse = ", "), "\n")

# =============================================================================
# SAVE COMBINED SUMMARIES
# =============================================================================

saveRDS(list(summary_df = summary_df, plot_est = plot_est, gelman_list = gelman_list),
        here("data/summary/truncated_normal_summaries.rds"))
cat("Saved: data/summary/truncated_normal_summaries.rds\n")

# =============================================================================
# EXTRACT PREDICTOR EFFECTS (beta parameters)
# =============================================================================

# Mirror the pattern from 04_tidyEffectSizes.r
beta_effects <- summary_df %>%
  filter(grepl("^beta", rowname)) %>%
  mutate(beta = as.character(beta))

cat("Beta effect rows:", nrow(beta_effects), "\n")

saveRDS(beta_effects, here("data/summary/truncated_normal_predictor_effects.rds"))
cat("Saved: data/summary/truncated_normal_predictor_effects.rds\n")
