# Summarize MCMC output from all single-taxon truncated normal models
# Pipeline order: run 02_combineModelChains_truncNorm.r first to merge per-chain RDS into
# combined sample files; then run this script. (This script can still combine or wrap
# chain-only outputs when no combined file exists, but production flow is 02 then 03.)
# Run from analysis/model_analysis/:
#   Rscript 02_combineModelChains_truncNorm.r
#   Rscript 03_summarizeModelOutputs_truncNorm.r

source("../../source.R")
library(dplyr)
library(purrr)
library(data.table)
library(coda)

# =============================================================================
# FIND AND ORGANIZE MODEL FILES
# =============================================================================

base_path <- here("data/model_outputs/truncated_normal")

# Find all sample files (both combined and per-chain)
all_files <- list.files(base_path, pattern = "^samples_.*_with_legacy_covariate.*\\.rds$",
                        recursive = TRUE, full.names = TRUE)

# Separate combined files and chain files
chain_files <- all_files[grepl("_chain[0-9]+\\.rds$", all_files)]
combined_files <- all_files[!grepl("_chain[0-9]+\\.rds$", all_files)]

# Deduplicate: prefer shorter paths (top-level over subdirectory copies)
combined_files <- combined_files[order(nchar(combined_files))]
combined_files <- combined_files[!duplicated(basename(combined_files))]

cat("Combined sample files found:", length(combined_files), "\n")
cat("Chain files found:", length(chain_files), "\n")

# Group chain files by model ID basename (strip _chainN suffix and path)
# This ensures chain1 in env_cycl/ and chain2 in env_cycl/ascomycota/ get grouped together
chain_basenames <- gsub("_chain[0-9]+\\.rds$", ".rds", basename(chain_files))
chain_groups <- split(chain_files, chain_basenames)

# Identify models that have chain files but no combined file
combined_basenames <- basename(combined_files)
models_needing_combination <- list()
for (model_basename in names(chain_groups)) {
  if (!model_basename %in% combined_basenames) {
    models_needing_combination[[model_basename]] <- chain_groups[[model_basename]]
  }
}

cat("Models needing chain combination:", length(models_needing_combination), "\n")

# =============================================================================
# COMBINE CHAINS FOR MODELS WITHOUT COMBINED FILES
# =============================================================================

if (length(models_needing_combination) > 0) {
  for (model_basename in names(models_needing_combination)) {
    chain_paths <- models_needing_combination[[model_basename]]
    # Deduplicate chain paths by basename (same chain in top-level and subdirectory)
    chain_paths <- chain_paths[order(nchar(chain_paths))]
    chain_paths <- chain_paths[!duplicated(basename(chain_paths))]

    if (length(chain_paths) < 2) {
      cat("  Skipping", model_basename, "- only", length(chain_paths), "chain(s)\n")
      next
    }

    cat("  Combining", length(chain_paths), "chains for", model_basename, "\n")

    combined <- tryCatch(
      microbialForecast::combine_chains(chain_paths),
      error = function(e) { cat("    ERROR:", e$message, "\n"); NULL }
    )

    if (!is.null(combined)) {
      # Save to the directory of the shortest chain path
      save_path <- file.path(dirname(chain_paths[1]), model_basename)
      tryCatch({
        saveRDS(combined, save_path)
        cat("    Saved combined file:", save_path, "\n")
        combined_files <- c(combined_files, save_path)
      }, error = function(e) {
        cat("    ERROR saving:", e$message, "\n")
      })
    }
  }
}

# Also handle models with only a single chain: create a combined-format file
# so summarize_beta_model sees a clean model_id without _chainN suffix.
for (model_basename in names(chain_groups)) {
  if (!model_basename %in% basename(combined_files)) {
    chain_paths <- chain_groups[[model_basename]]
    chain_paths <- chain_paths[order(nchar(chain_paths))]
    chain_paths <- chain_paths[!duplicated(basename(chain_paths))]
    if (length(chain_paths) == 1) {
      cat("  Wrapping single chain as combined file:", model_basename, "\n")
      save_path <- file.path(dirname(chain_paths[1]), model_basename)
      tryCatch({
        chain_data <- readRDS(chain_paths[1])
        # Wrap single chain samples in mcmc.list for consistency
        if (is.list(chain_data) && "samples" %in% names(chain_data)) {
          samp <- chain_data$samples
          if (is.matrix(samp)) samp <- coda::mcmc(samp)
          if (!coda::is.mcmc.list(samp)) samp <- coda::as.mcmc.list(list(samp))

          s2 <- chain_data$samples2
          if (!is.null(s2) && is.matrix(s2)) {
            s2_mcmc <- coda::mcmc(s2)
            attr(s2_mcmc, "mcpar") <- attr(samp[[1]], "mcpar")
            s2 <- coda::as.mcmc.list(list(s2_mcmc))
          } else {
            s2 <- list()
          }

          wrapped <- list(samples = samp, samples2 = s2,
                          metadata = chain_data$metadata)
          saveRDS(wrapped, save_path)
          combined_files <- c(combined_files, save_path)
          cat("    Saved:", save_path, "\n")
        }
      }, error = function(e) {
        cat("    ERROR wrapping single chain:", e$message, "\n")
      })
    }
  }
}

# Final file list: all combined files (including newly created ones)
# Deduplicate by model ID basename
file.list <- combined_files
file.list <- file.list[order(nchar(file.list))]
file.list <- file.list[!duplicated(basename(file.list))]

cat("\nFinal file list for summarization:", length(file.list), "\n")
for (f in file.list) cat("  ", basename(f), "\n")

# =============================================================================
# SUMMARIZE EACH MODEL FILE
# =============================================================================

results <- lapply(file.list, function(f) {
  cat("Processing:", basename(f), "\n")
  tryCatch(
    microbialForecast::summarize_beta_model(f, save_summary = TRUE,
                                            overwrite = TRUE, drop_other = TRUE),
    error = function(e) { cat("  ERROR:", e$message, "\n"); NULL }
  )
})

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

# Gelman diagnostics are in element 4 for summarize_beta_model output
gelman_list <- map_df(file_summaries, ~ {
  # summarize_beta_model returns list(summary_df, pred.means, pred.quantiles, gd)
  # Element 4 is gelman diagnostics
  if (length(.x) >= 4 && is.data.frame(.x[[4]])) .x[[4]] else data.frame()
})

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

# Extract beta effects from summary_df (populated by summarize_beta_model)
beta_effects <- summary_df %>% filter(grepl("^beta", rowname))

cat("Beta rows found in summary_df:", nrow(beta_effects), "\n")

if (nrow(beta_effects) > 0) {
  # Ensure beta column is character for downstream joins
  if ("beta" %in% colnames(beta_effects)) {
    beta_effects <- beta_effects %>% mutate(beta = as.character(beta))
  }
} else {
  # Fallback: extract directly from chain files if summarize_beta_model
  # did not produce beta rows (should not happen with the fix, but kept as safety)
  cat("No beta rows in summary_df; extracting from chain files directly\n")
  # Use the same covariate key mapping as summarize_beta_model (via all_covariates_key)
  # so that the primary and fallback paths produce consistent labels.
  # NOTE: For env_cycl, beta indices 1-8 are sin,cos,temp,mois,pH,pC,relEM,LAI
  # in the model, but all_covariates_key maps "1"="Temperature","7"="sin","8"="cos".
  # This is a pre-existing mislabeling that is consistent across all pipelines
  # (cloglog, CLR, truncnorm). Do NOT fix here -- fix globally if fixing at all.
  cov_key_env_cycl <- c("Temperature", "Moisture", "pH", "pC",
                         "Ectomycorrhizal trees", "LAI", "sin", "cos")
  cov_key_env_cov <- c("Temperature", "Moisture", "pH", "pC",
                        "Ectomycorrhizal trees", "LAI")
  cov_key_cycl_only <- c("sin", "cos")

  chain_files_all <- list.files(base_path, pattern = "^samples_.*_chain[0-9]+\\.rds$",
                            recursive = TRUE, full.names = TRUE)
  # Deduplicate by basename
  chain_files_all <- chain_files_all[order(nchar(chain_files_all))]
  chain_files_all <- chain_files_all[!duplicated(basename(chain_files_all))]

  beta_list <- list()
  for (cf in chain_files_all) {
    s <- tryCatch(readRDS(cf), error = function(e) NULL)
    if (is.null(s)) next

    # Handle both old and new chain formats
    if (is.list(s) && "samples" %in% names(s)) {
      samp <- s$samples
      meta <- s$metadata
    } else {
      next
    }

    if (is.null(samp) || is.null(meta)) next
    if (!is.matrix(samp)) samp <- as.matrix(samp)

    beta_idx <- grep("^beta\\[", colnames(samp))
    if (length(beta_idx) == 0) next

    mn <- meta$model_name
    cov_key <- if (mn == "env_cycl") cov_key_env_cycl
               else if (mn == "env_cov") cov_key_env_cov
               else cov_key_cycl_only
    if (length(beta_idx) != length(cov_key)) {
      cat("  WARNING: beta count mismatch in", basename(cf),
          "(", length(beta_idx), "vs", length(cov_key), ")\n")
      next
    }

    for (i in seq_along(beta_idx)) {
      vals <- samp[, beta_idx[i]]
      beta_list[[length(beta_list) + 1]] <- data.frame(
        Mean = mean(vals), SD = sd(vals), Median = median(vals),
        rowname = paste0("beta[", i, "]"),
        beta = cov_key[i], taxon = meta$species,
        significant = ifelse(quantile(vals, 0.025) > 0 | quantile(vals, 0.975) < 0, 1, 0),
        effSize = abs(mean(vals)),
        siteID = NA_character_, rank = meta$rank.name,
        model_name = mn,
        group = ifelse(grepl("bac", meta$rank.name), "16S", "ITS"),
        rank_only = sub("_.*", "", meta$rank.name),
        time_period = paste0(meta$min.date, "_", meta$max.date),
        fcast_type = "Taxonomic",
        pretty_group = ifelse(grepl("bac", meta$rank.name), "Bacteria", "Fungi"),
        model_id = meta$model_id,
        chain_file = basename(cf),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(beta_list) > 0) {
    beta_effects <- do.call(rbind, beta_list)
    # Average across chains for the same model
    beta_effects <- beta_effects %>%
      group_by(taxon, model_name, beta, time_period, rank, group, rank_only,
               fcast_type, pretty_group) %>%
      summarize(Mean = mean(Mean), SD = mean(SD), Median = mean(Median),
                significant = max(significant), effSize = mean(effSize),
                .groups = "drop") %>%
      mutate(rowname = "beta", siteID = NA_character_)
  } else {
    beta_effects <- data.frame()
  }
}

cat("Beta effect rows:", nrow(beta_effects), "\n")

saveRDS(beta_effects, here("data/summary/truncated_normal_predictor_effects.rds"))
cat("Saved: data/summary/truncated_normal_predictor_effects.rds\n")
