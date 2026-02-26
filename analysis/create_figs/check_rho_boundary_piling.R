# Check production models for rho boundary piling (posterior pushing against 1.0).
# Iterates over raw MCMC sample files (summary files may omit rho). Flags groups where
# the 95th percentile of rho is > 0.95.

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)
here::i_am("analysis/create_figs/check_rho_boundary_piling.R")

library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(ggplot2)

# Directory containing production model outputs (driver uncertainty = main production)
model_dir <- here("data", "model_outputs")
# Subfolders to scan (omit truncated_normal if you only want beta/driver-uncertainty)
subdirs <- c("cloglog_beta_driver_uncertainty", "cloglog_beta_fixed_drivers")
pattern <- "samples_.*\\.rds$"

sample_files <- character(0)
for (d in subdirs) {
  full_path <- file.path(model_dir, d)
  if (!dir.exists(full_path)) next
  f <- list.files(full_path, pattern = pattern, recursive = TRUE, full.names = TRUE)
  sample_files <- c(sample_files, f)
}

# Optional: restrict to one subdir via env, e.g. RHO_CHECK_DIR=cloglog_beta_driver_uncertainty
check_dir <- Sys.getenv("RHO_CHECK_DIR", "")
if (nzchar(check_dir)) {
  sample_files <- sample_files[grepl(check_dir, sample_files, fixed = TRUE)]
}

# Default: chain1 only (one file per group). Set FULL=1 to process all chains.
if (!nzchar(Sys.getenv("FULL", ""))) {
  sample_files <- sample_files[grepl("_chain1\\.rds$", sample_files)]
  message("Using chain1 only (", length(sample_files), " files). Set FULL=1 for all chains.")
}

message(sprintf("Found %d sample files to process.", length(sample_files)))
if (length(sample_files) == 0) {
  message("No files found. Check model_dir and subdirs, or set RHO_CHECK_DIR to a subfolder name.")
  quit(save = "no", status = 0)
}

# Extract rho stats from a single file. Reads once; uses only the rho column to limit memory.
extract_rho_stats <- function(filepath) {
  tryCatch({
    obj <- readRDS(filepath)
    samples <- if (is.list(obj) && "samples" %in% names(obj)) obj$samples else obj
    rm(obj)
    if (inherits(samples, "mcmc.list")) {
      if (!"rho" %in% colnames(samples[[1L]])) return(NULL)
      rho_draws <- unlist(lapply(samples, function(x) x[, "rho"]), use.names = FALSE)
    } else if (is.matrix(samples) || is.data.frame(samples)) {
      if (!"rho" %in% colnames(samples)) return(NULL)
      rho_draws <- c(samples[, "rho"])
    } else {
      return(NULL)
    }
    rm(samples)
    if (length(rho_draws) == 0L) return(NULL)

    path_parts <- str_split(filepath, .Platform$file.sep)[[1]]
    if (length(path_parts) < 3) return(NULL)
    filename <- path_parts[length(path_parts)]
    species <- path_parts[length(path_parts) - 1]
    model_name <- path_parts[length(path_parts) - 2]
    output_subdir <- path_parts[length(path_parts) - 3]

    tibble(
      output_subdir = output_subdir,
      model = model_name,
      species = species,
      file = filename,
      rho_mean = mean(rho_draws),
      rho_median = median(rho_draws),
      rho_95th = quantile(rho_draws, 0.95),
      rho_max = max(rho_draws),
      p_greater_09 = mean(rho_draws > 0.90),
      n_draws = length(rho_draws)
    )
  }, error = function(e) {
    warning(sprintf("Error reading %s: %s", filepath, e$message))
    return(NULL)
  })
}

message("Extracting rho statistics...")
rho_list <- map(sample_files, extract_rho_stats)
rho_summary_df <- bind_rows(rho_list[!sapply(rho_list, is.null)])

if (nrow(rho_summary_df) == 0) {
  message("No files contained a 'rho' parameter. Exiting.")
  quit(save = "no", status = 0)
}

# One row per file (multiple per group if multiple chains). Add boundary flag per file.
rho_summary_df <- rho_summary_df %>%
  mutate(boundary_flag = if_else(rho_95th > 0.95, "Piling against 1", "Safely below 1"))

# Per-group summary: flag group if any chain piles (max rho_95th across chains > 0.95)
group_summary <- rho_summary_df %>%
  group_by(output_subdir, model, species) %>%
  summarise(
    rho_mean = mean(rho_mean),
    rho_95th_max = max(rho_95th),
    rho_max_over_chains = max(rho_max),
    n_chains = n(),
    boundary_flag = if_else(any(rho_95th > 0.95), "Piling against 1", "Safely below 1"),
    .groups = "drop"
  )

message("\n--- Boundary piling (per file) ---")
print(table(rho_summary_df$boundary_flag))

message("\n--- Boundary piling (per group: any chain > 0.95) ---")
print(table(group_summary$boundary_flag))

out_path <- here("data", "summary", "rho_boundary_diagnostics.rds")
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(
  by_file = rho_summary_df,
  by_group = group_summary
), out_path)
message(sprintf("\nSaved diagnostics to: %s", out_path))

message("\nTop 10 groups by rho 95th percentile:")
group_summary %>%
  arrange(desc(rho_95th_max)) %>%
  select(species, model, rho_mean, rho_95th_max, rho_max_over_chains, boundary_flag) %>%
  head(10) %>%
  print()

message("\nGroups flagged as piling (if any):")
piling <- group_summary %>% filter(boundary_flag == "Piling against 1")
if (nrow(piling) > 0) {
  print(piling %>% select(output_subdir, model, species, rho_95th_max, rho_max_over_chains))
} else {
  message("None.")
}

# --- Single figure: all rho results (95th percentile per group, boundary at 0.95) ---
fig_dir <- here("figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

fig_rho <- ggplot(group_summary, aes(x = rho_95th_max)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 0.95, linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(
    x = "95th percentile of rho (per group)",
    y = "Number of groups",
    title = "Boundary piling: groups with rho 95th > 0.95 flagged"
  ) +
  theme_minimal(base_size = 12)
ggsave(file.path(fig_dir, "rho_boundary_piling.png"), fig_rho, width = 5, height = 4, dpi = 150)
message(sprintf("Saved: %s", file.path(fig_dir, "rho_boundary_piling.png")))
