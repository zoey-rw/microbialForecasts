# Covariate effect comparison: CLR vs Cloglog vs Dirichlet — env_cycl model
# Ascomycota only (all three models comparable for this taxon)
# Two-panel figure: (A) covariate effects, (B) observed vs estimated
# Colorblind-safe (Okabe-Ito) palette; shape + color double-encoding

source("source.R")
pacman::p_load(ggplot2, dplyr, stringr, patchwork)

# ---- Helpers ----

clean_beta <- function(x) ifelse(grepl("Ecto", as.character(x)), "EM trees", as.character(x))

# Human-readable covariate labels (display order = top to bottom on y-axis)
beta_labels <- c(
  "sin"         = "Seasonal (sin)",
  "cos"         = "Seasonal (cos)",
  "LAI"         = "Leaf area index",
  "EM trees"    = "EM tree cover",
  "pC"          = "Soil carbon",
  "pH"          = "Soil pH",
  "Moisture"    = "Soil moisture",
  "Temperature" = "Temperature"
)
beta_display_order <- rev(beta_labels)   # bottom = Temperature on y-axis

# Okabe-Ito palette — safe for deuteranopia, protanopia, tritanopia
model_colors <- c(
  "Cloglog"        = "#D55E00",   # vermillion
  "CLR"            = "#0072B2",   # deep blue
  "Trunc. normal"  = "#CC79A7",   # pink
  "Dirichlet"      = "#009E73"    # bluish green
)
model_shapes <- c(
  "Cloglog"        = 16,   # filled circle
  "CLR"            = 17,   # filled triangle
  "Trunc. normal"  = 18,   # filled diamond
  "Dirichlet"      = 15    # filled square
)

# ---- Extract betas directly from chain files (correct NIMBLE ordering) ----

# env_cycl NIMBLE ordering (same across all four model types):
# beta[1]=sin, beta[2]=cos, beta[3]=temp, beta[4]=mois,
# beta[5]=pH, beta[6]=pC, beta[7]=relEM, beta[8]=LAI
env_cycl_beta_names <- c("sin", "cos", "Temperature", "Moisture",
                         "pH", "pC", "EM trees", "LAI")

# Extract beta means/SDs from chain files, averaging across chains
extract_betas_from_chains <- function(chain_files, model_label,
                                      beta_names = env_cycl_beta_names,
                                      beta_pattern = "^beta\\[",
                                      species_idx = NULL) {
  chain_stats <- lapply(chain_files, function(f) {
    s <- readRDS(f)
    samp <- s$samples
    bcols <- grep(beta_pattern, colnames(samp), value = TRUE)
    # For Dirichlet: filter to species_idx (e.g. ascomycota = 1)
    if (!is.null(species_idx)) {
      nums <- regmatches(bcols, gregexpr("[0-9]+", bcols))
      spp <- as.integer(sapply(nums, `[`, 1))
      bcols <- bcols[spp == species_idx]
    }
    data.frame(param = bcols,
               mean  = colMeans(samp[, bcols, drop = FALSE]),
               sd    = apply(samp[, bcols, drop = FALSE], 2, sd))
  })
  all_stats <- do.call(rbind, chain_stats)
  multi <- aggregate(cbind(mean, sd) ~ param, data = all_stats, FUN = function(x) x)
  # Pooled: mean of means, sqrt(mean(var_within) + var(means))
  result <- data.frame(
    param = unique(all_stats$param),
    Mean  = tapply(all_stats$mean, all_stats$param, mean),
    SD    = tapply(all_stats$mean, all_stats$param, function(m) {
      sds <- all_stats$sd[all_stats$param == all_stats$param[1]]  # placeholder
      sqrt(mean(tapply(all_stats$sd[all_stats$param == unique(all_stats$param)[1]],
                       seq_along(sds), function(x) x^2)) + var(m))
    })
  )
  # Simpler: recompute properly
  result <- do.call(rbind, lapply(unique(all_stats$param), function(p) {
    rows <- all_stats[all_stats$param == p, ]
    data.frame(param = p,
               Mean = mean(rows$mean),
               SD = sqrt(mean(rows$sd^2) + var(rows$mean)))
  }))
  # Parse covariate index and assign names
  nums <- regmatches(result$param, gregexpr("[0-9]+", result$param))
  if (!is.null(species_idx)) {
    cov_idx <- as.integer(sapply(nums, `[`, 2))  # beta[species, covariate]
  } else {
    cov_idx <- as.integer(sapply(nums, `[`, 1))  # beta[covariate]
  }
  result$beta <- beta_names[cov_idx]
  result %>%
    filter(!is.na(beta)) %>%
    transmute(beta,
              beta_label = beta_labels[beta],
              Mean, lo = Mean - 1.96 * SD, hi = Mean + 1.96 * SD,
              model_type = model_label)
}

# Cloglog: all env_cycl ascomycota chain files
clog_chains <- list.files(here("data/model_outputs/cloglog_beta_driver_uncertainty/env_cycl/ascomycota"),
                          pattern = "^samples_env_cycl.*chain[0-9]\\.rds$", full.names = TRUE)
clog_f <- extract_betas_from_chains(clog_chains, "Cloglog")

# CLR: non-duplicate chain files
clr_chains <- list.files(here("data/model_outputs/CLR_regression/env_cycl/ascomycota"),
                         pattern = "^samples_env_cycl_ascomycota.*chain[0-9]\\.rds$", full.names = TRUE)
clr_chains <- clr_chains[!grepl("_clr_chain", clr_chains)]
clr_f <- extract_betas_from_chains(clr_chains, "CLR")

# Truncated normal: latest chain files (new warm-started + old)
tn_chains <- list.files(here("data/model_outputs/truncated_normal/env_cycl/ascomycota"),
                        pattern = "^samples_env_cycl_ascomycota_20130601.*chain[0-9]\\.rds$", full.names = TRUE)
tn_f <- extract_betas_from_chains(tn_chains, "Trunc. normal")

# Dirichlet: read from pipeline summary (03_summarizeModelOutputs_dirichlet.r)
dir_summary <- readRDS(here("data/summary/dirichlet_regression_summaries.rds"))
dir_betas <- dir_summary$summary_df %>%
  filter(taxon == "ascomycota", !is.na(beta), beta != "UNKNOWN") %>%
  mutate(beta = clean_beta(beta),
         beta_label = beta_labels[beta]) %>%
  filter(!is.na(beta_label)) %>%
  transmute(beta, beta_label,
            Mean, lo = Mean - 1.96 * SD, hi = Mean + 1.96 * SD,
            model_type = "Dirichlet")
dir_f <- dir_betas

forest_df <- bind_rows(clog_f, clr_f, tn_f, dir_f) %>%
  mutate(
    beta_label = factor(beta_label, levels = beta_display_order),
    model_type = factor(model_type, levels = names(model_colors))
  )

# ---- Panel A: Covariate effects ----

fig_a <- ggplot(forest_df, aes(x = Mean, y = beta_label,
                                color = model_type, shape = model_type)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey60", linewidth = 0.4) +
  geom_linerange(aes(xmin = lo, xmax = hi),
                 linewidth = 0.85,
                 position  = position_dodge(width = 0.6)) +
  geom_point(size     = 3.2,
             position = position_dodge(width = 0.6)) +
  scale_color_manual(values = model_colors, name = "Model") +
  scale_shape_manual(values = model_shapes, name = "Model") +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  theme_bw(base_size = 13) +
  theme(
    axis.title.y      = element_blank(),
    axis.text.y       = element_text(size = 12),
    legend.position   = "bottom",
    legend.title      = element_text(size = 11, face = "bold"),
    legend.text       = element_text(size = 11),
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3)
  ) +
  xlab("Posterior mean \u00b1 1.96 SD")

# ---- Panel B: Observed vs Estimated (calibration) ----
# Extract plot-level predicted vs observed from chain files directly

# Helper: extract plot_mu means from samples2, match to truth
extract_plot_ove <- function(chain_files, model_data, model_label, var_prefix = "plot_mu") {
  samples2_list <- lapply(chain_files, function(f) {
    s <- readRDS(f)
    if (!is.null(s$samples2)) s$samples2 else NULL
  })
  samples2_list <- samples2_list[!sapply(samples2_list, is.null)]
  if (length(samples2_list) == 0) return(data.frame())

  s2 <- do.call(rbind, samples2_list)
  mu_cols <- grep(paste0("^", var_prefix), colnames(s2))
  if (length(mu_cols) == 0) return(data.frame())

  mu_means <- colMeans(s2[, mu_cols, drop = FALSE])

  # Parse plot_mu[p, t] indices
  nums <- regmatches(names(mu_means), gregexpr("[0-9]+", names(mu_means)))
  plot_idx <- as.integer(sapply(nums, `[`, 1))
  time_idx <- as.integer(sapply(nums, `[`, 2))

  # Match to observed truth
  truth <- model_data$y
  plot_num <- model_data$plot_num
  timepoint <- model_data$timepoint

  est_df <- data.frame(plot_idx = plot_idx, time_idx = time_idx, estimated = as.numeric(mu_means))

  # Build truth lookup from core observations
  if (is.matrix(truth)) {
    obs_val <- truth[, 1]
  } else {
    obs_val <- truth
  }
  truth_df <- data.frame(plot_idx = plot_num, time_idx = timepoint, observed = obs_val) %>%
    group_by(plot_idx, time_idx) %>%
    summarize(observed = mean(observed, na.rm = TRUE), .groups = "drop")

  merged <- inner_join(est_df, truth_df, by = c("plot_idx", "time_idx")) %>%
    filter(!is.na(estimated), !is.na(observed)) %>%
    mutate(model_type = model_label)

  return(merged[, c("estimated", "observed", "model_type")])
}

# Cloglog: extract from chain samples2 for fair comparison (same method as other models)
clog_chain1 <- list.files(here("data/model_outputs/cloglog_beta_driver_uncertainty/env_cycl/ascomycota"),
                          pattern = "samples_env_cycl.*chain1.rds", full.names = TRUE)[1]
if (!is.na(clog_chain1) && file.exists(clog_chain1)) {
  clog_s1 <- readRDS(clog_chain1)
  clog_ove <- extract_plot_ove(clog_chain1, clog_s1$metadata$model_data, "Cloglog", var_prefix = "plot_mu")
  rm(clog_s1); gc()
} else {
  # Fallback to summary
  plot_est <- readRDS(here("data/summary/plot_estimates.rds"))
  clog_ove <- plot_est %>%
    filter(species == "ascomycota", model_name == "env_cycl") %>%
    group_by(plot_num, timepoint) %>%
    summarize(estimated = mean(Mean, na.rm = TRUE),
              observed = mean(truth, na.rm = TRUE), .groups = "drop") %>%
    mutate(model_type = "Cloglog")
}

# CLR: extract from single chain samples2 (mu = CLR scale) — too large for all chains
clr_chain1 <- here("data/model_outputs/CLR_regression/env_cycl/ascomycota",
                    "samples_env_cycl_ascomycota_20130601_20180101_with_legacy_covariate_clr_chain1.rds")
if (file.exists(clr_chain1)) {
  clr_s1 <- readRDS(clr_chain1)
  clr_ove <- extract_plot_ove(clr_chain1, clr_s1$metadata$model_data, "CLR", var_prefix = "mu")
  rm(clr_s1); gc()
} else {
  clr_ove <- data.frame()
}

# Truncated normal: extract from chain samples2 (plot_mu = proportion scale)
tn_chains <- list.files(here("data/model_outputs/truncated_normal/env_cycl/ascomycota"),
                        pattern = "samples_env_cycl.*20130601.*chain[0-9].rds", full.names = TRUE)
if (length(tn_chains) > 0) {
  tn_s1 <- readRDS(tn_chains[1])
  tn_ove <- extract_plot_ove(tn_chains[1], tn_s1$metadata$model_data, "Trunc. normal", var_prefix = "plot_mu")
  rm(tn_s1); gc()
} else {
  tn_ove <- data.frame()
}

# Dirichlet: use pipeline summary for obs vs estimated (Panel B)
dir_plot_est <- dir_summary$plot_est
if (!is.null(dir_plot_est) && nrow(dir_plot_est) > 0) {
  dir_ove <- dir_plot_est %>%
    filter(taxon == "ascomycota",
           !is.na(Mean), !is.na(truth)) %>%
    transmute(estimated = Mean, observed = truth,
              model_type = "Dirichlet")
} else {
  dir_ove <- data.frame()
}

# Combine all plot x time means (no subsampling)
ove_parts <- list()
if (nrow(clog_ove) > 0) ove_parts <- c(ove_parts, list(clog_ove))
if (nrow(clr_ove) > 0)  ove_parts <- c(ove_parts, list(clr_ove))
if (nrow(tn_ove) > 0)   ove_parts <- c(ove_parts, list(tn_ove))
if (nrow(dir_ove) > 0)  ove_parts <- c(ove_parts, list(dir_ove))

ove_df <- bind_rows(ove_parts) %>%
  filter(!is.na(estimated), !is.na(observed)) %>%
  mutate(model_type = factor(model_type, levels = names(model_colors)))

# Compute R-squared per model
rsq <- ove_df %>%
  group_by(model_type) %>%
  summarize(
    rsq = cor(observed, estimated, use = "complete.obs")^2,
    rmse = sqrt(mean((observed - estimated)^2, na.rm = TRUE)),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("R\u00b2 = %.3f\nRMSE = %.3f", rsq, rmse))

fig_b <- ggplot(ove_df, aes(x = observed, y = estimated, color = model_type)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.15, size = 0.8) +
  geom_text(data = rsq, aes(label = label),
            x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3,
            size = 3.2, show.legend = FALSE) +
  facet_wrap(~model_type, scales = "free") +
  scale_color_manual(values = model_colors, guide = "none") +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey95"),
    strip.text       = element_text(size = 11),
    panel.grid.minor = element_blank()
  ) +
  xlab("Observed") +
  ylab("Estimated")

# ---- Combine panels ----

combined <- fig_a / fig_b +
  plot_layout(heights = c(2, 1.2)) +
  plot_annotation(tag_levels = "A")

# ---- Save ----

outpath <- here("figures/fig_compare_CLR_betareg.png")
ggsave(outpath, combined, width = 8, height = 9, dpi = 300)
cat("Saved:", outpath, "\n")
