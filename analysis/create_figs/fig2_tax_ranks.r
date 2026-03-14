source("source.R")

library(ggallin)
library(deeptime)  # devtools::install_github("willgearty/deeptime")
library(tagger)    # devtools::install_github("eliocamp/tagger")
library(data.table)

tryCatch(mem.maxVSize(Inf), error = function(e) invisible(NULL))

# ── Data loading ───────────────────────────────────────────────────────────────
scores_list <- readRDS(here("data", "summary/scoring_metrics_plsr2.rds"))
sum.all     <- readRDS(here("data", "summary/predictor_effects.rds"))

parquet_file <- here("data/summary/parquet/all_hindcasts_plsr2.parquet")
rds_file     <- here("data/summary/all_hindcasts_plsr2.rds")
if (file.exists(parquet_file)) {
  if (requireNamespace("arrow",       quietly = TRUE)) {
    hindcast_in <- arrow::read_parquet(parquet_file)
  } else if (requireNamespace("nanoparquet", quietly = TRUE)) {
    hindcast_in <- nanoparquet::read_parquet(parquet_file)
  } else {
    hindcast_in <- readRDS(rds_file)
  }
} else if (file.exists(rds_file)) {
  hindcast_in <- readRDS(rds_file)
} else {
  stop("Neither Parquet nor RDS hindcast files found!")
}

converged      <- scores_list$converged_strict_list
converged_base <- gsub("_beta_regression$", "", converged)

# ── Taxonomy helpers ───────────────────────────────────────────────────────────
fg_names        <- microbialForecast:::keep_fg_names
rank_spec_names <- microbialForecast:::rank_spec_names
all_bacteria    <- unlist(rank_spec_names[grepl("_bac$", names(rank_spec_names))])
all_fungi       <- unlist(rank_spec_names[grepl("_fun$", names(rank_spec_names))])

create_rank_grouping <- function(data) {
  if ("rank_name" %in% colnames(data)) {
    data$pretty_name <- case_when(
      data$rank_name %in% c("genus_bac",  "genus_fun")  ~ "Genus",
      data$rank_name %in% c("family_bac", "family_fun") ~ "Family",
      data$rank_name %in% c("order_bac",  "order_fun")  ~ "Order",
      data$rank_name %in% c("class_bac",  "class_fun")  ~ "Class",
      data$rank_name %in% c("phylum_bac", "phylum_fun") ~ "Phylum",
      TRUE ~ "Functional group"
    )
    data$pretty_name <- factor(data$pretty_name,
                               levels = c("Genus", "Family", "Order",
                                          "Class", "Phylum", "Functional group"),
                               ordered = TRUE)
  }
  return(data)
}

# ── Scoring data ───────────────────────────────────────────────────────────────
hindcast_filter <- scores_list$scoring_metrics_long %>%
  filter(model_id %in% converged_base,
         metric %in% c("CRPS_truncated", "RSQ", "RSQ.1", "RMSE.norm", "RMSE.iqr"),
         site_prediction == "New time (observed site)") %>%
  create_rank_grouping()

calibration_filter <- scores_list$calibration_metrics_long %>%
  filter(model_id %in% converged_base,
         model_name == "cycl_only",
         metric %in% c("CRPS_truncated", "RSQ", "RSQ.1", "RMSE.norm"))

rmse_values <- hindcast_filter %>% filter(metric == "RMSE.norm")
crps_values <- hindcast_filter %>% filter(metric == "CRPS_truncated")
rsq_values  <- hindcast_filter %>% filter(metric == "RSQ")

# Tukey compact-letter displays — each computed once
tukey_rsq_rank <- rsq_values %>%
  group_by(pretty_group) %>%
  filter(n() >= 2) %>%
  reframe(tukey(x = pretty_name, y = score)) %>%
  rename(pretty_name = x) %>%
  mutate(metric = "RSQ")

tukey_rmse_rank <- rmse_values %>%
  group_by(pretty_group) %>%
  filter(n() >= 2) %>%
  reframe(tukey(x = pretty_name, y = score)) %>%
  rename(pretty_name = x) %>%
  mutate(metric = "RMSE.norm")

tukey_crps_rank <- crps_values %>%
  group_by(pretty_group) %>%
  filter(n() >= 2) %>%
  reframe(tukey(x = pretty_name, y = score)) %>%
  rename(pretty_name = x) %>%
  mutate(metric = "mean_crps")

# Rename metrics for combined plotting
rmse_values$metric <- "RMSE.norm"
crps_values$metric <- "mean_crps"
plotting_df_rank_scores <- bind_rows(rmse_values, crps_values)
plotting_tukey          <- bind_rows(tukey_crps_rank, tukey_rmse_rank)

# ── Shared aesthetics (Okabe-Ito color-blind safe) ─────────────────────────────
# Bacteria/Fungi kingdom colors — consistent with fig3_f_b_seasonality.r
kingdom_colors <- c(Bacteria = "#E69F00", Fungi = "#0072B2")
kingdom_shapes <- c(Bacteria = 16,        Fungi = 17)

# Model-type colors
model_colors <- c(env_cov   = "#56B4E9",
                  cycl_only = "#009E73",
                  env_cycl  = "#D55E00")

# Two-group pairwise comparison list for stat_compare_means
comparisons_kingdom <- list(c("Bacteria", "Fungi"))

base_theme <- theme_bw(base_size = 14) +
  theme(
    strip.background   = element_rect(fill = "grey92", color = NA),
    strip.text         = element_text(face = "bold"),
    axis.text.x        = element_text(angle = 30, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  )

# ── Panel A: nRMSE by taxonomic rank ──────────────────────────────────────────
rmse_rank_plot <- ggplot(rmse_values,
                         aes(x = pretty_name, y = score,
                             fill = pretty_group, color = pretty_group)) +
  geom_violin(alpha = 0.45, quantiles = 0.5, show.legend = FALSE) +
  geom_point(shape = 21, fill = "white", size = 2,
             position = position_jitter(width = 0.15, height = 0),
             alpha = 0.35, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6, color = "grey40") +
  geom_text(data = tukey_rmse_rank,
            aes(x = pretty_name, y = tot, label = Letters_Tukey),
            color = "black", size = 5, inherit.aes = FALSE) +
  facet_wrap(~pretty_group, scales = "free_y", nrow = 1) +
  scale_color_manual(values = kingdom_colors) +
  scale_fill_manual(values  = kingdom_colors) +
  scale_y_continuous(trans = pseudolog10_trans) +
  labs(x = NULL, y = "Relative forecast error (nRMSE)") +
  base_theme +
  theme(plot.margin = margin(0.3, 0.5, 0.3, 0.3, "cm"))

# ── Panel B: nRMSE by kingdom ──────────────────────────────────────────────────
rmse_fb_plot <- ggplot(rmse_values,
                       aes(x = pretty_group, y = score,
                           fill = pretty_group, color = pretty_group)) +
  geom_violin(alpha = 0.45, quantiles = 0.5, show.legend = FALSE) +
  geom_point(shape = 21, fill = "white", size = 2,
             position = position_jitter(width = 0.1, height = 0),
             alpha = 0.35, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6, color = "grey40") +
  stat_compare_means(comparisons = comparisons_kingdom,
                     method = "wilcox.test", label = "p.signif", size = 5) +
  scale_color_manual(values = kingdom_colors) +
  scale_fill_manual(values  = kingdom_colors) +
  scale_y_continuous(trans = pseudolog10_trans) +
  labs(x = NULL, y = "Relative forecast error (nRMSE)") +
  base_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.margin = margin(0.5, 0.3, 0.3, 0.3, "cm"))

# ── Panel D: new-site transferability by kingdom ───────────────────────────────
skill_score_species_data <- scores_list$skill_score_species %>%
  filter(model_id %in% converged_base)

newsite_fb_plot <- ggplot(skill_score_species_data,
                          aes(x = pretty_group, y = skill_score * 100,
                              fill = pretty_group, color = pretty_group)) +
  geom_violin(alpha = 0.45, quantiles = 0.5, show.legend = FALSE) +
  geom_point(shape = 21, fill = "white", size = 2,
             position = position_jitter(width = 0.1, height = 0),
             alpha = 0.35, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.7) +
  stat_compare_means(comparisons = comparisons_kingdom,
                     method = "wilcox.test", label = "p.signif", size = 5) +
  scale_color_manual(values = kingdom_colors) +
  scale_fill_manual(values  = kingdom_colors) +
  coord_cartesian(ylim = c(-100, 100)) +
  labs(x = NULL,
       y = "Transferability to new sites\n(% change in CRPS, clipped ±100%)") +
  base_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm"))

# ── CRPS rank/kingdom plots (supplementary / diagnostic) ──────────────────────
crps_rank_plot <- ggplot(plotting_df_rank_scores %>% filter(metric == "mean_crps"),
                         aes(x = pretty_name, y = score,
                             fill = pretty_group, color = pretty_group)) +
  geom_violin(alpha = 0.45, quantiles = 0.5, show.legend = FALSE) +
  geom_point(shape = 21, fill = "white", size = 2,
             position = position_jitter(width = 0.15, height = 0),
             alpha = 0.35, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6, color = "grey40") +
  geom_text(data = plotting_tukey %>% filter(metric == "mean_crps"),
            aes(x = pretty_name, y = tot, label = Letters_Tukey),
            color = "black", size = 5, inherit.aes = FALSE) +
  facet_wrap(~pretty_group, scales = "free_y", nrow = 1) +
  scale_color_manual(values = kingdom_colors) +
  scale_fill_manual(values  = kingdom_colors) +
  scale_y_continuous(trans = pseudolog10_trans) +
  labs(x = NULL, y = "Absolute forecast error (CRPS)") +
  base_theme +
  theme(plot.margin = margin(0.3, 0.5, 0.3, 0.3, "cm"))

crps_fb_plot <- ggplot(plotting_df_rank_scores %>% filter(metric == "mean_crps"),
                       aes(x = pretty_group, y = score,
                           fill = pretty_group, color = pretty_group)) +
  geom_violin(alpha = 0.45, quantiles = 0.5, show.legend = FALSE) +
  geom_point(shape = 21, fill = "white", size = 2,
             position = position_jitter(width = 0.1, height = 0),
             alpha = 0.35, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6, color = "grey40") +
  stat_compare_means(comparisons = comparisons_kingdom,
                     method = "wilcox.test", label = "p.signif", size = 5) +
  scale_color_manual(values = kingdom_colors) +
  scale_fill_manual(values  = kingdom_colors) +
  scale_y_continuous(trans = pseudolog10_trans) +
  labs(x = NULL, y = "Absolute forecast error (CRPS)") +
  base_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.margin = margin(0.5, 0.3, 0.3, 0.3, "cm"))

# ── Forecast horizon ───────────────────────────────────────────────────────────
fcast_horizon_in   <- readRDS(here("data", "summary/fcast_horizon_df.rds"))
fcast_horizon_data <- fcast_horizon_in[[3]] %>%
  filter(model_id %in% converged_base,
         metric == "rsq",
         horizon_parameter == "rsq_fcast_horizon") %>%
  select(species, model_name, model_id, value, pretty_group) %>%
  rename(forecast_horizon = value) %>%
  filter(!is.na(forecast_horizon))

# Vectorized kingdom classification
fcast_horizon_dt <- as.data.table(fcast_horizon_data)
fcast_horizon_dt[, pretty_group := fcase(
  species %in% all_fungi,    "Fungi",
  species %in% all_bacteria, "Bacteria",
  default = "Bacteria"
)]
fcast_horizon_data <- as.data.frame(fcast_horizon_dt)

if (length(unique(fcast_horizon_data$pretty_group)) < 2)
  cat("WARNING: forecast horizon data is missing one kingdom —",
      "check fcast_horizon_df.rds\n")

# ── Panel C: forecast horizon by kingdom ──────────────────────────────────────
horizon_plot_f <- ggplot(fcast_horizon_data,
                         aes(x = pretty_group, y = forecast_horizon,
                             fill = pretty_group, color = pretty_group)) +
  geom_violin(alpha = 0.45, quantiles = 0.5, show.legend = FALSE) +
  geom_point(shape = 21, fill = "white", size = 2,
             position = position_jitter(width = 0.1, height = 0),
             alpha = 0.45, show.legend = FALSE) +
  stat_compare_means(comparisons = comparisons_kingdom,
                     method = "wilcox.test", label = "p.signif", size = 5) +
  scale_color_manual(values = kingdom_colors) +
  scale_fill_manual(values  = kingdom_colors) +
  labs(x = NULL, y = "Forecast horizon\n(months since last observation)") +
  base_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm"))

# ── Forecast horizon summary stats ────────────────────────────────────────────
cat("\n=== FORECAST HORIZON ANALYSIS BY KINGDOM ===\n")
horizon_summary <- fcast_horizon_data %>%
  group_by(pretty_group) %>%
  summarise(n      = n(),
            mean   = round(mean(forecast_horizon,   na.rm = TRUE), 2),
            median = round(median(forecast_horizon, na.rm = TRUE), 2),
            sd     = round(sd(forecast_horizon,     na.rm = TRUE), 2),
            min    = min(forecast_horizon, na.rm = TRUE),
            max    = max(forecast_horizon, na.rm = TRUE),
            .groups = "drop")
print(horizon_summary)

bact_h <- fcast_horizon_data$forecast_horizon[fcast_horizon_data$pretty_group == "Bacteria"]
fung_h <- fcast_horizon_data$forecast_horizon[fcast_horizon_data$pretty_group == "Fungi"]
if (length(bact_h) > 1 && length(fung_h) > 1) {
  ht <- wilcox.test(bact_h, fung_h)
  cat("Wilcoxon rank-sum (Bacteria vs Fungi): W =", ht$statistic,
      ", p =", round(ht$p.value, 4), "\n")
}
max_h <- sum(fcast_horizon_data$forecast_horizon == 20, na.rm = TRUE)
cat("Models at max horizon (20 months):", max_h, "of", nrow(fcast_horizon_data),
    "(", round(100 * max_h / nrow(fcast_horizon_data), 1), "%)\n")

# ── Figure assembly ────────────────────────────────────────────────────────────
out_dir <- here("figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Tags applied here at assembly time — never baked into the plot objects above.
# This lets the same plot object be reused in different composites without
# conflicting or duplicated tags.
tag_theme <- theme(plot.tag = element_text(size = 14, face = "bold"))

# ── Figure 2 (main): 2 × 2 layout ─────────────────────────────────────────────
# Row 1: error by taxonomic rank — nRMSE (A) and CRPS (B) side-by-side.
#   Each panel is faceted Bacteria | Fungi so kingdom differences are visible
#   within the data, eliminating the need for a separate summary column.
# Row 2: kingdom-level outcome summaries — forecast horizon (C) and
#   new-site transferability (D) — same x-axis orientation as Row 1 facets.
fig2_main <- ggarrange2(
  ggarrange2(
    rmse_rank_plot + labs(tag = "A") + tag_theme,
    crps_rank_plot + labs(tag = "B") + tag_theme,
    nrow = 1),
  ggarrange2(
    horizon_plot_f  + labs(tag = "C") + tag_theme,
    newsite_fb_plot + labs(tag = "D") + tag_theme,
    nrow = 1),
  nrow = 2, heights = c(1.3, 1))

fig2_main
ggsave(file.path(out_dir, "figure_2_main.png"), fig2_main,
       width = 13, height = 11, dpi = 300, units = "in")
cat("Saved: figure_2_main.png\n")

# ── Supplementary S1: kingdom-level nRMSE and CRPS ────────────────────────────
# Shows the explicit Bacteria vs. Fungi Wilcoxon test that is implicit in the
# faceted rank plots above. Two panels, equal width and height.
fig_kingdom_supp <- ggarrange2(
  rmse_fb_plot + labs(tag = "A") + tag_theme,
  crps_fb_plot + labs(tag = "B") + tag_theme,
  nrow = 1)
ggsave(file.path(out_dir, "supp_figure_2_kingdom.png"), fig_kingdom_supp,
       width = 8, height = 5, dpi = 300, units = "in")
cat("Saved: supp_figure_2_kingdom.png\n")

# ── Supplementary S2: all metrics × groups × model types ──────────────────────
metric.labs.full <- c(
  RMSE.norm      = "Relative forecast\nerror (nRMSE)",
  mean_crps      = "Absolute forecast\nerror (CRPS)",
  CRPS_truncated = "Absolute forecast\nerror (CRPS)",
  RSQ            = "R-squared",
  RSQ.1          = "R-squared\nrelative to 1:1 line"
)

rank_models_data <- scores_list$scoring_metrics_long %>%
  filter(model_id %in% converged_base,
         metric %in% c("CRPS_truncated", "RMSE.norm", "RSQ", "RSQ.1"),
         site_prediction == "New time (observed site)") %>%
  create_rank_grouping()

rank_models <- ggplot(rank_models_data,
                      aes(x = pretty_name, y = score, color = model_name)) +
  geom_violin(quantiles = 0.5, show.legend = FALSE) +
  geom_point(size = 2,
             position = position_jitterdodge(jitter.width = 0.2),
             alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6, color = "grey40") +
  facet_grid(metric ~ pretty_group, scales = "free_y",
             labeller = labeller(metric = metric.labs.full)) +
  scale_color_manual(values = model_colors,
                     labels = model.labs,
                     name   = "Model predictors") +
  scale_y_continuous(trans = pseudolog10_trans) +
  labs(x = NULL, y = "Forecast metric") +
  base_theme +
  theme(plot.margin     = margin(1, 2, 1, 1, "cm"),
        legend.position = "right") +
  tag_facets() +
  theme(tagger.panel.tag.text       = element_text(size = 14),
        tagger.panel.tag.background = element_rect(fill = "white", color = "white"))

rank_models
ggsave(file.path(out_dir, "supp_rank_model_scores.png"), rank_models,
       width = 10, height = 20, dpi = 300, units = "in")
cat("Saved: supp_rank_model_scores.png\n")

# ============================================================================
# RESIDUALIZED nRMSE: remove abundance confound
# nRMSE = RMSE / abundance, so abundance mechanically drives nRMSE.
# Fit log(nRMSE) ~ log(abundance) per kingdom, plot residuals by rank.
# ============================================================================

resid_data <- scores_list$scoring_metrics %>%
  filter(model_id %in% converged_base,
         site_prediction == "New time (observed site)",
         !is.na(RMSE), !is.na(RMSE.norm), RMSE > 0, RMSE.norm > 0) %>%
  mutate(mean_abundance = pmax(RMSE / RMSE.norm, 0.001)) %>%
  create_rank_grouping()

# Fit log-log model per kingdom, extract residuals
resid_data <- resid_data %>%
  group_by(pretty_group) %>%
  mutate(log_nrmse  = log(pmin(RMSE.norm, 5)),
         log_abund  = log(mean_abundance),
         nrmse_resid = resid(lm(log_nrmse ~ log_abund))) %>%
  ungroup()

# Tukey tests on residuals
tukey_resid_rank <- resid_data %>%
  group_by(pretty_group) %>%
  filter(n() >= 2) %>%
  reframe(tukey(x = pretty_name, y = nrmse_resid)) %>%
  rename(pretty_name = x)

# Panel A: by-rank violin of residuals
resid_rank_plot <- ggplot(resid_data,
                          aes(x = pretty_name, y = nrmse_resid,
                              fill = pretty_group, color = pretty_group)) +
  geom_violin(alpha = 0.45, quantiles = 0.5, show.legend = FALSE) +
  geom_point(shape = 21, fill = "white", size = 2,
             position = position_jitter(width = 0.15, height = 0),
             alpha = 0.35, show.legend = FALSE) +
  facet_wrap(~pretty_group, scales = "free_y", nrow = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6, color = "grey40") +
  geom_text(data = tukey_resid_rank,
            aes(x = pretty_name, y = tot, label = Letters_Tukey),
            color = "black", size = 5, inherit.aes = FALSE) +
  scale_color_manual(values = kingdom_colors) +
  scale_fill_manual(values  = kingdom_colors) +
  labs(x = NULL,
       y = "Abundance-adjusted forecast error\n(nRMSE residuals, log scale)") +
  base_theme +
  theme(plot.margin = margin(0.3, 0.5, 0.3, 0.3, "cm"))

# Panel B: by-kingdom violin of residuals
resid_fb_plot <- ggplot(resid_data,
                        aes(x = pretty_group, y = nrmse_resid,
                            fill = pretty_group, color = pretty_group)) +
  geom_violin(alpha = 0.45, quantiles = 0.5, show.legend = FALSE) +
  geom_point(shape = 21, fill = "white", size = 2,
             position = position_jitter(width = 0.1, height = 0),
             alpha = 0.35, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6, color = "grey40") +
  stat_compare_means(comparisons = comparisons_kingdom,
                     method = "wilcox.test", label = "p.signif", size = 5) +
  scale_color_manual(values = kingdom_colors) +
  scale_fill_manual(values  = kingdom_colors) +
  labs(x = NULL,
       y = "Abundance-adjusted\nnRMSE residuals") +
  base_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.margin = margin(0.5, 0.3, 0.3, 0.3, "cm"))

# Panel C: scatter — nRMSE vs abundance, colored by rank
resid_scatter <- ggplot(resid_data,
                        aes(x = mean_abundance, y = pmin(RMSE.norm, 5),
                            color = pretty_name)) +
  geom_point(alpha = 0.5, size = 2.5) +
  geom_smooth(aes(group = pretty_group), method = "lm",
              formula = y ~ log(x), se = FALSE, color = "black", linewidth = 0.6) +
  facet_wrap(~pretty_group) +
  scale_x_sqrt() +
  coord_cartesian(ylim = c(0, 5)) +
  labs(x = "Mean observed abundance", y = "nRMSE (capped at 5)",
       color = "Taxonomic rank") +
  base_theme +
  theme(legend.position = "bottom")

fig_resid <- ggarrange2(
  ggarrange2(resid_rank_plot, resid_fb_plot, nrow = 1, widths = c(2, 1)),
  resid_scatter,
  nrow = 2, heights = c(1.2, 1))

ggsave(file.path(out_dir, "figure_2_nrmse_residuals.png"), fig_resid,
       width = 12, height = 14, dpi = 300, units = "in")
cat("Saved: figure_2_nrmse_residuals.png\n")

# Print diagnostic summary
cat("\n-- nRMSE residualization diagnostic --\n")
resid_summary <- resid_data %>%
  group_by(pretty_group, pretty_name) %>%
  summarise(n              = n(),
            mean_abundance = round(mean(mean_abundance, na.rm = TRUE), 4),
            mean_raw_nrmse = round(mean(pmin(RMSE.norm, 5), na.rm = TRUE), 3),
            mean_resid     = round(mean(nrmse_resid, na.rm = TRUE), 3),
            .groups = "drop")
print(as.data.frame(resid_summary))
