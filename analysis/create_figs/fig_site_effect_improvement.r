# Site effect improvement: CRPS improvement from modeled vs random site effects
# For each taxon, computes (random CRPS - modeled CRPS) / random CRPS.
# Positive = modeled site effect is more accurate than random.
#
# Panel A: Distribution of improvement (%) across taxa, faceted by model type
# Panel B: Overall mean ± SD improvement by model type and kingdom
# Panel C: Prediction vs truth scatter (env_cycl, unobserved sites)

source("source.R")
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ── Data ──────────────────────────────────────────────────────────────────────
scores_list    <- readRDS(here("data/summary/scoring_metrics_plsr2.rds"))
converged_base <- gsub("_beta_regression$", "", scores_list$converged_list)

model_labels <- c(cycl_only = "Cycl. only", env_cov = "Env. only",
                  env_cycl  = "Env. + Cycl.")

scores_wide <- scores_list$scoring_metrics %>%
  filter(model_id %in% converged_base,
         site_prediction %in% c("New time x site (modeled effect)",
                                "New time x site (random effect)"),
         !is.na(mean_crps_sample),
         is.finite(mean_crps_sample),
         !is.na(pretty_group)) %>%
  group_by(model_id, pretty_group, model_name, site_prediction) %>%
  summarise(crps = mean(mean_crps_sample, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = site_prediction, values_from = crps) %>%
  rename(modeled = `New time x site (modeled effect)`,
         random  = `New time x site (random effect)`) %>%
  filter(!is.na(modeled), !is.na(random)) %>%
  mutate(
    improvement_pct = (random - modeled) / random * 100,
    improved        = improvement_pct > 0,
    model_label     = recode(model_name, !!!model_labels)
  )

cat("Taxa included:", nrow(scores_wide), "\n")

# ── Shared aesthetics ─────────────────────────────────────────────────────────
kingdom_colors   <- c(Bacteria = "#E69F00", Fungi = "#0072B2")
improve_colors   <- c("TRUE" = "#0072B2", "FALSE" = "#D55E00")
model_colors     <- c("Cycl. only" = "#0072B2", "Env. only" = "#D55E00",
                      "Env. + Cycl." = "#009E73")

base_theme <- theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# ── Panel A: Histogram of CRPS improvement % by model type × kingdom ─────────
panel_a <- ggplot(scores_wide,
                  aes(x = improvement_pct, fill = improved)) +
  geom_histogram(bins = 25, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8,
             color = "grey20") +
  facet_grid(pretty_group ~ model_label, scales = "free_y") +
  scale_fill_manual(values = improve_colors,
                    labels = c("TRUE" = "Improved", "FALSE" = "Worsened"),
                    name   = NULL) +
  labs(x = "CRPS improvement (%)   [positive = modeled site effect is better]",
       y = "Number of taxa") +
  base_theme +
  theme(legend.position = "top")

# ── Panel B: Mean ± SD improvement by model type (per kingdom) ───────────────
improvement_summary <- scores_wide %>%
  group_by(pretty_group, model_name, model_label) %>%
  summarise(
    mean_pct = mean(improvement_pct, na.rm = TRUE),
    sd_pct   = sd(improvement_pct,   na.rm = TRUE),
    pct_improved = 100 * mean(improved, na.rm = TRUE),
    n        = n(),
    .groups  = "drop"
  )

panel_b <- ggplot(improvement_summary,
                  aes(x = model_label, y = mean_pct, color = pretty_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_pointrange(aes(ymin = mean_pct - sd_pct,
                      ymax = mean_pct + sd_pct),
                  position = position_dodge(width = 0.4),
                  linewidth = 0.8, size = 0.7) +
  scale_color_manual(values = kingdom_colors, name = "Kingdom") +
  labs(x = "Model type",
       y = "Mean CRPS improvement (%)") +
  base_theme +
  theme(legend.position = "right")

# ── Panel C: Prediction vs truth scatter (env_cycl, unobserved sites) ─────────
hindcasts <- readRDS(here("data/summary/all_hindcasts_plsr2.rds"))

scatter_df <- hindcasts %>%
  filter(model_name == "env_cycl",
         site_prediction %in% c("New time x site (modeled effect)",
                                "New time x site (random effect)"),
         !is.na(truth), !is.na(mean),
         is.finite(truth), is.finite(mean)) %>%
  mutate(
    site_prediction = recode(site_prediction,
      "New time x site (modeled effect)" = "Modeled site effect",
      "New time x site (random effect)"  = "Random site effect"
    )
  ) %>%
  sample_frac(0.3)   # downsample for speed; remove if you want all points

panel_c <- ggplot(scatter_df,
                  aes(x = truth, y = mean, color = pretty_group)) +
  geom_point(alpha = 0.15, size = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey30", linewidth = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  facet_grid(pretty_group ~ site_prediction) +
  scale_color_manual(values = kingdom_colors, guide = "none") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Observed relative abundance",
       y = "Predicted relative abundance") +
  base_theme

# ── Combine and save ──────────────────────────────────────────────────────────
fig_improvement <- (panel_a) / (panel_b | panel_c) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold", size = 14))
  ) +
  plot_layout(heights = c(1.2, 1))

out_dir <- here("data", "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "fig_site_effect_improvement.png"), fig_improvement,
       width = 13, height = 13, dpi = 200)
cat("Saved: data/figures/fig_site_effect_improvement.pdf / .png\n")

# ── Summary stats ─────────────────────────────────────────────────────────────
cat("\n=== Improvement by model type ===\n")
print(as.data.frame(improvement_summary), digits = 3)
