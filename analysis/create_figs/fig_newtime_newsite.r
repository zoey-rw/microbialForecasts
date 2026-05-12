# Site prediction types: observed, modeled, and random site effects
# Central question: do PLSR-modeled site effects improve new-site forecasts
# relative to random/mean effects? Panels show RSQ and CRPS across prediction types.
#
# Site prediction types:
#   "New time (observed site)"          — hindcast at a site seen during calibration
#   "New time x site (modeled effect)"  — new site, PLSR-predicted site effect
#   "New time x site (random effect)"   — new site, mean (random) site effect

source("source.R")
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(patchwork)

# ── Data ──────────────────────────────────────────────────────────────────────
scores_list    <- readRDS(here("data/summary/scoring_metrics_plsr2.rds"))
converged_base <- gsub("_beta_regression$", "", scores_list$converged_list)

site_pred_order <- c(
  "New time (observed site)",
  "New time x site (modeled effect)",
  "New time x site (random effect)"
)
site_pred_labels <- c(
  "New time\n(observed site)",
  "New time × site\n(modeled effect)",
  "New time × site\n(random effect)"
)

scores_df <- scores_list$scoring_metrics %>%
  filter(model_id %in% converged_base,
         model_name == "env_cycl",
         site_prediction %in% site_pred_order,
         !is.na(pretty_group)) %>%
  mutate(
    site_prediction = factor(site_prediction,
                             levels = site_pred_order,
                             labels = site_pred_labels)
  )

# ── Shared aesthetics ─────────────────────────────────────────────────────────
# kingdom_colors comes from source.R

site_colors <- c(
  "New time\n(observed site)"       = "#009E73",
  "New time × site\n(modeled effect)" = "#0072B2",
  "New time × site\n(random effect)"  = "#CC79A7"
)

base_theme <- theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(size = 9)
  )

comparisons_mods <- list(
  c("New time × site\n(modeled effect)", "New time × site\n(random effect)")
)

# ── Panel A: RSQ by prediction type ──────────────────────────────────────────
panel_a <- ggplot(scores_df,
                  aes(x = site_prediction, y = RSQ, fill = site_prediction)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.55) +
  geom_point(aes(color = site_prediction),
             position = position_jitter(width = 0.15, height = 0),
             alpha = 0.25, size = 1) +
  stat_compare_means(comparisons = comparisons_mods,
                     method = "wilcox.test", label = "p.signif", size = 4) +
  facet_wrap(~pretty_group) +
  scale_fill_manual(values  = site_colors, guide = "none") +
  scale_color_manual(values = site_colors, guide = "none") +
  labs(x = NULL, y = expression(R^2)) +
  base_theme

# ── Panel B: CRPS by prediction type ─────────────────────────────────────────
panel_b <- ggplot(scores_df %>% filter(!is.na(mean_crps_sample)),
                  aes(x = site_prediction, y = mean_crps_sample,
                      fill = site_prediction)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.55) +
  geom_point(aes(color = site_prediction),
             position = position_jitter(width = 0.15, height = 0),
             alpha = 0.25, size = 1) +
  stat_compare_means(comparisons = comparisons_mods,
                     method = "wilcox.test", label = "p.signif", size = 4) +
  facet_wrap(~pretty_group) +
  scale_fill_manual(values  = site_colors, guide = "none") +
  scale_color_manual(values = site_colors, guide = "none") +
  labs(x = NULL, y = "Probabilistic forecast error (CRPS)") +
  base_theme

# ── Panel C: Paired slope graph — random → modeled per taxon (env_cycl) ──────
# Each line = one taxon; blue = improved (modeled < random CRPS), red = worsened
paired <- scores_df %>%
  filter(site_prediction %in% c("New time × site\n(modeled effect)",
                                "New time × site\n(random effect)"),
         !is.na(mean_crps_sample)) %>%
  pivot_wider(id_cols     = c(model_id, pretty_group),
              names_from  = site_prediction,
              values_from = mean_crps_sample,
              values_fn   = mean) %>%
  rename(modeled = `New time × site\n(modeled effect)`,
         random  = `New time × site\n(random effect)`) %>%
  filter(!is.na(modeled), !is.na(random)) %>%
  mutate(improved = modeled < random)

panel_c <- paired %>%
  pivot_longer(cols = c(random, modeled),
               names_to  = "method",
               values_to = "crps") %>%
  mutate(method = factor(method,
                         levels = c("random", "modeled"),
                         labels = c("Random\nsite effect", "Modeled site\neffect (PLSR)"))) %>%
  ggplot(aes(x = method, y = crps, group = model_id)) +
  geom_line(aes(color = improved), alpha = 0.3, linewidth = 0.5) +
  geom_point(aes(color = improved), alpha = 0.45, size = 1.5) +
  facet_wrap(~pretty_group) +
  scale_color_manual(values = c("TRUE" = "#0072B2", "FALSE" = "#D55E00"),
                     labels = c("TRUE" = "Improved", "FALSE" = "Worsened"),
                     name = "CRPS outcome") +
  scale_y_log10() +
  labs(x = NULL, y = "Mean CRPS (log scale)") +
  base_theme +
  theme(legend.position = "right")

# ── Combine and save ──────────────────────────────────────────────────────────
# Add percentage improved annotation
pct_improved <- paired %>%
  group_by(pretty_group) %>%
  summarise(pct = round(100 * mean(improved, na.rm = TRUE), 1), .groups = "drop")
cat("Percent of taxa improved by modeled site effects:\n")
print(pct_improved)

fig_newsite <- panel_a / panel_b / panel_c +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold", size = 14))
  )

out_dir <- here("figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "fig_newtime_newsite.png"), fig_newsite,
       width = 10, height = 13, dpi = 200)
cat("Saved: figures/fig_newtime_newsite.png\n")
