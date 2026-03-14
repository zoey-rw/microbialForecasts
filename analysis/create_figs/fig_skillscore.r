# Skill score comparison: modeled vs random site effects across model types
# Central question: across cycl_only / env_cov / env_cycl models, does knowing
# a site's environmental context (PLSR-modeled site effect) improve new-site
# forecasts relative to assigning the average (random) site effect?
#
# Panel A: Raw CRPS by site_prediction type × model type (violin)
# Panel B: Skill score by model type (modeled effect, relative to observed-site baseline)
# Panel C: Skill score — modeled vs random effects head-to-head (env_cycl only)

source("source.R")
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(patchwork)

# ── Data ──────────────────────────────────────────────────────────────────────
scores_list    <- readRDS(here("data/summary/scoring_metrics_plsr2.rds"))
converged_base <- gsub("_beta_regression$", "", scores_list$converged_list)

# Model type labels
model_labels <- c(cycl_only = "Cycl. only", env_cov = "Env. only",
                  env_cycl  = "Env. + Cycl.")

site_pred_order <- c(
  "New time (observed site)",
  "New time x site (modeled effect)",
  "New time x site (random effect)"
)
site_pred_short <- c(
  "New time (observed site)"           = "Observed\nsite",
  "New time x site (modeled effect)"   = "Modeled\nsite effect",
  "New time x site (random effect)"    = "Random\nsite effect"
)

site_scores <- scores_list$scoring_metrics %>%
  filter(model_id %in% converged_base,
         site_prediction %in% site_pred_order,
         !is.na(pretty_group),
         !is.na(mean_crps_sample),
         is.finite(mean_crps_sample)) %>%
  mutate(
    site_prediction = recode(site_prediction, !!!site_pred_short),
    model_label     = recode(model_name, !!!model_labels)
  )

# ── Skill scores ──────────────────────────────────────────────────────────────
skill_df <- site_scores %>%
  select(model_id, pretty_group, model_name, model_label, site_prediction,
         mean_crps_sample) %>%
  group_by(model_id, pretty_group, model_name, model_label, site_prediction) %>%
  summarise(crps = mean(mean_crps_sample, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = site_prediction, values_from = crps) %>%
  filter(!is.na(`Observed\nsite`),
         !is.na(`Modeled\nsite effect`),
         !is.na(`Random\nsite effect`)) %>%
  mutate(
    skill_modeled = 1 - (`Modeled\nsite effect` / `Observed\nsite`),
    skill_random  = 1 - (`Random\nsite effect`  / `Observed\nsite`)
  )

# ── Shared aesthetics ─────────────────────────────────────────────────────────
kingdom_colors <- c(Bacteria = "#E69F00", Fungi = "#0072B2")
model_colors   <- c("Cycl. only" = "#0072B2", "Env. only" = "#D55E00",
                    "Env. + Cycl." = "#009E73")
site_colors    <- c("Observed\nsite"        = "#009E73",
                    "Modeled\nsite effect"  = "#0072B2",
                    "Random\nsite effect"   = "#CC79A7")

base_theme <- theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(size = 9)
  )

# ── Panel A: Raw CRPS by site_prediction × model type ────────────────────────
panel_a <- ggplot(site_scores,
                  aes(x = site_prediction, y = mean_crps_sample,
                      fill = site_prediction)) +
  geom_violin(quantiles = 0.5, alpha = 0.6, show.legend = FALSE) +
  geom_point(aes(color = site_prediction),
             position = position_jitter(width = 0.1, height = 0),
             alpha = 0.15, size = 0.8, show.legend = FALSE) +
  facet_grid(model_label ~ pretty_group, scales = "free_y") +
  scale_fill_manual(values  = site_colors) +
  scale_color_manual(values = site_colors) +
  scale_y_log10() +
  labs(x = NULL, y = "Forecast error (CRPS, log scale)") +
  base_theme +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 8))

# ── Panel B: Skill score by model type (modeled effect) ──────────────────────
# Positive skill = modeled new-site CRPS < observed-site CRPS baseline
tryCatch({
  stat_b <- skill_df %>%
    group_by(pretty_group) %>%
    rstatix::tukey_hsd(skill_modeled ~ model_name) %>%
    rstatix::add_y_position(step.increase = 0.4) %>%
    filter(p.adj < 0.1)
}, error = function(e) {
  stat_b <<- data.frame()
  cat("Skill Tukey failed:", e$message, "\n")
})

panel_b <- ggplot(skill_df,
                  aes(x = model_label, y = skill_modeled,
                      color = model_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_violin(quantiles = 0.5, alpha = 0.5, show.legend = FALSE) +
  geom_point(position = position_jitter(width = 0.12, height = 0),
             alpha = 0.3, size = 1.5, show.legend = FALSE) +
  facet_wrap(~pretty_group) +
  scale_color_manual(values = model_colors) +
  labs(x = NULL,
       y = "Skill score (modeled site effect\nvs. observed-site baseline)") +
  base_theme

# ── Panel C: Modeled vs random skill head-to-head (env_cycl) ─────────────────
skill_long <- skill_df %>%
  filter(model_name == "env_cycl") %>%
  pivot_longer(cols = c(skill_modeled, skill_random),
               names_to  = "method",
               values_to = "skill") %>%
  mutate(method = factor(method,
                         levels = c("skill_random", "skill_modeled"),
                         labels = c("Random\nsite effect", "Modeled site\neffect (PLSR)")))

panel_c <- ggplot(skill_long,
                  aes(x = method, y = skill, color = pretty_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_violin(quantiles = 0.5, aes(fill = pretty_group),
              alpha = 0.4, show.legend = FALSE) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7),
             alpha = 0.35, size = 1.5) +
  stat_compare_means(aes(group = method),
                     method = "wilcox.test", label = "p.signif", size = 4,
                     hide.ns = TRUE) +
  facet_wrap(~pretty_group) +
  scale_color_manual(values = kingdom_colors, guide = "none") +
  scale_fill_manual(values  = kingdom_colors, guide = "none") +
  labs(x = NULL,
       y = "Skill score (% decrease in CRPS\nvs. observed-site baseline)") +
  base_theme

# ── Combine and save ──────────────────────────────────────────────────────────
fig_skill <- panel_a / (panel_b | panel_c) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold", size = 14))
  ) +
  plot_layout(heights = c(1.8, 1))

out_dir <- here("figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "fig_skillscore.pdf"), fig_skill,
       width = 12, height = 14)
ggsave(file.path(out_dir, "fig_skillscore.png"), fig_skill,
       width = 12, height = 14, dpi = 200)
cat("Saved: figures/fig_skillscore.pdf / .png\n")

# ── Summary statistics ────────────────────────────────────────────────────────
cat("\n=== Skill score summary (env_cycl) ===\n")
skill_df %>%
  filter(model_name == "env_cycl") %>%
  group_by(pretty_group) %>%
  summarise(
    n             = n(),
    pct_improved  = round(100 * mean(skill_modeled > skill_random, na.rm = TRUE), 1),
    mean_skill_mod = round(mean(skill_modeled, na.rm = TRUE), 3),
    mean_skill_rnd = round(mean(skill_random,  na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  print()
