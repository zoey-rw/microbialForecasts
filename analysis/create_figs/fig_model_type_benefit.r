# Does adding environmental drivers improve microbial forecasts, or is the
# seasonal cycle alone sufficient? This figure compares the three predictor
# sets (seasonality only, environment only, environment + seasonality) head to
# head per taxon, rather than as overlapping violins.
#
# Interpretation (see manuscript Figs 4, 5, and the site-effect PLSR figure):
# the hindcasts supply the env models with the OBSERVED forecast-period drivers
# (only measurement uncertainty is propagated; run_hindcast.r L764), so env
# models are not penalised for unknown future drivers. They still do not beat
# the seasonality-only model because (i) the dynamic drivers are collinear with
# season (temperature-season r = 0.94, Fig 5), and (ii) the static drivers
# (pH, C) are absorbed by the site random effects. The harmonic already encodes
# the forecastable environmental signal; extra covariates add estimation
# variance. This is a parsimony result, NOT "environment is irrelevant".
source("source.R")

library(tidyr)

scores_list <- readRDS(here("data", "summary/scoring_metrics_plsr2.rds"))
converged   <- gsub("_beta_regression$", "", scores_list$converged_strict_list)

# One row per (taxon, model type) on the held-out forecast period.
sm <- scores_list$scoring_metrics %>%
  filter(model_id %in% converged,
         site_prediction == "New time (observed site)") %>%
  mutate(pretty_name = case_when(
    grepl("genus",  rank_name) ~ "Genus",
    grepl("family", rank_name) ~ "Family",
    grepl("order",  rank_name) ~ "Order",
    grepl("class",  rank_name) ~ "Class",
    grepl("phylum", rank_name) ~ "Phylum",
    TRUE                       ~ "Functional group"))

base_theme <- theme_bw(base_size = 14) +
  theme(strip.background = element_rect(fill = "grey92", color = NA),
        strip.text       = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# ── Paired deltas: each taxon contributes one cycl_only and one env_cycl score ──
wide <- sm %>%
  select(species, pretty_group, pretty_name, model_name, RSQ, CRPS_truncated) %>%
  pivot_wider(names_from = model_name,
              values_from = c(RSQ, CRPS_truncated))

# ΔR²  = env_cycl − cycl_only  (positive = environment helps)
# ΔCRPS% = relative change in CRPS (negative = environment helps, CRPS lower=better)
delta <- wide %>%
  filter(!is.na(RSQ_cycl_only), !is.na(RSQ_env_cycl),
         !is.na(CRPS_truncated_cycl_only), !is.na(CRPS_truncated_env_cycl)) %>%
  mutate(d_rsq  = RSQ_env_cycl - RSQ_cycl_only,
         d_crps = 100 * (CRPS_truncated_env_cycl - CRPS_truncated_cycl_only) /
                          CRPS_truncated_cycl_only)

# Per-kingdom annotation: % of taxa improved + paired Wilcoxon p.
annot <- delta %>%
  group_by(pretty_group) %>%
  summarise(
    n            = n(),
    pct_rsq_up   = round(100 * mean(d_rsq > 0)),
    p_rsq        = wilcox.test(RSQ_env_cycl, RSQ_cycl_only, paired = TRUE)$p.value,
    pct_crps_dn  = round(100 * mean(d_crps < 0)),
    p_crps       = wilcox.test(CRPS_truncated_env_cycl, CRPS_truncated_cycl_only,
                               paired = TRUE)$p.value,
    .groups = "drop")

psig <- function(p) ifelse(p < .001, "p<0.001", paste0("p=", signif(p, 2)))

# ── Panel A: ΔR² ───────────────────────────────────────────────────────────────
panelA <- ggplot(delta, aes(pretty_group, d_rsq, fill = pretty_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_violin(alpha = 0.4, color = NA, show.legend = FALSE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.8, show.legend = FALSE) +
  geom_jitter(width = 0.08, height = 0, alpha = 0.3, size = 1.4, show.legend = FALSE) +
  geom_text(data = annot, aes(x = pretty_group, y = Inf,
            label = paste0(pct_rsq_up, "% improved\n", psig(p_rsq))),
            vjust = 1.3, size = 3.6, inherit.aes = FALSE) +
  scale_fill_manual(values = kingdom_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  labs(x = NULL, y = expression(Delta*R^2*"  (env + seasonality  −  seasonality)"),
       subtitle = "Above 0 = environment improves variance explained") +
  base_theme

# ── Panel B: ΔCRPS% ─────────────────────────────────────────────────────────────
panelB <- ggplot(delta, aes(pretty_group, d_crps, fill = pretty_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_violin(alpha = 0.4, color = NA, show.legend = FALSE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.8, show.legend = FALSE) +
  geom_jitter(width = 0.08, height = 0, alpha = 0.3, size = 1.4, show.legend = FALSE) +
  geom_text(data = annot, aes(x = pretty_group, y = Inf,
            label = paste0(pct_crps_dn, "% improved\n", psig(p_crps))),
            vjust = 1.3, size = 3.6, inherit.aes = FALSE) +
  scale_fill_manual(values = kingdom_colors) +
  coord_cartesian(ylim = c(-60, 80)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  labs(x = NULL, y = "ΔCRPS (%)  (env + seasonality vs seasonality)",
       subtitle = "Below 0 = environment lowers forecast error") +
  base_theme

# ── Panel C: win rate — best of three models per taxon (lowest CRPS) ────────────
# Fairness fix: a "best of three" comparison is only valid for taxa where all
# three models converged. Otherwise the seasonality-only model wins by default
# whenever the (harder-to-converge) environmental models simply failed to
# converge for that taxon — inflating its apparent advantage. The convergence
# asymmetry is real evidence of parsimony but is reported separately, below,
# rather than smuggled into the win rate.
species_all3 <- sm %>%
  filter(!is.na(CRPS_truncated)) %>%
  distinct(species, model_name) %>%
  count(species) %>%
  filter(n == 3) %>%
  pull(species)

win <- sm %>%
  filter(!is.na(CRPS_truncated), species %in% species_all3) %>%
  group_by(species, pretty_group) %>%
  slice_min(CRPS_truncated, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  count(pretty_group, model_name) %>%
  group_by(pretty_group) %>%
  mutate(pct = 100 * n / sum(n), n_taxa = sum(n)) %>%
  ungroup() %>%
  mutate(model_name = factor(model_name,
                             levels = c("cycl_only", "env_cov", "env_cycl")))

# Convergence counts per model type (reported, not used as wins).
conv_counts <- sm %>%
  filter(!is.na(CRPS_truncated)) %>%
  distinct(species, pretty_group, model_name) %>%
  count(pretty_group, model_name)

panelC <- ggplot(win, aes(pretty_group, pct, fill = model_name)) +
  geom_col(width = 0.65, color = "white") +
  geom_text(aes(label = ifelse(pct >= 6, paste0(round(pct), "%"), "")),
            position = position_stack(vjust = 0.5), size = 3.8, color = "white",
            fontface = "bold") +
  scale_fill_manual(values = model_colors, labels = model.labs,
                    name = "Best model\n(lowest CRPS)") +
  labs(x = NULL, y = "Taxa best forecast by each model (%)",
       subtitle = paste0(
         "Among taxa with all three models converged, seasonality-only ",
         "wins most often\n(", win$n_taxa[win$pretty_group == "Bacteria"][1],
         " bacterial, ", win$n_taxa[win$pretty_group == "Fungi"][1],
         " fungal taxa)")) +
  base_theme

out_dir <- here("figures")
fig <- ggpubr::ggarrange(
  ggpubr::ggarrange(panelA, panelB, ncol = 2, labels = c("a", "b")),
  panelC, nrow = 2, labels = c("", "c"), heights = c(1, 0.95))

ggsave(file.path(out_dir, "fig_model_type_benefit.png"), fig,
       width = 11, height = 10, dpi = 300, units = "in")
cat("Saved: figures/fig_model_type_benefit.png\n")

# Console summary for the manuscript.
cat("\n=== Model-type comparison (strict-converged, observed-site forecast) ===\n")
print(as.data.frame(annot))
cat("\nWin rate (best of 3 by CRPS, taxa with all three converged):\n")
print(as.data.frame(win))
cat("\nConvergence counts per model type (parsimony evidence, reported separately):\n")
print(as.data.frame(conv_counts))
