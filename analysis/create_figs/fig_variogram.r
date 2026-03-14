# Spatial autocorrelation of site effects: variogram p-value analysis
# Shows that PLSR-modeled site effects capture most spatial autocorrelation.
# Panel A: % of taxa with significant spatial autocorrelation (raw vs PLSR residuals)
# Panel B: ECDF of variogram p-values vs uniform reference

source("source.R")
library(ggpubr)
library(spgs)

# ── Data loading ─────────────────────────────────────────────────────────────
variograms <- readRDS(here("data/summary/site_effect_variograms.rds"))
scores_list <- readRDS(here("data/summary/scoring_metrics_plsr2.rds"))
converged_strict <- scores_list$converged_strict_list

sig_results <- variograms[[1]] %>%
  mutate(
    model_name = ifelse(
      is.na(model_name),
      case_when(
        grepl("^env_cycl", model_id) ~ "env_cycl",
        grepl("^env_cov", model_id) ~ "env_cov",
        grepl("^cycl_only", model_id) ~ "cycl_only",
        TRUE ~ NA_character_
      ),
      model_name
    )
  )

sig_results_long <- sig_results %>%
  pivot_longer(cols = 5:6) %>%
  mutate(value = as.numeric(value)) %>%
  filter(model_id %in% converged_strict)

cat("Using", length(unique(sig_results_long$model_id)), "converged models\n")

# ── Uniformity tests ────────────────────────────────────────────────────────
# Test whether variogram p-values follow a uniform distribution
# (uniform = no spatial autocorrelation pattern)
for (stage in c("site effect", "site effect residuals")) {
  vals <- sig_results_long %>% filter(name == stage) %>% pull(value)
  if (length(vals) >= 2) {
    cat("\nUniformity test for", stage, ":\n")
    print(chisq.unif.test(vals))
  }
}

for (model_type in c("cycl_only", "env_cov", "env_cycl")) {
  vals <- sig_results_long %>%
    filter(name == "site effect residuals", model_name == model_type) %>%
    pull(value)
  if (length(vals) >= 2) {
    cat("\nResidual uniformity test for", model_type, ":\n")
    print(chisq.unif.test(vals))
  }
}

# ── Example variogram plot ──────────────────────────────────────────────────
if (length(variograms) >= 2 && length(variograms[[2]]) >= 12) {
  example_plot <- variograms[[2]][[12]]
  if (inherits(example_plot, "recordedplot")) {
    png(here("figures", "variogram_example.png"), width = 600, height = 800)
    replayPlot(example_plot)
    dev.off()
    cat("Generated variogram_example.png\n")
  }
}

# ── Display settings ────────────────────────────────────────────────────────
model_labs <- c(cycl_only = "Seasonal only",
                env_cov   = "Environment only",
                env_cycl  = "Environment + seasonal")
stage_labs <- c("site effect"           = "Raw site effects",
                "site effect residuals" = "After PLSR modeling")
stage_colors <- c("Raw site effects"     = "#D55E00",
                  "After PLSR modeling"  = "#0072B2")

sig_results_long <- sig_results_long %>%
  mutate(
    model_name  = factor(model_name, levels = c("cycl_only", "env_cov", "env_cycl")),
    stage_label = factor(stage_labs[name], levels = stage_labs),
    significant = value < 0.05
  )

# ── Panel A: % of taxa with significant spatial autocorrelation ─────────────
pct_sig <- sig_results_long %>%
  group_by(model_name, stage_label) %>%
  summarise(
    n_total = n(),
    n_sig   = sum(significant, na.rm = TRUE),
    pct_sig = 100 * n_sig / n_total,
    .groups = "drop"
  )

pA <- ggplot(pct_sig,
             aes(x = model_name, y = pct_sig, fill = stage_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_text(aes(label = paste0(round(pct_sig, 0), "%")),
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 4) +
  scale_x_discrete(labels = model_labs) +
  scale_fill_manual(values = stage_colors, name = NULL) +
  coord_cartesian(ylim = c(0, max(pct_sig$pct_sig) * 1.15)) +
  labs(x = NULL,
       y = "Taxa with significant\nspatial autocorrelation (%)",
       tag = "A") +
  theme_bw(base_size = 14) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        legend.position    = "top",
        axis.text.x        = element_text(angle = 20, hjust = 1))

# ── Panel B: ECDF of p-values vs uniform reference ─────────────────────────
pB <- ggplot(sig_results_long,
             aes(x = value, color = stage_label)) +
  stat_ecdf(linewidth = 0.8, pad = FALSE) +
  geom_abline(slope = 1, intercept = 0, color = "grey50",
              linetype = "dotted", linewidth = 0.6) +
  facet_wrap(~model_name, labeller = labeller(model_name = model_labs)) +
  scale_color_manual(values = stage_colors, name = NULL) +
  annotate("text", x = 0.65, y = 0.15, label = "uniform\nreference",
           color = "grey50", size = 3.5, fontface = "italic") +
  labs(x = "Variogram p-value",
       y = "Cumulative proportion of taxa",
       tag = "B") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_rect(fill = "grey92", color = NA),
        strip.text       = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position  = "top")

# ── Combine and save ────────────────────────────────────────────────────────
fig <- ggarrange(pA, pB, nrow = 1, widths = c(1, 1.5),
                 common.legend = TRUE, legend = "top")

ggsave(here("figures", "fig_variogram_summary.png"), fig,
       width = 12, height = 5.5, dpi = 300)
ggsave(here("figures", "fig_variogram_summary.pdf"), fig,
       width = 12, height = 5.5, dpi = 300)

cat("Saved: figures/fig_variogram_summary.png\n")
cat("Saved: figures/fig_variogram_summary.pdf\n")

cat("\n-- Spatial autocorrelation summary --\n")
print(as.data.frame(pct_sig))
