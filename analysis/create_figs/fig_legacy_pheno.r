# Microbial phenology strength figure
# Shows the prevalence and magnitude of seasonal patterns in microbial communities.
# Compares cycl_only (pure seasonality model) vs env_cycl (seasonality + environment)
# to show how much of the cyclic signal persists after accounting for env predictors.

library(tidyverse)
library(ggpubr)
library(ggrepel)

source("source.R")

# ── Data loading ─────────────────────────────────────────────────────────────
phenophase_in <- readRDS(here("data/clean/pheno_group_peak_phenophases.rds"))

# Element 1: one row per model_id — modal phenophase, amplitude, significance flags
mode_data <- phenophase_in[[1]]

# ── Display settings ─────────────────────────────────────────────────────────
pheno_levels <- c("greenup", "peak", "greendown", "dormancy")
pheno_labels <- c("Green-up", "Peak", "Senescence", "Dormancy")
pheno_colors <- c(
  "Green-up"   = "#009E73",
  "Peak"       = "#E69F00",
  "Senescence" = "#D55E00",
  "Dormancy"   = "#56B4E9"
)

model_labels <- c(
  cycl_only = "Seasonality only",
  env_cycl  = "Seasonality + environment"
)
fcast_labels <- c(Functional = "Functional groups", Taxonomic = "Taxonomic groups")

# ── Panel A: proportion of taxa with significant seasonality ──────────────────
# Compare cycl_only vs env_cycl: does adding env covariates change apparent seasonality?
pct_sig <- mode_data %>%
  filter(model_name %in% c("cycl_only", "env_cycl")) %>%
  group_by(model_name, fcast_type, pretty_group) %>%
  summarise(
    total   = n(),
    n_sig   = sum(significant_sin == 1 | significant_cos == 1, na.rm = TRUE),
    pct_sig = n_sig / total * 100,
    .groups = "drop"
  ) %>%
  mutate(
    model_name = recode(model_name, !!!model_labels),
    fcast_type = recode(fcast_type, !!!fcast_labels)
  )

pA <- ggplot(pct_sig,
             aes(x = pretty_group, y = pct_sig, fill = model_name)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_text(aes(label = paste0("n=", n_sig, "/", total)),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 3, color = "grey30") +
  facet_wrap(~fcast_type) +
  scale_fill_manual(
    values = c("Seasonality only" = "#0072B2",
               "Seasonality + environment" = "#009E73"),
    name   = "Model"
  ) +
  scale_y_continuous(limits = c(0, 115),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "% of taxa with significant seasonality") +
  theme_bw(base_size = 12) +
  theme(
    strip.background   = element_rect(fill = "grey92", color = NA),
    strip.text         = element_text(face = "bold"),
    legend.position    = "top",
    panel.grid.major.x = element_blank()
  )

# ── Panel B: amplitude distribution for significantly seasonal taxa ────────────
# Compare amplitude between cycl_only and env_cycl for the same taxa
# (taxa that are significant in either model)
amp_both <- mode_data %>%
  filter(model_name %in% c("cycl_only", "env_cycl"),
         significant_sin == 1 | significant_cos == 1,
         !is.na(amplitude)) %>%
  mutate(
    model_name = recode(model_name, !!!model_labels),
    fcast_type = recode(fcast_type, !!!fcast_labels)
  )

pB <- ggplot(amp_both,
             aes(x = pretty_group, y = amplitude, fill = model_name)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6,
               position = position_dodge(width = 0.75), width = 0.65) +
  geom_point(aes(color = model_name),
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
             alpha = 0.4, size = 1.5) +
  facet_wrap(~fcast_type) +
  scale_fill_manual(
    values = c("Seasonality only" = "#0072B2",
               "Seasonality + environment" = "#009E73"),
    name   = "Model"
  ) +
  scale_color_manual(
    values = c("Seasonality only" = "#0072B2",
               "Seasonality + environment" = "#009E73"),
    guide  = "none"
  ) +
  labs(x = NULL, y = "Seasonal amplitude") +
  theme_bw(base_size = 12) +
  theme(
    strip.background   = element_rect(fill = "grey92", color = NA),
    strip.text         = element_text(face = "bold"),
    legend.position    = "none",   # shared via pA
    panel.grid.major.x = element_blank()
  )

# ── Panel C: ranked functional group amplitudes (cycl_only) ──────────────────
# Lollipop chart of all functional groups, colored by peak phenophase
# Shows WHICH functional groups are most seasonal and WHEN they peak
fg_ranked <- mode_data %>%
  filter(
    model_name  == "cycl_only",
    fcast_type  == "Functional",
    !is.na(amplitude)
  ) %>%
  arrange(desc(amplitude)) %>%
  mutate(
    label      = tools::toTitleCase(gsub("_", " ", taxon)),
    sig        = significant_sin == 1 | significant_cos == 1,
    peak_phase = factor(sampling_season, levels = pheno_levels, labels = pheno_labels)
  )

pC <- ggplot(fg_ranked,
             aes(x = amplitude,
                 y = reorder(label, amplitude),
                 color = peak_phase,
                 alpha = sig)) +
  geom_segment(aes(xend = 0, yend = reorder(label, amplitude)),
               color = "grey75", linewidth = 0.5) +
  geom_point(aes(shape = pretty_group), size = 3.5) +
  scale_color_manual(values = pheno_colors, name = "Peak phenophase",
                     na.value = "grey60") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35),
                     guide  = "none") +
  scale_shape_manual(values = c(Bacteria = 16, Fungi = 17), name = "Kingdom") +
  labs(x = "Seasonal amplitude (cycl_only)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    legend.position    = "right"
  )

# ── Combine and save ──────────────────────────────────────────────────────────
fig_pheno <- ggarrange(
  ggarrange(pA, pB, nrow = 2, labels = c("A", "B"),
            common.legend = TRUE, legend = "top"),
  pC,
  ncol    = 2,
  labels  = c("", "C"),
  widths  = c(1, 0.9)
)

out_dir <- here("figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "fig_legacy_pheno.png"),
       fig_pheno, width = 12, height = 8, dpi = 200)

cat("Saved: figures/fig_legacy_pheno.png\n")
