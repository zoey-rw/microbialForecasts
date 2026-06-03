# Peak phenophase figure
# Shows which plant phenological season microbial taxa peak in
# Uses cycl_only model (seasonality-only, no env predictors) for clearest signal
# Filters to significantly seasonal taxa (significant sin or cos component)

library(tidyverse)
library(ggpubr)
library(ggrepel)

source("source.R")

# ── Data loading ─────────────────────────────────────────────────────────────
phenophase_in <- readRDS(here("data/clean/pheno_group_peak_phenophases.rds"))

# Per-model x phenophase mean modeled abundance, derived from the full monthly
# view (element 6). This is the foundation for all three panels: Panel B uses
# the full shape; Panels A and C derive each model's "peak phenophase" as the
# phenophase with the highest mean modeled abundance (mean-monthly approach),
# rather than the modal annual-argmax (element 1). The mean-monthly definition
# accounts for phenophase duration and matches the seasonal-niche
# interpretation in the main text.
seasonality_mode_all <- phenophase_in[[6]] %>%
  group_by(model_id, sampling_season, fcast_type, pretty_group, rank_only,
           model_name, taxon, amplitude, significant_sin, significant_cos) %>%
  summarise(mean_abun = mean(mean_modeled_abun, na.rm = TRUE),
            n = n(),
            .groups = "drop")

# One row per model: the phenophase with highest mean modeled abundance.
seasonality_mode_max <- seasonality_mode_all %>%
  group_by(model_id, fcast_type, pretty_group, rank_only, model_name, taxon,
           amplitude, significant_sin, significant_cos) %>%
  slice_max(mean_abun, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(model_id, fcast_type, pretty_group, rank_only, model_name, taxon,
         amplitude, significant_sin, significant_cos,
         sampling_season, mean_abun)

# ── Display settings ─────────────────────────────────────────────────────────
pheno_levels <- c("greenup", "peak", "greendown", "dormancy")
pheno_labels <- c("Green-up", "Peak", "Senescence", "Dormancy")
pheno_colors <- c(
  "Green-up"   = "#009E73",
  "Peak"       = "#E69F00",
  "Senescence" = "#D55E00",
  "Dormancy"   = "#56B4E9"
)
fcast_labels <- c(Functional = "Functional groups", Taxonomic = "Taxonomic groups")

# ── Filtering ─────────────────────────────────────────────────────────────────
sig_max <- seasonality_mode_max %>%
  filter(model_name == "cycl_only",
         significant_sin == 1 | significant_cos == 1) %>%
  mutate(
    sampling_season = factor(sampling_season, levels = pheno_levels, labels = pheno_labels),
    fcast_type      = recode(fcast_type, !!!fcast_labels)
  )

sig_all <- seasonality_mode_all %>%
  filter(model_name == "cycl_only",
         significant_sin == 1 | significant_cos == 1) %>%
  mutate(
    sampling_season = factor(sampling_season, levels = pheno_levels, labels = pheno_labels),
    fcast_type      = recode(fcast_type, !!!fcast_labels)
  )

# ── Panel A: proportion of taxa whose mean-monthly peak falls in each phase ──
prop_data <- sig_max %>%
  count(fcast_type, pretty_group, sampling_season) %>%
  group_by(fcast_type, pretty_group) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

pA <- ggplot(prop_data,
             aes(x = pretty_group, y = prop, fill = sampling_season)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  facet_wrap(~fcast_type) +
  scale_fill_manual(values = pheno_colors, name = "Peak phenophase") +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "Proportion of taxa") +
  theme_bw(base_size = 12) +
  theme(
    legend.position    = "right",
    strip.background   = element_rect(fill = "grey92", color = NA),
    strip.text         = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )

# ── Panel B: per-taxon normalized seasonal abundance profiles ─────────────────
# Each line = one taxon; y = mean relative abundance normalized 0–1 within taxon
# Shows the shape of the seasonal signal for each significantly seasonal taxon
profiles <- sig_all %>%
  group_by(model_id) %>%
  mutate(abun_norm = (mean_abun - min(mean_abun, na.rm = TRUE)) /
                     (max(mean_abun, na.rm = TRUE) - min(mean_abun, na.rm = TRUE))) %>%
  ungroup() %>%
  filter(is.finite(abun_norm))

# Summarise to median + IQR per group for a ribbon overlay
profile_ribbon <- profiles %>%
  group_by(fcast_type, pretty_group, sampling_season) %>%
  summarise(
    med  = median(abun_norm, na.rm = TRUE),
    q25  = quantile(abun_norm, 0.25, na.rm = TRUE),
    q75  = quantile(abun_norm, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

pB <- ggplot() +
  # individual taxon lines
  geom_line(data = profiles,
            aes(x = sampling_season, y = abun_norm,
                group = model_id, color = pretty_group),
            alpha = 0.25, linewidth = 0.4) +
  # median ribbon
  geom_ribbon(data = profile_ribbon,
              aes(x = sampling_season, ymin = q25, ymax = q75,
                  group = pretty_group, fill = pretty_group),
              alpha = 0.25) +
  geom_line(data = profile_ribbon,
            aes(x = sampling_season, y = med,
                group = pretty_group, color = pretty_group),
            linewidth = 1.2) +
  facet_grid(fcast_type ~ pretty_group) +
  scale_color_brewer(palette = "Set1", guide = "none") +
  scale_fill_brewer(palette  = "Set1", guide = "none") +
  labs(x = "Phenophase", y = "Normalized abundance (0–1 per taxon)") +
  theme_bw(base_size = 12) +
  theme(
    strip.background   = element_rect(fill = "grey92", color = NA),
    strip.text         = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )

# ── Panel C: labeled dot plot of functional groups by amplitude ───────────────
# One point per functional group taxon, x = amplitude, color = peak phenophase
# Shows WHICH functional groups are most seasonal and WHEN they peak
fg_amp <- sig_max %>%
  filter(fcast_type == "Functional groups") %>%
  arrange(pretty_group, desc(amplitude)) %>%
  mutate(
    label = gsub("_", " ", taxon),
    label = tools::toTitleCase(label)
  )

pC <- ggplot(fg_amp,
             aes(x = amplitude, y = reorder(label, amplitude),
                 color = sampling_season, shape = pretty_group)) +
  geom_segment(aes(xend = 0, yend = reorder(label, amplitude)),
               color = "grey75", linewidth = 0.5) +
  geom_point(size = 3.5) +
  scale_color_manual(values = pheno_colors, name = "Peak phenophase") +
  scale_shape_manual(values = c(Bacteria = 16, Fungi = 17), name = NULL) +
  labs(
    x = "Seasonal amplitude",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    legend.position    = "right"
  )

# ── Combine and save ──────────────────────────────────────────────────────────
fig_peak <- ggarrange(
  ggarrange(pA, pC, ncol = 2, labels = c("A", "C"), widths = c(1, 1.1)),
  pB,
  nrow    = 2,
  labels  = c("", "B"),
  heights = c(1, 1)
)

out_dir <- here("figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "figS13_phenophase_diagram.png"),
       fig_peak, width = 11, height = 9, dpi = 200)

cat("Saved: figures/figS13_phenophase_diagram.png\n")
