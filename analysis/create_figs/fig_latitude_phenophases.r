# Latitude × phenophase figure
# Shows how microbial abundance patterns across plant phenological phases
# differ with site latitude — testing whether microbes track plant phenology
# similarly at high vs low latitudes.

library(tidyverse)
library(ggpubr)

source("source.R")

# ── Data loading ─────────────────────────────────────────────────────────────
phenophase_in <- readRDS(here("data/clean/pheno_group_peak_phenophases.rds"))
site_descr    <- readRDS(here("data/clean/site_effect_predictors.rds"))

# Element 4 has per-observation abundance + all metadata
abun_data <- phenophase_in[[4]]   # siteID, dates, mean_modeled_abun, sampling_season,
                                   # fcast_type, pretty_group, model_name, taxon, amplitude,
                                   # significant_sin, significant_cos

# Element 1 has per-model modal phenophase + amplitude
mode_data <- phenophase_in[[1]]   # one row per model_id

# ── Latitude categories ───────────────────────────────────────────────────────
lat_df <- site_descr %>%
  select(siteID, latitude, MAT, MAP) %>%
  mutate(
    latitude_category = case_when(
      latitude > 44  ~ "High (>44\u00b0N)",
      latitude < 31  ~ "Low (<31\u00b0N)",
      TRUE           ~ "Mid (31-44\u00b0N)"
    ),
    latitude_category = factor(
      latitude_category,
      levels = c("Low (<31\u00b0N)", "Mid (31-44\u00b0N)", "High (>44\u00b0N)")
    )
  )

# ── Display settings ──────────────────────────────────────────────────────────
pheno_levels <- c("greenup", "peak", "greendown", "dormancy")
pheno_labels <- c("Green-up", "Peak", "Senescence", "Dormancy")
pheno_colors <- c(
  "Green-up"    = "#009E73",
  "Peak"        = "#E69F00",
  "Senescence"  = "#D55E00",
  "Dormancy"    = "#56B4E9"
)

lat_colors <- c(
  "Low (<31\u00b0N)"  = "#D55E00",
  "Mid (31-44\u00b0N)" = "#E69F00",
  "High (>44\u00b0N)" = "#56B4E9"
)

# ── Panel A: seasonal abundance profiles for key fungal FGs × latitude ────────
focal_fg <- c("ectomycorrhizal", "plant_pathogen", "saprotroph", "animal_pathogen")
fg_labels <- c(
  ectomycorrhizal = "Ectomycorrhizal",
  plant_pathogen  = "Plant pathogen",
  saprotroph      = "Saprotroph",
  animal_pathogen = "Animal pathogen"
)

pA_data <- abun_data %>%
  filter(
    model_name  == "cycl_only",
    fcast_type  == "Functional",
    pretty_group == "Fungi",
    taxon       %in% focal_fg,
    !is.na(mean_modeled_abun),
    !is.na(sampling_season)
  ) %>%
  left_join(lat_df, by = "siteID") %>%
  filter(!is.na(latitude_category)) %>%
  mutate(
    sampling_season = factor(sampling_season, levels = pheno_levels, labels = pheno_labels),
    taxon           = recode(taxon, !!!fg_labels)
  ) %>%
  group_by(taxon, latitude_category, sampling_season) %>%
  summarise(
    mean_abun = mean(mean_modeled_abun, na.rm = TRUE),
    se_abun   = sd(mean_modeled_abun,   na.rm = TRUE) / sqrt(sum(!is.na(mean_modeled_abun))),
    .groups   = "drop"
  )

pA <- ggplot(pA_data,
             aes(x = sampling_season, y = mean_abun,
                 color = latitude_category, group = latitude_category)) +
  geom_line(linewidth = 0.8) +
  geom_errorbar(aes(ymin = mean_abun - se_abun, ymax = mean_abun + se_abun),
                width = 0.15, linewidth = 0.5) +
  geom_point(size = 2.5) +
  facet_wrap(~taxon, scales = "free_y", nrow = 1) +
  scale_color_manual(values = lat_colors, name = "Latitude") +
  labs(
    x = "Plant phenophase",
    y = "Mean modeled abundance"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background  = element_rect(fill = "grey92", color = NA),
    strip.text        = element_text(face = "bold"),
    axis.text.x       = element_text(angle = 30, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.position   = "right"
  )

# ── Panel B: same for bacterial functional groups ─────────────────────────────
focal_bac <- c("cellulolytic", "n_fixation", "denitrification", "nitrification")
bac_labels <- c(
  cellulolytic   = "Cellulolytic",
  n_fixation     = "N-fixation",
  denitrification = "Denitrification",
  nitrification  = "Nitrification"
)
bac_present <- intersect(focal_bac, unique(abun_data$taxon))

if (length(bac_present) >= 2) {
  pB_data <- abun_data %>%
    filter(
      model_name   == "cycl_only",
      fcast_type   == "Functional",
      pretty_group == "Bacteria",
      taxon        %in% bac_present,
      !is.na(mean_modeled_abun),
      !is.na(sampling_season)
    ) %>%
    left_join(lat_df, by = "siteID") %>%
    filter(!is.na(latitude_category)) %>%
    mutate(
      sampling_season = factor(sampling_season, levels = pheno_levels, labels = pheno_labels),
      taxon           = recode(taxon, !!!bac_labels)
    ) %>%
    group_by(taxon, latitude_category, sampling_season) %>%
    summarise(
      mean_abun = mean(mean_modeled_abun, na.rm = TRUE),
      se_abun   = sd(mean_modeled_abun,   na.rm = TRUE) / sqrt(sum(!is.na(mean_modeled_abun))),
      .groups   = "drop"
    )

  pB <- ggplot(pB_data,
               aes(x = sampling_season, y = mean_abun,
                   color = latitude_category, group = latitude_category)) +
    geom_line(linewidth = 0.8) +
    geom_errorbar(aes(ymin = mean_abun - se_abun, ymax = mean_abun + se_abun),
                  width = 0.15, linewidth = 0.5) +
    geom_point(size = 2.5) +
    facet_wrap(~taxon, scales = "free_y", nrow = 1) +
    scale_color_manual(values = lat_colors, name = "Latitude") +
    labs(
      x = "Plant phenophase",
      y = "Mean modeled abundance"
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.background   = element_rect(fill = "grey92", color = NA),
      strip.text         = element_text(face = "bold"),
      axis.text.x        = element_text(angle = 30, hjust = 1),
      panel.grid.major.x = element_blank(),
      legend.position    = "bottom",
      legend.key.size    = unit(0.8, "lines")
    )
} else {
  pB <- NULL
}

# ── Panel C: seasonal CV by latitude — Bacteria and Fungi combined ────────────
# For each site × significantly seasonal taxon, compute seasonal CV.
# Bacteria and Fungi are shown together (colored), faceted by model type.
# Kruskal-Wallis p-value tests whether CV differs across latitude bands.
sig_models <- mode_data %>%
  filter(model_name == "env_cycl",
         significant_sin == 1 | significant_cos == 1) %>%
  pull(model_id)

site_cv <- abun_data %>%
  filter(model_id %in% sig_models,
         !is.na(mean_modeled_abun), !is.na(sampling_season)) %>%
  group_by(model_id, fcast_type, pretty_group, siteID) %>%
  summarise(
    seasonal_cv = sd(mean_modeled_abun, na.rm = TRUE) /
                  mean(mean_modeled_abun, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(lat_df, by = "siteID") %>%
  filter(!is.na(latitude_category), is.finite(seasonal_cv))

kingdom_colors <- c(Bacteria = "#E69F00", Fungi = "#0072B2")

# Compute KW p-values per fcast_type × pretty_group stratum
kw_annot <- site_cv %>%
  group_by(fcast_type, pretty_group) %>%
  summarise(
    p_val = kruskal.test(seasonal_cv ~ latitude_category)$p.value,
    n     = n(),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(pretty_group, " KW p=", format.pval(p_val, digits = 2, eps = 0.001)),
    # Place annotation at top of each facet, centered on mid latitude
    latitude_category = factor("Mid (31-44\u00b0N)",
                               levels = levels(site_cv$latitude_category)),
    seasonal_cv = max(site_cv$seasonal_cv, na.rm = TRUE) *
                  ifelse(pretty_group == "Bacteria", 0.98, 0.88)
  )

pC <- ggplot(site_cv,
             aes(x = latitude_category, y = seasonal_cv, fill = pretty_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.6,
               position = position_dodge(width = 0.75)) +
  geom_point(aes(color = pretty_group),
             position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
             alpha = 0.15, size = 0.8) +
  geom_text(data = kw_annot,
            aes(x = latitude_category, y = seasonal_cv, label = label,
                color = pretty_group),
            inherit.aes = FALSE, size = 3, hjust = 0.5, show.legend = FALSE) +
  facet_wrap(~fcast_type) +
  scale_fill_manual(values  = kingdom_colors, name = "Kingdom") +
  scale_color_manual(values = kingdom_colors, guide = "none") +
  labs(
    x = "Site latitude",
    y = "Seasonal CV\n(SD / mean across phenophases)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background   = element_rect(fill = "grey92", color = NA),
    strip.text         = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    legend.position    = "right"
  )

# ── Panel D: functional group lollipop (from fig_peak_phenophase Panel C) ─────
# All functional groups ranked by cycl_only seasonal amplitude,
# colored by peak phenophase, faded if not significant.
pheno_levels <- c("greenup", "peak", "greendown", "dormancy")
pheno_labels <- c("Green-up", "Peak", "Senescence", "Dormancy")
pheno_colors <- c(
  "Green-up"   = "#009E73",
  "Peak"       = "#E69F00",
  "Senescence" = "#D55E00",
  "Dormancy"   = "#56B4E9"
)

fg_ranked <- mode_data %>%
  filter(model_name == "cycl_only", fcast_type == "Functional",
         significant_sin == 1 | significant_cos == 1,
         !is.na(amplitude)) %>%
  arrange(desc(amplitude)) %>%
  mutate(
    label      = tools::toTitleCase(gsub("_", " ", taxon)),
    peak_phase = factor(sampling_season, levels = pheno_levels, labels = pheno_labels)
  )

pD <- ggplot(fg_ranked,
             aes(x = amplitude,
                 y = reorder(label, amplitude),
                 color = peak_phase)) +
  geom_segment(aes(xend = 0, yend = reorder(label, amplitude)),
               color = "grey75", linewidth = 0.5) +
  geom_point(aes(shape = pretty_group), size = 3) +
  scale_color_manual(values = pheno_colors, name = "Peak phenophase",
                     na.value = "grey60") +
  scale_shape_manual(values = c(Bacteria = 16, Fungi = 17), name = "Kingdom") +
  labs(x = "Seasonal amplitude", y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    legend.position    = "right",
    axis.text.y        = element_text(size = 9)
  )

# ── Combine and save ──────────────────────────────────────────────────────────
pB_noleg <- if (!is.null(pB)) pB + theme(legend.position = "none") else NULL

top_row <- if (!is.null(pB_noleg)) {
  ggarrange(pA, pB_noleg, nrow = 2, labels = c("A", "B"))
} else {
  pA
}

bottom_row <- ggarrange(pC, pD, ncol = 2, labels = c("C", "D"), widths = c(1, 1))

fig_lat <- ggarrange(
  top_row, bottom_row,
  nrow = 2, heights = c(1.6, 1)
)

out_dir <- here("figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "fig_latitude_phenophases.png"),
       fig_lat, width = 15, height = 12, dpi = 200)

cat("Saved: figures/fig_latitude_phenophases.png\n")
