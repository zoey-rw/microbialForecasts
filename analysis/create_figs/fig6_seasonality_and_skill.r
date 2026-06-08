# Figure 6: How seasonality maps onto forecast skill across phenology, guild,
# and geography.
#
#   Panel A — Modeled seasonal abundance profiles for two fungal functional
#             groups (Ectomycorrhizal, Saprotroph) across plant phenophases,
#             split by latitude band. Reuses the focal-FG view from
#             figS10_error_by_latitude Panel A. env_cycl model.
#   Panel B — Relative seasonal abundance (min-max scaled within guild) across
#             plant phenophases for the four converged fungal guilds.
#             Reuses the fig_fg_complex view. env_cycl model.
#   Panel C — Seasonal CV (realized seasonal magnitude) of each significantly
#             seasonal functional group, colored by peak phenophase. Moved here
#             from the retired standalone phenophase figure. env_cycl model.
#   Panel D — Proportion of modeled groups (taxonomic + functional) with
#             skilled forecasts at each NEON site, plotted against site
#             latitude. Reuses Panel F from fig_horizon_by_seasonality.

source("source.R")
suppressMessages({
  library(tidyverse)
  library(data.table)
  library(ggpubr)
  library(ggrepel)
  library(patchwork)
  library(ggtext)
})

# Guilds featured in panels A-B; their labels are bolded in panel C so the
# reader can locate them within the full set of functional groups.
featured_guilds <- c("ectomycorrhizal", "saprotroph",
                     "plant_pathogen", "nitrification")

# ── Shared theme + display constants ─────────────────────────────────────────
# Visual language harmonized across all three panels: same base theme, same
# typography, same latitude palette. Every place latitude appears, it is
# encoded redundantly (color + linetype + point shape) so the figure remains
# readable for colorblind viewers and in greyscale.
BASE_SIZE <- 12
base_theme <- theme_bw(base_size = BASE_SIZE) +
  theme(
    strip.background   = element_rect(fill = "grey92", color = NA),
    strip.text         = element_text(face = "bold", size = BASE_SIZE),
    axis.title         = element_text(size = BASE_SIZE),
    panel.grid.minor   = element_blank(),
    legend.title       = element_text(size = BASE_SIZE - 1, face = "bold"),
    legend.text        = element_text(size = BASE_SIZE - 1),
    plot.tag           = element_text(face = "bold", size = BASE_SIZE + 1)
  )

# Start the phenological year at dormancy so both panels A and B read
# left-to-right in the same temporal order.
pheno_levels <- c("dormancy", "greenup", "peak", "greendown")
pheno_labels <- c("Dormancy", "Green-up", "Peak", "Senescence")

# Latitude encoding — three redundant channels.
# Colors are Wong's vermillion / amber / sky-blue (colorblind-safe).
# Linetypes give a second channel that survives greyscale printing.
# Point shapes use the three most-distinguishable ggplot shapes (16/17/15).
lat_levels <- c("Low (<31°N)", "Mid (31-44°N)", "High (>44°N)")
lat_colors    <- c("Low (<31°N)" = "#D55E00",
                   "Mid (31-44°N)" = "#E69F00",
                   "High (>44°N)" = "#0072B2")
lat_linetypes <- c("Low (<31°N)" = "solid",
                   "Mid (31-44°N)" = "longdash",
                   "High (>44°N)" = "dotted")
lat_shapes    <- c("Low (<31°N)" = 16,   # filled circle
                   "Mid (31-44°N)" = 17, # filled triangle
                   "High (>44°N)" = 15)  # filled square

# =============================================================================
# Panel A — focal FG seasonal profiles by latitude (from figS10 Panel A)
# =============================================================================
phenophase_in <- readRDS(here("data/clean/pheno_group_peak_phenophases.rds"))
site_descr    <- readRDS(here("data/clean/site_effect_predictors.rds"))
abun_data     <- phenophase_in[[6]]

lat_df <- site_descr %>%
  select(siteID, latitude) %>%
  mutate(
    latitude_category = case_when(
      latitude > 44 ~ "High (>44°N)",
      latitude < 31 ~ "Low (<31°N)",
      TRUE          ~ "Mid (31-44°N)"
    ),
    latitude_category = factor(latitude_category, levels = lat_levels)
  )

# Helper to build a band-aggregated phenophase profile for one or more focal
# taxa. One line per latitude band, with color + linetype + point shape all
# mapped to latitude so the encoding survives colorblindness and greyscale.
build_band_panel <- function(focal_taxa, taxon_labels) {
  panel_data <- abun_data %>%
    filter(model_name == "env_cycl",
           fcast_type == "Functional",
           taxon %in% focal_taxa,
           !is.na(mean_modeled_abun),
           !is.na(sampling_season)) %>%
    left_join(lat_df, by = "siteID") %>%
    filter(!is.na(latitude_category)) %>%
    mutate(
      sampling_season = factor(sampling_season, levels = pheno_levels,
                               labels = pheno_labels),
      taxon = factor(recode(taxon, !!!taxon_labels),
                     levels = unname(taxon_labels))
    ) %>%
    group_by(taxon, latitude_category, sampling_season) %>%
    summarise(
      mean_abun = mean(mean_modeled_abun, na.rm = TRUE),
      se_abun   = sd(mean_modeled_abun, na.rm = TRUE) /
                  sqrt(sum(!is.na(mean_modeled_abun))),
      .groups   = "drop"
    )

  ggplot(panel_data,
         aes(x = sampling_season, y = mean_abun,
             color    = latitude_category,
             linetype = latitude_category,
             shape    = latitude_category,
             group    = latitude_category)) +
    geom_errorbar(aes(ymin = mean_abun - se_abun, ymax = mean_abun + se_abun),
                  width = 0.15, linewidth = 0.5, linetype = "solid",
                  show.legend = FALSE) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.6, fill = "white") +
    facet_wrap(~taxon, scales = "free_y", nrow = 1) +
    scale_color_manual(values = lat_colors,    name = "Latitude") +
    scale_linetype_manual(values = lat_linetypes, name = "Latitude") +
    scale_shape_manual(values = lat_shapes,      name = "Latitude") +
    guides(
      color    = guide_legend(
        override.aes = list(linetype = unname(lat_linetypes),
                            shape    = unname(lat_shapes),
                            size     = 3)),
      linetype = "none",
      shape    = "none") +
    labs(x = "Plant phenophase",
         y = "Mean relative abundance (modeled)") +
    base_theme +
    theme(axis.text.x        = element_text(angle = 30, hjust = 1),
          panel.grid.major.x = element_blank(),
          legend.position    = "right",
          legend.key.width   = unit(1.4, "lines"))
}

pA <- build_band_panel(
  focal_taxa = c("ectomycorrhizal", "saprotroph"),
  taxon_labels = c(ectomycorrhizal = "Ectomycorrhizal",
                   saprotroph = "Saprotroph")
)

# =============================================================================
# Panel B — fungal guild phenology (from fig_fg_complex)
# =============================================================================
# Panel B pairs a dormancy-peaked fungal guild (plant pathogen) with a
# greenup-peaked bacterial guild (nitrification). Same per-site + band-median
# overlay so the spread of sites within each latitude band is visible.
pB <- build_band_panel(
  focal_taxa = c("plant_pathogen", "nitrification"),
  taxon_labels = c(plant_pathogen = "Plant pathogen",
                   nitrification = "Nitrification")
)

# =============================================================================
# Panel C — Seasonal magnitude (CV) by functional group
# =============================================================================
# Magnitude = realized seasonal CV of env_cycl modeled abundance: per site,
# SD/mean across the month-of-year climatology (replicate years collapsed to a
# 12-point seasonal cycle), then the median across sites. Collapsing years
# isolates the repeatable seasonal swing from interannual/measurement scatter.
# This replaces the sin/cos harmonic amplitude, which under env_cycl is inflated
# by collinearity with the (also-seasonal) temperature/moisture predictors and
# misranks how seasonal a guild really is.
#
# This panel intentionally encodes magnitude only, NOT a "peak phenophase":
# a single peak phase is not well defined for most guilds. Pooling modeled
# abundance across sites, the four phenophase means are nearly tied for all but
# the most strongly seasonal guilds, and a per-site argmax agrees with the modal
# phase at >=50% of sites for only 8 of 33 groups (e.g. nitrification's peak
# falls in green-up at just 32% of sites). Coloring every guild by an argmax
# phase therefore asserts a continental phenological niche that the data do not
# support and contradicts the per-band curves in panels A-B. Where a guild does
# have a robust, consistent peak it is shown directly in A-B.
seasonal_cv_tax <- abun_data %>%
  filter(model_name == "env_cycl", !is.na(mean_modeled_abun)) %>%
  mutate(moy = lubridate::month(dates)) %>%
  group_by(model_id, siteID, moy) %>%
  summarise(moy_mean = mean(mean_modeled_abun, na.rm = TRUE), .groups = "drop") %>%
  group_by(model_id, siteID) %>%
  summarise(cv = sd(moy_mean, na.rm = TRUE) / mean(moy_mean, na.rm = TRUE),
            .groups = "drop") %>%
  filter(is.finite(cv)) %>%
  group_by(model_id) %>%
  summarise(seasonal_cv = median(cv, na.rm = TRUE), .groups = "drop")

# Every converged functional group (kingdom lookup). No significance filter:
# the metric is realized CV, independent of env_cycl harmonic significance.
fg_meta <- abun_data %>%
  filter(model_name == "env_cycl", fcast_type == "Functional") %>%
  distinct(model_id, taxon, pretty_group)

fg_cv <- seasonal_cv_tax %>%
  inner_join(fg_meta, by = "model_id") %>%
  filter(is.finite(seasonal_cv)) %>%
  mutate(
    label = tools::toTitleCase(gsub("_", " ", taxon)),
    # Bold the guilds featured in panels A-B (rendered via ggtext markdown).
    label_md = ifelse(taxon %in% featured_guilds,
                      paste0("**", label, "**"), label)
  )

pCV <- ggplot(fg_cv,
              aes(x = seasonal_cv, y = reorder(label_md, seasonal_cv),
                  shape = pretty_group)) +
  geom_segment(aes(xend = 0, yend = reorder(label_md, seasonal_cv)),
               color = "grey75", linewidth = 0.5) +
  geom_point(size = 2.8, color = "grey20") +
  scale_shape_manual(values = c(Bacteria = 16, Fungi = 17), name = NULL) +
  labs(x = "Seasonal CV (SD / mean of monthly climatology)", y = NULL) +
  base_theme +
  theme(axis.text.y        = ggtext::element_markdown(),
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.major.x = element_blank(),
        legend.position    = "right")

# =============================================================================
# Panel D — Proportion of skilled groups per site vs latitude (horizon Panel F)
# =============================================================================
fi <- readRDS(here("data/summary/fcast_horizon_input.rds"))
model_mean <- as.data.table(fi[[4]])
null_site  <- as.data.table(fi[[3]])
if ("null_type" %in% names(null_site)) {
  null_site <- null_site[null_type == "site_mean"]
}
nc <- grep("^null_", names(null_site), value = TRUE)
if (length(nc) > 0) setnames(null_site, nc, gsub("^null_", "", nc))

null_rsq <- null_site[model_name == "env_cycl",
                      .(null_RSQ = mean(RSQ.1, na.rm = TRUE)), by = model_id]

site_data <- model_mean[model_name == "env_cycl"]
site_data <- merge(site_data, null_rsq, by = "model_id", all.x = TRUE)

site_hz <- site_data[RSQ.1 > null_RSQ & is.finite(RSQ.1),
                     .(site_horizon = max(months_since_obs, na.rm = TRUE)),
                     by = .(model_id, species, siteID, pretty_group)]
all_combos <- unique(site_data[, .(model_id, species, siteID, pretty_group)])
site_hz <- merge(all_combos, site_hz,
                 by = c("model_id", "species", "siteID", "pretty_group"),
                 all.x = TRUE)
site_hz[is.na(site_horizon), site_horizon := 0]

lat_only <- readRDS(here("data/clean/site_effect_predictors.rds"))[
  , c("siteID", "latitude")]
site_hz <- merge(site_hz, lat_only, by = "siteID", all.x = TRUE)

site_avg <- site_hz[, .(prop_skilled = mean(site_horizon > 0)),
                    by = .(siteID, latitude)] %>%
  as.data.frame() %>%
  mutate(latitude_category = case_when(
    latitude > 44 ~ "High (>44°N)",
    latitude < 31 ~ "Low (<31°N)",
    TRUE          ~ "Mid (31-44°N)"),
    latitude_category = factor(latitude_category, levels = lat_levels))

pC <- ggplot(site_avg,
             aes(x = latitude, y = prop_skilled,
                 color = latitude_category,
                 shape = latitude_category)) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, linewidth = 0.8,
              color = "grey30", fill = "grey80") +
  geom_point(size = 2.6, alpha = 0.85) +
  geom_text_repel(aes(label = siteID), size = 2.6, alpha = 0.8,
                  color = "grey25", max.overlaps = 12,
                  segment.color = "grey60", segment.size = 0.25,
                  show.legend = FALSE) +
  scale_color_manual(values = lat_colors, name = "Latitude") +
  scale_shape_manual(values = lat_shapes, name = "Latitude") +
  guides(
    color = guide_legend(
      override.aes = list(linetype = unname(lat_linetypes),
                          shape    = unname(lat_shapes),
                          size     = 3)),
    shape = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Latitude (°N)",
       y = "Proportion of groups with\nskilled forecasts") +
  base_theme +
  theme(legend.position    = "right",
        legend.key.width   = unit(1.4, "lines"),
        panel.grid.major.x = element_blank())

tryCatch({
  cor_prop <- cor.test(site_avg$latitude, site_avg$prop_skilled)
  prop_label <- paste0("italic(r) == ", round(cor_prop$estimate, 2),
                       "*\", \"~italic(p) == ",
                       round(cor_prop$p.value, 3))
  pC <- pC + annotate("text", x = max(site_avg$latitude) - 2, y = 0.97,
                      label = prop_label, parse = TRUE,
                      size = 3.8, hjust = 1)
}, error = function(e) cat("Latitude correlation failed:", e$message, "\n"))

# =============================================================================
# Compose: A and B span the top two rows (each has internal facets); the bottom
# row pairs the per-guild seasonal-CV panel (C, wider) with the skill-vs-latitude
# scatter (D).
# =============================================================================
# Latitude is mapped to color + linetype/shape in A/B/D. ggplot's automatic
# guide-merging fails across the patchwork, so the Latitude legend is rendered
# once on Panel B and suppressed on A and D; Panel C carries its own
# peak-phenophase legend. A flat wrap_plots design keeps panel tagging reliable.
pA_noleg <- pA + theme(legend.position = "none")
pD_noleg <- pC + theme(legend.position = "none")

design <- "
AAAAA
BBBBB
CCCDD
"
fig6 <- wrap_plots(A = pA_noleg, B = pB, C = pCV, D = pD_noleg,
                   design = design) +
  plot_layout(heights = c(0.85, 0.85, 1.35)) +
  plot_annotation(tag_levels = "A")

out_dir <- here("figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
ggsave(file.path(out_dir, "fig6_seasonality_and_skill.png"),
       fig6, width = 12, height = 15, dpi = 300)

cat("Saved: figures/fig6_seasonality_and_skill.png\n")
