# How seasonality maps onto forecast skill across guild and geography.
#
#   Panel A — Seasonal CV (realized seasonal magnitude) of every functional
#             group, shaped by kingdom. Featured guilds (Panel B) are bolded.
#   Panel B — Modeled seasonal abundance profiles for four example functional
#             groups (Ectomycorrhizal, Saprotroph, Plant pathogen, Cellulolytic)
#             across plant phenophases, split by latitude band. One stacked
#             column of small multiples. env_cycl model.
#   Panel C — Proportion of modeled groups (taxonomic + functional) with
#             skilled forecasts at each NEON site, plotted against site latitude.
#
# Additional output: fig6_FG_phenophase_survey.png — the Panel B view drawn for
# every converged functional group, so candidate guilds for Panel B can be
# compared side by side.

source("source.R")
suppressMessages({
  library(tidyverse)
  library(data.table)
  library(ggpubr)
  library(ggrepel)
  library(patchwork)
  library(ggtext)
})

# ── Which guilds are featured as example profiles (Panel B) ──────────────────
# The four example profiles, drawn top-to-bottom in this order: two fungal
# guilds with strong latitude structure, then a contrasting pair — a
# dormancy-peaked fungal guild (plant pathogen) and a growing-season-peaked
# bacterial guild (cellulolytic). See the survey figure to compare alternatives.
featured_guilds <- c("ectomycorrhizal", "saprotroph",
                     "plant_pathogen", "cellulolytic")

# Canonical display names for every functional group, from the package. These
# encode the experimental basis of each group (e.g. "Glycerol-enriched",
# "Chloramphenicol-resistant", "Acidic stress-tolerant") rather than the terse
# internal keys, so the figure is interpretable without the methods text.
pretty_names_vec <- unlist(microbialForecast:::pretty_names)
guild_labels <- pretty_names_vec[featured_guilds]   # labels for the examples

# ── Shared theme + display constants ─────────────────────────────────────────
# Visual language harmonized across all panels: same base theme, same
# typography, same latitude palette. Latitude is encoded redundantly
# (color + linetype + point shape) so the figure stays readable for
# colorblind viewers and in greyscale.
BASE_SIZE <- 20
base_theme <- theme_bw(base_size = BASE_SIZE) +
  theme(
    strip.background   = element_rect(fill = "grey92", color = NA),
    strip.text         = element_text(face = "bold", size = BASE_SIZE),
    axis.title         = element_text(size = BASE_SIZE),
    panel.grid.minor   = element_blank(),
    legend.title       = element_text(size = BASE_SIZE - 1, face = "bold"),
    legend.text        = element_text(size = BASE_SIZE - 1),
    legend.background  = element_blank(),
    legend.key         = element_blank(),
    plot.tag           = element_text(face = "bold", size = BASE_SIZE + 2)
  )

# Start the phenological year at dormancy so panels read left-to-right in the
# same temporal order.
pheno_levels <- c("dormancy", "greenup", "peak", "greendown")
pheno_labels <- c("Dormancy", "Green-up", "Peak", "Senescence")

# Latitude encoding — three redundant channels.
# Colors are Wong's vermillion / amber / sky-blue (colorblind-safe).
# Linetypes give a second channel that survives greyscale printing.
# Open point shapes (white fill + colored ring) keep markers legible where
# they sit on top of the connecting lines.
lat_levels    <- c("Low (<31°N)", "Mid (31-44°N)", "High (>44°N)")
lat_colors    <- c("Low (<31°N)"   = "#D55E00",
                   "Mid (31-44°N)" = "#E69F00",
                   "High (>44°N)"  = "#0072B2")
lat_linetypes <- c("Low (<31°N)"   = "solid",
                   "Mid (31-44°N)" = "longdash",
                   "High (>44°N)"  = "dotted")
lat_shapes    <- c("Low (<31°N)"   = 21,   # circle  (open, fillable)
                   "Mid (31-44°N)" = 24,   # triangle
                   "High (>44°N)"  = 22)   # square

# Classify a latitude into the three display bands.
lat_band <- function(latitude) {
  factor(case_when(latitude > 44 ~ "High (>44°N)",
                   latitude < 31 ~ "Low (<31°N)",
                   TRUE          ~ "Mid (31-44°N)"),
         levels = lat_levels)
}

# =============================================================================
# Shared data + band-aggregated phenophase profiles
# =============================================================================
phenophase_in <- readRDS(here("data/clean/pheno_group_peak_phenophases.rds"))
site_descr    <- readRDS(here("data/clean/site_effect_predictors.rds"))
abun_data     <- phenophase_in[[6]]

lat_df <- site_descr %>%
  select(siteID, latitude) %>%
  mutate(latitude_category = lat_band(latitude))

# Mean ± SE modeled abundance per (taxon × latitude band × phenophase) for the
# env_cycl functional-group models. Used by Panels A/B and the survey.
fg_profiles <- abun_data %>%
  filter(model_name == "env_cycl",
         fcast_type == "Functional",
         !is.na(mean_modeled_abun),
         !is.na(sampling_season)) %>%
  left_join(lat_df, by = "siteID") %>%
  filter(!is.na(latitude_category)) %>%
  mutate(sampling_season = factor(sampling_season, levels = pheno_levels,
                                  labels = pheno_labels)) %>%
  group_by(taxon, pretty_group, latitude_category, sampling_season) %>%
  summarise(mean_abun = mean(mean_modeled_abun, na.rm = TRUE),
            se_abun   = sd(mean_modeled_abun, na.rm = TRUE) /
                        sqrt(sum(!is.na(mean_modeled_abun))),
            .groups   = "drop")

# Shared geoms / scales for any latitude-band phenophase profile, so Panels A,
# B and the survey are visually identical.
profile_layers <- function() {
  list(
    geom_errorbar(aes(ymin = mean_abun - se_abun, ymax = mean_abun + se_abun),
                  width = 0.16, linewidth = 0.8, linetype = "solid",
                  alpha = 0.9, show.legend = FALSE),
    geom_line(linewidth = 1.6),
    geom_point(size = 4.8, stroke = 1.6, fill = "white"),
    scale_color_manual(values = lat_colors,    name = "Latitude"),
    scale_linetype_manual(values = lat_linetypes, name = "Latitude"),
    scale_shape_manual(values = lat_shapes,    name = "Latitude"),
    guides(
      color    = guide_legend(
        override.aes = list(linetype = unname(lat_linetypes),
                            shape    = unname(lat_shapes),
                            fill     = "white", stroke = 1.6, size = 4.5,
                            linewidth = 1.2)),
      linetype = "none", shape = "none"),
    labs(x = "Plant phenophase", y = "Mean relative abundance (modeled)"),
    base_theme,
    theme(axis.text.x        = element_text(angle = 30, hjust = 1),
          panel.grid.major.x = element_blank())
  )
}

# Build a stacked (single-column) phenophase panel for a set of focal taxa.
# The example profiles form one continuous column of small multiples (Panel B);
# the survey reuses the same builder across many columns.
build_band_panel <- function(focal_taxa, ncol = 1) {
  fg_profiles %>%
    filter(taxon %in% focal_taxa) %>%
    mutate(taxon = factor(recode(taxon, !!!guild_labels),
                          levels = unname(guild_labels[focal_taxa]))) %>%
    ggplot(aes(x = sampling_season, y = mean_abun,
               color = latitude_category, linetype = latitude_category,
               shape = latitude_category, group = latitude_category)) +
    profile_layers() +
    facet_wrap(~taxon, scales = "free_y", ncol = ncol) +
    theme(legend.position  = "bottom",
          legend.key.width = unit(1.8, "lines"))
}

# Panel B: the four example guilds in one stacked column.
pExamples <- build_band_panel(featured_guilds)

# =============================================================================
# Survey: the Panel B view for every functional group (decision aid)
# =============================================================================
# Order facets by realized seasonal swing (largest first) so the guilds with
# the clearest phenophase structure — the best candidates for Panel B —
# appear at the top-left.
swing_order <- fg_profiles %>%
  group_by(taxon) %>%
  summarise(swing = mean(tapply(mean_abun, latitude_category,
                                function(v) (max(v) - min(v)) / mean(v))),
            .groups = "drop") %>%
  arrange(desc(swing))

survey_df <- fg_profiles %>%
  mutate(
    facet_label = paste0(pretty_names_vec[taxon], " (", pretty_group, ")"),
    facet_label = factor(facet_label,
                         levels = paste0(
                           pretty_names_vec[swing_order$taxon], " (",
                           fg_profiles$pretty_group[
                             match(swing_order$taxon, fg_profiles$taxon)], ")")),
    featured = taxon %in% featured_guilds
  )

survey_fig <- ggplot(survey_df,
                     aes(x = sampling_season, y = mean_abun,
                         color = latitude_category, linetype = latitude_category,
                         shape = latitude_category, group = latitude_category)) +
  profile_layers() +
  facet_wrap(~facet_label, scales = "free_y", ncol = 5) +
  theme(legend.position  = "bottom",
        legend.key.width = unit(1.6, "lines"),
        strip.text       = element_text(size = 10, face = "bold"),
        axis.text        = element_text(size = 8)) +
  labs(title = "Seasonal abundance profile by latitude band — all functional groups",
       subtitle = "Ordered by seasonal swing (largest first); featured guilds appear as the example profiles")

out_dir <- here("figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
ggsave(file.path(out_dir, "fig6_FG_phenophase_survey.png"),
       survey_fig, width = 16, height = 16, dpi = 200)
cat("Saved: figures/fig6_FG_phenophase_survey.png\n")

# =============================================================================
# Panel A — Seasonal magnitude (CV) by functional group
# =============================================================================
# Magnitude = realized seasonal CV of env_cycl modeled abundance: per site,
# SD/mean across the month-of-year climatology (replicate years collapsed to a
# 12-point seasonal cycle), then the median across sites. Collapsing years
# isolates the repeatable seasonal swing from interannual/measurement scatter.
# CV (rather than the sin/cos harmonic amplitude) avoids the collinearity
# between the env_cycl harmonics and the temperature/moisture predictors.
#
# This panel intentionally encodes magnitude only, NOT a "peak phenophase":
# a single peak phase is not well defined for most guilds. Pooling modeled
# abundance across sites, the four phenophase means are nearly tied for all but
# the most strongly seasonal guilds, and a per-site argmax agrees with the modal
# phase at >=50% of sites for only 8 of 33 groups. Where a guild does have a
# robust, consistent peak it is shown directly in the example profiles (B).
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

# Show only the most-seasonal groups; the survey figure has the full set, and
# the "top 20" framing belongs in the caption (the total guild count differs
# from the count used elsewhere in the manuscript, so it is left off the panel).
N_TOP <- 20
fg_cv <- seasonal_cv_tax %>%
  inner_join(fg_meta, by = "model_id") %>%
  filter(is.finite(seasonal_cv)) %>%
  slice_max(seasonal_cv, n = N_TOP) %>%
  mutate(
    label = pretty_names_vec[taxon],
    # Bold the guilds shown as example profiles (rendered via ggtext markdown).
    label_md = ifelse(taxon %in% featured_guilds,
                      paste0("**", label, "**"), label)
  )

# Kingdom is encoded by color (the Bacteria/Fungi palette used throughout the
# manuscript) plus a redundant point shape. This legend is specific to Panel A
# and is drawn inside the panel — it must NOT be merged with the latitude legend
# that serves Panels B and C.
pCV <- ggplot(fg_cv,
              aes(x = seasonal_cv, y = reorder(label_md, seasonal_cv),
                  color = pretty_group, shape = pretty_group)) +
  geom_segment(aes(xend = 0, yend = reorder(label_md, seasonal_cv)),
               color = "grey75", linewidth = 0.9) +
  geom_point(size = 5.2) +
  scale_color_manual(values = kingdom_colors, name = NULL) +
  scale_shape_manual(values = c(Bacteria = 16, Fungi = 17), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  guides(color = guide_legend(override.aes = list(size = 5.5))) +
  labs(x = "Seasonal CV (SD / mean of monthly climatology)", y = NULL) +
  base_theme +
  theme(axis.text.y        = ggtext::element_markdown(size = BASE_SIZE - 2),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
        panel.grid.major.x = element_blank(),
        legend.position      = c(0.98, 0.04),
        legend.justification = c(1, 0),
        legend.background    = element_rect(fill = "white", color = "grey80"),
        legend.margin        = margin(4, 7, 4, 7))

# =============================================================================
# Panel D — Proportion of skilled groups per site vs latitude
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
  mutate(latitude_category = lat_band(latitude))

pD <- ggplot(site_avg,
             aes(x = latitude, y = prop_skilled)) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, linewidth = 1.4,
              color = "grey30", fill = "grey80") +
  geom_point(aes(color = latitude_category, shape = latitude_category),
             size = 4.6, stroke = 1.6, fill = "white") +
  geom_text_repel(aes(label = siteID), size = 4.0, alpha = 0.85,
                  color = "grey25", max.overlaps = 16,
                  segment.color = "grey65", segment.size = 0.3,
                  show.legend = FALSE) +
  scale_color_manual(values = lat_colors, name = "Latitude") +
  scale_shape_manual(values = lat_shapes, name = "Latitude") +
  # This panel carries the shared Latitude legend (also describes Panel B). The
  # override draws the linetype + open marker used in the profile panels.
  guides(
    color = guide_legend(override.aes = list(linetype = unname(lat_linetypes),
                                             shape = unname(lat_shapes),
                                             fill = "white", stroke = 1.1,
                                             size = 4)),
    shape = "none") +
  # No site exceeds ~0.67; cap just above the data so the panel isn't mostly
  # empty white space above the points.
  scale_y_continuous(limits = c(0, 0.75), breaks = seq(0, 0.75, 0.25)) +
  labs(x = "Latitude (°N)",
       y = "Proportion of groups with\nskilled forecasts") +
  base_theme +
  theme(panel.grid.major.x = element_blank(),
        legend.position    = "bottom",
        legend.key.width   = unit(1.8, "lines"))

tryCatch({
  cor_prop <- cor.test(site_avg$latitude, site_avg$prop_skilled)
  prop_label <- paste0("italic(r) == ", round(cor_prop$estimate, 2),
                       "*\", \"~italic(p) == ",
                       round(cor_prop$p.value, 3))
  pD <- pD + annotate("text", x = max(site_avg$latitude) - 2, y = 0.73,
                      label = prop_label, parse = TRUE,
                      size = 6, hjust = 1)
}, error = function(e) cat("Latitude correlation failed:", e$message, "\n"))

# =============================================================================
# Compose: panels read in order. A (seasonal-CV ranking) sits top-left;
# B (the four example guild profiles, one stacked column) sits top-right; C
# (skill vs latitude) spans the full width along the bottom. Legends are NOT
# collected: Panel A carries its own kingdom legend inside the panel, while the
# Latitude legend (shared by B and C) is drawn once beneath C. free() releases
# C from column-alignment with A so it spans the true full width instead of
# being pushed right to line up with A's long guild labels.
# =============================================================================
design <- "
AAAAABBBB
AAAAABBBB
AAAAABBBB
AAAAABBBB
AAAAABBBB
AAAAABBBB
CCCCCCCCC
CCCCCCCCC
CCCCCCCCC
"
fig <- wrap_plots(
  A = pCV,
  B = pExamples + theme(legend.position = "none"),
  C = free(pD),
  design = design
) + plot_annotation(tag_levels = "A")

ggsave(file.path(out_dir, "seasonality_and_skill.png"),
       fig, width = 14, height = 18, dpi = 300)

cat("Saved: figures/seasonality_and_skill.png\n")
