# Latitude × phenophase figure
# Shows how microbial abundance patterns across plant phenological phases
# differ with site latitude — testing whether microbes track plant phenology
# similarly at high vs low latitudes.

library(tidyverse)
library(ggpubr)
library(lme4)
library(lmerTest)   # Satterthwaite p-values for the two-way mixed model in panel B

source("source.R")

# ── Data loading ─────────────────────────────────────────────────────────────
phenophase_in <- readRDS(here("data/clean/pheno_group_peak_phenophases.rds"))
site_descr    <- readRDS(here("data/clean/site_effect_predictors.rds"))

# Element 6 is the full monthly view: every model x site x month modeled
# estimate with its assigned phenophase. We use this for seasonal-profile
# panels so the means reflect modeled abundance across all months in each
# phenophase, not just the per-site-year peak month (which is element 4).
abun_data <- phenophase_in[[6]]

# Element 1 has per-model modal phenophase + amplitude
mode_data <- phenophase_in[[1]]   # one row per model_id

# ── Latitude categories ───────────────────────────────────────────────────────
lat_df <- site_descr %>%
  select(siteID, latitude, MAT, MAP) %>%
  mutate(
    latitude_category = case_when(
      latitude > 44  ~ "High (>44°N)",
      latitude < 31  ~ "Low (<31°N)",
      TRUE           ~ "Mid (31-44°N)"
    ),
    latitude_category = factor(
      latitude_category,
      levels = c("Low (<31°N)", "Mid (31-44°N)", "High (>44°N)")
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
  "Low (<31°N)"  = "#D55E00",
  "Mid (31-44°N)" = "#E69F00",
  "High (>44°N)" = "#56B4E9"
)

# ── Panel A: seasonal relative-abundance profiles for focal FGs × latitude ────
# Three functional groups with the clearest latitude contrasts, spanning both
# kingdoms: Ectomycorrhizal & Saprotroph (Fungi) and N-fixation (Bacteria).
# n_fixation only exists among bacteria and the two fungal names are unique to
# fungi, so filtering on taxon alone unambiguously selects the right kingdom.
focal_fg <- c("ectomycorrhizal", "saprotroph", "n_fixation")
fg_labels <- c(
  ectomycorrhizal = "Ectomycorrhizal",
  saprotroph      = "Saprotroph",
  n_fixation      = "N-fixation"
)

pA_data <- abun_data %>%
  filter(
    # env_cycl: dormancy predictions are anchored to year-round soil sensors
    # rather than extrapolated from the Apr-Oct-dominated NEON sampling window.
    model_name  == "env_cycl",
    fcast_type  == "Functional",
    taxon       %in% focal_fg,
    !is.na(mean_modeled_abun),
    !is.na(sampling_season)
  ) %>%
  left_join(lat_df, by = "siteID") %>%
  filter(!is.na(latitude_category)) %>%
  mutate(
    sampling_season = factor(sampling_season, levels = pheno_levels, labels = pheno_labels),
    taxon           = factor(recode(taxon, !!!fg_labels), levels = unname(fg_labels))
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
    y = "Mean relative abundance (modeled)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background  = element_rect(fill = "grey92", color = NA),
    strip.text        = element_text(face = "bold"),
    axis.text.x       = element_text(angle = 30, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.position   = "right"
  )

# ── Panel B: latitude × kingdom interaction in seasonal variability ───────────
# Seasonal CV (SD / mean of modeled abundance across phenophases) measures how
# strongly a taxon's abundance swings over the season. We test whether that
# variability depends on latitude, on kingdom, and -- the key question -- whether
# the fungi-vs-bacteria difference itself changes across latitude bands.
#
# Unit of analysis is one CV per taxon per latitude band (median across the
# sites it occupies in that band). This avoids the pseudoreplication of the
# previous per-taxon-by-site version, where the same taxon counted once per site
# and inflated the sample size into automatic significance.
sig_models <- mode_data %>%
  filter(model_name == "env_cycl",
         significant_sin == 1 | significant_cos == 1) %>%
  pull(model_id)

site_cv <- abun_data %>%
  filter(model_id %in% sig_models,
         !is.na(mean_modeled_abun), !is.na(sampling_season)) %>%
  group_by(model_id, pretty_group, siteID) %>%
  summarise(
    seasonal_cv = sd(mean_modeled_abun, na.rm = TRUE) /
                  mean(mean_modeled_abun, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(lat_df, by = "siteID") %>%
  filter(!is.na(latitude_category), is.finite(seasonal_cv))

# Collapse to one value per taxon per latitude band
taxon_cv <- site_cv %>%
  group_by(model_id, pretty_group, latitude_category) %>%
  summarise(cv = median(seasonal_cv, na.rm = TRUE), .groups = "drop")

# kingdom_colors comes from source.R

# Two-way mixed model: tests latitude, kingdom, and their interaction, with
# taxon as a random intercept so repeated bands of the same taxon aren't treated
# as independent.
cv_mod <- lmer(cv ~ latitude_category * pretty_group + (1 | model_id),
               data = taxon_cv)
cv_aov <- anova(cv_mod)   # Type III, Satterthwaite df via lmerTest

pfmt <- function(p) ifelse(p < 0.001, "p < 0.001",
                           paste0("p = ", formatC(p, format = "f", digits = 3)))
stat_lab <- paste0(
  "Latitude: ", pfmt(cv_aov["latitude_category", "Pr(>F)"]), "\n",
  "Kingdom: ", pfmt(cv_aov["pretty_group", "Pr(>F)"]), "\n",
  "Latitude × Kingdom: ",
  pfmt(cv_aov["latitude_category:pretty_group", "Pr(>F)"])
)

# Clip the y-axis so the bulk of the distribution is legible; only a couple of
# taxon-band points (both high-CV fungi) fall above the cap.
y_cap   <- 0.82
n_above <- sum(taxon_cv$cv > y_cap)
cap_txt <- if (n_above == 0) NULL else paste0(
  n_above, " high-CV fungal ", if (n_above == 1) "taxon" else "taxa",
  " above the axis cap"
)

pB <- ggplot(taxon_cv,
             aes(x = latitude_category, y = cv,
                 fill = pretty_group, color = pretty_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.45, width = 0.6,
               position = position_dodge(width = 0.75)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12,
                                             dodge.width = 0.75),
             alpha = 0.35, size = 1) +
  # Connect kingdom means across latitude so the interaction is visible: the
  # fungi-bacteria gap is wide at low latitude and narrows toward the poles.
  stat_summary(aes(group = pretty_group), fun = mean, geom = "line",
               linewidth = 0.9, position = position_dodge(width = 0.75),
               show.legend = FALSE) +
  stat_summary(aes(group = pretty_group), fun = mean, geom = "point",
               shape = 18, size = 3, position = position_dodge(width = 0.75),
               show.legend = FALSE) +
  annotate("text", x = 0.55, y = y_cap, label = stat_lab,
           hjust = 0, vjust = 1, size = 3.3, lineheight = 0.95) +
  scale_fill_manual(values  = kingdom_colors, name = NULL) +
  scale_color_manual(values = kingdom_colors, name = NULL) +
  coord_cartesian(ylim = c(0, y_cap)) +
  labs(
    x = "Site latitude",
    y = "Seasonal CV per taxon\n(SD / mean across phenophases)",
    caption = cap_txt
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position    = "right",
    plot.caption       = element_text(size = 8, color = "grey40")
  )

# ── Combine and save ──────────────────────────────────────────────────────────
# A: seasonal relative-abundance profiles (3 focal FGs) across latitude.
# B: latitude x kingdom interaction in seasonal variability.
# Panel B is centered at partial width so a single 3-group panel isn't stretched.
pB_row <- ggarrange(NULL, pB, NULL, ncol = 3, widths = c(0.18, 0.64, 0.18))

fig_lat <- ggarrange(
  pA, pB_row,
  nrow = 2, heights = c(1, 0.95), labels = c("A", "B")
)

out_dir <- here("figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "figS10_error_by_latitude.png"),
       fig_lat, width = 13, height = 9, dpi = 200)

cat("Saved: figures/figS10_error_by_latitude.png\n")
