# Main figure (fig3): Precision parameter (A), Proportion significant predictors (B),
# Fungal guild phenology (C, formerly D).
# Forecast error by functional group evidence source is now saved as a separate
# supplementary figure (figS_fg_evidence_source) rather than a panel of fig3.
# Source scripts: compare_core_sd_rho.r, fig5_eff_size.r, fig_compareFunctionalCategories.r

source("source.R")
library(ggallin)
library(rstatix)
library(ggpubr)
library(ggh4x)
library(ggrepel)
library(patchwork)
library(data.table)
library(lubridate)

# ── Shared constants ─────────────────────────────────────────────────────────
# kingdom_colors and fcast_type_colors come from source.R.
# Bacteria/Fungi = orange/blue (Wong); Taxonomic/Functional = green/pink (Wong).
BASE_SIZE <- 12
BACT_COLOR  <- kingdom_colors[["Bacteria"]]
FUNGI_COLOR <- kingdom_colors[["Fungi"]]
TAX_COLOR   <- fcast_type_colors[["Taxonomic"]]
FUNC_COLOR  <- fcast_type_colors[["Functional"]]

base_theme <- theme_bw(base_size = BASE_SIZE) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold", size = BASE_SIZE),
    axis.title       = element_text(size = BASE_SIZE),
    panel.grid.minor = element_blank()
  )

# =============================================================================
# Panel A: Precision parameter estimates
# =============================================================================
rho_core_in <- readRDS(here("data", "summary/rho_core_sd_effects.rds")) %>%
  filter(model_name != "all_covariates") %>%
  select(-any_of("pretty_name")) %>%
  mutate(model_id = gsub("_beta_regression$", "", model_id))

driver_uncertainty_pattern <- "20130601_20180101_with_legacy_covariate"
if (nrow(rho_core_in) > 0 && sum(grepl(driver_uncertainty_pattern, rho_core_in$model_id)) < nrow(rho_core_in)) {
  rho_core_in <- rho_core_in %>% filter(grepl(driver_uncertainty_pattern, model_id))
}

in_list <- readRDS(here("data/summary/fcast_horizon_input.rds"))
fcast_horizon_null_site <- in_list[[3]]
if ("model_id" %in% colnames(fcast_horizon_null_site)) {
  fcast_horizon_null_site$model_id <- gsub("_beta_regression$", "", fcast_horizon_null_site$model_id)
}
abundance_by_model <- fcast_horizon_null_site %>%
  group_by(model_id) %>%
  summarise(abundance = mean(abundance, na.rm = TRUE), .groups = "drop")
rho_core_in <- merge(rho_core_in, abundance_by_model, by = "model_id", all.x = TRUE)

precision_data <- rho_core_in %>%
  filter(rowname == "precision") %>%
  mutate(adj_sd = ifelse(is.na(abundance) | abundance == 0, Mean, Mean / abundance))

precision_plot_data <- precision_data %>%
  filter(model_name == "env_cycl") %>%
  mutate(fcast_type = recode(fcast_type, "functional" = "Functional", "taxon" = "Taxonomic"))

precision_stats <- data.frame()
if (nrow(precision_plot_data) >= 4) {
  group_counts_precision <- precision_plot_data %>%
    group_by(fcast_type) %>%
    summarise(n = n(), .groups = "drop")
  if (all(group_counts_precision$n >= 2)) {
    precision_stats <- precision_plot_data %>%
      t_test(adj_sd ~ fcast_type) %>%
      add_significance() %>%
      add_xy_position(x = "fcast_type", dodge = 0.8)
  }
}

# Cap extreme outliers at 99th percentile for cleaner display
cap_val <- quantile(precision_plot_data$adj_sd, 0.99, na.rm = TRUE)
precision_plot_data <- precision_plot_data %>%
  mutate(adj_sd_plot = pmin(adj_sd, cap_val))

pA <- ggplot(precision_plot_data, aes(x = fcast_type, y = adj_sd_plot, fill = fcast_type)) +
  geom_violin(alpha = 0.45, trim = FALSE, draw_quantiles = 0.5, show.legend = FALSE) +
  geom_point(shape = 21, fill = "white", size = 2,
             position = position_jitter(width = 0.1, height = 0),
             alpha = 0.35, show.legend = FALSE) +
  scale_y_log10(breaks = c(10, 100, 1000, 10000, 100000),
                labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                limits = c(10, 2e6)) +
  scale_fill_manual(values = c("Functional" = FUNC_COLOR, "Taxonomic" = TAX_COLOR)) +
  labs(x = NULL, y = "Core variability") +
  base_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

if (nrow(precision_stats) > 0) {
  bracket_y <- 500000
  pA <- pA +
    annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y, linewidth = 0.4) +
    annotate("segment", x = 1, xend = 1, y = bracket_y, yend = bracket_y * 0.7, linewidth = 0.4) +
    annotate("segment", x = 2, xend = 2, y = bracket_y, yend = bracket_y * 0.7, linewidth = 0.4) +
    annotate("text", x = 1.5, y = bracket_y * 1.5, label = precision_stats$p.signif,
             size = 4.5, hjust = 0.5)
}

# =============================================================================
# Panel B: Proportion of significant predictors (bacteria, env_cycl)
# =============================================================================
sum_all <- readRDS(here("data/summary/predictor_effects.rds"))
seasonal_amplitude_in <- readRDS(here("data/summary/seasonal_amplitude.rds"))
cycl_only_vals_scores <- seasonal_amplitude_in[[6]] %>%
  filter(model_name == "cycl_only") %>%
  mutate(cycl_amplitude = amplitude) %>%
  select(-any_of(c("sin", "cos", "max", "amplitude_orig"))) %>%
  pivot_longer(cols = cycl_amplitude, values_to = "effSize", names_to = "beta")
env_cycl_vals_scores <- seasonal_amplitude_in[[6]] %>%
  filter(model_name == "env_cycl") %>%
  mutate(residual_amplitude = amplitude) %>%
  select(-any_of(c("sin", "cos", "max", "amplitude_orig"))) %>%
  pivot_longer(cols = residual_amplitude, values_to = "effSize", names_to = "beta")

df_cal_fg_tax <- sum_all %>%
  filter(time_period == "2013-06_2018-01") %>%
  filter(!beta %in% c("sin", "cos"))
df_cal_fg_tax <- rbindlist(list(df_cal_fg_tax, env_cycl_vals_scores, cycl_only_vals_scores), fill = TRUE)
df_cal_fg_tax <- df_cal_fg_tax %>% filter(time_period != "2015-11_2018-01")
df_cal_fg_tax$fcast_type <- recode(df_cal_fg_tax$fcast_type,
  "functional" = "Functional", "taxon" = "Taxonomic",
  "Functional" = "Functional", "Taxonomic" = "Taxonomic"
)
df_cal_fg_tax$beta_pretty <- recode(df_cal_fg_tax$beta,
  "residual_amplitude" = "Seasonality", "cycl_amplitude" = "Seasonality",
  "pC" = "% Carbon", "LAI" = "Leaf area\nindex",
  "Ectomycorrhizal\ntrees" = "EM trees"
)
df_cal_fg_tax$beta_pretty <- factor(df_cal_fg_tax$beta_pretty,
  levels = c("Seasonality", "EM trees", "Leaf area\nindex",
             "% Carbon", "pH", "Temperature", "Moisture")
)

df_cal_fg_tax_fixed <- df_cal_fg_tax %>%
  mutate(significant = case_when(
    is.na(significant) & beta %in% c("residual_amplitude", "cycl_amplitude") ~
      as.numeric(significant_sin == 1 | significant_cos == 1),
    TRUE ~ significant
  ))

# Deduplicate on model_id + beta_pretty to get one row per model per predictor.
# Use max(significant) so if any duplicate row is significant, the model counts.
df_cal_fg_tax_sig <- df_cal_fg_tax_fixed %>%
  filter(model_name == "env_cycl") %>%
  group_by(fcast_type, model_id, beta_pretty) %>%
  summarise(significant = max(significant, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    significant = ifelse(is.infinite(significant), NA_real_, significant),
    fcast_type = factor(fcast_type, levels = c("Functional", "Taxonomic"))
  )

sig_summary <- df_cal_fg_tax_sig %>%
  group_by(fcast_type, beta_pretty) %>%
  summarise(
    total = n(),
    significant = sum(significant, na.rm = TRUE),
    sig_rate = significant / total,
    .groups = "drop"
  )

sig_test_results <- df_cal_fg_tax_sig %>%
  group_by(beta_pretty) %>%
  summarise(
    functional_sig = sum(significant[fcast_type == "Functional"]),
    functional_total = sum(fcast_type == "Functional"),
    taxonomic_sig = sum(significant[fcast_type == "Taxonomic"]),
    taxonomic_total = sum(fcast_type == "Taxonomic"),
    .groups = "drop"
  ) %>%
  filter(functional_total >= 2 & taxonomic_total >= 2) %>%
  rowwise() %>%
  mutate(
    p_value = tryCatch(
      prop.test(c(functional_sig, taxonomic_sig), c(functional_total, taxonomic_total))$p.value,
      error = function(e) NA_real_
    ),
    sig_label = case_when(
      is.na(p_value) ~ "",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  ungroup()

# Grouped bar chart — predictors on x, grouped by forecast type
pB <- ggplot(sig_summary, aes(x = beta_pretty, y = sig_rate, fill = fcast_type)) +
  geom_col(alpha = 0.7, position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = paste0(significant, "/", total), group = fcast_type),
            position = position_dodge(width = 0.8), vjust = -0.4, size = 2.8) +
  geom_text(data = sig_test_results %>% filter(sig_label != ""),
            aes(x = beta_pretty, y = 1.05, label = sig_label),
            size = 4, color = "black", inherit.aes = FALSE) +
  scale_fill_manual(values = c("Functional" = FUNC_COLOR, "Taxonomic" = TAX_COLOR),
                    name = "Forecast type") +
  scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.25),
                     labels = c("0", "0.25", "0.5", "0.75", "1")) +
  labs(x = NULL, y = "Proportion significant") +
  base_theme +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = BASE_SIZE - 1),
        legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"))

# =============================================================================
# Supplementary figure: Forecast error by functional group evidence source
# (formerly panel C of fig3; now saved standalone as figS_fg_evidence_source)
# =============================================================================
scores_list <- readRDS(here("data", "summary/scoring_metrics_plsr2.rds"))

functional_taxa <- c(
  "cellulose_complex", "acetogen_anaerobic", "assim_nitrate_reduction", "assim_nitrite_reduction",
  "benomyl_antibiotic", "cellobiose_complex", "chitin_complex", "chitinolytic", "copiotroph",
  "dissim_nitrate_reduction", "dissim_nitrite_reduction", "erythromycin_antibiotic",
  "gentamycin_antibiotic", "glucose_simple", "glycerol_simple", "heat_stress", "herbicide_stress",
  "lignolytic", "n_fixation", "oligotroph", "animal_pathogen", "lichenized",
  "streptomycin_antibiotic", "sucrose_complex", "talaromyces",
  "saprotroph", "ectomycorrhizal", "plant_pathogen", "endophyte"
)

# Determine the taxon column name
taxon_col <- intersect(c("taxon", "rank_name", "species"),
                       colnames(scores_list$scoring_metrics_long))[1]

# Filter to env_cycl model only to avoid duplicate points per taxon
fg_source_data <- scores_list$scoring_metrics_long %>%
  filter(.data[[taxon_col]] %in% functional_taxa,
         metric == "RMSE.norm",
         site_prediction == "New time (observed site)",
         model_name == "env_cycl") %>%
  mutate(score = pmax(score, 0)) %>%
  distinct()

fg_source_data$fg_source <- assign_fg_sources(fg_source_data[[taxon_col]])
fg_source_data <- fg_source_data %>% filter(!is.na(fg_source))

# Fungi functional groups come from FUNGuild, not literature review
fg_source_data$fg_source <- ifelse(
  fg_source_data$pretty_group == "Fungi",
  "Scientific consensus (FUNGuild)",
  fg_source_data$fg_source
)

local_pretty_names <- c(
  "assim_nitrite_reduction" = "Assim. nitrite red.",
  "dissim_nitrite_reduction" = "Dissim. nitrite red.",
  "assim_nitrate_reduction" = "Assim. nitrate red.",
  "n_fixation" = "N fixers",
  "dissim_nitrate_reduction" = "Dissim. nitrate red.",
  "chitinolytic" = "Chitin degraders",
  "lignolytic" = "Lignin degraders",
  "copiotroph" = "Copiotrophs",
  "oligotroph" = "Oligotrophs",
  "benomyl_antibiotic" = "Benomyl-res.",
  "glucose_simple" = "Glucose-enr.",
  "streptomycin_antibiotic" = "Streptomycin-res.",
  "sucrose_complex" = "Sucrose-enr.",
  "acetogen_anaerobic" = "Acetogen anaerobic",
  "erythromycin_antibiotic" = "Erythromycin-res.",
  "gentamycin_antibiotic" = "Gentamycin-res.",
  "glycerol_simple" = "Glycerol-enr.",
  "cellobiose_complex" = "Cellobiose-enr.",
  "cellulose_complex" = "Cellulose-enr.",
  "chitin_complex" = "Chitin-enr.",
  "herbicide_stress" = "Herbicide stress-tol.",
  "heat_stress" = "Heat stress-tol.",
  "lichenized" = "Lichenized fungi",
  "animal_pathogen" = "Animal pathogens",
  "saprotroph" = "Saprotrophs",
  "ectomycorrhizal" = "Ectomycorrhizae",
  "plant_pathogen" = "Plant pathogens",
  "endophyte" = "Endophytes"
)
fg_source_data$pretty_fg <- recode(fg_source_data[[taxon_col]], !!!local_pretty_names)

# Shorten source labels for x-axis
fg_source_data$fg_source <- recode(fg_source_data$fg_source,
  "Experimental enrichment" = "Experimental\nenrichment",
  "Literature review" = "Literature\nreview",
  "Literature review + genomic pathway" = "Lit. review +\ngenomic pathway",
  "Scientific consensus (FUNGuild)" = "Scientific consensus\n(FUNGuild)"
)

# Tukey test on fg_source
stat_pvalue_fg_source <- fg_source_data %>%
  rstatix::tukey_hsd(score ~ fg_source)

# One label per unique functional group
fg_label_data <- fg_source_data %>%
  distinct(pretty_fg, fg_source, pretty_group, .keep_all = TRUE)

# Identify significant Tukey comparisons
sig_tukey <- stat_pvalue_fg_source %>% filter(p.adj < 0.05)

# Build manual bracket data for log-scale compatibility
fg_levels <- levels(factor(fg_source_data$fg_source))
x_map <- setNames(seq_along(fg_levels), fg_levels)
max_score <- max(fg_source_data$score, na.rm = TRUE)

bracket_annotations <- list()
if (nrow(sig_tukey) > 0) {
  sig_tukey <- sig_tukey %>%
    filter(grepl("Scientific consensus", group1) | grepl("Scientific consensus", group2))

  bracket_y_start <- max_score * 1.8
  step_mult <- 1.5

  for (i in seq_len(nrow(sig_tukey))) {
    g1 <- sig_tukey$group1[i]
    g2 <- sig_tukey$group2[i]
    x1 <- x_map[g1]
    x2 <- x_map[g2]
    y_bar <- bracket_y_start * step_mult^(i - 1)
    y_tick <- y_bar * 0.85
    lab <- sig_tukey$p.adj.signif[i]

    bracket_annotations <- c(bracket_annotations, list(
      annotate("segment", x = x1, xend = x2, y = y_bar, yend = y_bar, linewidth = 0.35),
      annotate("segment", x = x1, xend = x1, y = y_bar, yend = y_tick, linewidth = 0.35),
      annotate("segment", x = x2, xend = x2, y = y_bar, yend = y_tick, linewidth = 0.35),
      annotate("text", x = (x1 + x2) / 2, y = y_bar * 1.1, label = lab,
               size = 3.5, hjust = 0.5, vjust = 0)
    ))
  }
  y_ceiling <- bracket_y_start * step_mult^nrow(sig_tukey) * 1.3
} else {
  y_ceiling <- max_score * 3
}

p_fg_source <- ggplot(fg_source_data, aes(x = fg_source, y = as.numeric(score),
                                  color = pretty_group)) +
  geom_point(size = 2.5, alpha = 0.5,
             position = position_jitter(width = 0.15, height = 0, seed = 42)) +
  geom_text_repel(data = fg_label_data,
                  aes(label = pretty_fg),
                  size = 2.5, max.overlaps = 25,
                  box.padding = 0.35, point.padding = 0.2, min.segment.length = 0.2,
                  show.legend = FALSE, seed = 42) +
  scale_color_manual(values = c("Bacteria" = BACT_COLOR, "Fungi" = FUNGI_COLOR),
                     name = NULL) +
  scale_y_log10(limits = c(min(fg_source_data$score, na.rm = TRUE) * 0.8, y_ceiling)) +
  labs(x = NULL, y = "Forecast error (nRMSE, log scale)") +
  bracket_annotations +
  base_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = BASE_SIZE - 1),
        legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"))

# =============================================================================
# Panel C (formerly D): Fungal guild seasonal phenology aligned to plant phenophase
# =============================================================================
converged <- scores_list$converged_list

# Element [[6]] contains every model x site x month modeled estimate with its
# assigned phenophase. We need this (not [[4]]) so the per-phenophase means
# reflect modeled abundance across all months in each phenophase, not just the
# single per-site-year peak month.
pheno_data <- readRDS(here("data/clean/pheno_group_peak_phenophases.rds"))[[6]]

fungal_guilds <- c("saprotroph", "ectomycorrhizal", "plant_pathogen",
                   "animal_pathogen")

# Restrict to env_cycl so dormancy predictions are anchored to year-round soil
# temperature and moisture sensors, not pure sinusoidal extrapolation. NEON
# cores are sampled mostly Apr-Oct, so cycl_only's dormancy predictions are
# largely extrapolated and over-weight the long dormancy window.
guild_data <- pheno_data %>%
  filter(taxon %in% fungal_guilds,
         model_id %in% converged,
         model_name == "env_cycl") %>%
  mutate(pretty_name = recode(taxon, !!!microbialForecast:::pretty_names))

# Aggregate: mean abundance per guild x phenophase, min-max scaled
guild_pheno <- guild_data %>%
  group_by(taxon, pretty_name, sampling_season) %>%
  summarise(mean_abun = mean(mean_modeled_abun, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  group_by(taxon) %>%
  mutate(scaled = (mean_abun - min(mean_abun)) / (max(mean_abun) - min(mean_abun))) %>%
  ungroup()

season_labels <- c(dormancy = "Dormancy", greenup = "Green-up",
                   peak = "Peak", greendown = "Senescence")
guild_pheno <- guild_pheno %>%
  mutate(season_label = factor(season_labels[as.character(sampling_season)],
                               levels = season_labels))

# Colors and linetypes for accessibility (labels placed directly on lines).
# Avoids kingdom orange/blue (#E69F00, #0072B2) which are reserved for Bacteria/Fungi.
guild_colors <- c(
  "Saprotrophs"         = "#56B4E9",  # sky blue
  "Ectomycorrhizae"     = "#009E73",  # green
  "Plant pathogens"     = "#D55E00",  # vermillion
  "Animal pathogens"    = "#CC79A7"   # pink
)
guild_linetypes <- c(
  "Saprotrophs"         = "solid",
  "Ectomycorrhizae"     = "dashed",
  "Plant pathogens"     = "dotdash",
  "Animal pathogens"    = "dotted"
)

# Phenophase background shading.
# Peak uses Wong yellow (#F0E442) to avoid colliding with the reserved kingdom orange.
phenophase_fills <- c(
  "Dormancy"   = "grey85",
  "Green-up"   = "#009E73",
  "Peak"       = "#F0E442",
  "Senescence" = "#D55E00"
)

# Labels at right end of each line (Senescence) using ggrepel for auto-separation
label_data <- guild_pheno %>%
  filter(season_label == "Senescence")

pD <- ggplot(guild_pheno, aes(x = season_label, y = scaled,
                               color = pretty_name, linetype = pretty_name,
                               group = pretty_name)) +
  # Phenophase background
  geom_rect(data = data.frame(
    season_label = factor(season_labels, levels = season_labels),
    fill_label = names(phenophase_fills)
  ),
  aes(xmin = as.numeric(season_label) - 0.5,
      xmax = as.numeric(season_label) + 0.5,
      ymin = -Inf, ymax = Inf, fill = fill_label),
  inherit.aes = FALSE, alpha = 0.12) +
  scale_fill_manual(values = phenophase_fills, guide = "none") +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  # Labels at line endpoints, repelled to avoid overlap
  geom_text_repel(data = label_data,
                  aes(label = pretty_name),
                  size = 3.2, fontface = "bold",
                  direction = "y", hjust = 0, nudge_x = 0.15,
                  segment.size = 0.3, segment.color = "grey50",
                  show.legend = FALSE, seed = 42) +
  scale_color_manual(values = guild_colors, guide = "none") +
  scale_linetype_manual(values = guild_linetypes, guide = "none") +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.25))) +
  scale_y_continuous(labels = scales::percent_format(),
                     breaks = seq(0, 1, 0.25)) +
  labs(x = "Plant phenophase",
       y = "Relative seasonal abundance\n(min-max scaled within guild)") +
  base_theme +
  theme(axis.text.x = element_text(size = BASE_SIZE))

# =============================================================================
# Assemble main figure: A (violin) full width on top, B (bars) + C (phenology)
# below. tag_levels = "A" letters the panels A, B, C in patchwork order.
# =============================================================================
bottom_row <- pB + pD + plot_layout(widths = c(2, 1.4))
combined <- pA / bottom_row +
  plot_layout(heights = c(0.85, 1.15)) +
  plot_annotation(tag_levels = "A")

ggsave(here("figures", "fig3_functional_group_error.pdf"), combined,
       width = 13, height = 11, dpi = 300)
ggsave(here("figures", "fig3_functional_group_error.png"), combined,
       width = 13, height = 11, dpi = 300)

cat("Saved: figures/fig3_functional_group_error.pdf\n")
cat("Saved: figures/fig3_functional_group_error.png\n")

# =============================================================================
# Supplementary figure: forecast error by functional group evidence source
# =============================================================================
ggsave(here("figures", "figS_fg_evidence_source.pdf"), p_fg_source,
       width = 9, height = 7, dpi = 300)
ggsave(here("figures", "figS_fg_evidence_source.png"), p_fg_source,
       width = 9, height = 7, dpi = 300)

cat("Saved: figures/figS_fg_evidence_source.pdf\n")
cat("Saved: figures/figS_fg_evidence_source.png\n")
