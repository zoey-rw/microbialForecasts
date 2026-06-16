# Main figure: Precision parameter (A), Proportion significant predictors (B).
# Forecast error by functional group evidence source is saved as a separate
# supplementary figure (fg_evidence_source.png).
# Fungal guild phenology is in seasonality_and_skill.r.
# Related scripts: compare_core_sd_rho.r, eff_size.r, fig_compareFunctionalCategories.r

source("source.R")
library(ggallin)
library(rstatix)
library(ggpubr)
library(ggh4x)
library(ggrepel)
library(ggsignif)
library(patchwork)
library(data.table)
library(lubridate)

# ── Shared constants ─────────────────────────────────────────────────────────
# kingdom_colors and fcast_type_colors come from source.R.
# Bacteria/Fungi = orange/blue (Wong); Taxonomic/Functional = green/pink (Wong).
BASE_SIZE <- 12
SIG_SIZE  <- 5.5   # shared significance-asterisk text size across panels
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

# Trim the violins to the data range so their smoothed tails don't run into
# the significance bracket, then place the bracket above all plotted data
# (points + trimmed violin) with clear headroom for the asterisks.
data_max  <- max(precision_plot_data$adj_sd_plot, na.rm = TRUE)
bracket_y <- data_max * 2.5    # bracket bar sits above the highest data
y_ceiling <- data_max * 12     # headroom for the asterisks above the bar

pA <- ggplot(precision_plot_data, aes(x = fcast_type, y = adj_sd_plot, fill = fcast_type)) +
  geom_violin(alpha = 0.45, trim = TRUE, draw_quantiles = 0.5, show.legend = FALSE) +
  geom_point(shape = 21, fill = "white", size = 2,
             position = position_jitter(width = 0.1, height = 0),
             alpha = 0.35, show.legend = FALSE) +
  scale_y_log10(breaks = c(10, 100, 1000, 10000, 100000),
                labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                limits = c(10, y_ceiling)) +
  scale_fill_manual(values = c("Functional" = FUNC_COLOR, "Taxonomic" = TAX_COLOR)) +
  labs(x = NULL, y = "Core variability") +
  base_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

if (nrow(precision_stats) > 0) {
  pA <- pA +
    annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y, linewidth = 0.7) +
    annotate("segment", x = 1, xend = 1, y = bracket_y, yend = bracket_y * 0.65, linewidth = 0.7) +
    annotate("segment", x = 2, xend = 2, y = bracket_y, yend = bracket_y * 0.65, linewidth = 0.7) +
    annotate("text", x = 1.5, y = bracket_y * 2.4, label = precision_stats$p.signif,
             size = SIG_SIZE, fontface = "bold", hjust = 0.5)
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
            size = SIG_SIZE, fontface = "bold", color = "black", inherit.aes = FALSE) +
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
# =============================================================================
scores_list <- readRDS(here("data", "summary/scoring_metrics_plsr2.rds"))

# Use the canonical functional-group list so the Tukey test below sees every
# converged functional group (taxonomic groups, substrate enrichments, and
# bacterial N-cyclers).
functional_taxa <- microbialForecast:::keep_fg_names

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

# Split the comparison into two questions, each with the appropriate test:
#  (1) Kingdom-level: does FUNGuild (the sole fungal assignment method) differ
#      from bacterial FGs as a whole? Pooled Bacteria vs Fungi Wilcoxon, since
#      the fungal side has only one category.
#  (2) Within bacteria: does the evidence type used to assign bacterial FGs
#      predict forecast accuracy? Tukey HSD across the 3 bacterial categories.
# Running them jointly inflates the multiple-comparison correction across
# unrelated questions and conflates the kingdom effect (fungi are systematically
# harder to forecast — see the forecast-error figure) with the evidence-type effect.

stat_pvalue_fg_source <- fg_source_data %>%
  filter(pretty_group == "Bacteria") %>%
  rstatix::tukey_hsd(score ~ fg_source)

king_test <- wilcox.test(score ~ pretty_group, data = fg_source_data)
king_meds <- fg_source_data %>%
  group_by(pretty_group) %>%
  summarise(med = median(score, na.rm = TRUE), n = n(), .groups = "drop")

sig_tukey <- stat_pvalue_fg_source %>% filter(p.adj < 0.05)

fg_levels <- levels(factor(fg_source_data$fg_source))
max_score <- max(fg_source_data$score, na.rm = TRUE)
log_max <- log10(max_score)

# Build inset text panel showing both the kingdom-level Wilcoxon and the
# within-bacteria Tukey results. Use scientific notation for very small
# p-values and an explicit "<" comparison so the direction of each significant
# difference (which side has lower nRMSE, i.e. better forecasts) is visible
# without having to decode signed estimates or compare medians visually.
fmt_p <- function(p) {
  if (p < 0.001) sprintf("p = %.1e", p) else sprintf("p = %.3f", p)
}

# Direction of each Tukey comparison: positive estimate means group2_mean
# exceeds group1_mean (so group1 has lower nRMSE, i.e. better forecasts).
sig_tukey <- sig_tukey %>%
  mutate(
    better = ifelse(estimate > 0, group1, group2),
    worse  = ifelse(estimate > 0, group2, group1)
  )

bact_med <- king_meds$med[king_meds$pretty_group == "Bacteria"]
fungi_med <- king_meds$med[king_meds$pretty_group == "Fungi"]
king_better <- ifelse(bact_med < fungi_med, "Bacteria", "Fungi/FUNGuild")
king_worse  <- ifelse(bact_med < fungi_med, "Fungi/FUNGuild", "Bacteria")

# Build two compact inset strings — one per test — so each take-home result
# stands alone and can be rendered at a readable font size. The Wilcoxon
# inset will sit above the FUNGuild column (the kingdom contrast involves the
# fungal side); the Tukey inset will sit above the bacterial columns.
abbr_cat <- function(x) {
  x <- gsub("Experimental enrichment",               "Exp. enr.",   x)
  x <- gsub("Literature review \\+ genomic pathway", "Lit.+pathway", x)
  x <- gsub("Literature review",                     "Lit. review", x)
  x
}

king_text <- sprintf(
  "Bacteria vs Fungi (Wilcoxon):\n%s (n=%d) < %s (n=%d), %s",
  ifelse(king_better == "Fungi/FUNGuild", "Fungi", king_better),
  king_meds$n[king_meds$pretty_group == "Bacteria"],
  ifelse(king_worse == "Fungi/FUNGuild", "Fungi", king_worse),
  king_meds$n[king_meds$pretty_group == "Fungi"],
  fmt_p(king_test$p.value)
)

if (nrow(sig_tukey) > 0) {
  bact_lines <- paste(
    sprintf("%s < %s, %s",
            abbr_cat(gsub("\n", " ", sig_tukey$better)),
            abbr_cat(gsub("\n", " ", sig_tukey$worse)),
            vapply(sig_tukey$p.adj, fmt_p, character(1))),
    collapse = "\n")
  bact_text <- paste0("Within bacteria (Tukey HSD):\n", bact_lines)
} else {
  bact_text <- NA_character_
}

# Tighten the y-axis to the actual data range (max score ~1.5) instead of
# letting the inset push the ceiling out to 100.
y_ceiling <- 10^(log_max + 0.4)

# Precompute the jittered x positions so the same value drives both the points
# and the repelled labels — using position_jitter() separately on each layer
# would assign independent random offsets and the leader lines would point to
# the wrong dot. A continuous x scale lets us share x_jit across layers and
# place the inset label at an arbitrary numeric coordinate.
fg_levels_order <- fg_levels
fg_source_data <- fg_source_data %>%
  mutate(x_factor = factor(fg_source, levels = fg_levels_order),
         x_int    = as.numeric(x_factor))
set.seed(42)
fg_source_data$x_jit <- fg_source_data$x_int +
  runif(nrow(fg_source_data), -0.15, 0.15)

# Carry the jitter through to the label-anchor table: one row per group, so
# distinct() collapses each group's many rows to the first observation's x_jit
# (and matching score), which the label then anchors to via a leader line.
fg_label_data <- fg_source_data %>%
  distinct(pretty_fg, fg_source, pretty_group, .keep_all = TRUE) %>%
  mutate(pretty_fg = gsub("_antibiotic$", "-res.", pretty_fg),
         pretty_fg = gsub("_simple$",     "-enr.", pretty_fg),
         pretty_fg = gsub("_complex$",    "-enr.", pretty_fg),
         pretty_fg = gsub("_stress$",     " stress", pretty_fg),
         pretty_fg = gsub("_",            " ", pretty_fg),
         pretty_fg = sub("^(.)", "\\U\\1", pretty_fg, perl = TRUE))

# Label only the 3 best and 3 worst groups per evidence-type column so the
# Experimental-enrichment cluster (22 groups) doesn't drown the others, while
# the bottom-of-rank and top-of-rank exemplars in each category are still
# named. Categories with <= 6 groups get every group labelled.
fg_label_extremes <- fg_label_data %>%
  group_by(fg_source) %>%
  arrange(score, .by_group = TRUE) %>%
  mutate(rk_low = row_number(), rk_high = n() - rk_low + 1) %>%
  filter(rk_low <= 3 | rk_high <= 3) %>%
  ungroup() %>%
  select(-rk_low, -rk_high)

# Place two insets at the top, each near the data it summarises:
#   - within-bacteria Tukey above the bacterial columns (top-left)
#   - kingdom Wilcoxon above the FUNGuild column (top-right)
# Both sit at the same y above all labels (Ectomycorrhizae at y~1.7,
# Nitrification at y~0.75, Pyruvate-enr. at y~0.4), so the bigger font size
# can't crowd any leader lines.
inset_y <- 10^(log_max + 0.2)
label_aes <- list(size = 3.5, color = "grey15", fill = "white",
                  label.r = unit(0.15, "lines"),
                  label.padding = unit(0.45, "lines"))
bracket_annotations <- list(
  do.call(annotate, c(list("label", x = 0.55, y = inset_y,
                           label = bact_text, hjust = 0, vjust = 1),
                      label_aes)),
  do.call(annotate, c(list("label",
                           x = length(fg_levels_order) + 0.45, y = inset_y,
                           label = king_text, hjust = 1, vjust = 1),
                      label_aes))
)
if (is.na(bact_text)) bracket_annotations <- bracket_annotations[2]

p_fg_source <- ggplot(fg_source_data, aes(x = x_jit, y = as.numeric(score),
                                  color = pretty_group)) +
  geom_point(size = 2.5, alpha = 0.5) +
  geom_text_repel(data = fg_label_extremes,
                  aes(label = pretty_fg),
                  size = 3.2, max.overlaps = Inf,
                  box.padding = 0.5, point.padding = 0.3,
                  min.segment.length = 0,
                  force = 3, force_pull = 0.7,
                  max.iter = 30000, max.time = 2,
                  segment.size = 0.25, segment.alpha = 0.7,
                  segment.color = "grey50",
                  show.legend = FALSE, seed = 42) +
  scale_color_manual(values = c("Bacteria" = BACT_COLOR, "Fungi" = FUNGI_COLOR),
                     name = NULL) +
  scale_x_continuous(breaks = seq_along(fg_levels_order),
                     labels = fg_levels_order,
                     limits = c(0.5, length(fg_levels_order) + 0.5),
                     expand = expansion(0)) +
  scale_y_log10(limits = c(min(fg_source_data$score, na.rm = TRUE) * 0.8, y_ceiling)) +
  labs(x = NULL, y = "Forecast error (nRMSE, log scale)") +
  bracket_annotations +
  base_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = BASE_SIZE - 1),
        legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"))

# =============================================================================
# Assemble main figure: A (precision violin) and B (predictor significance) side
# by side. Fungal guild phenology lives in seasonality_and_skill.r.
# tag_levels = "A" letters the panels A, B in patchwork order.
# =============================================================================
combined <- pA + pB +
  plot_layout(widths = c(1, 1.6)) +
  plot_annotation(tag_levels = "A")

ggsave(here("figures", "functional_group_error.png"), combined,
       width = 13, height = 6, dpi = 300)

cat("Saved: figures/functional_group_error.png\n")

# =============================================================================
# Supplementary figure: forecast error by functional group evidence source
# =============================================================================
ggsave(here("figures", "fg_evidence_source.png"), p_fg_source,
       width = 9, height = 7, dpi = 300)

cat("Saved: figures/fg_evidence_source.png\n")
