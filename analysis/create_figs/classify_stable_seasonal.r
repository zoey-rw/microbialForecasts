# Classify each group based on seasonality and degree of temporal persistence
# using lag-1 (rho) parameter and seasonal amplitude from env_cycl models.
# Since models use a bounded prior (rho ~ dbeta(1,1)), rho is constrained to (0,1).
# We classify groups as high vs. low persistence instead of stable vs. chaotic.

source("source.R")
library(dplyr)
library(here)

# Optional: restrict to converged models
scores_list <- readRDS(here("data", "summary/scoring_metrics_plsr2.rds"))
converged <- scores_list$converged_list
converged_base <- gsub("_(combined|beta_regression)$", "", converged)

# Load rho (lag-1) estimates: one row per model_id per model_name
rho_in <- readRDS(here("data", "summary/rho_core_sd_effects.rds")) %>%
  filter(model_name != "all_covariates", rowname == "rho") %>%
  mutate(model_id_base = gsub("_(combined|beta_regression)$", "", model_id))

# Assert one row per base ID — silently dropping duplicates with slice() would
# hide upstream errors (suffix-stripping collapsing two genuine models,
# accidental double-combine, etc.). Better to fail loudly.
assert_unique <- function(df, key, where) {
  dups <- df %>% count(.data[[key]]) %>% filter(n > 1)
  if (nrow(dups) > 0) {
    stop("Duplicate ", key, " in ", where, " (", nrow(dups),
         " keys); investigate before continuing. First few:\n",
         paste(utils::capture.output(print(head(dups, 5))), collapse = "\n"))
  }
  df
}

# Use env_cycl so both rho and seasonal terms come from the same model
rho_env_cycl <- rho_in %>%
  filter(model_name == "env_cycl") %>%
  filter(model_id_base %in% converged_base) %>%
  transmute(
    model_id = model_id_base,
    taxon,
    pretty_group,
    rank,
    rank_only,
    fcast_type,
    rho_mean = Mean,
    rho_sd = SD,
    abs_rho = abs(Mean)
  ) %>%
  assert_unique("model_id", "rho_env_cycl")

# Load seasonal amplitude
seasonal_in <- readRDS(here("data/summary/seasonal_amplitude.rds"))
seas_vals <- seasonal_in[[6]]

# Normalize model_id and filter to same period and converged
seas_env_cycl <- seas_vals %>%
  filter(model_name == "env_cycl") %>%
  mutate(model_id_base = gsub("_(combined|beta_regression)$", "", model_id)) %>%
  filter(model_id_base %in% converged_base)

# Standardize column names
if (!"significant_sin" %in% names(seas_env_cycl) && "Mean_sin" %in% names(seas_env_cycl)) {
  seas_env_cycl <- seas_env_cycl %>% rename(sin = Mean_sin, cos = Mean_cos)
}
if ("significant_sin" %in% names(seas_env_cycl)) {
  seas_env_cycl <- seas_env_cycl %>%
    mutate(significant_seasonal = (significant_sin == 1 | significant_cos == 1))
} else {
  seas_env_cycl <- seas_env_cycl %>% mutate(significant_seasonal = NA)
}

seas_env_cycl <- seas_env_cycl %>%
  select(model_id_base, amplitude, significant_seasonal) %>%
  assert_unique("model_id_base", "seas_env_cycl") %>%
  rename(model_id = model_id_base)

# Merge rho and seasonal amplitude on model_id
class_df <- rho_env_cycl %>%
  inner_join(seas_env_cycl, by = "model_id")

# Thresholds for the quadrant classification. The persistence threshold is
# 0.1 (not 0.3) because env_cycl posterior ρ rarely exceeds 0.12; a higher
# cutoff would leave the "high persistence" quadrant empty and waste the
# visual contrast in the figure.
AMPLITUDE_THRESHOLD <- 0.05
PERSISTENCE_THRESHOLD <- 0.1

class_df <- class_df %>%
  mutate(
    is_seasonal = (amplitude > AMPLITUDE_THRESHOLD) | (significant_seasonal %in% TRUE),
    is_high_persistence = abs_rho >= PERSISTENCE_THRESHOLD,
    classification = case_when(
      is_seasonal &  is_high_persistence ~ "seasonal_high_persistence",
      is_seasonal & !is_high_persistence ~ "seasonal_low_persistence",
      !is_seasonal &  is_high_persistence ~ "non_seasonal_high_persistence",
      !is_seasonal & !is_high_persistence ~ "non_seasonal_low_persistence",
      TRUE ~ NA_character_
    )
  )

# Summary: counts and proportions
summary_counts <- class_df %>%
  count(classification, name = "n") %>%
  mutate(prop = n / sum(n))
total_n <- nrow(class_df)

cat("\n=== DEGREE OF PERSISTENCE (rho >=", PERSISTENCE_THRESHOLD, "vs <", PERSISTENCE_THRESHOLD, ") ===\n")
cat("High Persistence (rho >= ", PERSISTENCE_THRESHOLD, "): ", sum(class_df$is_high_persistence, na.rm = TRUE), " of ", total_n, "\n", sep = "")
cat("Low Persistence (rho < ", PERSISTENCE_THRESHOLD, "): ", sum(!class_df$is_high_persistence, na.rm = TRUE), " of ", total_n, "\n", sep = "")
cat("Rho range:", round(min(class_df$rho_mean, na.rm = TRUE), 4), "to", round(max(class_df$rho_mean, na.rm = TRUE), 4), "\n")

cat("\n=== SEASONAL vs NON-SEASONAL (amplitude >", AMPLITUDE_THRESHOLD, "or sig terms) ===\n")
cat("Seasonal:", sum(class_df$is_seasonal, na.rm = TRUE), "of", total_n, "\n")
cat("Non-seasonal:", sum(!class_df$is_seasonal, na.rm = TRUE), "of", total_n, "\n")

cat("\n=== FOUR-WAY CLASSIFICATION (env_cycl, converged) ===\n")
print(summary_counts)
cat("\nProportions:\n")
print(summary_counts %>% mutate(prop = round(prop, 3)))

# By kingdom
by_kingdom <- class_df %>%
  count(pretty_group, classification, name = "n") %>%
  group_by(pretty_group) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
cat("\nBy kingdom (pretty_group):\n")
print(as.data.frame(by_kingdom))

# Save
out_path <- here("data", "summary/persistence_seasonal_classification.rds")
saveRDS(list(
  classification = class_df,
  summary_counts = summary_counts,
  by_kingdom = by_kingdom,
  amplitude_threshold = AMPLITUDE_THRESHOLD,
  persistence_threshold = PERSISTENCE_THRESHOLD,
  n_total = total_n
), out_path)
cat("\nSaved classification and summary to", out_path, "\n")

# ── Figure 1: Scatter plot of rho vs amplitude with marginal densities ───────
library(ggplot2)
library(cowplot)
library(ggrepel)
library(ggExtra)
library(tidyr)

# kingdom_colors comes from source.R

fig_out_dir <- here("figures")
if (!dir.exists(fig_out_dir)) dir.create(fig_out_dir, recursive = TRUE)

# Build env_cycl-equivalent dataset for any model that has seasonal amplitude.
# env_cov has rho but no sin/cos terms, so only env_cycl and cycl_only qualify.
build_class <- function(model_n) {
  r <- rho_in %>%
    filter(model_name == model_n, model_id_base %in% converged_base) %>%
    transmute(model_id = model_id_base, taxon, pretty_group, rank, rank_only,
              fcast_type, rho_mean = Mean, abs_rho = abs(Mean),
              model_name = model_n) %>%
    assert_unique("model_id", paste0("rho (", model_n, ")"))
  s <- seas_vals %>%
    filter(model_name == model_n) %>%
    mutate(model_id_base = gsub("_(combined|beta_regression)$", "", model_id)) %>%
    filter(model_id_base %in% converged_base)
  if (!"significant_sin" %in% names(s) && "Mean_sin" %in% names(s)) {
    s <- s %>% rename(sin = Mean_sin, cos = Mean_cos)
  }
  if ("significant_sin" %in% names(s)) {
    s <- s %>% mutate(significant_seasonal = (significant_sin == 1 | significant_cos == 1))
  } else {
    s <- s %>% mutate(significant_seasonal = NA)
  }
  s <- s %>% select(model_id_base, amplitude, significant_seasonal) %>%
    assert_unique("model_id_base", paste0("seas (", model_n, ")")) %>%
    rename(model_id = model_id_base)
  inner_join(r, s, by = "model_id") %>%
    mutate(
      is_seasonal = (amplitude > AMPLITUDE_THRESHOLD) | (significant_seasonal %in% TRUE),
      is_high_persistence = abs_rho >= PERSISTENCE_THRESHOLD,
      classification = case_when(
        is_seasonal &  is_high_persistence ~ "seasonal_high_persistence",
        is_seasonal & !is_high_persistence ~ "seasonal_low_persistence",
        !is_seasonal &  is_high_persistence ~ "non_seasonal_high_persistence",
        !is_seasonal & !is_high_persistence ~ "non_seasonal_low_persistence",
        TRUE ~ NA_character_
      )
    )
}

# Helper: produce a scatter with marginal kernel densities for one model.
# Labels n_per_quadrant most-extreme taxa per occupied quadrant.
make_marginal_scatter <- function(df, model_label, n_per_quadrant = 3) {
  df <- df %>%
    group_by(classification) %>%
    mutate(extreme_score = case_when(
      classification == "seasonal_high_persistence"     ~ amplitude + abs_rho,
      classification == "seasonal_low_persistence"     ~ amplitude + (1 - abs_rho),
      classification == "non_seasonal_high_persistence" ~ (1 - amplitude) + abs_rho,
      classification == "non_seasonal_low_persistence" ~ (1 - amplitude) + (1 - abs_rho),
      TRUE ~ 0
    )) %>%
    arrange(desc(extreme_score)) %>%
    mutate(rank_in_quad = row_number(),
           label_taxon = ifelse(rank_in_quad <= n_per_quadrant, taxon, NA_character_)) %>%
    ungroup() %>%
    select(-extreme_score, -rank_in_quad)

  # Breakdown line: collapse quadrant counts into a single subtitle so they
  # never compete with data points for in-panel real estate. Compact labels
  # (S-Lo / NS-Hi) keep the line short enough to fit alongside ggMarginal.
  pretty_class <- c(
    seasonal_high_persistence     = "S-Hi",
    seasonal_low_persistence      = "S-Lo",
    non_seasonal_high_persistence = "NS-Hi",
    non_seasonal_low_persistence  = "NS-Lo"
  )
  breakdown <- df %>%
    count(classification, name = "n") %>%
    mutate(prop = n / sum(n),
           pretty = paste0(pretty_class[classification], " ",
                           n, " (", round(prop * 100), "%)")) %>%
    arrange(desc(n)) %>%
    pull(pretty) %>%
    paste(collapse = " · ")
  breakdown <- paste0(
    "S = seasonal, NS = non-seasonal, Hi/Lo = high/low ρ",
    "\n", breakdown)

  # Data-driven axis limits so the data fills the panel
  x_max <- max(df$amplitude, na.rm = TRUE) * 1.05
  y_max <- max(df$abs_rho,   na.rm = TRUE) * 1.12
  threshold_visible <- PERSISTENCE_THRESHOLD <= y_max

  # Test whether Bacteria vs Fungi distributions differ on each axis.
  # Annotate both results so a non-significant axis is just as visible as a
  # significant one.
  fmt_p <- function(p) {
    if (is.na(p)) return("n.s.")
    star <- if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "n.s."
    ptext <- if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)
    paste0(star, " (", ptext, ")")
  }
  kingdoms <- unique(df$pretty_group)
  if (length(kingdoms) == 2) {
    p_rho <- suppressWarnings(wilcox.test(abs_rho   ~ pretty_group, data = df))$p.value
    p_amp <- suppressWarnings(wilcox.test(amplitude ~ pretty_group, data = df))$p.value
    test_label <- paste0("Bacteria vs Fungi (Wilcoxon):\n",
                         "ρ:  ", fmt_p(p_rho), "\n",
                         "amplitude:  ", fmt_p(p_amp))
  } else {
    p_rho <- NA_real_; p_amp <- NA_real_; test_label <- NULL
  }

  base <- ggplot(df, aes(x = amplitude, y = abs_rho,
                         color = pretty_group, fill = pretty_group)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_text_repel(aes(label = label_taxon),
                    size = 2.8, fontface = "italic",
                    max.overlaps = Inf, na.rm = TRUE,
                    force = 4, force_pull = 0.5,
                    segment.color = "gray55", segment.size = 0.3,
                    box.padding = 0.45, point.padding = 0.25,
                    min.segment.length = 0.1, show.legend = FALSE) +
    geom_vline(xintercept = AMPLITUDE_THRESHOLD, linetype = "dashed", color = "gray45")

  if (threshold_visible) {
    base <- base +
      geom_hline(yintercept = PERSISTENCE_THRESHOLD, linetype = "dashed", color = "gray45") +
      annotate("text", x = x_max * 0.99, y = PERSISTENCE_THRESHOLD,
               label = "  high persistence  ", size = 2.8,
               fontface = "italic", color = "gray45",
               hjust = 1, vjust = -0.4) +
      annotate("text", x = x_max * 0.99, y = PERSISTENCE_THRESHOLD,
               label = "  low persistence  ", size = 2.8,
               fontface = "italic", color = "gray45",
               hjust = 1, vjust = 1.4)
  }

  if (!is.null(test_label)) {
    base <- base +
      annotate("label", x = x_max * 0.98, y = y_max * 0.97,
               label = test_label,
               hjust = 1, vjust = 1, size = 2.7, color = "gray15",
               fill = "white", label.size = 0.3, label.r = unit(0.1, "lines"),
               lineheight = 1.1)
  }

  base <- base +
    scale_color_manual(values = kingdom_colors, name = NULL) +
    scale_fill_manual(values = kingdom_colors, name = NULL) +
    coord_cartesian(xlim = c(0, x_max), ylim = c(0, y_max), clip = "off") +
    labs(x = "Seasonal amplitude",
         y = expression("Temporal persistence (" * rho * ")"),
         title = paste0("Model: ", model_label),
         subtitle = breakdown) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top",
          legend.justification = "right",
          legend.background = element_blank(),
          legend.margin = margin(t = 0, b = 0),
          plot.title = element_text(size = 11, color = "gray15"),
          plot.subtitle = element_text(size = 9, color = "gray40",
                                       margin = margin(b = 6)),
          plot.title.position = "plot",
          panel.grid.minor = element_blank())

  list(
    plot  = ggMarginal(base, type = "density",
                       groupColour = TRUE, groupFill = TRUE, alpha = 0.35),
    tests = list(p_rho = p_rho, p_amplitude = p_amp,
                 n_bacteria = sum(df$pretty_group == "Bacteria"),
                 n_fungi    = sum(df$pretty_group == "Fungi"))
  )
}

# Build the two datasets and render
class_env_cycl  <- build_class("env_cycl")
class_cycl_only <- build_class("cycl_only")

scatter_env_cycl <- make_marginal_scatter(
  class_env_cycl, "env_cycl (environmental + seasonal)")

# Manuscript Figure S19. cycl_only data is still built above so the by-model
# diagnostic figure below can use it; we just don't publish it as a separate
# figure here.
ggsave(file.path(fig_out_dir, "figS19_persistence_seasonality_env_cycl.png"),
       scatter_env_cycl$plot, width = 8, height = 6.5, dpi = 200)
cat("Saved: figures/figS19_persistence_seasonality_env_cycl.png\n")

report_tests <- function(scatter, label) {
  t <- scatter$tests
  fmt <- function(p) if (is.na(p)) "n.s." else if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)
  cat(sprintf("[%s] n_bacteria=%d  n_fungi=%d  ρ %s  amplitude %s\n",
              label, t$n_bacteria, t$n_fungi, fmt(t$p_rho), fmt(t$p_amplitude)))
}
cat("\nValues to cross-check against manuscript.tex Fig S19 caption:\n")
report_tests(scatter_env_cycl, "S19 env_cycl")
cat("\nQuadrant breakdown (env_cycl, persistence threshold = ",
    PERSISTENCE_THRESHOLD, "):\n", sep = "")
print(class_env_cycl %>% count(classification, name = "n") %>%
      mutate(pct = round(100 * n / sum(n))))

# Keep class_df pointing at env_cycl for downstream code (rank panel, exemplars)
class_df <- class_env_cycl

# ── Figure 1b: Consistency across model types ────────────────────────────────

class_all <- bind_rows(class_env_cycl, class_cycl_only) %>%
  mutate(model_name = factor(model_name, levels = c("env_cycl", "cycl_only"),
                             labels = c("env_cycl (environmental + seasonal)",
                                        "cycl_only (seasonal only)")))

fig_by_model <- ggplot(class_all, aes(x = amplitude, y = abs_rho, color = pretty_group)) +
  geom_point(size = 1.8, alpha = 0.6) +
  geom_hline(yintercept = PERSISTENCE_THRESHOLD, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = AMPLITUDE_THRESHOLD, linetype = "dashed", color = "gray50") +
  facet_wrap(~ model_name) +
  scale_color_manual(values = kingdom_colors, name = NULL) +
  labs(x = "Seasonal amplitude",
       y = expression("Temporal persistence (" * rho * ")"),
       title = "Persistence vs. seasonality across model structures") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave(file.path(fig_out_dir, "fig_persistence_seasonality_by_model.png"), fig_by_model,
       width = 9, height = 5, dpi = 200)
cat("Saved: figures/fig_persistence_seasonality_by_model.png\n")

# Within-taxon agreement: same taxon estimated under both models
agree_wide <- class_all %>%
  mutate(model_short = ifelse(grepl("env_cycl", model_name), "env_cycl", "cycl_only")) %>%
  select(model_short, taxon, pretty_group, abs_rho, amplitude) %>%
  pivot_wider(names_from = model_short, values_from = c(abs_rho, amplitude))
agree_cor <- agree_wide %>%
  summarise(
    rho_cor   = suppressWarnings(cor(abs_rho_env_cycl, abs_rho_cycl_only, use = "complete.obs")),
    amp_cor   = suppressWarnings(cor(amplitude_env_cycl, amplitude_cycl_only, use = "complete.obs")),
    rho_cor_b = suppressWarnings(cor(abs_rho_env_cycl[pretty_group == "Bacteria"],
                                     abs_rho_cycl_only[pretty_group == "Bacteria"], use = "complete.obs")),
    rho_cor_f = suppressWarnings(cor(abs_rho_env_cycl[pretty_group == "Fungi"],
                                     abs_rho_cycl_only[pretty_group == "Fungi"], use = "complete.obs")),
    amp_cor_b = suppressWarnings(cor(amplitude_env_cycl[pretty_group == "Bacteria"],
                                     amplitude_cycl_only[pretty_group == "Bacteria"], use = "complete.obs")),
    amp_cor_f = suppressWarnings(cor(amplitude_env_cycl[pretty_group == "Fungi"],
                                     amplitude_cycl_only[pretty_group == "Fungi"], use = "complete.obs"))
  )
cat("\nWithin-taxon Pearson correlations between env_cycl and cycl_only estimates:\n")
print(round(as.data.frame(agree_cor), 3))

# ── Figure 1c: Consistency across taxonomic rank and functional groups ───────
class_df_ranks <- class_df %>%
  mutate(
    rank_lab = factor(
      rank_only,
      levels = c("phylum", "class", "order", "family", "genus", "functional"),
      labels = c("Phylum", "Class", "Order", "Family", "Genus", "Functional group")
    )
  )

fig_by_rank <- ggplot(class_df_ranks, aes(x = amplitude, y = abs_rho, color = pretty_group)) +
  geom_point(size = 1.8, alpha = 0.65) +
  geom_hline(yintercept = PERSISTENCE_THRESHOLD, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = AMPLITUDE_THRESHOLD, linetype = "dashed", color = "gray50") +
  facet_wrap(~ rank_lab, ncol = 3) +
  scale_color_manual(values = kingdom_colors, name = NULL) +
  labs(x = "Seasonal amplitude",
       y = expression("Temporal persistence (" * rho * ")"),
       title = "Persistence vs. seasonality across taxonomic ranks (env_cycl)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave(file.path(fig_out_dir, "fig_persistence_seasonality_by_rank.png"), fig_by_rank,
       width = 9, height = 6.5, dpi = 200)
cat("Saved: figures/fig_persistence_seasonality_by_rank.png\n")

# Consistency tests: do distributions differ by rank / kingdom?
cat("\nKruskal-Wallis tests across rank (env_cycl):\n")
cat("  amplitude ~ rank_only: p =",
    signif(kruskal.test(amplitude ~ rank_only, data = class_df_ranks)$p.value, 3), "\n")
cat("  abs_rho   ~ rank_only: p =",
    signif(kruskal.test(abs_rho   ~ rank_only, data = class_df_ranks)$p.value, 3), "\n")
cat("Wilcoxon tests across kingdom (env_cycl):\n")
cat("  amplitude (Bacteria vs Fungi): p =",
    signif(wilcox.test(amplitude ~ pretty_group, data = class_df_ranks)$p.value, 3), "\n")
cat("  abs_rho   (Bacteria vs Fungi): p =",
    signif(wilcox.test(abs_rho   ~ pretty_group, data = class_df_ranks)$p.value, 3), "\n")

# Quadrant proportions by rank — does the dominant quadrant shift across ranks?
quad_by_rank <- class_df_ranks %>%
  count(rank_lab, classification) %>%
  group_by(rank_lab) %>%
  mutate(prop = round(n / sum(n), 3)) %>%
  ungroup()
cat("\nClassification proportions by rank (env_cycl):\n")
print(as.data.frame(quad_by_rank))

# ── Figure 2: Example hindcast panels for each quadrant ──────────────────────
cat("\nGenerating example hindcast panels for each classification quadrant...\n")

# Data-driven exemplars: one most-extreme taxon per occupied quadrant. Empty
# quadrants (e.g. non_seasonal/high-persistence under threshold 0.1) are skipped
# entirely so panels can't end up with NA labels for missing taxa.
quadrant_pretty <- c(
  seasonal_high_persistence     = "Seasonal,\nhigh persistence",
  seasonal_low_persistence      = "Seasonal,\nlow persistence",
  non_seasonal_high_persistence = "Non-seasonal,\nhigh persistence",
  non_seasonal_low_persistence  = "Non-seasonal,\nlow persistence"
)
exemplars <- class_df %>%
  group_by(classification) %>%
  mutate(extreme_score = case_when(
    classification == "seasonal_high_persistence"     ~ amplitude + abs_rho,
    classification == "seasonal_low_persistence"     ~ amplitude + (1 - abs_rho),
    classification == "non_seasonal_high_persistence" ~ (1 - amplitude) + abs_rho,
    classification == "non_seasonal_low_persistence" ~ (1 - amplitude) + (1 - abs_rho),
    TRUE ~ 0
  )) %>%
  slice_max(extreme_score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(species = taxon, pretty_group, classification,
            rho_mean, amplitude,
            quadrant = unname(quadrant_pretty[classification]))

cat("Exemplar taxa per occupied quadrant:\n")
print(exemplars %>% select(classification, species, pretty_group,
                           rho_mean, amplitude))

# For each exemplar, find an observed-site plot with the most truth values
# (better visualization than picking an arbitrary plot).
cat("Loading hindcast data from per-site files...\n")
hindcast_dir <- here("data/hindcasts/driver_uncertainty")

pick_best_plot <- function(species) {
  pat <- paste0("^hindcasts_env_cycl_", species,
                "_20130601_20180101_with_legacy_covariate_.*_observed\\.rds$")
  files <- list.files(hindcast_dir, pattern = pat, full.names = TRUE)
  if (length(files) == 0) return(NULL)
  best_df <- NULL; best_plot <- NA_character_; best_n <- 0L
  for (f in files) {
    df <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(df) || !"truth" %in% names(df)) next
    counts <- df %>% filter(!is.na(truth)) %>% count(plotID)
    if (nrow(counts) == 0) next
    top <- counts %>% slice_max(n, n = 1, with_ties = FALSE)
    if (top$n > best_n) {
      best_n   <- top$n
      best_plot <- top$plotID
      best_df  <- df %>% filter(plotID == top$plotID)
    }
  }
  if (is.null(best_df)) return(NULL)
  list(plotID = best_plot, hindcast = best_df, n_obs = best_n)
}

hind_list <- list()
for (i in seq_len(nrow(exemplars))) {
  sp <- exemplars$species[i]
  picked <- pick_best_plot(sp)
  if (is.null(picked)) {
    cat("  No hindcast files found for", sp, "- skipping panel\n")
    next
  }
  cat("  ", sp, "@", picked$plotID,
      "(", picked$n_obs, "observations)\n")
  hind_list[[length(hind_list) + 1L]] <- picked$hindcast %>%
    mutate(panel_species = sp, panel_plotID = picked$plotID)
}
if (length(hind_list) == 0) stop("No exemplar hindcast files found.")

hind_sub <- bind_rows(hind_list) %>%
  inner_join(exemplars %>% select(species, quadrant),
             by = c("panel_species" = "species"))

# Hindcast files may carry NA pretty_group on some rows; rebuild from class_df
# using the canonical (non-NA) per-taxon mapping.
pg_lookup <- class_df %>%
  filter(!is.na(pretty_group)) %>%
  select(taxon, pretty_group) %>% distinct()
hind_sub <- hind_sub %>%
  select(-pretty_group) %>%
  left_join(pg_lookup, by = c("panel_species" = "taxon"))

# Trim to start at first observation per panel
first_obs <- hind_sub %>%
  filter(!is.na(truth)) %>%
  group_by(panel_species, panel_plotID) %>%
  summarise(first_date = min(dates), .groups = "drop")

hind_sub <- hind_sub %>%
  left_join(first_obs, by = c("panel_species", "panel_plotID")) %>%
  filter(dates >= first_date) %>%
  select(-first_date)

# Attach the rho/amp values used in the panel header
hind_sub <- hind_sub %>%
  left_join(class_df %>% select(taxon, rho_mean, amplitude),
            by = c("panel_species" = "taxon"))

# Panel label with taxon name, rho, amplitude. Skip the kingdom in parens
# when pretty_group is NA (which only happens if class_df is missing it).
hind_sub <- hind_sub %>%
  mutate(
    kingdom_tag = ifelse(is.na(pretty_group), "",
                         paste0(" (", pretty_group, ")")),
    panel_label = paste0(quadrant, "\n",
                         panel_species, kingdom_tag, " @ ", panel_plotID, "\n",
                         "rho=", round(rho_mean, 2),
                         ", amp=", round(amplitude, 2))
  )

fig_quadrants <- ggplot(hind_sub, aes(x = dates)) +
  geom_ribbon(data = ~ filter(.x, fcast_period == "calibration"),
              aes(ymin = lo, ymax = hi, fill = pretty_group), alpha = 0.25) +
  geom_ribbon(data = ~ filter(.x, fcast_period == "calibration"),
              aes(ymin = lo_25, ymax = hi_75, fill = pretty_group), alpha = 0.4) +
  geom_ribbon(data = ~ filter(.x, fcast_period == "hindcast"),
              aes(ymin = lo, ymax = hi, fill = pretty_group), alpha = 0.12) +
  geom_ribbon(data = ~ filter(.x, fcast_period == "hindcast"),
              aes(ymin = lo_25, ymax = hi_75, fill = pretty_group), alpha = 0.2) +
  geom_line(data = ~ filter(.x, fcast_period == "calibration"),
            aes(y = med, color = pretty_group), linewidth = 0.7) +
  geom_line(data = ~ filter(.x, fcast_period == "hindcast"),
            aes(y = med, color = pretty_group), linewidth = 0.5, alpha = 0.6) +
  geom_point(aes(y = truth), color = "black", size = 1.3, alpha = 0.8) +
  facet_wrap(~ panel_label, scales = "free", ncol = 2) +
  scale_fill_manual(values = kingdom_colors) +
  scale_color_manual(values = kingdom_colors) +
  labs(x = "Date", y = "Relative abundance") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8, lineheight = 1.1),
    panel.spacing = unit(0.5, "cm"),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(file.path(fig_out_dir, "fig_classification_examples.png"), fig_quadrants,
       width = 9, height = 7, dpi = 200)
cat("Saved: figures/fig_classification_examples.png\n")
