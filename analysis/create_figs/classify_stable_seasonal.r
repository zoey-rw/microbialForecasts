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
  group_by(model_id) %>%
  slice(1L) %>%
  ungroup()

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
  group_by(model_id_base) %>%
  slice(1L) %>%
  ungroup() %>%
  rename(model_id = model_id_base)

# Merge rho and seasonal amplitude on model_id
class_df <- rho_env_cycl %>%
  inner_join(seas_env_cycl, by = "model_id")

# Thresholds
AMPLITUDE_THRESHOLD <- 0.05
PERSISTENCE_THRESHOLD <- 0.3

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

# ── Figure 1: Scatter plot of rho vs amplitude with quadrant thresholds ──────
library(ggplot2)
library(cowplot)
library(ggrepel)

kingdom_colors <- c("Bacteria" = "#4DAF4A", "Fungi" = "#FF7F00")

fig_out_dir <- here("data", "figures")
if (!dir.exists(fig_out_dir)) dir.create(fig_out_dir, recursive = TRUE)

# Label the most extreme point in each quadrant
class_df <- class_df %>%
  group_by(classification) %>%
  mutate(is_extreme = case_when(
    classification == "seasonal_high_persistence" ~ amplitude + abs_rho,
    classification == "seasonal_low_persistence" ~ amplitude + (1 - abs_rho),
    classification == "non_seasonal_high_persistence" ~ (1 - amplitude) + abs_rho,
    classification == "non_seasonal_low_persistence" ~ (1 - amplitude) + (1 - abs_rho),
    TRUE ~ 0
  )) %>%
  mutate(label_taxon = ifelse(is_extreme == max(is_extreme), taxon, NA_character_)) %>%
  ungroup() %>%
  select(-is_extreme)

# Quadrant count labels
quad_labels <- summary_counts %>%
  mutate(
    x = ifelse(grepl("non_seasonal", classification), 0.01, max(class_df$amplitude, na.rm = TRUE) * 0.9),
    y = ifelse(grepl("high", classification), max(class_df$abs_rho, na.rm = TRUE) * 0.95, 0.01),
    label = paste0("n=", n, " (", round(prop * 100), "%)")
  )

fig_scatter <- ggplot(class_df, aes(x = amplitude, y = abs_rho)) +
  geom_point(aes(color = pretty_group), size = 2, alpha = 0.6) +
  geom_text_repel(aes(label = label_taxon), size = 2.5, max.overlaps = 10,
                  na.rm = TRUE, segment.color = "gray60", segment.size = 0.3) +
  geom_hline(yintercept = PERSISTENCE_THRESHOLD, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = AMPLITUDE_THRESHOLD, linetype = "dashed", color = "gray40") +
  geom_text(data = quad_labels, aes(x = x, y = y, label = label),
            size = 3, color = "gray30", hjust = 0, vjust = 0) +
  scale_color_manual(values = kingdom_colors, name = "Kingdom") +
  annotate("text", x = max(class_df$amplitude, na.rm = TRUE) * 0.5,
           y = max(class_df$abs_rho, na.rm = TRUE) * 0.97,
           label = "High persistence", size = 3, fontface = "italic", color = "gray50") +
  annotate("text", x = max(class_df$amplitude, na.rm = TRUE) * 0.5,
           y = 0.01,
           label = "Low persistence", size = 3, fontface = "italic", color = "gray50") +
  labs(x = "Seasonal amplitude", y = expression("Temporal persistence (" * rho * ")")) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.85, 0.85),
        legend.background = element_rect(fill = "white", color = "gray80"),
        panel.grid.minor = element_blank())

ggsave(file.path(fig_out_dir, "fig_persistence_seasonality_scatter.png"), fig_scatter,
       width = 7, height = 5.5, dpi = 200)
cat("Saved: data/figures/fig_persistence_seasonality_scatter.png\n")

# ── Figure 2: Example hindcast panels for each quadrant ──────────────────────
cat("\nGenerating example hindcast panels for each classification quadrant...\n")

# Representative taxa (most extreme in each quadrant) and best-calibrated plots
exemplars <- data.frame(
  species   = c("ectomycorrhizal", "wd2101.soil.group", "chaetosphaeriales", "orbiliomycetes"),
  plotID    = c("OSBS_001",        "SCBI_006",          "ORNL_049",          "DSNY_041"),
  quadrant  = c("Seasonal,\nlow persistence",
                "Non-seasonal,\nlow persistence",
                "Seasonal,\nhigh persistence",
                "Non-seasonal,\nhigh persistence"),
  stringsAsFactors = FALSE
)

# Load hindcast data from fresh per-site files (re-run by 06_hindcast_observed.r)
cat("Loading hindcast data from per-site files...\n")
hindcast_dir <- here("data/hindcasts/driver_uncertainty")

# Map exemplars to their site hindcast files
exemplars$siteID <- gsub("_.*", "", exemplars$plotID)
hind_list <- list()
for (i in seq_len(nrow(exemplars))) {
  fname <- file.path(hindcast_dir,
    paste0("hindcasts_env_cycl_", exemplars$species[i],
           "_20130601_20180101_with_legacy_covariate_",
           exemplars$siteID[i], "_observed.rds"))
  if (file.exists(fname)) {
    df <- readRDS(fname) %>%
      filter(plotID == exemplars$plotID[i])
    hind_list[[i]] <- df
    cat("  Loaded", basename(fname), ":", nrow(df), "rows\n")
  } else {
    # Fall back to the big summary file
    cat("  Not found:", basename(fname), "- falling back to all_hindcasts_plsr2.rds\n")
    if (!exists("hindcast_all")) hindcast_all <- readRDS(here("data/summary/all_hindcasts_plsr2.rds"))
    hind_list[[i]] <- hindcast_all %>%
      filter(species == exemplars$species[i], model_name == "env_cycl",
             plotID == exemplars$plotID[i])
  }
}
hind_sub <- bind_rows(hind_list) %>%
  inner_join(exemplars %>% select(species, plotID, quadrant), by = c("species", "plotID"))
if (exists("hindcast_all")) { rm(hindcast_all); gc() }

# Fill pretty_group from classification data (new hindcast files have NA for some rows)
pg_lookup <- class_df %>% select(taxon, pretty_group) %>% distinct()
hind_sub <- hind_sub %>%
  select(-pretty_group) %>%
  left_join(pg_lookup, by = c("species" = "taxon"))

# Trim to start at first observation per panel
first_obs <- hind_sub %>%
  filter(!is.na(truth)) %>%
  group_by(species, plotID) %>%
  summarise(first_date = min(dates), .groups = "drop")

hind_sub <- hind_sub %>%
  left_join(first_obs, by = c("species", "plotID")) %>%
  filter(dates >= first_date) %>%
  select(-first_date)

# Add rho and amplitude values for subtitle
hind_sub <- hind_sub %>%
  left_join(class_df %>% select(taxon, rho_mean, amplitude),
            by = c("species" = "taxon"))

# Panel label with taxon name, rho, amplitude
hind_sub <- hind_sub %>%
  mutate(
    panel_label = paste0(quadrant, "\n",
                         species, " (", pretty_group, ")\n",
                         "rho=", round(rho_mean, 2),
                         ", amp=", round(amplitude, 2)),
    # Order quadrants logically: seasonal first, then by persistence
    quadrant_f = factor(quadrant,
                        levels = c("Seasonal,\nhigh persistence",
                                   "Seasonal,\nlow persistence",
                                   "Non-seasonal,\nhigh persistence",
                                   "Non-seasonal,\nlow persistence"))
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
cat("Saved: data/figures/fig_classification_examples.png\n")
