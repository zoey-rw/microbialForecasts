# Fig: What predicts forecast accuracy?
# Multi-panel forest plot showing standardized effects of organism traits
# and environmental sensitivities on forecast error (nRMSE),
# with separate models for fungi and bacteria.

source("source.R")
pacman::p_load(ggplot2, dplyr, tidyr, cowplot, lme4)
options(na.action = "na.omit")

# =============================================================================
# 1. Load and prepare data
# =============================================================================

scores_list <- readRDS(here("data/summary/scoring_metrics_plsr2.rds"))
converged <- unique(scores_list$scoring_metrics$model_id)
converged_base <- gsub("_(combined|beta_regression)$", "", converged)

# Scoring metrics (restrict to observed-site hindcasts)
scores_df <- scores_list$scoring_metrics %>%
  filter(model_id %in% converged,
         site_prediction == "New time (observed site)") %>%
  select(model_id, species, pretty_group, model_name, RMSE.norm, RSQ, mean_crps_sample)

# Temporal autocorrelation (rho) and core-level precision
rho_core <- readRDS(here("data/summary/rho_core_sd_effects.rds")) %>%
  filter(time_period == "20130601_20180101",
         model_name != "all_covariates") %>%
  mutate(model_id = gsub("_(combined|beta_regression)$", "", model_id)) %>%
  filter(model_id %in% converged_base) %>%
  select(model_id, taxon, rowname, Mean) %>%
  pivot_wider(names_from = "rowname", values_from = "Mean", values_fn = mean)

# Spatial autocorrelation (Moran's I)
moran_df <- readRDS("data/clean/moran_stat.rds") %>%
  select(taxon, mean_morans)

# Spatial heterogeneity (CV across sites)
cv_data <- scores_list$cv_metric_scaled %>%
  filter(model_id %in% converged,
         cv_type == "mean_per_site_cv",
         !is.na(cv)) %>%
  select(model_id, species, cv) %>%
  distinct()

# Seasonal amplitude
seas_amplitude <- readRDS(here("data/summary/seasonal_amplitude.rds"))[[6]] %>%
  filter(time_period == "20130601_20180101") %>%
  mutate(model_id = gsub("_(combined|beta_regression)$", "", model_id)) %>%
  filter(model_id %in% converged_base) %>%
  select(model_id, taxon, model_name, amplitude)

# Environmental predictor effect sizes
beta_wide <- readRDS(here("data/summary/predictor_effects.rds")) %>%
  filter(time_period == "20130601_20180101") %>%
  mutate(model_id = gsub("_(combined|beta_regression)$", "", model_id)) %>%
  filter(model_id %in% converged_base) %>%
  select(model_id, beta, effSize) %>%
  group_by(model_id, beta) %>%
  summarise(effSize = mean(effSize, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = "beta", values_from = "effSize")

# =============================================================================
# 2. Merge into a single analysis data frame
# =============================================================================

# Start from scores, merge rho/precision by species = taxon
df <- scores_df %>%
  left_join(rho_core, by = c("model_id", "species" = "taxon")) %>%
  left_join(moran_df, by = c("species" = "taxon")) %>%
  left_join(seas_amplitude %>% select(-model_name),
            by = c("model_id", "species" = "taxon")) %>%
  left_join(beta_wide, by = "model_id") %>%
  left_join(cv_data, by = c("model_id", "species")) %>%
  filter(!is.na(RMSE.norm), !is.na(pretty_group))

# Unlist any list-columns from the beta pivot
list_cols <- names(df)[sapply(df, is.list)]
for (col in list_cols) {
  df[[col]] <- as.numeric(sapply(df[[col]], function(x) if (is.list(x)) x[[1]] else x))
}

cat("Merged dataset:", nrow(df), "rows,", length(unique(df$species)), "taxa\n")
cat("Group counts:", paste(names(table(df$pretty_group)), table(df$pretty_group), sep = "=", collapse = ", "), "\n")
cat("Moran's I available:", sum(!is.na(df$mean_morans)), "/", nrow(df), "\n")
cat("CV available:", sum(!is.na(df$cv)), "/", nrow(df), "\n")
cat("Precision available:", sum(!is.na(df$precision)), "/", nrow(df), "\n")

# =============================================================================
# 3. Scale predictors within group x model_name
# =============================================================================

# Rename columns for cleaner plotting
df_scaled <- df %>%
  group_by(pretty_group, model_name) %>%
  mutate(
    nRMSE              = as.numeric(scale(RMSE.norm)),
    `Temporal memory`  = as.numeric(scale(rho)),
    `Seasonal amplitude` = as.numeric(scale(amplitude)),
    `Core variability` = as.numeric(scale(1 / precision)),
    `Temporal variability` = as.numeric(scale(cv)),
    `Spatial autocorrelation` = as.numeric(scale(mean_morans)),
    Temperature        = as.numeric(scale(Temperature)),
    Moisture           = as.numeric(scale(Moisture)),
    pH                 = as.numeric(scale(pH)),
    `Percent C`        = as.numeric(scale(pC)),
    `EcM trees`        = as.numeric(scale(`Ectomycorrhizal\ntrees`)),
    LAI                = as.numeric(scale(LAI))
  ) %>%
  ungroup()

# =============================================================================
# 4. Fit models: separate for fungi and bacteria
# =============================================================================

# Predictor formula (organism traits + environmental sensitivities)
# Using species as random intercept to account for multiple model_names per taxon
base_formula <- nRMSE ~
  `Temporal memory` + `Seasonal amplitude` + `Core variability` +
  `Temporal variability` + `Spatial autocorrelation` +
  Temperature + Moisture + pH + `Percent C` + `EcM trees` + LAI

# Try mixed model; fall back to lm if insufficient grouping
fit_model <- function(data, group_label) {
  n_species <- length(unique(data$species))
  cat(group_label, ":", nrow(data), "obs,", n_species, "taxa\n")

  if (n_species >= 5) {
    # Mixed model with random intercept per species
    tryCatch({
      mod <- lmer(update(base_formula, . ~ . + (1 | species)), data = data,
                  control = lmerControl(optimizer = "bobyqa"))
      cat("  -> mixed model (lmer) fit successfully\n")
      return(mod)
    }, error = function(e) {
      cat("  -> lmer failed:", e$message, "; falling back to lm\n")
    })
  }

  # Fallback: plain OLS
  mod <- lm(base_formula, data = data)
  cat("  -> lm fit successfully\n")
  return(mod)
}

fungi_mod <- fit_model(df_scaled %>% filter(pretty_group == "Fungi"), "Fungi")
bact_mod  <- fit_model(df_scaled %>% filter(pretty_group == "Bacteria"), "Bacteria")

# =============================================================================
# 5. Extract coefficients into a tidy data frame
# =============================================================================

extract_coefs <- function(mod, group_label) {
  if (inherits(mod, "lmerMod")) {
    cc <- summary(mod)$coefficients
    ci <- confint(mod, method = "Wald", parm = "beta_")
    est <- data.frame(
      term      = rownames(cc),
      estimate  = cc[, "Estimate"],
      se        = cc[, "Std. Error"],
      ci_lo     = ci[, 1],
      ci_hi     = ci[, 2],
      pvalue    = 2 * pnorm(-abs(cc[, "t value"])),  # approx p from t
      stringsAsFactors = FALSE
    )
  } else {
    cc <- summary(mod)$coefficients
    ci <- confint(mod)
    est <- data.frame(
      term     = rownames(cc),
      estimate = cc[, "Estimate"],
      se       = cc[, "Std. Error"],
      ci_lo    = ci[, 1],
      ci_hi    = ci[, 2],
      pvalue   = cc[, "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
  }
  est$group <- group_label
  est <- est %>%
    filter(term != "(Intercept)") %>%
    mutate(term = gsub("^`|`$", "", term))  # strip backticks from coefficient names
  return(est)
}

coef_df <- bind_rows(
  extract_coefs(fungi_mod, "Fungi"),
  extract_coefs(bact_mod, "Bacteria")
)

# Assign predictor categories for visual grouping
trait_terms <- c("Temporal memory", "Seasonal amplitude",
                 "Core variability", "Temporal variability",
                 "Spatial autocorrelation")
env_terms   <- c("Temperature", "Moisture", "pH", "Percent C",
                  "EcM trees", "LAI")

coef_df$category <- ifelse(coef_df$term %in% trait_terms,
                           "Organism traits",
                           "Environmental sensitivities")

# Significance indicators
coef_df$sig_label <- case_when(
  coef_df$pvalue < 0.001 ~ "***",
  coef_df$pvalue < 0.01  ~ "**",
  coef_df$pvalue < 0.05  ~ "*",
  TRUE                   ~ ""
)

# Order terms: traits on top, env below; within each, order by mean |estimate|
term_order <- coef_df %>%
  group_by(term) %>%
  summarise(mean_abs = mean(abs(estimate), na.rm = TRUE), .groups = "drop") %>%
  left_join(distinct(coef_df, term, category), by = "term") %>%
  arrange(factor(category, levels = c("Environmental sensitivities",
                                      "Organism traits")),
          mean_abs) %>%
  pull(term)

coef_df$term <- factor(coef_df$term, levels = term_order)

# =============================================================================
# 6. Plot: overlaid forest plot with category grouping
# =============================================================================

# Color palette
group_colors <- c("Fungi" = "#0072B2", "Bacteria" = "#E69F00")

# Dodge width for overlaying fungi vs bacteria
pd <- position_dodge(width = 0.6)

# Create a significance column for mapping
coef_df$significant <- coef_df$pvalue < 0.05

p_forest <- ggplot(coef_df, aes(x = term, y = estimate, group = group, color = group)) +
  # Zero reference line
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50", linetype = "dashed") +
  # 95% CI segments
  geom_linerange(aes(ymin = ci_lo, ymax = ci_hi), linewidth = 0.6,
                 position = pd, alpha = 0.8) +
  # Point estimates: filled = p<0.05, open = n.s.
  geom_point(aes(shape = significant), size = 3, position = pd,
             stroke = 0.8) +
  # Significance stars (offset to right of CI)
  geom_text(aes(y = pmax(ci_hi, 0) + 0.08, label = sig_label),
            position = pd, size = 3, show.legend = FALSE, vjust = 0.5,
            fontface = "bold") +
  # Category separation
  facet_grid(category ~ ., scales = "free_y", space = "free_y",
             switch = "y") +
  coord_flip() +
  scale_color_manual(values = group_colors, name = NULL) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "p < 0.05", "FALSE" = "n.s."),
                     name = NULL) +
  labs(x = NULL,
       y = "Standardized effect on forecast error (nRMSE)") +
  theme_bw(base_size = 13) +
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey93", color = NA),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 10,
                                     margin = margin(r = 5)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.8, "lines"),
    legend.position = "top",
    legend.text = element_text(size = 11),
    legend.spacing.x = unit(0.5, "cm"),
    legend.margin = margin(b = -5),
    plot.margin = margin(8, 12, 5, 5),
    axis.text.y = element_text(size = 11)
  ) +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2, override.aes = list(size = 2.5)))

# Sample sizes for caption
n_fungi <- sum(df_scaled$pretty_group == "Fungi" & !is.na(df_scaled$nRMSE))
n_bact  <- sum(df_scaled$pretty_group == "Bacteria" & !is.na(df_scaled$nRMSE))
cat("Sample sizes  --  Fungi:", n_fungi, "  Bacteria:", n_bact, "\n")

# =============================================================================
# 7. Supplementary panel: scatter of top predictors vs nRMSE
# =============================================================================

# Identify the 4 predictors with largest absolute pooled estimates
top4 <- coef_df %>%
  group_by(term) %>%
  summarise(mean_abs = mean(abs(estimate)), .groups = "drop") %>%
  slice_max(mean_abs, n = 4) %>%
  pull(term) %>% as.character()

# Build scatter data manually to avoid any_of() issues with special characters
scatter_rows <- list()
for (pred in top4) {
  if (pred %in% colnames(df_scaled)) {
    scatter_rows[[pred]] <- df_scaled %>%
      ungroup() %>%
      transmute(species, pretty_group, nRMSE,
                predictor = pred, value = .data[[pred]]) %>%
      filter(!is.na(value), !is.na(nRMSE))
  }
}
scatter_df <- bind_rows(scatter_rows)
cat("Scatter data:", nrow(scatter_df), "rows for", length(scatter_rows), "predictors\n")

p_scatter <- ggplot(scatter_df, aes(x = value, y = nRMSE, color = pretty_group)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  facet_wrap(~ predictor, scales = "free_x", nrow = 1) +
  scale_color_manual(values = group_colors) +
  labs(x = "Scaled predictor value", y = "Scaled nRMSE", color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey92", color = NA),
        panel.grid.minor = element_blank())

# =============================================================================
# 8. Combine and save
# =============================================================================

combined <- plot_grid(
  p_forest, p_scatter,
  ncol = 1, rel_heights = c(3, 1.2),
  labels = c("A", "B"), label_size = 16
)

out_dir <- here("data", "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "fig_patchiness_predictability.png"), combined,
       width = 9, height = 10, dpi = 200)

cat("Saved: data/figures/fig_patchiness_predictability.png\n")

# =============================================================================
# 9. Print model summaries
# =============================================================================

cat("\n--- Fungi model ---\n")
print(summary(fungi_mod))
cat("\n--- Bacteria model ---\n")
print(summary(bact_mod))

cat("\n=== DONE ===\n")
