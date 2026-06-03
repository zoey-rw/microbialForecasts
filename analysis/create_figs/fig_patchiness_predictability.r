# Fig 5: What predicts per-taxon forecast skill?
# Forest plot of standardized regression coefficients from a model of forecast
# R-squared on organism traits and predictor sensitivities. Uses only env_cycl
# (seasonality + environment) models that converged at Rhat<1.1. Temperature
# sensitivity and seasonal amplitude are grouped as "Seasonal forcing" to show
# the temperature-explained and residual-seasonality decomposition.

source("source.R")
pacman::p_load(ggplot2, dplyr, tidyr, cowplot, lme4)
options(na.action = "na.omit")

# =============================================================================
# 1. Load and prepare data
# =============================================================================

scores_list <- readRDS(here("data/summary/scoring_metrics_plsr2.rds"))
# converged_list uses "_beta_regression" suffix; scoring_metrics does not.
converged <- unique(gsub("_(combined|beta_regression)$", "",
                         scores_list$converged_list))
cat("Converged models (Rhat<1.1):", length(converged), "\n")

# Restrict to env_cycl models at the observed-site time horizon.
scores_df <- scores_list$scoring_metrics %>%
  filter(model_id %in% converged,
         site_prediction == "New time (observed site)",
         model_name == "env_cycl") %>%
  select(model_id, species, pretty_group, model_name, RMSE.norm, RSQ,
         mean_crps_sample)

# Temporal autocorrelation (rho) and core-level precision
rho_core <- readRDS(here("data/summary/rho_core_sd_effects.rds")) %>%
  filter(time_period == "2013-06_2018-01",
         model_name == "env_cycl") %>%
  mutate(model_id = gsub("_(combined|beta_regression)$", "", model_id)) %>%
  filter(model_id %in% converged) %>%
  select(model_id, taxon, rowname, Mean) %>%
  pivot_wider(names_from = "rowname", values_from = "Mean", values_fn = mean)

# Spatial autocorrelation in raw relative abundances (Moran's I)
moran_df <- readRDS("data/clean/moran_stat.rds") %>%
  select(taxon, mean_morans)

# Month-to-month CV (averaged across sites within a plot)
cv_data <- scores_list$cv_metric_scaled %>%
  filter(model_id %in% converged,
         cv_type == "mean_per_site_cv",
         !is.na(cv)) %>%
  select(model_id, species, cv) %>%
  distinct()

# Seasonal amplitude from the sin/cos terms of the env_cycl fits
seas_amplitude <- readRDS(here("data/summary/seasonal_amplitude.rds"))[[6]] %>%
  filter(time_period == "2013-06_2018-01",
         model_name == "env_cycl") %>%
  mutate(model_id = gsub("_(combined|beta_regression)$", "", model_id)) %>%
  filter(model_id %in% converged) %>%
  select(model_id, taxon, amplitude)

# Standardized effect sizes of each environmental predictor on relative
# abundance (the predictor "sensitivity" measures)
beta_wide <- readRDS(here("data/summary/predictor_effects.rds")) %>%
  filter(time_period == "2013-06_2018-01") %>%
  mutate(model_id = gsub("_(combined|beta_regression)$", "", model_id)) %>%
  filter(model_id %in% converged) %>%
  select(model_id, beta, effSize) %>%
  group_by(model_id, beta) %>%
  summarise(effSize = mean(effSize, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = "beta", values_from = "effSize")

# =============================================================================
# 2. Merge into one analysis data frame
# =============================================================================

df <- scores_df %>%
  left_join(rho_core,   by = c("model_id", "species" = "taxon")) %>%
  left_join(moran_df,   by = c("species" = "taxon")) %>%
  left_join(seas_amplitude, by = c("model_id", "species" = "taxon")) %>%
  left_join(beta_wide,  by = "model_id") %>%
  left_join(cv_data,    by = c("model_id", "species")) %>%
  filter(!is.na(RSQ), !is.na(pretty_group))

# Some columns can come through as list columns from the pivot — flatten them.
list_cols <- names(df)[sapply(df, is.list)]
for (col in list_cols) {
  df[[col]] <- as.numeric(sapply(df[[col]], function(x) if (is.list(x)) x[[1]] else x))
}

cat("env_cycl converged taxa:", nrow(df), "rows,",
    length(unique(df$species)), "taxa\n")
cat("Bacteria:", sum(df$pretty_group == "Bacteria"),
    " Fungi:", sum(df$pretty_group == "Fungi"), "\n")

# =============================================================================
# 3. Scale predictors within group
# =============================================================================

# Standardize predictors globally across both kingdoms (not within kingdom),
# so coefficients in the pooled regression are comparable.
df_scaled <- df %>%
  mutate(
    RSQ_scaled                  = as.numeric(scale(RSQ)),
    `Temporal memory`           = as.numeric(scale(rho)),
    `Seasonal amplitude`        = as.numeric(scale(amplitude)),
    `Core variability`          = as.numeric(scale(1 / precision)),
    `Temporal variability`      = as.numeric(scale(cv)),
    `Spatial autocorrelation`   = as.numeric(scale(mean_morans)),
    Temperature                 = as.numeric(scale(Temperature)),
    Moisture                    = as.numeric(scale(Moisture)),
    pH                          = as.numeric(scale(pH)),
    `Percent C`                 = as.numeric(scale(pC)),
    `EM trees`                 = as.numeric(scale(`Ectomycorrhizal\ntrees`)),
    LAI                         = as.numeric(scale(LAI)),
    Kingdom                     = factor(pretty_group, levels = c("Bacteria","Fungi"))
  )

# =============================================================================
# 4. Fit one pooled model with kingdom as a covariate
# =============================================================================
# Pooling bacteria + fungi (n = ~148) and controlling for kingdom-level
# intercept differences. Kingdom coefficient = mean shift in scaled forecast
# R^2 for fungi relative to bacteria.

base_formula <- RSQ_scaled ~ Kingdom +
  `Temporal memory` + `Seasonal amplitude` + `Core variability` +
  `Temporal variability` + `Spatial autocorrelation` +
  Temperature + Moisture + pH + `Percent C` + `EM trees` + LAI

pooled_mod <- lm(base_formula, data = df_scaled)
cat("Pooled model: n =", nobs(pooled_mod),
    "  R^2 =", round(summary(pooled_mod)$r.squared, 3), "\n")

# =============================================================================
# 5. Extract coefficients into a tidy data frame
# =============================================================================

extract_coefs <- function(mod) {
  cc <- summary(mod)$coefficients
  ci <- confint(mod)
  shared <- intersect(rownames(cc), rownames(ci))
  est <- data.frame(
    term     = shared,
    estimate = cc[shared, "Estimate"],
    se       = cc[shared, "Std. Error"],
    ci_lo    = ci[shared, 1],
    ci_hi    = ci[shared, 2],
    pvalue   = cc[shared, "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
  est %>%
    filter(term != "(Intercept)") %>%
    mutate(term = gsub("^`|`$", "", term),
           term = gsub("^KingdomFungi$", "Fungi vs. Bacteria", term))
}

coef_df <- extract_coefs(pooled_mod)

# =============================================================================
# 6. Assign categories
# =============================================================================
# Three groups so the temperature / residual-seasonality decomposition is
# visually obvious: the two paired seasonal predictors sit together.

# Two top-level categories:
#   - "Organism traits" = variables that exist independent of the forecasting
#     model (group identity + raw-data summaries).
#   - "Model-estimated\nparameters" = everything that comes out of the env_cycl
#     model fit. Within this category, the paired Temperature / Seasonal
#     amplitude predictors will appear together at the top because they have
#     the largest absolute effect sizes.
trait_terms <- c("Fungi vs. Bacteria",
                 "Spatial autocorrelation", "Temporal variability")
model_terms <- c("Temperature", "Seasonal amplitude",
                 "Temporal memory", "Core variability",
                 "Moisture", "pH", "Percent C", "EM trees", "LAI")
kingdom_term <- "Fungi vs. Bacteria"

coef_df$category <- case_when(
  coef_df$term %in% trait_terms ~ "Organism traits",
  coef_df$term %in% model_terms ~ "Model-estimated\nparameters",
  TRUE                          ~ NA_character_
)

coef_df$sig_label <- case_when(
  coef_df$pvalue < 0.001 ~ "***",
  coef_df$pvalue < 0.01  ~ "**",
  coef_df$pvalue < 0.05  ~ "*",
  TRUE                   ~ ""
)

# Order: within each category, order by mean absolute estimate (smaller at
# bottom, larger at top). Across categories, top-to-bottom = Seasonal forcing,
# Other environmental sensitivities, Organism traits.
category_order <- c("Organism traits", "Model-estimated\nparameters")

term_order <- coef_df %>%
  group_by(term) %>%
  summarise(mean_abs = mean(abs(estimate), na.rm = TRUE), .groups = "drop") %>%
  left_join(distinct(coef_df, term, category), by = "term") %>%
  arrange(factor(category, levels = category_order), mean_abs) %>%
  pull(term)

coef_df$term     <- factor(coef_df$term, levels = term_order)
coef_df$category <- factor(coef_df$category, levels = category_order)

# =============================================================================
# 7. Forest plot
# =============================================================================

point_color <- "#444444"  # neutral dark grey for pooled coefficients
group_colors <- c("Fungi" = "#0072B2", "Bacteria" = "#E69F00")  # for panel B

coef_df$significant <- coef_df$pvalue < 0.05

p_forest <- ggplot(coef_df, aes(x = term, y = estimate)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50",
             linetype = "dashed") +
  geom_linerange(aes(ymin = ci_lo, ymax = ci_hi), linewidth = 0.7,
                 color = point_color, alpha = 0.85) +
  geom_point(aes(shape = significant), size = 3.2, stroke = 0.9,
             color = point_color) +
  geom_text(aes(y = pmax(ci_hi, 0) + 0.08, label = sig_label),
            size = 3.5, show.legend = FALSE, vjust = 0.5, fontface = "bold",
            color = point_color) +
  facet_grid(category ~ ., scales = "free_y", space = "free_y",
             switch = "y") +
  coord_flip() +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "p < 0.05", "FALSE" = "n.s."),
                     name = NULL) +
  labs(x = NULL,
       y = expression("Standardized effect on forecast R"^2)) +
  theme_bw(base_size = 13) +
  theme(
    strip.placement   = "outside",
    strip.background  = element_rect(fill = "grey93", color = NA),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 10,
                                     margin = margin(r = 5)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.spacing      = unit(0.8, "lines"),
    legend.position    = "top",
    legend.text        = element_text(size = 11),
    legend.margin      = margin(b = -5),
    plot.margin        = margin(8, 12, 5, 5),
    axis.text.y        = element_text(size = 11)
  ) +
  guides(shape = guide_legend(override.aes = list(size = 2.8)))

cat("Pooled n =", nobs(pooled_mod),
    " Pooled model R^2 =", round(summary(pooled_mod)$r.squared, 3), "\n")

# =============================================================================
# 8. Scatter panel for the four largest-effect predictors
# =============================================================================

top4 <- coef_df %>%
  filter(term != kingdom_term) %>%   # Kingdom is a factor — not a scatter axis
  group_by(term) %>%
  summarise(mean_abs = mean(abs(estimate)), .groups = "drop") %>%
  slice_max(mean_abs, n = 4) %>%
  pull(term) %>% as.character()

scatter_rows <- list()
for (pred in top4) {
  if (pred %in% colnames(df_scaled)) {
    scatter_rows[[pred]] <- df_scaled %>%
      ungroup() %>%
      transmute(species, pretty_group, RSQ_scaled,
                predictor = pred, value = .data[[pred]]) %>%
      filter(!is.na(value), !is.na(RSQ_scaled))
  }
}
scatter_df <- bind_rows(scatter_rows)
scatter_df$predictor <- factor(scatter_df$predictor, levels = top4)

p_scatter <- ggplot(scatter_df,
                    aes(x = value, y = RSQ_scaled)) +
  geom_point(alpha = 0.4, size = 1, color = point_color) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8,
              color = point_color, fill = "grey60") +
  facet_wrap(~ predictor, scales = "free_x", nrow = 1) +
  labs(x = "Scaled predictor value",
       y = expression("Scaled forecast R"^2)) +
  theme_bw(base_size = 12) +
  theme(legend.position    = "none",
        strip.background   = element_rect(fill = "grey92", color = NA),
        panel.grid.minor   = element_blank())

# =============================================================================
# 9. Combine and save
# =============================================================================

combined <- plot_grid(
  p_forest, p_scatter,
  ncol = 1, rel_heights = c(3, 1.2),
  labels = c("A", "B"), label_size = 16
)

out_dir <- here("figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "fig_patchiness_predictability.png"), combined,
       width = 9, height = 10, dpi = 200)
cat("Saved: figures/fig_patchiness_predictability.png\n")

# =============================================================================
# 10. Print model summaries for sanity
# =============================================================================

cat("\n--- Pooled model ---\n")
print(summary(pooled_mod))

cat("\n=== DONE ===\n")
