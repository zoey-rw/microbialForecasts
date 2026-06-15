# Visualize effect size estimates (beta covariates) from all models
source("source.R")
library(stringr)
library(gridExtra)
library(ggpubr)
library(rstatix)
library(ggh4x)

# ── Shared constants ─────────────────────────────────────────────────────────
# kingdom_colors and fcast_type_colors come from source.R.
BASE_SIZE <- 14
func_tax_colors <- fcast_type_colors

base_theme <- theme_bw(base_size = BASE_SIZE) +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold", size = BASE_SIZE),
    axis.title       = element_text(size = BASE_SIZE),
    panel.grid.minor = element_blank()
  )

# Helper: significance label from p-value
sig_label_from_p <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ ""
  )
}

# ── Data loading ─────────────────────────────────────────────────────────────
converged        <- readRDS(here("data/summary/weak_converged_taxa_list.rds"))
converged_strict <- readRDS(here("data/summary/converged_taxa_list.rds"))
sum.all          <- readRDS(here("data/summary/predictor_effects.rds"))

seasonal_amplitude_in <- readRDS(here("data/summary/seasonal_amplitude.rds"))
cycl_only_vals_scores <- seasonal_amplitude_in[[6]] %>%
  filter(model_name == "cycl_only") %>%
  mutate(cycl_amplitude = amplitude) %>%
  select(-c(sin, cos, max, amplitude_orig, max)) %>%
  pivot_longer(cols = cycl_amplitude, values_to = "effSize", names_to = "beta")

env_cycl_vals_scores <- seasonal_amplitude_in[[6]] %>%
  filter(model_name == "env_cycl") %>%
  mutate(residual_amplitude = amplitude) %>%
  select(-c(sin, cos, max, amplitude_orig, max)) %>%
  pivot_longer(cols = residual_amplitude, values_to = "effSize", names_to = "beta")

# ── Combine and clean ────────────────────────────────────────────────────────
df_cal_fg_tax <- sum.all %>%
  filter(time_period == "2013-06_2018-01") %>%
  filter(!beta %in% c("sin", "cos"))

df_cal_fg_tax <- rbindlist(list(df_cal_fg_tax, env_cycl_vals_scores, cycl_only_vals_scores), fill = TRUE)

df_cal_fg_tax <- df_cal_fg_tax %>%
  filter(time_period != "2015-11_2018-01")

df_cal_fg_tax$fcast_type <- recode(df_cal_fg_tax$fcast_type,
  "functional" = "Functional", "taxon" = "Taxonomic",
  "Functional" = "Functional", "Taxonomic" = "Taxonomic"
)

# ── Facet labels (shared) ────────────────────────────────────────────────────
df_cal_fg_tax$beta_pretty <- recode(df_cal_fg_tax$beta,
  "residual_amplitude" = "Seasonality",
  "cycl_amplitude"     = "Seasonality",
  "pC"                 = "% Carbon",
  "LAI"                = "Leaf area\nindex",
  "Ectomycorrhizal\ntrees" = "EM trees"
)

df_cal_fg_tax$beta_pretty <- factor(df_cal_fg_tax$beta_pretty,
  levels = c("Seasonality", "EM trees", "Leaf area\nindex",
             "% Carbon", "pH", "Temperature", "Moisture")
)

beta.labs  <- setNames(levels(df_cal_fg_tax$beta_pretty), levels(df_cal_fg_tax$beta_pretty))
group.labs <- c("Functional" = "Functional", "Taxonomic" = "Taxonomic")
model.labs <- c("env_cycl" = "Environmental + Cyclic", "cycl_only" = "Cyclic Only", "env_cov" = "Environmental Only")
label_fn2  <- as_labeller(c(model.labs, beta.labs, group.labs))

# ── Fix significance for seasonality rows (define once) ──────────────────────
df_cal_fg_tax_fixed <- df_cal_fg_tax %>%
  mutate(significant = case_when(
    is.na(significant) & beta %in% c("residual_amplitude", "cycl_amplitude") ~
      as.numeric(significant_sin == 1 | significant_cos == 1),
    TRUE ~ significant
  ))

# ══════════════════════════════════════════════════════════════════════════════
# Bacteria vs Fungi effect sizes (env_cycl, faceted by fcast_type x predictor)
# ══════════════════════════════════════════════════════════════════════════════

# Statistical tests (env_cycl only, groups with sufficient data)
groups_with_sufficient_data <- df_cal_fg_tax %>%
  filter(model_name == "env_cycl") %>%
  group_by(fcast_type, model_name, beta_pretty, pretty_group) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(fcast_type, model_name, beta_pretty) %>%
  summarise(n_groups = n_distinct(pretty_group), min_n = min(n), .groups = "drop") %>%
  filter(n_groups >= 2, min_n >= 2)

df_cal_fg_tax_for_tests <- df_cal_fg_tax %>%
  filter(model_name == "env_cycl") %>%
  semi_join(groups_with_sufficient_data, by = c("fcast_type", "model_name", "beta_pretty"))

cat("Groups with sufficient data for statistical tests:", nrow(groups_with_sufficient_data), "\n")

stat.test <- tryCatch({
  if (nrow(df_cal_fg_tax_for_tests) > 0) {
    df_cal_fg_tax_for_tests %>%
      group_by(fcast_type, beta_pretty) %>%
      t_test(effSize ~ pretty_group) %>%
      adjust_pvalue(method = "fdr") %>%
      add_significance() %>%
      add_xy_position(x = "pretty_group")
  } else {
    data.frame()
  }
}, error = function(e) {
  cat("Warning: Statistical tests failed:", e$message, "\n")
  data.frame()
})

b_vs_f_fcast_type_plot <- ggplot(
  data = df_cal_fg_tax %>% filter(model_name == "env_cycl"),
  aes(x = pretty_group, color = pretty_group, y = effSize)
) +
  geom_violin(alpha = 0.45, trim = FALSE, draw_quantiles = 0.5,
              show.legend = FALSE, color = 1) +
  geom_point(shape = 21, fill = "white", size = 3,
             position = position_jitter(width = 0.1, height = 0),
             alpha = 0.35, show.legend = FALSE) +
  xlab(NULL) +
  ylab("Absolute effect size") +
  facet_nested(fcast_type ~ beta_pretty, labeller = label_fn2) +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
  scale_color_manual(values = kingdom_colors) +
  scale_y_log10() +
  {if (nrow(stat.test) > 0) stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = TRUE, size = 4) else NULL}

ggsave(here("figures", "effsize_f_b.png"), b_vs_f_fcast_type_plot,
       width = 14, height = 7, dpi = 300)
cat("Saved: figures/effsize_f_b.png\n")

# ══════════════════════════════════════════════════════════════════════════════
# Taxonomic vs Functional effect sizes (Bacteria only, env_cycl)
# ══════════════════════════════════════════════════════════════════════════════

df_cal_fg_tax_v2 <- df_cal_fg_tax %>%
  filter(pretty_group == "Bacteria", model_name == "env_cycl") %>%
  mutate(fcast_type = factor(fcast_type, levels = c("Functional", "Taxonomic")))

groups_with_sufficient_data_v2 <- df_cal_fg_tax_v2 %>%
  group_by(beta_pretty, fcast_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(beta_pretty) %>%
  summarise(n_groups = n_distinct(fcast_type), min_n = min(n), .groups = "drop") %>%
  filter(n_groups >= 2, min_n >= 2)

df_cal_fg_tax_for_tests_v2 <- df_cal_fg_tax_v2 %>%
  semi_join(groups_with_sufficient_data_v2, by = "beta_pretty")

stat.test.v2 <- tryCatch({
  if (nrow(df_cal_fg_tax_for_tests_v2) > 0) {
    df_cal_fg_tax_for_tests_v2 %>%
      group_by(beta_pretty) %>%
      t_test(effSize ~ fcast_type) %>%
      adjust_pvalue(method = "fdr") %>%
      add_significance() %>%
      add_xy_position(x = "fcast_type")
  } else {
    data.frame()
  }
}, error = function(e) {
  cat("Warning: Statistical tests failed:", e$message, "\n")
  data.frame()
})

tax_vs_func_plot <- ggplot(
  data = df_cal_fg_tax_v2,
  aes(x = fcast_type, color = fcast_type, y = effSize)
) +
  geom_violin(alpha = 0.45, trim = FALSE, draw_quantiles = 0.5,
              show.legend = FALSE, color = 1) +
  geom_point(shape = 21, fill = "white", size = 2.5,
             position = position_jitter(width = 0.1, height = 0),
             alpha = 0.35, show.legend = FALSE) +
  xlab("Forecast Type") +
  ylab("Absolute effect size") +
  facet_nested(~ beta_pretty, labeller = label_fn2) +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
  scale_color_manual(values = func_tax_colors) +
  scale_y_log10() +
  {if (nrow(stat.test.v2) > 0) stat_pvalue_manual(stat.test.v2, label = "p.adj.signif", hide.ns = TRUE, size = 4) else NULL}

ggsave(here("figures", "effsize_tax_vs_func.png"), tax_vs_func_plot,
       width = 14, height = 5, dpi = 300)
cat("Saved: figures/effsize_tax_vs_func.png\n")

# ══════════════════════════════════════════════════════════════════════════════
# Significance rates — Functional vs Taxonomic (Bacteria only, env_cycl)
# ══════════════════════════════════════════════════════════════════════════════

df_cal_fg_tax_sig <- df_cal_fg_tax_fixed %>%
  filter(pretty_group == "Bacteria", model_name == "env_cycl") %>%
  distinct(fcast_type, model_id, pretty_group, beta_pretty, significant) %>%
  mutate(fcast_type = factor(fcast_type, levels = c("Functional", "Taxonomic")))

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
    functional_sig   = sum(significant[fcast_type == "Functional"]),
    functional_total = sum(fcast_type == "Functional"),
    taxonomic_sig    = sum(significant[fcast_type == "Taxonomic"]),
    taxonomic_total  = sum(fcast_type == "Taxonomic"),
    .groups = "drop"
  ) %>%
  filter(functional_total >= 2 & taxonomic_total >= 2) %>%
  rowwise() %>%
  mutate(
    p_value   = tryCatch(prop.test(c(functional_sig, taxonomic_sig),
                                   c(functional_total, taxonomic_total))$p.value,
                         error = function(e) NA_real_),
    sig_label = sig_label_from_p(p_value)
  ) %>%
  ungroup()

sig_plot <- ggplot(sig_summary, aes(x = fcast_type, y = sig_rate, fill = fcast_type)) +
  geom_col(alpha = 0.7, show.legend = FALSE, width = 0.65) +
  geom_text(aes(label = paste0(significant, "/", total)), vjust = -0.5, size = 3.5) +
  geom_text(data = sig_test_results,
            aes(x = 1.5, y = 0.9, label = sig_label),
            size = 5, color = "black", inherit.aes = FALSE) +
  xlab("Forecast Type") +
  ylab("Proportion of Significant Predictors") +
  facet_nested(~ beta_pretty, labeller = label_fn2) +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
  scale_fill_manual(values = func_tax_colors) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))

ggsave(here("figures", "significance_tax_vs_func.png"), sig_plot,
       width = 12, height = 5, dpi = 300)
cat("Saved: figures/significance_tax_vs_func.png\n")

# ══════════════════════════════════════════════════════════════════════════════
# Significance rates — Bacteria vs Fungi (env_cycl)
# ══════════════════════════════════════════════════════════════════════════════

df_cal_fg_tax_sig_bf <- df_cal_fg_tax_fixed %>%
  filter(model_name == "env_cycl") %>%
  distinct(fcast_type, model_id, pretty_group, beta_pretty, significant) %>%
  mutate(fcast_type = factor(fcast_type, levels = c("Functional", "Taxonomic")))

sig_summary_bf <- df_cal_fg_tax_sig_bf %>%
  group_by(fcast_type, beta_pretty, pretty_group) %>%
  summarise(
    total = n(),
    significant = sum(significant, na.rm = TRUE),
    sig_rate = significant / total,
    .groups = "drop"
  )

sig_test_results_bf <- df_cal_fg_tax_sig_bf %>%
  group_by(fcast_type, beta_pretty) %>%
  summarise(
    bac_sig   = sum(significant[pretty_group == "Bacteria"], na.rm = TRUE),
    bac_total = sum(pretty_group == "Bacteria"),
    fun_sig   = sum(significant[pretty_group == "Fungi"], na.rm = TRUE),
    fun_total = sum(pretty_group == "Fungi"),
    .groups = "drop"
  ) %>%
  filter(bac_total >= 2 & fun_total >= 2) %>%
  rowwise() %>%
  mutate(
    p_value   = tryCatch(prop.test(c(bac_sig, fun_sig), c(bac_total, fun_total))$p.value,
                         error = function(e) NA_real_),
    sig_label = sig_label_from_p(p_value)
  ) %>%
  ungroup()

sig_bf_plot <- ggplot(sig_summary_bf,
                      aes(x = pretty_group, y = sig_rate, fill = pretty_group)) +
  geom_col(alpha = 0.7, show.legend = FALSE, width = 0.65) +
  geom_text(aes(label = paste0(significant, "/", total)), vjust = -0.5, size = 3.5) +
  geom_text(data = sig_test_results_bf,
            aes(x = 1.5, y = 0.9, label = sig_label),
            size = 5, color = "black", inherit.aes = FALSE) +
  xlab(NULL) +
  ylab("Proportion of Significant Predictors") +
  facet_nested(fcast_type ~ beta_pretty, labeller = label_fn2) +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
  scale_fill_manual(values = kingdom_colors) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))

ggsave(here("figures", "parameter_violin.png"), sig_bf_plot,
       width = 14, height = 5, dpi = 300)
cat("Saved: figures/parameter_violin.png\n")

# ══════════════════════════════════════════════════════════════════════════════
# Proportion of groups with >=1 or >=2 significant predictors
# ══════════════════════════════════════════════════════════════════════════════

# Shared counts: per-group number of significant predictors (env_cycl)
group_sig_counts_overall <- df_cal_fg_tax_fixed %>%
  filter(model_name == "env_cycl") %>%
  distinct(fcast_type, model_id, pretty_group, beta_pretty, significant) %>%
  group_by(fcast_type, model_id, pretty_group) %>%
  summarise(
    n_predictors  = n(),
    n_significant = sum(significant, na.rm = TRUE),
    has_1plus_sig = n_significant >= 1,
    has_2plus_sig = n_significant >= 2,
    .groups = "drop"
  )

# Helper to build a bar plot for proportion of groups with N+ sig predictors
make_prop_plot <- function(data, threshold_col, ylab_text, fill_col, fill_values) {
  summary_df <- data %>%
    group_by(!!sym(fill_col)) %>%
    summarise(
      total_groups    = n(),
      groups_meeting  = sum(!!sym(threshold_col)),
      prop            = groups_meeting / total_groups,
      .groups = "drop"
    )

  # Prop test between two groups
  grp_vals <- split(data[[threshold_col]], data[[fill_col]])
  sig_lab <- ""
  if (length(grp_vals) == 2 && all(sapply(grp_vals, length) >= 2)) {
    p_val <- tryCatch(
      prop.test(sapply(grp_vals, sum), sapply(grp_vals, length))$p.value,
      error = function(e) NA_real_
    )
    sig_lab <- sig_label_from_p(p_val)
  }

  ggplot(summary_df, aes(x = !!sym(fill_col), y = prop, fill = !!sym(fill_col))) +
    geom_col(alpha = 0.7, show.legend = FALSE, width = 0.6) +
    geom_text(aes(label = paste0(groups_meeting, "/", total_groups)),
              vjust = -0.5, size = 4.5) +
    annotate("text", x = 1.5, y = 0.9, label = sig_lab, size = 6, color = "black") +
    xlab(NULL) +
    ylab(ylab_text) +
    base_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
    scale_fill_manual(values = fill_values) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))
}

# ── >=2 sig predictors, Functional vs Taxonomic ──────────────────────────────
group_sig_counts_ft <- group_sig_counts_overall %>%
  mutate(fcast_type = factor(fcast_type, levels = c("Functional", "Taxonomic")))

group_sig_plot <- make_prop_plot(
  group_sig_counts_ft, "has_2plus_sig",
  "Proportion with \u22652 significant predictors",
  "fcast_type", func_tax_colors
)

ggsave(here("figures", "groups_with_multiple_sig_predictors.png"), group_sig_plot,
       width = 4, height = 5, dpi = 300)
cat("Saved: figures/groups_with_multiple_sig_predictors.png\n")

# ── >=1 sig predictor, Bacteria vs Fungi ─────────────────────────────────────
group_sig_plot_1plus <- make_prop_plot(
  group_sig_counts_overall, "has_1plus_sig",
  "Proportion with \u22651 significant predictor",
  "pretty_group", kingdom_colors
)

ggsave(here("figures", "groups_with_1plus_sig_predictors_bac_fun.png"), group_sig_plot_1plus,
       width = 4, height = 5, dpi = 300)
cat("Saved: figures/groups_with_1plus_sig_predictors_bac_fun.png\n")

# ── >=2 sig predictors, Bacteria vs Fungi ────────────────────────────────────
group_sig_plot_2plus <- make_prop_plot(
  group_sig_counts_overall, "has_2plus_sig",
  "Proportion with \u22652 significant predictors",
  "pretty_group", kingdom_colors
)

ggsave(here("figures", "groups_with_2plus_sig_predictors_bac_fun.png"), group_sig_plot_2plus,
       width = 4, height = 5, dpi = 300)
cat("Saved: figures/groups_with_2plus_sig_predictors_bac_fun.png\n")
