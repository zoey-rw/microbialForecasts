# Covariate effect comparison: CLR vs Cloglog vs Dirichlet — env_cycl model
# Ascomycota only (all three models comparable for this taxon)
# Two-panel figure: (A) covariate effects, (B) observed vs estimated
# Colorblind-safe (Okabe-Ito) palette; shape + color double-encoding

source("source.R")
pacman::p_load(ggplot2, dplyr, stringr, patchwork)

# ---- Helpers ----

clean_beta <- function(x) ifelse(grepl("Ecto", as.character(x)), "EM trees", as.character(x))

# Human-readable covariate labels (display order = top to bottom on y-axis)
beta_labels <- c(
  "sin"         = "Seasonal (sin)",
  "cos"         = "Seasonal (cos)",
  "LAI"         = "Leaf area index",
  "EM trees"    = "EM tree cover",
  "pC"          = "Soil carbon",
  "pH"          = "Soil pH",
  "Moisture"    = "Soil moisture",
  "Temperature" = "Temperature"
)
beta_display_order <- rev(beta_labels)   # bottom = Temperature on y-axis

# Okabe-Ito palette — safe for deuteranopia, protanopia, tritanopia
model_colors <- c(
  "Cloglog"        = "#D55E00",   # vermillion
  "CLR"            = "#0072B2",   # deep blue
  "Trunc. normal"  = "#CC79A7",   # pink
  "Dirichlet"      = "#009E73"    # bluish green
)
model_shapes <- c(
  "Cloglog"        = 16,   # filled circle
  "CLR"            = 17,   # filled triangle
  "Trunc. normal"  = 18,   # filled diamond
  "Dirichlet"      = 15    # filled square
)

# ---- Load predictor-effect summaries ----

clr_eff      <- readRDS(here("data/summary/clr_predictor_effects.rds"))
cloglog_eff  <- readRDS(here("data/summary/predictor_effects.rds"))
truncnorm_eff <- readRDS(here("data/summary/truncated_normal_predictor_effects.rds"))
dirich_eff   <- readRDS(here("data/summary/dirichlet_predictor_effects.rds"))

mk_forest_df <- function(df, model_label) {
  df %>%
    filter(model_name == "env_cycl") %>%
    mutate(beta = clean_beta(as.character(beta))) %>%
    filter(beta %in% names(beta_labels)) %>%
    transmute(
      beta,
      beta_label = beta_labels[beta],
      Mean,
      lo         = Mean - 1.96 * SD,
      hi         = Mean + 1.96 * SD,
      model_type = model_label
    )
}

# All three models: ascomycota only
clr_f <- clr_eff %>%
  filter(model_name == "env_cycl",
         taxon       == "ascomycota",
         time_period == "20130601_20180101") %>%
  mutate(beta = as.character(beta)) %>%
  group_by(model_name, beta) %>% slice(1) %>% ungroup() %>%
  mk_forest_df("CLR")

clog_f <- cloglog_eff %>%
  filter(model_name == "env_cycl",
         taxon       == "ascomycota",
         time_period == "20130601_20180101") %>%
  mutate(beta = as.character(beta)) %>%
  mk_forest_df("Cloglog")

tn_f <- truncnorm_eff %>%
  filter(taxon == "ascomycota") %>%
  mutate(beta = clean_beta(as.character(beta))) %>%
  filter(beta %in% names(beta_labels)) %>%
  mk_forest_df("Trunc. normal")

dir_f <- dirich_eff %>%
  filter(model_name == "env_cycl", taxon == "ascomycota") %>%
  mutate(beta = clean_beta(as.character(beta))) %>%
  filter(beta %in% names(beta_labels)) %>%
  transmute(
    beta,
    beta_label = beta_labels[beta],
    Mean,
    lo         = Mean - 1.96 * SD,
    hi         = Mean + 1.96 * SD,
    model_type = "Dirichlet"
  )

forest_df <- bind_rows(clr_f, clog_f, tn_f, dir_f) %>%
  mutate(
    beta_label = factor(beta_label, levels = beta_display_order),
    model_type = factor(model_type, levels = names(model_colors))
  )

# ---- Panel A: Covariate effects ----

fig_a <- ggplot(forest_df, aes(x = Mean, y = beta_label,
                                color = model_type, shape = model_type)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey60", linewidth = 0.4) +
  geom_linerange(aes(xmin = lo, xmax = hi),
                 linewidth = 0.85,
                 position  = position_dodge(width = 0.6)) +
  geom_point(size     = 3.2,
             position = position_dodge(width = 0.6)) +
  scale_color_manual(values = model_colors, name = "Model") +
  scale_shape_manual(values = model_shapes, name = "Model") +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  theme_bw(base_size = 13) +
  theme(
    axis.title.y      = element_blank(),
    axis.text.y       = element_text(size = 12),
    legend.position   = "bottom",
    legend.title      = element_text(size = 11, face = "bold"),
    legend.text       = element_text(size = 11),
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3)
  ) +
  xlab("Posterior mean \u00b1 1.96 SD")

# ---- Panel B: Observed vs Estimated (calibration) ----

# Load cloglog plot estimates
plot_est <- readRDS(here("data/summary/plot_estimates.rds"))

# CLR plot estimates: use CLR-scale obs vs est if predictions are available
clr_sum <- readRDS(here("data/summary/clr_regression_summaries.rds"))
clr_pe  <- clr_sum$plot_est
clr_ove <- clr_pe %>%
  filter(species == "ascomycota",
         !is.na(model_name), model_name == "env_cycl",
         !is.na(Mean), Mean != 0, !is.na(truth)) %>%
  transmute(estimated = Mean, observed = truth, model_type = "CLR")
has_clr_panel <- nrow(clr_ove) > 10

# Dirichlet plot estimates from summary file
dirich_sum <- readRDS(here("data/summary/dirichlet_regression_summaries.rds"))
dirich_pe <- dirich_sum$plot_est

# Subset to env_cycl ascomycota for comparable obs vs est
clog_ove <- plot_est %>%
  filter(species == "ascomycota", model_name == "env_cycl") %>%
  transmute(estimated = Mean, observed = truth, model_type = "Cloglog")

dir_ove <- dirich_pe %>%
  filter(taxon == "ascomycota",
         !is.na(model_name), model_name == "env_cycl") %>%
  transmute(estimated = Mean, observed = truth, model_type = "Dirichlet")

# Sample to reasonable size for plotting
set.seed(42)
ove_parts <- list(
  clog_ove %>% sample_n(min(2000, n())),
  dir_ove  %>% sample_n(min(2000, n()))
)
if (has_clr_panel) {
  ove_parts <- c(ove_parts, list(clr_ove %>% sample_n(min(2000, n()))))
}
ove_df <- bind_rows(ove_parts) %>%
  filter(!is.na(estimated), !is.na(observed)) %>%
  mutate(model_type = factor(model_type, levels = names(model_colors)))

# Compute R-squared per model
rsq <- ove_df %>%
  group_by(model_type) %>%
  summarize(
    rsq = cor(observed, estimated, use = "complete.obs")^2,
    rmse = sqrt(mean((observed - estimated)^2, na.rm = TRUE)),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("R\u00b2 = %.3f\nRMSE = %.3f", rsq, rmse))

fig_b <- ggplot(ove_df, aes(x = observed, y = estimated, color = model_type)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.15, size = 0.8) +
  geom_text(data = rsq, aes(label = label),
            x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3,
            size = 3.2, show.legend = FALSE) +
  facet_wrap(~model_type, scales = if (has_clr_panel) "free" else "fixed") +
  scale_color_manual(values = model_colors, guide = "none") +
  { if (!has_clr_panel) coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) } +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey95"),
    strip.text       = element_text(size = 11),
    panel.grid.minor = element_blank()
  ) +
  xlab("Observed") +
  ylab("Estimated")

# ---- Combine panels ----

combined <- fig_a / fig_b +
  plot_layout(heights = c(2, 1.2)) +
  plot_annotation(tag_levels = "A")

# ---- Save ----

outpath <- here("figures/fig_compare_CLR_betareg.png")
ggsave(outpath, combined, width = 8, height = 9, dpi = 300)
cat("Saved:", outpath, "\n")
