# Covariate effect comparison: CLR vs Cloglog vs Dirichlet — env_cycl model
# Ascomycota phylum, 20130601–20180101 calibration period
# Colorblind-safe (Okabe-Ito) palette; shape + color double-encoding

source("source.R")
pacman::p_load(ggplot2, dplyr, stringr)

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
  "Cloglog"   = "#D55E00",   # vermillion
  "CLR"       = "#0072B2",   # deep blue
  "Dirichlet" = "#009E73"    # bluish green
)
model_shapes <- c(
  "Cloglog"   = 16,   # filled circle
  "CLR"       = 17,   # filled triangle
  "Dirichlet" = 15    # filled square
)

# ---- Load predictor-effect summaries ----

clr_eff     <- readRDS(here("data/summary/clr_predictor_effects.rds"))
cloglog_eff <- readRDS(here("data/summary/predictor_effects.rds"))
dirich_eff  <- readRDS(here("data/summary/dirichlet_predictor_effects.rds"))

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

forest_df <- bind_rows(clr_f, clog_f, dir_f) %>%
  mutate(
    beta_label = factor(beta_label, levels = beta_display_order),
    model_type = factor(model_type, levels = names(model_colors))
  )

# ---- Figure ----

fig <- ggplot(forest_df, aes(x = Mean, y = beta_label,
                              color = model_type, shape = model_type)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey60", linewidth = 0.4) +
  geom_linerange(aes(xmin = lo, xmax = hi),
                 linewidth = 0.85,
                 position  = position_dodge(width = 0.6)) +
  geom_point(size     = 3.2,
             position = position_dodge(width = 0.6)) +
  scale_color_manual(values = model_colors,
                     name   = "Model") +
  scale_shape_manual(values = model_shapes,
                     name   = "Model") +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  theme_classic(base_size = 13) +
  theme(
    axis.title.y      = element_blank(),
    axis.text.y       = element_text(size = 12),
    legend.position   = c(0.82, 0.18),
    legend.background = element_rect(fill = "white", color = "grey80",
                                     linewidth = 0.3),
    legend.key.height = unit(1.3, "lines"),
    legend.title      = element_text(size = 11, face = "bold"),
    legend.text       = element_text(size = 11),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3)
  ) +
  xlab("Posterior mean \u00b1 1.96 SD")

# ---- Caption (printed to console; paste into manuscript) ----

caption <- paste0(
  "Fig. X. Posterior mean covariate effects (\u00b11.96 SD) from three model types ",
  "fit to Ascomycota phylum relative abundance across NEON soil fungal communities ",
  "(calibration period 2013\u20132018; environmental + seasonal model structure). ",
  "Filled circles: cloglog beta-regression (univariate, proportion-scale response with driver uncertainty). ",
  "Filled triangles: CLR model (univariate, centered log-ratio response with driver uncertainty). ",
  "Filled squares: Dirichlet model (multivariate, all fungal phyla modeled jointly as a compositional vector). ",
  "Seasonal predictors represent monthly sine and cosine harmonics. ",
  "Temperature and soil moisture include driver uncertainty propagated from repeated measurements. ",
  "EM tree cover = relative basal area of ectomycorrhizal-associated tree species; ",
  "LAI = leaf area index. ",
  "Colors and shapes are redundantly encoded for colorblind accessibility. ",
  "Note: Dirichlet model results are preliminary (median Rhat = 1.18 at 10k iterations; ",
  "target < 1.1). Wider credible intervals partly reflect incomplete convergence."
)
cat("\n--- FIGURE CAPTION ---\n", caption, "\n\n")

# ---- Save ----

outpath <- here("figures/fig_compare_CLR_betareg.png")
ggsave(outpath, fig, width = 7, height = 5.5, dpi = 300)
cat("Saved:", outpath, "\n")
