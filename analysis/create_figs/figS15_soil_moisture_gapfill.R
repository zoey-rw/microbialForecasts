source("source.R")
pacman::p_load(ggplot2, dplyr, tidyr, cowplot, ggpubr)

# Soil moisture gap-fill, multi-source.
# A: monthly series per site, coloured by source (NEON, gap-filled with SCAN/SMAP/SMOS).
# B-D: calibration fits for each satellite product against NEON observations.
# Input monthly_soil_moisture.rds is built by
# data_construction/covariate_prep/soil_moisture/combine_soil_moisture_sources.r

mois_path <- here("data/clean/monthly_soil_moisture.rds")
if (!file.exists(mois_path)) stop("Missing data/clean/monthly_soil_moisture.rds")

# Restrict to the (site, month) cells the models actually received
# (all_predictor_data.rds$mois), so the figure shows only the model input.
ap_path <- here("data/clean/all_predictor_data.rds")
if (!file.exists(ap_path)) stop("Missing data/clean/all_predictor_data.rds")
mat <- readRDS(ap_path)$mois
ij  <- which(!is.na(mat), arr.ind = TRUE)
model_cells <- paste(rownames(mat)[ij[, 1]],
                     gsub("[^0-9]", "", colnames(mat)[ij[, 2]]))

m <- readRDS(mois_path) %>%
  filter(!is.na(moisture_out), !is.na(siteID)) %>%
  mutate(date = as.Date(date),
         source = factor(source, levels = c("NEON", "SCAN", "SMAP", "SMOS")),
         .cell = paste(siteID, format(date, "%Y%m"))) %>%
  filter(.cell %in% model_cells) %>%
  select(-.cell)

source_colors <- c(NEON = "#D55E00", SCAN = "#009E73",
                   SMAP = "#56B4E9", SMOS = "#CC79A7")

site_panel_theme <- theme_bw(base_size = 9) +
  theme(strip.background = element_rect(fill = "grey92", color = NA),
        strip.text       = element_text(face = "bold", size = 8),
        axis.text        = element_text(size = 7),
        axis.text.x      = element_text(angle = 30, hjust = 1),
        panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        legend.title     = element_blank())

pA <- ggplot(m, aes(x = date, y = moisture_out, color = source)) +
  geom_point(size = 0.7, alpha = 0.85) +
  geom_errorbar(aes(ymin = low, ymax = hi), alpha = 0.25, linewidth = 0.25) +
  facet_wrap(~siteID, ncol = 5, scales = "free_y") +
  scale_color_manual(values = source_colors) +
  labs(x = "date", y = "Soil moisture",
       title = NULL) +
  site_panel_theme +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

calibration_panel <- function(df, est_col, label) {
  sub <- df %>%
    select(neon_mean, est = all_of(est_col)) %>%
    filter(!is.na(neon_mean), !is.na(est))
  if (nrow(sub) < 3) {
    return(ggplot() + theme_void() +
             annotate("text", x = 0.5, y = 0.5,
                      label = paste("No overlap for", label)))
  }
  fit <- lm(est ~ neon_mean, data = sub)
  r2  <- signif(summary(fit)$adj.r.squared, 4)
  ggplot(sub, aes(x = neon_mean, y = est)) +
    geom_point(alpha = 0.55, size = 1.1, color = "black") +
    stat_smooth(method = "lm", formula = y ~ x, color = "#D55E00",
                se = TRUE, linewidth = 0.7) +
    labs(title    = "Calibration fits",
         subtitle = paste0("Adj R2 = ", r2),
         x        = "Observed values from NEON",
         y        = paste("Predicted values from", label)) +
    theme_bw(base_size = 10) +
    theme(plot.subtitle = element_text(color = "grey25"),
          panel.grid.minor = element_blank())
}

pB <- calibration_panel(m, "SCAN_est", "SCAN")
pC <- calibration_panel(m, "SMAP_est", "SMAP")
pD <- calibration_panel(m, "SMOS_est", "SMOS")

right_col <- plot_grid(pB, pC, pD, ncol = 1, labels = c("B", "C", "D"),
                       label_size = 14, label_fontface = "bold",
                       align = "v")
full <- plot_grid(pA + theme(plot.margin = margin(5, 5, 5, 5)),
                  right_col,
                  ncol = 2, rel_widths = c(2, 1),
                  labels = c("A", ""), label_size = 14, label_fontface = "bold")

out_dir <- here("figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
ggsave(file.path(out_dir, "figS15_soil_moisture_gapfill.png"), full,
       width = 14, height = 9, dpi = 200, bg = "white")
cat("Saved: figures/figS15_soil_moisture_gapfill.png\n")
