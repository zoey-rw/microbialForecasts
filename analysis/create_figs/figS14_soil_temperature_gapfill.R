source("source.R")
pacman::p_load(ggplot2, dplyr, tidyr, cowplot, ggpubr)

# Soil temperature gap-fill (sibling of figS15_soil_moisture_gapfill.R).
# A: monthly series per site, coloured by source (NEON, gap-filled with DAYMET).
# B: DAYMET calibration fit against NEON observations.
# Input monthly_soil_temperature.rds is built by
# data_construction/covariate_prep/soil_temperature/combine_soil_temperature_sources.r

temp_path <- here("data/clean/monthly_soil_temperature.rds")
if (!file.exists(temp_path)) stop("Missing data/clean/monthly_soil_temperature.rds")

# Restrict to the (site, month) cells the models actually received
# (all_predictor_data.rds$temp), so the figure shows only the model input.
ap_path <- here("data/clean/all_predictor_data.rds")
if (!file.exists(ap_path)) stop("Missing data/clean/all_predictor_data.rds")
mat <- readRDS(ap_path)$temp
ij  <- which(!is.na(mat), arr.ind = TRUE)
model_cells <- paste(rownames(mat)[ij[, 1]],
                     gsub("[^0-9]", "", colnames(mat)[ij[, 2]]))

d <- readRDS(temp_path) %>%
  filter(!is.na(temperature_out), !is.na(siteID)) %>%
  mutate(date = as.Date(date),
         source = factor(source, levels = c("NEON", "DAYMET")),
         .cell = paste(siteID, format(date, "%Y%m"))) %>%
  filter(.cell %in% model_cells) %>%
  select(-.cell)

source_colors <- c(NEON = "#D55E00", DAYMET = "#009E73")

site_panel_theme <- theme_bw(base_size = 9) +
  theme(strip.background = element_rect(fill = "grey92", color = NA),
        strip.text       = element_text(face = "bold", size = 8),
        axis.text        = element_text(size = 7),
        axis.text.x      = element_text(angle = 30, hjust = 1),
        panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        legend.title     = element_blank())

pA <- ggplot(d, aes(x = date, y = temperature_out, color = source)) +
  geom_point(size = 0.7, alpha = 0.85) +
  geom_errorbar(aes(ymin = low, ymax = hi), alpha = 0.25, linewidth = 0.25) +
  facet_wrap(~siteID, ncol = 5, scales = "free_y") +
  scale_color_manual(values = source_colors) +
  labs(x = "date", y = "Soil temperature (°C)", title = NULL) +
  site_panel_theme +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

cal <- d %>%
  select(neon = NEON_temp_mean, est = daymet_est) %>%
  filter(!is.na(neon), !is.na(est))
fit <- lm(est ~ neon, data = cal)
r2  <- signif(summary(fit)$adj.r.squared, 4)
pB <- ggplot(cal, aes(x = neon, y = est)) +
  geom_point(alpha = 0.5, size = 1.1, color = "black") +
  stat_smooth(method = "lm", formula = y ~ x, color = "#D55E00",
              se = TRUE, linewidth = 0.7) +
  labs(title    = "Calibration fits",
       subtitle = paste0("Adj R2 = ", r2),
       x        = "Observed values from NEON",
       y        = "Predicted values from DAYMET") +
  theme_bw(base_size = 10) +
  theme(plot.subtitle = element_text(color = "grey25"),
        panel.grid.minor = element_blank())

full <- plot_grid(pA + theme(plot.margin = margin(5, 5, 5, 5)),
                  pB,
                  ncol = 2, rel_widths = c(2, 1),
                  labels = c("A", "B"), label_size = 14, label_fontface = "bold")

out_dir <- here("figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
ggsave(file.path(out_dir, "soil_temperature_gapfill.png"), full,
       width = 14, height = 9, dpi = 200, bg = "white")
cat("Saved: figures/soil_temperature_gapfill.png\n")
