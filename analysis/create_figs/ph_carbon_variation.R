source("source.R")
pacman::p_load(ggplot2, dplyr, tidyr, cowplot)

# Variation across sampling plots and timepoints, within 3 example NEON sites,
# for observed soil pH (A) and percent carbon (B).
#
# Source: NEON DP1.10086.001 (Soil physical and chemical properties, periodic),
# specifically sls_soilData. Each point is one soil core measurement; the
# production model collapses these to one plot-level mean (treated as constant
# over time) — this figure illustrates the within-plot / across-timepoint
# variability that the constant-mean assumption sets aside.

core_path <- here("data/clean/soilCore_raw_measurements.rds")
if (!file.exists(core_path)) {
  stop("Missing data/clean/soilCore_raw_measurements.rds. ",
       "Extract it from the NEON sls_soilData CSV (raw lives on HPC/hard drive); ",
       "keep columns siteID, plotID, collectDate, soilInCaClpH, organicCPercent.")
}

core <- readRDS(core_path)

# Use BART (northeast forest), NIWO (alpine), TALL (southeastern forest):
# three sites with enough plots and core-level measurements to make the
# boxplots informative.
example_sites <- c("BART", "NIWO", "TALL")

ph_sub <- core %>%
  filter(siteID %in% example_sites, !is.na(soilInCaClpH)) %>%
  mutate(plotID = factor(plotID))

pc_sub <- core %>%
  filter(siteID %in% example_sites, !is.na(organicCPercent)) %>%
  mutate(plotID = factor(plotID))

# Shared aesthetic — keep boxes uncoloured (per the reference layout); guild /
# kingdom palettes don't apply here since we're showing raw chemistry only.
boxplot_theme <- theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey92", color = NA),
        strip.text       = element_text(face = "bold"),
        axis.text.x      = element_text(angle = 60, hjust = 1, size = 7),
        plot.title       = element_text(face = "bold", size = 11),
        plot.subtitle    = element_text(size = 9, color = "grey30"),
        panel.grid.minor = element_blank())

panel_pH <- ggplot(ph_sub, aes(x = plotID, y = soilInCaClpH)) +
  geom_boxplot(outlier.shape = NA, fill = "white", color = "grey25",
               linewidth = 0.35, width = 0.65) +
  geom_jitter(width = 0.18, height = 0, alpha = 0.5, size = 0.7,
              color = "black") +
  facet_grid(~siteID, scales = "free_x", space = "free_x") +
  labs(x = NULL,
       y = expression(paste("pH in CaCl"[2], ", soil core")),
       title    = "A  pH at 3 representative NEON sites",
       subtitle = "Plot-level means were treated as constant over time.") +
  boxplot_theme

panel_pC <- ggplot(pc_sub, aes(x = plotID, y = organicCPercent)) +
  geom_boxplot(outlier.shape = NA, fill = "white", color = "grey25",
               linewidth = 0.35, width = 0.65) +
  geom_jitter(width = 0.18, height = 0, alpha = 0.5, size = 0.7,
              color = "black") +
  facet_grid(~siteID, scales = "free_x", space = "free_x") +
  labs(x = "Plot ID", y = "%C, soil core",
       title    = "B  Percent carbon at 3 representative NEON sites",
       subtitle = "Plot-level means were treated as constant over time.") +
  boxplot_theme

p <- plot_grid(panel_pH, panel_pC, ncol = 1, align = "v", rel_heights = c(1, 1))

out_dir <- here("figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
ggsave(file.path(out_dir, "ph_carbon_variation.png"), p,
       width = 11, height = 8, dpi = 200, bg = "white")
cat("Saved: figures/ph_carbon_variation.png\n")
cat("Example sites:", paste(example_sites, collapse = ", "), "\n")
cat(sprintf("pH measurements: %d across %d plots; %%C measurements: %d across %d plots\n",
            nrow(ph_sub), length(unique(ph_sub$plotID)),
            nrow(pc_sub), length(unique(pc_sub$plotID))))
