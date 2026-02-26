#!/usr/bin/env Rscript
# Example hindcast figure: map + ribbon panels for cellulolytic and saprotroph
# at CPER (Central Plains) and BART (Bartlett Experimental Forest)
#
# To regenerate hindcasts for CPER, HARV, WOOD, BART and selected functional groups
# with diagnostic figures: from project root run
#   Rscript 06_hindcast_observed.r --taxon="cellulolytic|oligotroph|plant_pathogen|saprotroph|ectomycorrhizal|denitrification" --mname=env_cycl --sites=CPER,HARV,WOOD,BART --figs=true [--force=true]
# Then inspect figures/hindcast_diagnostics/observed_sites/driver_uncertainty/ and
# update best_plots / species / siteID filters below to match the chosen plot.

source("source.R")
library(ggplot2)
library(dplyr)
library(cowplot)
library(maps)

# ── 1. Load and filter data ─────────────────────────────────────────────────
cat("Loading hindcast data from per-site files (with fallback)...\n")
hindcast_dir <- here("data/hindcasts/driver_uncertainty")
best_plots <- c("CPER_004", "BART_042")
needed_taxa <- c("cellulolytic", "saprotroph")

# Try per-site hindcast files first
taxa_sites <- list(
  c("cellulolytic", "CPER"), c("cellulolytic", "BART"),
  c("saprotroph",   "CPER"), c("saprotroph",   "BART")
)
hind_list <- list()
for (ts in taxa_sites) {
  fname <- file.path(hindcast_dir,
    paste0("hindcasts_env_cycl_", ts[1],
           "_20130601_20180101_with_legacy_covariate_", ts[2], "_observed.rds"))
  if (file.exists(fname)) {
    hind_list[[length(hind_list) + 1]] <- readRDS(fname)
    cat("  Loaded", basename(fname), "\n")
  } else {
    cat("  Not found:", basename(fname), "\n")
  }
}
hindcast <- bind_rows(hind_list)

# Check if any needed plots are missing and fall back to summary file
found_combos <- hindcast %>%
  filter(plotID %in% best_plots) %>%
  distinct(species, plotID)
needed_combos <- expand.grid(species = needed_taxa, plotID = best_plots,
                             stringsAsFactors = FALSE)
missing <- anti_join(needed_combos, found_combos, by = c("species", "plotID"))

if (nrow(missing) > 0) {
  cat("  Missing from per-site files:", paste(missing$species, missing$plotID, collapse = ", "), "\n")
  cat("  Falling back to all_hindcasts_plsr2.rds for missing combos...\n")
  fallback <- readRDS(here("data/summary/all_hindcasts_plsr2.rds")) %>%
    filter(model_name == "env_cycl") %>%
    inner_join(missing, by = c("species", "plotID"))
  cat("  Got", nrow(fallback), "rows from fallback\n")
  hindcast <- bind_rows(hindcast, fallback)
  rm(fallback); gc()
}

# Fill any NA pretty_group values
hindcast <- hindcast %>%
  mutate(pretty_group = case_when(
    species == "cellulolytic" ~ "Bacteria",
    species == "saprotroph"  ~ "Fungi",
    TRUE ~ pretty_group
  ))

cat("Loaded", nrow(hindcast), "rows\n")

# ── 2. Use pre-selected well-calibrated plots ───────────────────────────────
# CPER_004: cellulolytic 75%/100% cal coverage, saprotroph 39%/83%
# BART_042: cellulolytic 100%/100%, saprotroph 67%/100% cal coverage
cat("Selected plots:", paste(best_plots, collapse = ", "), "\n")

# Trim predictions to start no earlier than first observation per panel
first_obs <- hindcast %>%
  filter(plotID %in% best_plots, !is.na(truth)) %>%
  group_by(plotID, species) %>%
  summarise(first_date = min(dates), .groups = "drop")

plot_data <- hindcast %>%
  filter(plotID %in% best_plots) %>%
  left_join(first_obs, by = c("plotID", "species")) %>%
  filter(dates >= first_date) %>%
  select(-first_date) %>%
  mutate(
    site_label = case_when(
      siteID == "CPER" ~ "Central Plains, CO",
      siteID == "BART" ~ "Bartlett Forest, NH"
    ),
    taxon_label = case_when(
      species == "cellulolytic" ~ "Bacteria (cellulolytic)",
      species == "saprotroph"  ~ "Fungi (saprotroph)"
    ),
    taxon_color = case_when(
      species == "cellulolytic" ~ "Bacteria",
      species == "saprotroph"  ~ "Fungi"
    )
  )

# ── 3. Panel A: NEON site map ───────────────────────────────────────────────
us_map <- map_data("state")

neon_sites <- data.frame(
  siteID = c("ABBY","BARR","BART","BLAN","BONA","CLBJ","CPER","DCFS","DEJU",
             "DELA","DSNY","GRSM","GUAN","HARV","HEAL","JERC","JORN","KONA",
             "KONZ","LAJA","LENO","MLBS","MOAB","NIWO","NOGP","OAES","ONAQ",
             "ORNL","OSBS","PUUM","RMNP","SCBI","SERC","SJER","SOAP","SRER",
             "STEI","STER","TALL","TEAK","TOOL","TREE","UKFS","UNDE","WOOD",
             "WREF","YELL"),
  lat = c(45.76, 71.28, 44.06, 39.06, 65.15, 33.40, 40.82, 47.16, 63.88,
          32.54, 28.13, 35.69, 18.11, 42.54, 63.88, 31.19, 32.59, 39.11,
          39.10, 18.02, 31.85, 37.38, 38.25, 40.05, 46.77, 35.41, 40.18,
          35.96, 29.69, 19.55, 40.28, 38.89, 38.89, 37.11, 37.03, 31.91,
          45.51, 40.46, 32.95, 37.01, 68.66, 45.49, 39.04, 46.23, 47.13,
          45.82, 44.95),
  lon = c(-122.33, -156.61, -71.29, -78.07, -147.50, -97.57, -104.75, -99.11,
          -145.75, -87.80, -81.44, -83.50, -66.87, -72.17, -149.21, -84.47,
          -106.84, -96.61, -96.56, -67.08, -88.16, -80.52, -109.39, -105.58,
          -100.92, -99.06, -112.45, -84.28, -81.99, -155.32, -105.55, -78.14,
          -76.56, -119.73, -119.26, -110.84, -89.59, -96.62, -87.39, -119.01,
          -149.37, -89.59, -95.19, -89.54, -99.24, -121.95, -110.59),
  stringsAsFactors = FALSE
)

# Continental US sites only
neon_conus <- neon_sites %>%
  filter(lon > -130, lon < -60, lat > 24, lat < 50)

# Calibration vs validation (held-out) sites
validation_sites <- c("ABBY","BARR","BONA","DEJU","HEAL","KONA","LAJA",
                       "LENO","MLBS","RMNP","SOAP","TOOL","WREF","YELL")
neon_conus <- neon_conus %>%
  mutate(site_type = ifelse(siteID %in% validation_sites, "Held-out", "Calibration"))

highlight <- neon_conus %>% filter(siteID %in% c("CPER", "BART"))

panel_map <- ggplot() +
  geom_polygon(data = us_map, aes(x = long, y = lat, group = group),
               fill = "gray95", color = "gray70", linewidth = 0.2) +
  geom_point(data = neon_conus %>% filter(site_type == "Calibration"),
             aes(x = lon, y = lat), shape = 16, color = "gray40", size = 1.8) +
  geom_point(data = neon_conus %>% filter(site_type == "Held-out"),
             aes(x = lon, y = lat), shape = 17, color = "gray40", size = 1.8) +
  geom_point(data = highlight, aes(x = lon, y = lat),
             shape = 16, color = "red", size = 3.5) +
  annotate("segment", x = highlight$lon[highlight$siteID == "CPER"] + 3,
           y = highlight$lat[highlight$siteID == "CPER"] - 2,
           xend = highlight$lon[highlight$siteID == "CPER"] + 0.3,
           yend = highlight$lat[highlight$siteID == "CPER"] - 0.3,
           arrow = arrow(length = unit(0.15, "cm")), color = "black") +
  annotate("text", x = highlight$lon[highlight$siteID == "CPER"] + 3,
           y = highlight$lat[highlight$siteID == "CPER"] - 2.5,
           label = "Central Plains\n(CPER)", size = 3, fontface = "bold") +
  annotate("segment", x = highlight$lon[highlight$siteID == "BART"] - 5,
           y = highlight$lat[highlight$siteID == "BART"] - 1,
           xend = highlight$lon[highlight$siteID == "BART"] - 0.3,
           yend = highlight$lat[highlight$siteID == "BART"] - 0.3,
           arrow = arrow(length = unit(0.15, "cm")), color = "black") +
  annotate("text", x = highlight$lon[highlight$siteID == "BART"] - 5,
           y = highlight$lat[highlight$siteID == "BART"] - 1.5,
           label = "Bartlett Forest\n(BART)", size = 3, fontface = "bold") +
  # Legend for site types
  annotate("point", x = -124, y = 27, shape = 16, color = "gray40", size = 1.8) +
  annotate("text", x = -122.5, y = 27, label = "Calibration (n=30)",
           size = 2.5, hjust = 0) +
  annotate("point", x = -124, y = 25.5, shape = 17, color = "gray40", size = 1.8) +
  annotate("text", x = -122.5, y = 25.5, label = "Held-out (n=14)",
           size = 2.5, hjust = 0) +
  coord_fixed(1.3, xlim = c(-125, -66), ylim = c(25, 50)) +
  theme_void() +
  theme(plot.margin = margin(5, 5, 5, 5))

# ── 4. Panels B-E: Hindcast ribbons ─────────────────────────────────────────
taxon_colors <- c("Bacteria" = "#4DAF4A", "Fungi" = "#FF7F00")

# Build separate layers for calibration vs hindcast with different alphas
panel_hindcasts <- ggplot(plot_data, aes(x = dates)) +
  # 95% ribbon - calibration
  geom_ribbon(data = ~ filter(.x, fcast_period == "calibration"),
              aes(ymin = lo, ymax = hi, fill = taxon_color), alpha = 0.3) +
  # 50% ribbon - calibration
  geom_ribbon(data = ~ filter(.x, fcast_period == "calibration"),
              aes(ymin = lo_25, ymax = hi_75, fill = taxon_color), alpha = 0.4) +
  # 95% ribbon - hindcast
  geom_ribbon(data = ~ filter(.x, fcast_period == "hindcast"),
              aes(ymin = lo, ymax = hi, fill = taxon_color), alpha = 0.15) +
  # 50% ribbon - hindcast
  geom_ribbon(data = ~ filter(.x, fcast_period == "hindcast"),
              aes(ymin = lo_25, ymax = hi_75, fill = taxon_color), alpha = 0.25) +
  # Median line - calibration
  geom_line(data = ~ filter(.x, fcast_period == "calibration"),
            aes(y = med, color = taxon_color), linewidth = 0.7) +
  # Median line - hindcast (lighter)
  geom_line(data = ~ filter(.x, fcast_period == "hindcast"),
            aes(y = med, color = taxon_color), linewidth = 0.5, alpha = 0.6) +
  # Observed points
  geom_point(aes(y = truth), color = "black", size = 1.2, alpha = 0.8) +
  # plotID label inside each panel (top-right, with background)
  geom_label(data = plot_data %>%
               group_by(taxon_label, site_label, plotID) %>%
               slice(1),
             aes(x = -Inf, y = Inf, label = plotID),
             vjust = 1.3, hjust = -0.1, size = 3, color = "gray20",
             fill = "white", alpha = 0.7, label.size = 0,
             label.padding = unit(0.15, "lines")) +
  facet_grid(rows = vars(taxon_label), cols = vars(site_label),
             scales = "free") +
  scale_fill_manual(values = taxon_colors) +
  scale_color_manual(values = taxon_colors) +
  labs(x = "Date", y = "Relative abundance") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 11, face = "bold"),
    panel.spacing = unit(0.4, "cm"),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

# ── 5. Combine panels ───────────────────────────────────────────────────────
fig <- plot_grid(
  panel_map, panel_hindcasts,
  ncol = 1, rel_heights = c(0.4, 0.6),
  labels = c("A", ""), label_size = 14
)

# Add B-E labels to the facets via a second label layer
# The facet_grid produces 4 panels: top-left=B, top-right=C, bottom-left=D, bottom-right=E
# We'll use draw_label for manual placement
fig_labeled <- ggdraw(fig) +
  draw_label("B", x = 0.02, y = 0.57, fontface = "bold", size = 14) +
  draw_label("C", x = 0.50, y = 0.57, fontface = "bold", size = 14) +
  draw_label("D", x = 0.02, y = 0.29, fontface = "bold", size = 14) +
  draw_label("E", x = 0.50, y = 0.29, fontface = "bold", size = 14)

# ── 6. Save ─────────────────────────────────────────────────────────────────
out_dir <- here("data", "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "fig_exampleHindcasts.png"), fig_labeled,
       width = 10, height = 10, dpi = 200)

cat("Saved: data/figures/fig_exampleHindcasts.png\n")
