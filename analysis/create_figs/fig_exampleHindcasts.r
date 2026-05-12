#!/usr/bin/env Rscript
# Example hindcast figure: map + ribbon panels for cellulolytic and saprotroph
# at BART (Bartlett Forest, NH), CPER (Central Plains, CO), and WOOD (Woodworth, ND)
#
# To regenerate hindcasts from project root:
#   Rscript 06_hindcast_observed.r --sequential=true --force=true --figs=true \
#     --taxon=cellulolytic --sites=BART,CPER,WOOD
#   Rscript 06_hindcast_observed.r --sequential=true --force=true --figs=true \
#     --taxon=saprotroph --sites=BART,CPER,WOOD

source("source.R")
library(ggplot2)
library(dplyr)
library(cowplot)
library(maps)

# ── 1. Load and filter data ─────────────────────────────────────────────────
cat("Loading hindcast data from per-site files...\n")
hindcast_dir <- here("data/hindcasts/driver_uncertainty")
best_plots <- c("BART_042", "CPER_004", "WOOD_001")
needed_taxa <- c("cellulolytic", "saprotroph")
site_ids <- c("BART", "CPER", "WOOD")

taxa_sites <- expand.grid(taxon = needed_taxa, site = site_ids, stringsAsFactors = FALSE)
hind_list <- list()
for (i in seq_len(nrow(taxa_sites))) {
  fname <- file.path(hindcast_dir,
    paste0("hindcasts_env_cycl_", taxa_sites$taxon[i],
           "_20130601_20180101_with_legacy_covariate_", taxa_sites$site[i], "_observed.rds"))
  if (file.exists(fname)) {
    hind_list[[length(hind_list) + 1]] <- readRDS(fname)
    cat("  Loaded", basename(fname), "\n")
  } else {
    cat("  Not found:", basename(fname), "\n")
  }
}
hindcast <- bind_rows(hind_list)
cat("Loaded", nrow(hindcast), "rows\n")

# ── 2. Prepare plot data ────────────────────────────────────────────────────
# Ensure fcast_period is set
if (!"fcast_period" %in% names(hindcast)) {
  hindcast$fcast_period <- ifelse(!is.na(hindcast$dateID) & hindcast$dateID > 201801,
                                  "hindcast", "calibration")
}

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
    site_label = factor(case_when(
      siteID == "BART" ~ "Bartlett Forest, NH",
      siteID == "CPER" ~ "Central Plains, CO",
      siteID == "WOOD" ~ "Woodworth, ND"
    ), levels = c("Bartlett Forest, NH", "Central Plains, CO", "Woodworth, ND")),
    taxon_label = case_when(
      species == "cellulolytic" ~ "Bacteria (cellulolytic)",
      species == "saprotroph"  ~ "Fungi (saprotroph)"
    ),
    taxon_color = case_when(
      species == "cellulolytic" ~ "Bacteria",
      species == "saprotroph"  ~ "Fungi"
    )
  )

cat("Selected plots:", paste(best_plots, collapse = ", "), "\n")
cat("Plot data:", nrow(plot_data), "rows\n")

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

# Calibration vs held-out sites
validation_sites <- c("ABBY","BARR","BONA","DEJU","HEAL","KONA","LAJA",
                       "LENO","MLBS","RMNP","SOAP","TOOL","WREF","YELL")
neon_conus <- neon_conus %>%
  mutate(site_type = ifelse(siteID %in% validation_sites, "Held-out", "Calibration"))

highlight <- neon_conus %>% filter(siteID %in% site_ids)

# Label positions: offset from site points to avoid overlap
label_pos <- data.frame(
  siteID = c("BART", "CPER", "WOOD"),
  label  = c("Bartlett Forest\n(BART)", "Central Plains\n(CPER)", "Woodworth\n(WOOD)"),
  nudge_x = c(3, 5, -6),
  nudge_y = c(3, -4, 2),
  stringsAsFactors = FALSE
) %>% left_join(highlight, by = "siteID")

panel_map <- ggplot() +
  geom_polygon(data = us_map, aes(x = long, y = lat, group = group),
               fill = "gray95", color = "gray70", linewidth = 0.2) +
  geom_point(data = neon_conus %>% filter(site_type == "Calibration"),
             aes(x = lon, y = lat, shape = site_type), color = "gray40", size = 1.8) +
  geom_point(data = neon_conus %>% filter(site_type == "Held-out"),
             aes(x = lon, y = lat, shape = site_type), color = "gray40", size = 1.8) +
  geom_point(data = highlight, aes(x = lon, y = lat),
             shape = 16, color = "red", size = 3.5) +
  geom_segment(data = label_pos,
               aes(x = lon + nudge_x + sign(nudge_x) * -0.5,
                   y = lat + nudge_y + sign(nudge_y) * -0.5,
                   xend = lon + sign(nudge_x) * 0.3,
                   yend = lat + sign(nudge_y) * 0.3),
               arrow = arrow(length = unit(0.12, "cm")), color = "black") +
  geom_label(data = label_pos,
             aes(x = lon + nudge_x, y = lat + nudge_y, label = label),
             size = 2.8, fontface = "bold", label.size = 0,
             fill = "white", alpha = 0.85, label.padding = unit(0.2, "lines")) +
  scale_shape_manual(values = c("Calibration" = 16, "Held-out" = 17),
                     name = NULL) +
  coord_fixed(1.3, xlim = c(-125, -66), ylim = c(25, 50)) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(5, 5, 5, 5),
        legend.position = c(0.15, 0.15),
        legend.background = element_rect(fill = "white", color = NA),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"))

# ── 4. Hindcast ribbon panels ──────────────────────────────────────────────
taxon_colors <- kingdom_colors  # from source.R

# Calibration boundary date for vertical line
cal_boundary <- as.Date("2018-01-01")

panel_hindcasts <- ggplot(plot_data, aes(x = dates)) +
  # 95% CI ribbon - calibration (darker)
  geom_ribbon(data = ~ filter(.x, fcast_period == "calibration"),
              aes(ymin = lo, ymax = hi, fill = taxon_color), alpha = 0.25) +
  # 50% CI ribbon - calibration
  geom_ribbon(data = ~ filter(.x, fcast_period == "calibration"),
              aes(ymin = lo_25, ymax = hi_75, fill = taxon_color), alpha = 0.35) +
  # 95% CI ribbon - hindcast (lighter)
  geom_ribbon(data = ~ filter(.x, fcast_period == "hindcast"),
              aes(ymin = lo, ymax = hi, fill = taxon_color), alpha = 0.12) +
  # 50% CI ribbon - hindcast
  geom_ribbon(data = ~ filter(.x, fcast_period == "hindcast"),
              aes(ymin = lo_25, ymax = hi_75, fill = taxon_color), alpha = 0.20) +
  # Median line - calibration
  geom_line(data = ~ filter(.x, fcast_period == "calibration"),
            aes(y = med, color = taxon_color), linewidth = 0.6) +
  # Median line - hindcast
  geom_line(data = ~ filter(.x, fcast_period == "hindcast"),
            aes(y = med, color = taxon_color), linewidth = 0.5, alpha = 0.6) +
  # Calibration/hindcast boundary
  geom_vline(xintercept = cal_boundary, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  # Observed points
  geom_point(aes(y = truth), color = "black", size = 1.2, alpha = 0.8) +
  # plotID label inside each panel
  geom_label(data = plot_data %>%
               group_by(taxon_label, site_label, plotID) %>%
               slice(1),
             aes(x = -Inf, y = Inf, label = plotID),
             vjust = 1.3, hjust = -0.1, size = 2.8, color = "gray30",
             fill = "white", alpha = 0.8, linewidth = 0,
             label.padding = unit(0.15, "lines")) +
  facet_grid(rows = vars(taxon_label), cols = vars(site_label),
             scales = "free") +
  scale_fill_manual(values = taxon_colors) +
  scale_color_manual(values = taxon_colors) +
  scale_x_date(date_labels = "%Y", date_breaks = "2 years") +
  labs(x = "Date", y = "Relative abundance") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 10, face = "bold"),
    strip.text.y = element_text(size = 9),
    panel.spacing = unit(0.4, "cm"),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

# ── 5. Combine panels ───────────────────────────────────────────────────────
fig <- plot_grid(
  panel_map, panel_hindcasts,
  ncol = 1, rel_heights = c(0.38, 0.62),
  labels = c("A", ""), label_size = 14
)

# Facet labels: 3 cols x 2 rows = B-G
fig_labeled <- ggdraw(fig) +
  theme(plot.background = element_rect(fill = "white", color = NA)) +
  draw_label("B", x = 0.02, y = 0.59, fontface = "bold", size = 14) +
  draw_label("C", x = 0.35, y = 0.59, fontface = "bold", size = 14) +
  draw_label("D", x = 0.66, y = 0.59, fontface = "bold", size = 14) +
  draw_label("E", x = 0.02, y = 0.30, fontface = "bold", size = 14) +
  draw_label("F", x = 0.35, y = 0.30, fontface = "bold", size = 14) +
  draw_label("G", x = 0.66, y = 0.30, fontface = "bold", size = 14)

# ── 6. Save ─────────────────────────────────────────────────────────────────
out_dir <- here("figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "fig_exampleHindcasts.png"), fig_labeled,
       width = 12, height = 10, dpi = 200)

cat("Saved: data/figures/fig_exampleHindcasts.png\n")
