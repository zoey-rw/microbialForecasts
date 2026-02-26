#!/usr/bin/env Rscript
# Comparison of env_cycl model uncertainty calibration across two species
source("source.R")
library(ggplot2)
library(dplyr)
library(lubridate)

# Read hindcast data (07_tidyHindcasts_v2 writes plsr2 from driver_uncertainty site files)
hindcast_path <- here("data/summary/all_hindcasts_plsr2.rds")
if (!file.exists(hindcast_path)) hindcast_path <- here("data/summary/all_hindcasts_raw.rds")
hindcast_in <- readRDS(hindcast_path)

# Filter for two env_cycl species (basidiomycota and rozellomycota); use median (med) for center line
hindcast <- hindcast_in %>%
  filter(species %in% c("basidiomycota", "rozellomycota"), model_name == "env_cycl") %>%
  mutate(
    timepoint = as.numeric(as.factor(dates)),
    month = lubridate::month(dates),
    month_label = lubridate::month(dates, label = T)
  )

cat("Available data points:", nrow(hindcast), "\n")
cat("Available species:", paste(unique(hindcast$species), collapse = ", "), "\n")

# Use BART plots for both species
select_plots <- c("BART_001", "BART_002")
cat("Using BART plots:", paste(select_plots, collapse=", "), "\n")

# Filter hindcasts for selected plots
select_hindcasts_select_plots <- hindcast %>%
  filter(plotID %in% select_plots)

# One env_cycl model per species so ribbons are not duplicated per date
model_ids_to_keep <- hindcast %>%
  filter(plotID %in% select_plots) %>%
  distinct(species, model_id) %>%
  group_by(species) %>%
  slice(1L) %>%
  pull(model_id)
select_hindcasts_select_plots <- select_hindcasts_select_plots %>%
  filter(model_id %in% model_ids_to_keep)

cat("Filtered data points (one env_cycl model per species):", nrow(select_hindcasts_select_plots), "\n")

# Create facet labels (both env_cycl)
tax.labs <- c("Basidiomycota (env_cycl)", "Rozellomycota (env_cycl)")
names(tax.labs) <- c("basidiomycota", "rozellomycota")

# Plot uses median (med) for center line
model_comparison <- ggplot(select_hindcasts_select_plots, 
                  aes(x = dates, y = med, group = plotID)) +
  facet_grid(rows = vars(species),
             cols = vars(plotID),
             drop = T,
             scales = "free",
             labeller = labeller(species = tax.labs)) +
  # Uncertainty bands for calibration period
  geom_ribbon(data = ~filter(.x, fcast_period == "calibration"),
              aes(x = dates, ymin = lo, ymax = hi), alpha = 0.2, fill = "steelblue") +
  geom_ribbon(data = ~filter(.x, fcast_period == "calibration"),
              aes(x = dates, ymin = lo_25, ymax = hi_75), alpha = 0.4, fill = "steelblue") +
  # Uncertainty bands for hindcast period
  geom_ribbon(data = ~filter(.x, fcast_period == "hindcast"),
              aes(x = dates, ymin = lo, ymax = hi), alpha = 0.1, fill = "coral") +
  geom_ribbon(data = ~filter(.x, fcast_period == "hindcast"),
              aes(x = dates, ymin = lo_25, ymax = hi_75), alpha = 0.2, fill = "coral") +
  # Median lines (calibration and hindcast)
  geom_line(data = ~filter(.x, fcast_period == "calibration"),
            aes(x = dates, y = med), alpha = 0.8, color = "steelblue", linewidth = 0.8) +
  geom_line(data = ~filter(.x, fcast_period == "hindcast"),
            aes(x = dates, y = med), alpha = 0.6, color = "coral", linewidth = 0.8, linetype = "dashed") +
  # Add vertical line to separate calibration and hindcast periods
  geom_vline(xintercept = as.Date("2018-01-01"), linetype = "dashed", color = "darkgray", linewidth = 0.5) +
  xlab("Date") +
  ylab("Relative Abundance") +
  theme_bw(base_size = 14) +
  theme(panel.spacing = unit(0.3, "cm"),
        legend.position = "none",
        strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5)) +
  labs(title = "Model Type Uncertainty Calibration Comparison",
       subtitle = "Steel blue: Calibration period, Coral: Hindcast period (env_cycl, median)",
       caption = "Vertical line indicates end of calibration period (2018-01-01)")

# Save the plot
png(here("figures", "hindcast_model_type_comparison.png"), width = 1200, height = 800)
print(model_comparison)
dev.off()

cat("Plot saved to figures/hindcast_model_type_comparison.png\n")

# Print summary
cat("\n=== MODEL TYPE UNCERTAINTY SUMMARY ===\n")
cat("Both panels: env_cycl (Environmental + Cyclical). Center line: median (med).\n")
