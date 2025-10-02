#!/usr/bin/env Rscript
# Hindcast examples script for 2013-2018 legacy models
source("source.R")
library(ggplot2)
library(dplyr)
library(lubridate)

# Read in hindcast data
hindcast_in <- readRDS(here("data/summary/all_hindcasts_raw_plsr2.rds"))

# Filter for 2013-2018 models and available data
hindcast <- hindcast_in %>%
  filter(time_period == "20130601_20180101" & taxon == "chytridiomycota") %>%
  mutate(
    # Add missing columns for plotting
    pred_mean = mean,
    fcast_period = ifelse(dates <= as.Date("2018-01-01"), "calibration", "hindcast"),
    timepoint = as.numeric(as.factor(dates)),
    month = lubridate::month(dates),
    month_label = lubridate::month(dates, label = T),
    # Add truth column (will be NA for now since it's not in the hindcast data)
    truth = NA_real_
  )

cat("Available data points:", nrow(hindcast), "\n")
cat("Available plots:", paste(unique(hindcast$plotID), collapse = ", "), "\n")

# Use available plots - check for BART first, then use DEJU
available_plots <- unique(hindcast$plotID)
bart_plots <- available_plots[grepl("BART", available_plots)]
if (length(bart_plots) > 0) {
  select_plots <- bart_plots[1:min(2, length(bart_plots))]
  cat("Using BART plots:", paste(select_plots, collapse=", "), "\n")
} else {
  select_plots <- c("DEJU_006", "DEJU_009")
  cat("BART not available, using DEJU plots:", paste(select_plots, collapse=", "), "\n")
}

# Filter hindcasts for selected plots
select_hindcasts_select_plots <- hindcast %>% 
  filter(plotID %in% select_plots)

cat("Filtered data points:", nrow(select_hindcasts_select_plots), "\n")

# Create facet labels
tax.labs <- c("Fungi (chytridiomycota)")
names(tax.labs) <- c("chytridiomycota")

# Create the main plot
examples <- ggplot(select_hindcasts_select_plots, 
                  aes(fill = species, x = dates, y = mean, group = plotID)) +
  facet_grid(rows = vars(species),
             cols = vars(plotID),
             drop = T,
             scales = "free",
             labeller = labeller(species = tax.labs)) +
  geom_ribbon(data = ~filter(.x, fcast_period == "calibration"),
              aes(x = dates, ymin = lo, ymax = hi), alpha = 0.35, fill = "blue") +
  geom_ribbon(data = ~filter(.x, fcast_period == "calibration"),
              aes(x = dates, ymin = lo_25, ymax = hi_75), alpha = .65, fill = "blue") +
  geom_line(data = ~filter(.x, fcast_period == "calibration"),
            alpha = 0.8, color = "blue") +
  # Add hindcast period (2018-2020)
  geom_ribbon(data = ~filter(.x, fcast_period == "hindcast"),
              aes(x = dates, ymin = lo, ymax = hi), alpha = 0.35, fill = "red") +
  geom_ribbon(data = ~filter(.x, fcast_period == "hindcast"),
              aes(x = dates, ymin = lo_25, ymax = hi_75), alpha = .65, fill = "red") +
  geom_line(data = ~filter(.x, fcast_period == "hindcast"),
            alpha = 0.8, color = "red") +
  # geom_point(aes(y = as.numeric(truth)), position = position_jitter()) +  # No truth column available
  xlab(NULL) + labs(fill = '') +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  theme_classic(base_size = 18) +
  theme(panel.spacing = unit(.1, "cm"),
        legend.position = "none",
        plot.margin = unit(c(.1, .2, 2, 0), "cm")) + 
  ylab("Abundance") +
  ggtitle("2013-2020 Chytridiomycota Fungi Hindcasts (Blue=Calibration, Red=Hindcast)")

# Display the plot
print(examples)

# Save the plot
png(here("figures", "hindcast_examples_original_script.png"), width = 1000, height = 600)
print(examples)
dev.off()

cat("Plot saved to figures/hindcast_examples_original_script.png\n")

