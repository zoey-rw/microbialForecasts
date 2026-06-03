source("source.R")

# Peak-month-by-rank, faceted by model type. cycl_only is the pure seasonal model;
# env_cycl carries sin/cos AND env covariates, so its sin/cos describe seasonality
# after env covariates are controlled for ("residual seasonal trend"). env_cov has
# no sin/cos and is excluded. Faceting prevents the two distributions from
# overlaying into a misleading bimodal appearance.
#
# Historically this script also produced figS17_calibration_periods.png and
# several other comparison plots for a workflow that fit separate models per
# calibration period. That workflow has been retired -- models are now fit on
# the full available calibration period (2013-06_2018-01) -- so those plots
# and figS17 have been removed.

if (!file.exists(here("data/summary/seasonal_amplitude.rds"))) {
  stop("seasonal_amplitude.rds not found")
}

seas_vals <- readRDS(here("data/summary/seasonal_amplitude.rds"))
cycl_calibration    <- seas_vals$cycl_only_vals
all_cov_calibration <- seas_vals$env_cycl_vals

if (is.data.frame(cycl_calibration) && nrow(cycl_calibration) > 0 &&
    is.data.frame(all_cov_calibration) && nrow(all_cov_calibration) > 0) {
  peak_df <- bind_rows(
    cycl_calibration    %>% mutate(panel = "cycl_only:\npeak seasonal trend"),
    all_cov_calibration %>% mutate(panel = "env_cycl:\npeak residual seasonal trend\n(after env covariates)")
  )
  p3 <- ggplot(peak_df, aes(x = max, y = rank, color = panel)) +
    geom_point(alpha = 0.6) +
    facet_wrap(~ panel, nrow = 1) +
    scale_x_continuous(breaks = seq(0, 12, 2), limits = c(0, 12)) +
    xlab("Month in which peak is observed") +
    ylab("rank") +
    theme_bw() +
    theme(legend.position = "none")
  png(here("figures","compare_calibration_periods_peak_by_model.png"),
      width = 1400, height = 700)
  print(p3)
  dev.off()
  cat("Saved: figures/compare_calibration_periods_peak_by_model.png\n")
} else {
  cat("cycl_only_vals or env_cycl_vals data not available\n")
}
