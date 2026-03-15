#!/usr/bin/env Rscript
# Diagnostic script to visualize forecast horizon calculations
# Shows pooled data points, loess fits, null baselines, and calculated horizons

library(here)
source("source.R")
library(data.table)
library(ggplot2)
library(gridExtra)

# Load data
cat("Loading horizon data...\n")
in_list <- readRDS(here("data/summary/fcast_horizon_input.rds"))
model_mean <- as.data.table(in_list[[4]])  # per model x site x months_since_obs
null_all <- as.data.table(in_list[[3]])

# Strip null_ prefix and filter to site_mean
null_site <- copy(null_all[null_type == "site_mean"])
nc <- grep("^null_", names(null_site), value = TRUE)
if (length(nc) > 0) setnames(null_site, nc, gsub("^null_", "", nc))

# Load horizon results
res <- readRDS(here("data/summary/fcast_horizon_df.rds"))
horizon_results <- as.data.table(res[[1]])

# Select example models: pick 6 with a range of horizons
# Exclude models where all metrics are NA (no valid predictions)
valid_models <- model_mean[is.finite(months_since_obs) & months_since_obs > 0,
                            .(has_data = any(!is.na(RSQ) | !is.na(mean_crps) | !is.na(RMSE.norm))),
                            by = model_id][has_data == TRUE]$model_id
horizon_results <- horizon_results[model_id %in% valid_models]
horizon_results <- horizon_results[order(rmse_fcast_horizon)]
n <- nrow(horizon_results)
example_idx <- unique(round(seq(1, n, length.out = 6)))
example_models <- horizon_results$model_id[example_idx]

cat("Creating diagnostic figures for", length(example_models), "example models...\n")

create_diagnostic_plot <- function(mid) {
  # Pooled data for this model
  pooled <- model_mean[model_id == mid & is.finite(months_since_obs) &
                         months_since_obs >= 0 & months_since_obs <= 20]
  null_row <- null_site[model_id == mid]

  if (nrow(pooled) < 4 || nrow(null_row) == 0) return(NULL)

  # Null baselines (median across sites)
  null_rsq <- median(null_row$RSQ, na.rm = TRUE)
  null_crps <- median(null_row$mean_crps, na.rm = TRUE)
  null_rmse <- median(null_row$RMSE.norm, na.rm = TRUE)

  # Floor RSQ at 0 for plotting
  pooled[, RSQ_plot := pmax(RSQ, 0)]

  # Get computed horizons
  hr <- horizon_results[model_id == mid]
  rsq_h <- if (nrow(hr) > 0) hr$rsq_fcast_horizon[1] else NA
  crps_h <- if (nrow(hr) > 0) hr$crps_fcast_horizon[1] else NA
  rmse_h <- if (nrow(hr) > 0) hr$rmse_fcast_horizon[1] else NA

  # Fit loess on pooled data
  xgrid <- seq(0, max(pooled$months_since_obs), by = 0.5)
  safe_predict <- function(y, x) {
    ok <- is.finite(y) & is.finite(x)
    if (sum(ok) < 4) return(rep(NA_real_, length(xgrid)))
    fit <- tryCatch(loess(y[ok] ~ x[ok], span = 0.75,
                          control = loess.control(surface = "direct")),
                    error = function(e) NULL)
    if (is.null(fit)) return(rep(NA_real_, length(xgrid)))
    predict(fit, xgrid)
  }

  pred_rsq <- safe_predict(pooled$RSQ_plot, pooled$months_since_obs)
  pred_crps <- safe_predict(pooled$mean_crps, pooled$months_since_obs)
  pred_rmse <- safe_predict(pooled$RMSE.norm, pooled$months_since_obs)

  # Build long data for facets
  pts <- rbind(
    pooled[, .(x = months_since_obs, y = RSQ_plot, metric = "RSQ (NSE)")],
    pooled[, .(x = months_since_obs, y = mean_crps, metric = "CRPS")],
    pooled[, .(x = months_since_obs, y = RMSE.norm, metric = "RMSE.norm")]
  )

  curves <- rbind(
    data.table(x = xgrid, y = pred_rsq, metric = "RSQ (NSE)"),
    data.table(x = xgrid, y = pred_crps, metric = "CRPS"),
    data.table(x = xgrid, y = pred_rmse, metric = "RMSE.norm")
  )

  nulls <- data.table(
    metric = c("RSQ (NSE)", "CRPS", "RMSE.norm"),
    y = c(null_rsq, null_crps, null_rmse)
  )

  horizons <- data.table(
    metric = c("RSQ (NSE)", "CRPS", "RMSE.norm"),
    x = c(rsq_h, crps_h, rmse_h)
  )
  horizons <- horizons[is.finite(x) & x < 20]

  species <- hr$species[1]
  model_name <- hr$model_name[1]

  ggplot() +
    geom_point(data = pts[is.finite(y)], aes(x = x, y = y),
               color = "grey40", size = 1.5, alpha = 0.5) +
    geom_line(data = curves[is.finite(y)], aes(x = x, y = y),
              color = "steelblue", linewidth = 1) +
    geom_hline(data = nulls, aes(yintercept = y),
               linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_vline(data = horizons, aes(xintercept = x),
               linetype = "dotted", color = "forestgreen", linewidth = 0.8) +
    facet_wrap(~ metric, scales = "free_y", ncol = 3) +
    theme_bw(base_size = 11) +
    labs(title = paste0(species, " (", model_name, ")"),
         subtitle = paste0("Horizons: RSQ=", round(rsq_h, 1),
                           "  CRPS=", round(crps_h, 1),
                           "  RMSE=", round(rmse_h, 1), " months"),
         x = "Months since last observation", y = NULL) +
    xlim(0, 20)
}

# Create plots
plots <- lapply(seq_along(example_models), function(i) {
  cat("  Plot", i, "of", length(example_models), ":", example_models[i], "\n")
  create_diagnostic_plot(example_models[i])
})
plots <- Filter(Negate(is.null), plots)

if (length(plots) > 0) {
  # Save individual plots
  for (i in seq_along(plots)) {
    fn <- paste0("diagnostic_horizon_", gsub("[^A-Za-z0-9]", "_", example_models[i]), ".png")
    ggsave(here("figures", fn), plots[[i]], width = 12, height = 3.5, dpi = 150)
  }

  # Combined
  combined <- do.call(grid.arrange, c(plots, ncol = 1))
  ggsave(here("figures", "diagnostic_horizon_combined.png"),
         combined, width = 12, height = 3.5 * length(plots), dpi = 150)
  cat("Saved: figures/diagnostic_horizon_combined.png\n")
}

cat("Diagnostic figure generation complete!\n")
