#!/usr/bin/env Rscript
# Diagnostic script to visualize forecast horizon calculations
# Shows data points, loess fits, null baselines, and calculated horizons

library(here)
source(here("source.R"))
library(data.table)
library(ggplot2)
library(dplyr)

# Ensure figures directory exists
if (!dir.exists(here("figures"))) {
  dir.create(here("figures"), recursive = TRUE)
}

# Load the horizon input data
cat("Loading horizon input data...\n")
in_list <- readRDS(here("data/summary/fcast_horizon_input.rds"))
to_plot <- in_list[[1]]
to_plot_null <- in_list[[2]]
fcast_horizon_null_site <- in_list[[3]]
fcast_horizon_model_mean <- in_list[[4]]

setDT(to_plot)
setDT(to_plot_null)
setDT(fcast_horizon_null_site)
setDT(fcast_horizon_model_mean)

# Load horizon results (if exists)
cat("Loading horizon results...\n")
if (file.exists(here("data/summary/fcast_horizon_results.rds"))) {
  fcast_horizon_results <- readRDS(here("data/summary/fcast_horizon_results.rds"))
  setDT(fcast_horizon_results)
} else {
  cat("Warning: fcast_horizon_results.rds not found, creating empty data.table\n")
  fcast_horizon_results <- data.table(model_id = character(), 
                                      rsq_fcast_horizon = numeric(),
                                      crps_fcast_horizon = numeric(),
                                      rmse_fcast_horizon = numeric())
}

# Select a few example models to visualize
model_id_list <- unique(to_plot$model_id)
example_models <- head(model_id_list, 6)  # First 6 models

cat("Creating diagnostic figures for", length(example_models), "example models...\n")

# Function to create diagnostic plot for a single model
create_diagnostic_plot <- function(model_id, to_plot, to_plot_null, fcast_horizon_results) {
  # Get data for this model
  single_tax <- to_plot[model_id == model_id]
  single_tax_null <- to_plot_null[model_id == model_id]
  
  # Get horizon results
  horizon_row <- fcast_horizon_results[model_id == model_id]
  
  if (nrow(single_tax) == 0 || nrow(single_tax_null) == 0) {
    cat("Skipping", model_id, "- insufficient data\n")
    return(NULL)
  }
  
  # Filter finite months_since_obs
  single_tax_finite <- single_tax[is.finite(months_since_obs) & months_since_obs <= 20]
  
  # Aggregate to grid (same as in main script)
  single_tax_grid <- single_tax_finite[, 
    .(RSQ = mean(RSQ, na.rm = TRUE),
      mean_crps = mean(mean_crps, na.rm = TRUE),
      RMSE.norm = mean(RMSE.norm, na.rm = TRUE)),
    by = months_since_obs]
  single_tax_grid <- single_tax_grid[is.finite(months_since_obs) & 
                                      months_since_obs >= 0 &
                                      (is.finite(RSQ) | is.finite(mean_crps) | is.finite(RMSE.norm))]
  
  # Get null baselines (handle missing mean_crps)
  if ("mean_crps" %in% names(single_tax_null)) {
    crps_col <- "mean_crps"
  } else if ("CRPS_truncated" %in% names(single_tax_null)) {
    crps_col <- "CRPS_truncated"
  } else {
    crps_col <- "CRPS"
  }
  
  null_collapse <- single_tax_null[, 
    .(RSQ_null = mean(RSQ, na.rm = TRUE),
      CRPS_null = mean(get(crps_col), na.rm = TRUE),
      RMSE_null = mean(RMSE.norm, na.rm = TRUE))]
  
  # Get horizons
  rsq_horizon <- if (nrow(horizon_row) > 0) horizon_row$rsq_fcast_horizon[1] else NA
  crps_horizon <- if (nrow(horizon_row) > 0) horizon_row$crps_fcast_horizon[1] else NA
  rmse_horizon <- if (nrow(horizon_row) > 0) horizon_row$rmse_fcast_horizon[1] else NA
  
  # Fit loess curves
  xgrid <- seq(0, 20, by = 0.1)
  
  safe_loess <- function(y, x) {
    valid <- is.finite(x) & is.finite(y)
    if (sum(valid) < 2) return(NULL)
    ux <- unique(x[valid])
    if (length(ux) >= 4 && sum(valid) >= 4 && sd(x[valid], na.rm=TRUE) > 0) {
      tryCatch({
        stats::loess(y[valid] ~ x[valid], span = 0.8, control = loess.control(surface = "direct"))
      }, error = function(e) {
        if (length(ux) >= 2 && sd(x[valid], na.rm=TRUE) > 0) {
          stats::lm(y[valid] ~ x[valid])
        } else {
          NULL
        }
      })
    } else if (length(ux) >= 2 && sd(x[valid], na.rm=TRUE) > 0) {
      stats::lm(y[valid] ~ x[valid])
    } else {
      NULL
    }
  }
  
  predict_on_grid <- function(fit, newx) {
    if (is.null(fit)) return(rep(NA_real_, length(newx)))
    tryCatch({
      as.numeric(stats::predict(fit, data.frame(x = newx)))
    }, error = function(e) {
      rep(NA_real_, length(newx))
    })
  }
  
  # Fit curves
  fit_rsq <- safe_loess(single_tax_grid$RSQ, single_tax_grid$months_since_obs)
  fit_crps <- safe_loess(single_tax_grid$mean_crps, single_tax_grid$months_since_obs)
  fit_rmse <- safe_loess(single_tax_grid$RMSE.norm, single_tax_grid$months_since_obs)
  
  # Predict on grid
  yhat_rsq <- predict_on_grid(fit_rsq, xgrid)
  yhat_crps <- predict_on_grid(fit_crps, xgrid)
  yhat_rmse <- predict_on_grid(fit_rmse, xgrid)
  
  # Create prediction data frames
  pred_rsq <- data.frame(x = xgrid, y = yhat_rsq, metric = "RSQ")
  pred_crps <- data.frame(x = xgrid, y = yhat_crps, metric = "CRPS")
  pred_rmse <- data.frame(x = xgrid, y = yhat_rmse, metric = "RMSE")
  pred_all <- rbind(pred_rsq, pred_crps, pred_rmse)
  
  # Prepare data for plotting
  grid_long <- rbind(
    single_tax_grid[, .(x = months_since_obs, y = RSQ, metric = "RSQ")],
    single_tax_grid[, .(x = months_since_obs, y = mean_crps, metric = "CRPS")],
    single_tax_grid[, .(x = months_since_obs, y = RMSE.norm, metric = "RMSE")]
  )
  
  # Null baselines
  null_df <- data.frame(
    metric = c("RSQ", "CRPS", "RMSE"),
    y = c(null_collapse$RSQ_null[1], null_collapse$CRPS_null[1], null_collapse$RMSE_null[1])
  )
  
  # Horizon lines
  horizon_df <- data.frame(
    metric = c("RSQ", "CRPS", "RMSE"),
    x = c(rsq_horizon, crps_horizon, rmse_horizon)
  )
  horizon_df <- horizon_df[is.finite(horizon_df$x), ]
  
  # Create plot
  p <- ggplot() +
    # Data points
    geom_point(data = grid_long, aes(x = x, y = y), color = "black", size = 2, alpha = 0.6) +
    # Loess fit
    geom_line(data = pred_all, aes(x = x, y = y), color = "blue", linewidth = 1, alpha = 0.7) +
    # Null baseline
    geom_hline(data = null_df, aes(yintercept = y), linetype = "dashed", color = "red", linewidth = 1) +
    # Horizon lines
    geom_vline(data = horizon_df, aes(xintercept = x), linetype = "dotted", color = "green", linewidth = 1) +
    # Facet by metric
    facet_wrap(~ metric, scales = "free_y", ncol = 3) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("Model:", model_id),
      subtitle = paste("Horizons: RSQ =", round(rsq_horizon, 1), 
                       ", CRPS =", round(crps_horizon, 1),
                       ", RMSE =", round(rmse_horizon, 1)),
      x = "Months since last observation",
      y = "Metric value"
    ) +
    xlim(0, 20) +
    theme(plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 9))
  
  return(p)
}

# Create plots for example models
plots_list <- list()
for (i in seq_along(example_models)) {
  model_id <- example_models[i]
  cat("Creating plot", i, "of", length(example_models), ":", model_id, "\n")
  p <- create_diagnostic_plot(model_id, to_plot, to_plot_null, fcast_horizon_results)
  if (!is.null(p)) {
    plots_list[[i]] <- p
  }
}

# Combine plots
if (length(plots_list) > 0) {
  cat("Saving diagnostic figures...\n")
  
  # Save individual plots
  for (i in seq_along(plots_list)) {
    model_id <- example_models[i]
    filename <- paste0("diagnostic_horizon_", gsub("[^A-Za-z0-9]", "_", model_id), ".png")
    filepath <- here("figures", filename)
    ggsave(filepath, plots_list[[i]], width = 14, height = 4, dpi = 150)
  }
  
  # Create a combined plot
  if (length(plots_list) <= 6) {
    library(gridExtra)
    combined <- do.call(grid.arrange, c(plots_list, ncol = 2))
    ggsave(here("figures", "diagnostic_horizon_combined.png"), 
           combined, width = 16, height = 4 * ceiling(length(plots_list) / 2), dpi = 150)
  }
  
  cat("Diagnostic figures saved to figures/ directory\n")
} else {
  cat("No plots created - check data availability\n")
}

# Create summary comparison plot
cat("Creating summary comparison plot...\n")
if (nrow(fcast_horizon_results) > 0 && "model_id" %in% names(fcast_horizon_results)) {
  # Extract model_name from model_id if needed
  if (!"model_name" %in% names(fcast_horizon_results)) {
    fcast_horizon_results[, model_name := gsub("_.*", "", model_id)]
  }
  if (!"pretty_group" %in% names(fcast_horizon_results)) {
    # Try to get from to_plot
    model_info <- unique(to_plot[, .(model_id, pretty_group, model_name)])
    fcast_horizon_results <- fcast_horizon_results[model_info, on = "model_id", nomatch = 0]
  }
  
  if (nrow(fcast_horizon_results) > 0 && "model_name" %in% names(fcast_horizon_results)) {
    summary_data <- fcast_horizon_results[, .(
      rsq_mean = mean(rsq_fcast_horizon, na.rm = TRUE),
      crps_mean = mean(crps_fcast_horizon, na.rm = TRUE),
      rmse_mean = mean(rmse_fcast_horizon, na.rm = TRUE),
      n_models = .N
    ), by = .(model_name, pretty_group)]
  } else {
    summary_data <- data.table()
  }
} else {
  summary_data <- data.table()
}

if (nrow(summary_data) > 0) {
  summary_long <- data.table::melt(summary_data, 
                       id.vars = c("model_name", "pretty_group", "n_models"),
                       measure.vars = c("rsq_mean", "crps_mean", "rmse_mean"),
                       variable.name = "metric",
                       value.name = "horizon")
  
  summary_long$metric <- factor(summary_long$metric, 
                                levels = c("rsq_mean", "crps_mean", "rmse_mean"),
                                labels = c("RSQ", "CRPS", "RMSE"))
  
  p_summary <- ggplot(summary_long, aes(x = model_name, y = horizon, fill = pretty_group)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ metric, scales = "free_y") +
    theme_minimal(base_size = 14) +
    labs(
      title = "Average Forecast Horizons by Model Type",
      x = "Model Type",
      y = "Horizon (months)",
      fill = "Kingdom"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(here("figures", "diagnostic_horizon_summary.png"), 
         p_summary, width = 12, height = 6, dpi = 150)
  
  cat("Summary plot saved\n")
}

cat("Diagnostic figure generation complete!\n")

