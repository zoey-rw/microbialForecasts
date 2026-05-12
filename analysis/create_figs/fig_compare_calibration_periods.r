source("source.R")

# Check if required data files exist
if (!file.exists(here("data/summary/seasonal_amplitude.rds"))) {
  stop("seasonal_amplitude.rds not found")
}

seas_vals = readRDS(here("data/summary/seasonal_amplitude.rds"))
# Named slots from 04_tidyEffectSizes.r: seas_vals_long, env_cycl_vals, cycl_only_vals,
# env_cycl_vals_refit, cycl_only_vals_refit, seas_vals (all rows, wide).
# all_seas_vals is the full wide-format table used by the per-taxon trend block below.
all_seas_vals = seas_vals$seas_vals

# Check if plot_seasonal_trend function exists
if (!exists("plot_seasonal_trend")) {
  cat("plot_seasonal_trend function not found - creating placeholder\n")
  plot_seasonal_trend <- function(input_coefficients, input_dates) {
    # Placeholder function
    ggplot() + 
      geom_text(aes(x = 0.5, y = 0.5, label = "plot_seasonal_trend function not available")) +
      theme_void() +
      ggtitle("Placeholder - function not available")
  }
}

taxon_name="ectomycorrhizal_ectomycorrhizal"
# input_dates = unique(select_hindcasts$dateID)  # select_hindcasts not defined
input_dates = c("201601", "201602", "201603", "201605",
                "201606", "201607", "201608", "201609", "201610", "201611", "201612")

# Check if the required data exists
if (length(all_seas_vals) > 0 && nrow(all_seas_vals) > 0) {
  taxon_row = all_seas_vals %>% filter(taxon==!!taxon_name & model_name=="all_covariates" & time_period == "2015-11_2018-01")
  if (nrow(taxon_row) > 0) {
    input_coefficients = data.frame(sin = taxon_row$sin,
                                   cos = taxon_row$cos)
    trend_201511_201801 = plot_seasonal_trend(input_coefficients, input_dates) + ggtitle("Estimated from 2015-11_2018-01")
  } else {
    cat("No data found for taxon:", taxon_name, "and time period 2015-11_2018-01\n")
  }

  taxon_row = all_seas_vals %>% filter(taxon==!!taxon_name & model_name=="all_covariates" & time_period == "2015-11_2020-01")
  if (nrow(taxon_row) > 0) {
    input_coefficients = data.frame(sin = taxon_row$sin,
                                   cos = taxon_row$cos)
    trend_201511_202001 = plot_seasonal_trend(input_coefficients, input_dates) + ggtitle("Estimated from 2015-11_2020-01")
  } else {
    cat("No data found for taxon:", taxon_name, "and time period 2015-11_2020-01\n")
  }

  taxon_row = all_seas_vals %>% filter(taxon==!!taxon_name & model_name=="all_covariates" & time_period == "20130601_20151101")
  if (nrow(taxon_row) > 0) {
    input_coefficients = data.frame(sin = taxon_row$sin,
                                   cos = taxon_row$cos)
    trend_20130601_20151101 = plot_seasonal_trend(input_coefficients, input_dates) + ggtitle("Estimated from 2013-06_2015-11")
  } else {
    cat("No data found for taxon:", taxon_name, "and time period 20130601_20151101\n")
  }
  
  # Only try to combine plots if we have the data
  if (exists("trend_20130601_20151101") && exists("trend_201511_201801") && exists("trend_201511_202001")) {
    tryCatch({
      combined_plot <- ggarrange(trend_20130601_20151101, trend_201511_201801, trend_201511_202001, nrow = 3)
      png(here("figures","compare_calibration_periods_seasonal_trends.png"), width = 1200, height = 1800)
      print(combined_plot)
      dev.off()
    }, error = function(e) {
      cat("Error combining plots:", e$message, "\n")
    })
  }
} else {
  cat("No seasonal amplitude data available\n")
}

# Load plot estimates from separate file (plot_est is intentionally empty in summaries)
if (file.exists(here("data/summary/plot_estimates.rds"))) {
  cat("Loading plot_estimates from separate file...\n")
  plot_est_all = readRDS(here("data/summary/plot_estimates.rds"))
  cat("Loaded", nrow(plot_est_all), "rows from plot_estimates.rds\n")
  
  # Check if model_name column exists
  if ("model_name" %in% names(plot_est_all)) {
    plot_est = plot_est_all %>% filter(model_name=="env_cycl")
    plot_est_cycl = plot_est_all %>% filter(model_name=="cycl_only")
    cat("Filtered to", nrow(plot_est), "env_cycl rows and", nrow(plot_est_cycl), "cycl_only rows\n")
  } else {
    cat("Warning: model_name column not found in plot_estimates\n")
    cat("Available columns:", paste(names(plot_est_all), collapse = ", "), "\n")
    plot_est = plot_est_all
    plot_est_cycl = plot_est_all
  }
} else if (file.exists(here("data/summary/logit_beta_fixed_priors_summaries.rds"))) {
  # Fallback to summaries file (though plot_est should be empty there)
  summaries = readRDS(here("data/summary/logit_beta_fixed_priors_summaries.rds"))
  if ("plot_est" %in% names(summaries) && nrow(summaries$plot_est) > 0) {
    if ("model_name" %in% names(summaries$plot_est)) {
      plot_est = summaries$plot_est %>% filter(model_name=="env_cycl")
      plot_est_cycl = summaries$plot_est %>% filter(model_name=="cycl_only")
    } else {
      plot_est = summaries$plot_est
      plot_est_cycl = summaries$plot_est
    }
  } else {
    cat("Warning: plot_est not available in summaries file\n")
    plot_est = data.frame()
    plot_est_cycl = data.frame()
  }
} else {
  cat("Warning: Neither plot_estimates.rds nor summaries file found\n")
  plot_est = data.frame()
  plot_est_cycl = data.frame()
}
  
  # Check if date_num column exists
  if ("date_num" %in% names(plot_est)) {
    plot_est <- plot_est[!(plot_est$time_period %in% c("2015-11_2018-01","2015-11_2020-01") & plot_est$date_num == 5),]
  }
  if ("date_num" %in% names(plot_est_cycl)) {
    plot_est_cycl <- plot_est_cycl[!(plot_est_cycl$time_period %in% c("2015-11_2018-01","2015-11_2020-01") & plot_est_cycl$date_num == 5),]
  }
  
  # Check if plotID column exists, otherwise use siteID
  plot_id_col <- if ("plotID" %in% names(plot_est_cycl)) "plotID" else "siteID"
  
  # Filter for available data
  available_taxa <- unique(plot_est_cycl$taxon)
  available_plots <- unique(plot_est_cycl[[plot_id_col]])
  
  cat("Available taxa:", paste(available_taxa, collapse = ", "), "\n")
  cat("Available plots/sites:", paste(head(available_plots, 5), collapse = ", "), "\n")
  
  # Only create plots if we have data
  if (nrow(plot_est_cycl) > 0 && length(available_taxa) > 0 && length(available_plots) > 0) {
    # Use first available taxon and plots
    first_taxon <- available_taxa[1]
    first_plots <- head(available_plots, 2)
    
    cat("Creating plot for taxon:", first_taxon, "and plots:", paste(first_plots, collapse = ", "), "\n")
    
    # Filter data for plotting
    plot_data <- plot_est_cycl %>%
      filter(taxon == first_taxon & .data[[plot_id_col]] %in% first_plots)
    
    if (nrow(plot_data) > 0) {
      # Create basic plot - use Mean instead of 50% if 50% doesn't exist
      y_col <- if ("50%" %in% names(plot_data)) "50%" else if ("Mean" %in% names(plot_data)) "Mean" else stop("No y column found")
      p <- ggplot(plot_data, aes(x = dates, y = .data[[y_col]], color = .data[[plot_id_col]])) +
        geom_line() +
        theme_bw() +
        ggtitle(paste("Seasonal trends for", first_taxon)) +
        ylab("Modeled abundance") +
        xlab("Date")
      
      png(here("figures","figS17_calibration_periods.png"), width = 1200, height = 800)
      print(p)
      dev.off()
      cat("Plot created successfully\n")
    } else {
      cat("No data available for plotting with current filters\n")
    }
    
    # Create plot with available data
    plot_data_2 <- plot_est_cycl %>% 
      filter(.data[[plot_id_col]] %in% first_plots,
             taxon == first_taxon &
             model_name == "cycl_only" &
             time_period %in% c("20130601_20151101","2015-11_2018-01","2015-11_2020-01"))
    
    if (nrow(plot_data_2) > 0) {
      # Use Mean instead of 50% if 50% doesn't exist
      y_col_2 <- if ("50%" %in% names(plot_data_2)) "50%" else if ("Mean" %in% names(plot_data_2)) "Mean" else stop("No y column found")
      p <- ggplot(plot_data_2,
             aes(fill=species, x = dates, y = .data[[y_col_2]], group=.data[[plot_id_col]])) +
        theme_bw() +
        scale_fill_brewer(palette = "Paired") +
        theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
              legend.position = "bottom", legend.title = element_text(NULL),
              plot.margin = unit(c(.2, .2, 2, .2), "cm")) + 
        ylab(NULL) +
        geom_point(aes(x = dates, y = as.numeric(truth))) + 
        xlab(NULL) + 
        labs(fill='') +
        facet_grid(~.data[[plot_id_col]])
      
      png(here("figures","compare_calibration_periods_cycl_only.png"), width = 1600, height = 800)
      print(p)
      dev.off()
      cat("Second plot created successfully\n")
    } else {
      cat("No data available for second plot\n")
    }
  } else {
    cat("No data available for plotting\n")
  }

# Peak-month-by-rank, faceted by model type. cycl_only is the pure seasonal model;
# env_cycl carries sin/cos AND env covariates, so its sin/cos describe seasonality
# after env covariates are controlled for ("residual seasonal trend"). env_cov has
# no sin/cos and is excluded. Faceting prevents the two distributions from
# overlaying into a misleading bimodal appearance.
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
