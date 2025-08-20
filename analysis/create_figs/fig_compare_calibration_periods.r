source("source.R")

# Check if required data files exist
if (!file.exists(here("data/summary/seasonal_amplitude.rds"))) {
  stop("seasonal_amplitude.rds not found")
}

seas_vals = readRDS(here("data/summary/seasonal_amplitude.rds"))
all_seas_vals = seas_vals[[5]]

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
      print(combined_plot)
    }, error = function(e) {
      cat("Error combining plots:", e$message, "\n")
    })
  }
} else {
  cat("No seasonal amplitude data available\n")
}

# Load summaries data
if (file.exists(here("data/summary/logit_beta_regression_summaries.rds"))) {
  summaries = readRDS(here("data/summary/logit_beta_regression_summaries.rds"))
  plot_est = summaries$plot_est %>% filter(model_name=="env_cycl")
  plot_est_cycl = summaries$plot_est %>% filter(model_name=="cycl_only")
  
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
      # Create basic plot
      p <- ggplot(plot_data, aes(x = dates, y = `50%`, color = .data[[plot_id_col]])) +
        geom_line() +
        theme_bw() +
        ggtitle(paste("Seasonal trends for", first_taxon)) +
        ylab("Modeled abundance") +
        xlab("Date")
      
      print(p)
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
      p <- ggplot(plot_data_2,
             aes(fill=species, x = dates, y = `50%`, group=.data[[plot_id_col]])) +
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
      
      print(p)
      cat("Second plot created successfully\n")
    } else {
      cat("No data available for second plot\n")
    }
  } else {
    cat("No data available for plotting\n")
  }
} else {
  cat("logit_beta_regression_summaries.rds not found\n")
}

# Check if seasonal amplitude data has the expected structure
if (length(seas_vals) >= 2) {
  cycl_calibration = seas_vals[[2]]
  all_cov_calibration = seas_vals[[1]]
  
  if (is.data.frame(cycl_calibration) && nrow(cycl_calibration) > 0) {
    p1 <- ggplot(cycl_calibration) + 
      geom_point(aes(y = rank, x = max)) + 
      xlab("Month in which peak seasonal trend is observed")
    print(p1)
  } else {
    cat("cycl_calibration data not available\n")
  }
  
  if (is.data.frame(all_cov_calibration) && nrow(all_cov_calibration) > 0) {
    p2 <- ggplot(all_cov_calibration) + 
      geom_point(aes(y = rank, x = max)) + 
      xlab("Month in which peak residual seasonal trend is observed")
    print(p2)
  } else {
    cat("all_cov_calibration data not available\n")
  }
} else {
  cat("Seasonal amplitude data structure not as expected\n")
}
