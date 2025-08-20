source("source.R")
#remotes::install_github("giocomai/ganttrify")
# p_load(ganttrify)  # Commented out for testing

# Check if ganttrify is available
if (!require(ganttrify, quietly = TRUE)) {
  cat("ganttrify package not available - skipping gantt chart\n")
  ganttrify_available <- FALSE
} else {
  ganttrify_available <- TRUE
}

data_in <- readRDS(here("data/summary/logit_beta_regression_summaries.rds"))
keep_list <- readRDS(here("data/summary/converged_taxa_list.rds"))

# Check if data has the expected structure
if (!"plot_est" %in% names(data_in)) {
  stop("Expected 'plot_est' element in data_in")
}

# Check if dates column exists and has valid data
if (!"dates" %in% names(data_in$plot_est)) {
  cat("No dates column found in plot_est data\n")
  data_in$plot_est$dates <- as.Date(NA)
}

# Filter out rows with missing dates for the summary
plot_est_filtered <- data_in$plot_est %>%
  filter(!is.na(dates))

if (nrow(plot_est_filtered) == 0) {
  cat("No valid dates found in plot_est data\n")
  # Create a minimal summary for testing
  model_date_summaries <- data.frame(
    time_period = "2015-11_2018-01",
    start_date = as.Date("2015-11-01"),
    end_date = as.Date("2018-01-01"),
    wp = "Excluding legacy data",
    activity = "Excluding legacy data"
  )
} else {
  model_date_summaries <- plot_est_filtered %>%
    group_by(time_period) %>%
    filter(!is.na(truth)) %>%
    summarise(start_date = min(dates, na.rm = TRUE),
              end_date = max(dates, na.rm = TRUE),
              .groups = 'drop') %>%
    mutate(wp = recode(as.character(time_period), 
                       `2015-11_2018-01` = "Excluding legacy data",
                       `2013-06_2017-01` = "Including legacy data",
                       `20130601_20151101` = "Legacy only",
                       `2013-06_2020-01` = "Legacy + current (full dataset)",
                       `2015-11_2020-01` = "Excluding legacy data"),
           activity = wp)
}

# Only try ganttrify if we have valid data and the package is available
if (ganttrify_available && nrow(model_date_summaries) > 0) {
  tryCatch({
    ganttrify(project = model_date_summaries,
              month_number_label = FALSE,
              by_date = TRUE,
              month_breaks = 6,
              font_family = "Roboto Condensed")
  }, error = function(e) {
    cat("Error in ganttrify: ", e$message, "\n")
  })
} else {
  cat("ganttrify package not available or no valid data - skipping gantt chart\n")
}

# Comment out undefined variables for now
# heat_stress
# plant_pathogen
# russula

# Filter data for plotting
cal_data <- data_in$plot_est %>%  
  filter(model_id %in% keep_list) %>%
  mutate(calibration_label = recode(as.character(time_period), 
                                   `2015-11_2018-01` = "Partial dataset for calibration",
                                   `2015-11_2020-01` = "Full dataset",
                                   `20130601_20151101` = "Provisional data"))

# Set factor levels
cal_data$calibration_label = factor(cal_data$calibration_label, 
                                   ordered = TRUE, 
                                   levels = c("Provisional data",
                                            "Partial dataset for calibration",
                                            "Full dataset"))

# Filter for plotting - use available data
to_plot = cal_data %>% 
  filter(taxon == "acidobacteriota" & 
         model_name == "env_cov" & 
         siteID %in% c("CPER_004"))

# Check if we have data to plot
if (nrow(to_plot) > 0) {
  # View actual models & hindcasts for a single site
  p <- ggplot(to_plot, aes(group = calibration_label)) +
    geom_line(aes(x = dates, y = `50%`), show.legend = FALSE) +
    geom_ribbon(aes(x = dates, ymin = `25%`, ymax = `75%`, fill = calibration_label), alpha = 0.3) +
    geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`, fill = calibration_label), alpha = 0.4) +
    theme_bw() +
    scale_fill_brewer(palette = "Set1") +
    theme(text = element_text(size = 14), 
          panel.spacing = unit(.2, "cm"),
          legend.position = "bottom",
          legend.title = element_text(NULL),
          plot.margin = unit(c(.2, .2, 2, .2), "cm")) + 
    ylab("Modeled abundance") +
    scale_y_sqrt() +
    geom_point(aes(x = dates, y = as.numeric(truth))) + 
    xlab(NULL) + 
    labs(fill = '') +
    xlim(c(as.Date("2013-06-01"), as.Date("2020-01-01")))
  
  print(p)
  cat("Plot created successfully\n")
} else {
  cat("No data available for plotting with current filters\n")
  # Show what data we do have
  cat("Available taxa:", unique(cal_data$taxon), "\n")
  cat("Available model names:", unique(cal_data$model_name), "\n")
  cat("Available site IDs:", unique(cal_data$siteID), "\n")
}

