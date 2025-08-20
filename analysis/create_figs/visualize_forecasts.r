library(ggforce)
source("source.R")

# view output of functional group forecasts.

# Check if forecast plotting data files exist
forecast_data_files <- c(
  "data/forecast_plotting_data_coh_multi_its.rds",
  "data/forecast_plotting_data_no_coh.rds", 
  "data/forecast_plotting_data_coh_16s.rds",
  "data/forecast_plotting_multi_no.rds"
)

# Check which files exist
existing_files <- forecast_data_files[file.exists(here(forecast_data_files))]

if (length(existing_files) == 0) {
  stop("No forecast plotting data files found. Required files:\n", 
       paste("  -", forecast_data_files, collapse = "\n"),
       "\n\nPlease regenerate these files from the model analysis pipeline.")
}

cat("Loading forecast plotting data from:", paste(existing_files, collapse = ", "), "\n")

# Load existing data files
dat_list <- list()
for (i in seq_along(existing_files)) {
  tryCatch({
    dat_list[[i]] <- readRDS(here(existing_files[i]))
  }, error = function(e) {
    stop("Error loading ", existing_files[i], ": ", e$message)
  })
}

# Combine data files
plot.allrank <- data.table::rbindlist(dat_list, fill = TRUE)

# Check if we have data
if (nrow(plot.allrank) == 0) {
  stop("No data found in forecast plotting files. Please check data generation.")
}

# Filter data
plot.allrank <- plot.allrank %>% filter(rank == "No_cohesion")
plot.subset <- plot.allrank

plot.name <- "CPER_002"
plot.name <- "OSBS_001"

# acetogen_anaerobic, alkaline_stress, denitrification, erythromycin_antibiotic, heat_stress, nitrification
# Ectomycorrhizal - somethign wrong with validation?
plot.allrank$pretty_rank <- ifelse(plot.allrank$rank == "No_cohesion", "No cohesion",
                                   ifelse(plot.allrank$rank == "coh_multi", "Fungal-bacterial cohesion", NA))

# Subset to group that cohesion mattered
plot.subset <- plot.allrank[plot.allrank$rank %in% c("No_cohesion","coh_multi") & 
                           plot.allrank$group.name %in% c("erythromycin_antibiotic", "nitrification"),]
plot.subset[plot.subset$group.name=="erythromycin_antibiotic",]$pretty_name <- "Erythromycin-resistant"

# Filter for plotting - use available data
to_plot = plot.allrank %>% 
  filter(rank == "No_cohesion" & 
         group.name == "erythromycin_antibiotic")

# Check if we have data to plot
if (nrow(to_plot) > 0) {
  # View actual models & hindcasts for a single site
  p <- ggplot(to_plot, aes(group = data.type)) +
    geom_line(aes(x = time, y = value), show.legend = FALSE) +
    geom_ribbon(aes(x = time, ymin = lo_2.5, ymax = hi_97.5, fill = data.type), alpha = 0.3) +
    theme_bw() +
    scale_fill_brewer(palette = "Set1") +
    theme(text = element_text(size = 14), 
          panel.spacing = unit(.2, "cm"),
          legend.position = "bottom",
          plot.margin = unit(c(.2, .2, 2, .2), "cm")) + 
    ylab("Abundance") +
    xlab(NULL) +
    labs(fill='') +
    ggtitle("Example Forecast Plot")
  
  print(p)
  cat("Example forecast plot created successfully\n")
} else {
  cat("No data available for plotting\n")
}

# Create plot subset for additional analysis
plot.subset <- plot.allrank %>% 
  filter(rank %in% c("No_cohesion") & 
         group.name %in% c("erythromycin_antibiotic"))

# Check if plot.subset has data
if (nrow(plot.subset) > 0) {
  # WITH validation points
  p1 <- ggplot(plot.subset, aes(x = time, y = value)) + 
    geom_ribbon(aes(ymin = lo_2.5, ymax = hi_97.5), alpha = .2) +
    geom_point(aes(color = data.type), size = 4) + 
    labs(x = "", y = "Abundance (CLR-transformed)", title = paste0("Microbial abundance at NEON plot: ", plot.name)) +
    theme_bw() +
    scale_fill_brewer(palette = "Set1") +
    theme(text = element_text(size = 14), 
          panel.spacing = unit(.2, "cm"),
          legend.position = "bottom",
          plot.margin = unit(c(.2, .2, 2, .2), "cm")) + 
    ylab("Abundance") +
    xlab(NULL) +
    labs(fill='') +
    ggtitle("Forecast with Validation Points")
  
  print(p1)
  cat("Forecast with validation points plot created successfully\n")
  
  ### without validation points
  plot.subset_no_val <- plot.subset
  plot.subset_no_val[plot.subset_no_val$data.type=="forecast",]$value <- NaN
  
  p2 <- ggplot(plot.subset_no_val, aes(x = time, y = value)) + 
    geom_ribbon(aes(ymin = lo_2.5, ymax = hi_97.5), alpha = .2) +
    geom_point(aes(color = data.type), size = 4) + 
    labs(x = "", y = "Abundance (CLR-transformed)", title = paste0("Microbial abundance at NEON plot: ", plot.name)) +
    theme_bw() +
    scale_fill_brewer(palette = "Set1") +
    theme(text = element_text(size = 14), 
          panel.spacing = unit(.2, "cm"),
          legend.position = "bottom",
          plot.margin = unit(c(.2, .2, 2, .2), "cm")) + 
    ylab("Abundance") +
    xlab(NULL) +
    labs(fill='') +
    ggtitle("Forecast without Validation Points")
  
  print(p2)
  cat("Forecast without validation points plot created successfully\n")
} else {
  cat("No data available for plot.subset analysis\n")
}

# Check if plot.allrank has the expected structure for plotting
if (is.list(plot.allrank) && length(plot.allrank) >= 5) {
  fg_plots <- plot.allrank[[1]]
  bac_phy_plots <- plot.allrank[[2]]
  fun_phy_plots <- plot.allrank[[3]]
  genus_plots <- c(plot.allrank[[4]], plot.allrank[[5]])
  
  # Create plots if data is available
  if (is.list(fg_plots) && length(fg_plots) > 0) {
    ggarrange(plotlist = fg_plots, nrow = 4, ncol = 2)
  }
  if (is.list(bac_phy_plots) && length(bac_phy_plots) > 0) {
    ggarrange(plotlist = bac_phy_plots, ncol = 2, nrow = 4)
  }
  if (is.list(fun_phy_plots) && length(fun_phy_plots) > 0) {
    ggarrange(plotlist = fun_phy_plots, ncol = 2, nrow = 4)
  }
  if (length(genus_plots) > 0) {
    ggarrange(plotlist = genus_plots, ncol = 2, nrow = 4)
  }
  
  # Modify plots if they exist
  if (is.list(fg_plots) && length(fg_plots) >= 16) {
    fg_plots[[1]] <- fg_plots[[1]] + scale_y_continuous(breaks=c(-1,-2)) + labs(y="")
    fg_plots[[8]] <- fg_plots[[8]] + scale_y_continuous(breaks=c(-1,-2), limits = c(-2.2, -.8)) + labs(y="")
    fg_plots[[16]] <- fg_plots[[16]] + labs(title = paste0("Plant pathogens at plot: CPER_001")) + labs(y="")
  }
  
  if (length(genus_plots) >= 17) {
    genus_plots[[1]] <- genus_plots[[1]] + labs(y="")
    genus_plots[[17]] <- genus_plots[[17]] + labs(y="")
  }
  
  # Create combined figure
  if (is.list(fg_plots) && length(fg_plots) >= 16 && length(genus_plots) >= 17) {
    plots_for_mike <- list(fg_plots[[1]], fg_plots[[8]], fg_plots[[16]], genus_plots[[1]], genus_plots[[17]])
    figure <- ggarrange(plotlist = plots_for_mike, nrow=5, ncol = 1)
    annotate_figure(figure, left = text_grob("Abundance (CLR-transformed)", size = 29, rot = 90))
  }
} else {
  cat("plot.allrank does not have expected structure for plotting. Data may need to be regenerated.\n")
}
