# Visualize hindcasts for best/worst of each group
pacman::p_load(scoringRules, reshape2, parallel, lubridate, data.table, ggforce, ggrepel)
source("source.R")
# library(sp)  # Not available
# library(rgdal)  # Not available for this R version
library(cowplot)
library(ggrepel)

# Check if hindcast data exists
if (!file.exists(here("data/summary/all_hindcasts.rds"))) {
  stop("all_hindcasts.rds not found. Please regenerate this file from the model analysis pipeline.")
}

hindcast_in <- readRDS(here("data/summary/all_hindcasts.rds"))

# Get top plots
not_na <- hindcast_in[!is.na(hindcast_in$truth),]
top_plots <- names(tail(sort(table(not_na$plotID)), 8))
hindcast_in <- hindcast_in[hindcast_in$plotID %in% top_plots,]

cat("Hindcast data loaded successfully\n")

# Check if we have data for plotting
if (nrow(hindcast_in) > 0) {
  # Create example plot
  p1 <- ggplot(hindcast_in %>% filter(siteID=="CPER" & taxon_name=="ectomycorrhizal")) +
    geom_line(aes(x = dates, y = med), show.legend = FALSE, na.rm = TRUE) +
    geom_ribbon(aes(x = dates, ymin = lo, ymax = hi), alpha=0.6, fill="blue", na.rm = TRUE) +
    geom_ribbon(aes(x = dates, ymin = lo_25, ymax = hi_75),fill="red", alpha=0.6, na.rm = TRUE) +
    theme_bw()+
    scale_fill_brewer(palette = "Paired") +
    theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
          legend.position = "bottom",legend.title = element_text(NULL),
          plot.margin = unit(c(.2, .2, 2, .2), "cm")) + 
    ylab(NULL) +
    facet_grid(plotID~model_name) +
    geom_point(aes(x = dates, y = as.numeric(truth))) + 
    xlab(NULL) + 
    labs(fill='')
  
  print(p1)
  cat("Single forecast inspection plot created successfully\n")
} else {
  cat("No data available for plotting\n")
}

