# Look at relationships between Pinaceae, soil fungi, and soil bacteria at NEON sites

# Load R packages, install if necessary.
#devtools::install_github("admahood/neondiversity") 
library(neondiversity) 
library(tidyverse)
library(zoo)
library(corrplot) 

##### 1. Plant data #####
# Download data from NEON sites (the ones I have microbial data for)
# plant_info <- download_plant_div(sites = "all") 
# saveRDS(plant_info, "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEON_plant_data.rds")
plant_info <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEON_plant_data.rds")

# Calculate a bunch of metrics for Poaceae. Documentation for this is here: https://github.com/admahood/neondiversity
# Takes a couple minutes
plant_stats <- get_diversity_info(plant_info, families = "Poaceae", spp = "Poa pratensis L.")

# Let's check out some of this plant data 
# Histogram of species richness
hist(plant_stats$nspp_total)
# Histogram of relative cover (RC) of grasses (Poaceae)
hist(plant_stats$rc_Poaceae)
# Histogram of relative cover (RC) of exotic species
hist(plant_stats$rel_cover_exotic)

# Remove non-numeric columns and all-zero columns, just so we can look at strong correlations
plant_numeric <- plant_stats %>% dplyr::select_if(function(col) is.numeric(col) & any(col != 0))
# Create correlation matrix
cor_mat <- cor(plant_numeric)
# Subset to strong correlations
strong <- abs(cor_mat[1,]) > .2
cor_mat_strong <- cor_mat[strong, strong]
# Visualize correlations
corrplot::corrplot(cor_mat_strong, method="circle", type = "lower")

# Based on the literature, we want to include metrics of plant species richness, cover of grasses (as a functional group),
# and invasive species status.

# There is no strong correlation between relative cover of Poaceae (grasses) and total number of species
# There is a weak correlation between relative cover of Poaceae (grasses) and relative cover of exotic species
# So let's subset accordingly:
plant_subset <- plant_stats %>% 
  dplyr::select(plotID, year, rc_Poaceae, nspp_total, rel_cover_exotic) %>% 
  mutate(siteID = substr(plotID, 1, 4)) %>% 
  mutate(year_date = as.Date(paste0(year, "-01-1"), "%Y-%m-%d"))

# Get all plotID x year combinations
poss_site_dates <- plant_subset %>% tidyr::expand(nesting(siteID, year_date))
poss_site_plots <- plant_subset %>% tidyr::expand(nesting(siteID, plotID))
poss_plot_dates <- merge(poss_site_dates, poss_site_plots)

# Let's gap-fill the data so that a plot retains its values until a new measurement is taken.
# If there are no measurements at the beginning of the dataset, 
# we fill using first measurement ("downup" .direction does this)
plant_subset <- merge(plant_subset, poss_plot_dates, all = T)
plant_subset_fill <- plant_subset  %>% 
  group_by(plotID) %>% padr::pad(interval = "year") %>% #tidyr::fill(.direction = "updown") %>% 
  mutate_all(funs(na.locf(., na.rm = FALSE))) %>% 
mutate_all(funs(na.locf(., na.rm = FALSE, fromLast = TRUE))) 
  
#na.locf(fromLast = TRUE)

# before gap-filling
dsny <- plant_subset[which(plant_subset$siteID == "DSNY"),]
ggplot(dsny) + geom_point(aes(x = year_date, y = nspp_total)) + facet_wrap(~plotID) 

# to mark which points were gap-filled, we'll identify those with incorrect years (since na.locf filled those too)
plant_subset_fill$real_year <- substr(plant_subset_fill$year_date, 1, 4)
plant_subset_fill$gap_filled <- ifelse(plant_subset_fill$real_year != plant_subset_fill$year, TRUE, FALSE)
dsny_fill <- plant_subset_fill[which(plant_subset_fill$siteID == "DSNY"),]
p <- ggplot(dsny_fill) + geom_point(aes(x = year_date, y = nspp_total, shape = gap_filled)) + 
  facet_wrap(~plotID) + 
  labs(title = "Plant species richness, gap-filled with subsequent/previous observations", 
      subtitle = "Example from one NEON site (DSNY)",
       x = "year", y = "Plant species richness", shape = "Gap-filled?") + 
  #scale_color_manual(values=c("#F8766D", "#00BA38"))
  scale_shape_manual(values = c(19, 4)) + theme (legend.position = c(0.7, 0.04), text = element_text(size = 16)) 
p
ggsave(p, filename = "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/figures/plant_diversity_gapfill.png", device = "png", width = 12, height = 12, units = "in")

saveRDS(plant_subset_fill, "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/annual_plant_data.rds")
