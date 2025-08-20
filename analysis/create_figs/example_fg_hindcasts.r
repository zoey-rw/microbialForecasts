# Create example figure with one hindcast from each type
library(tidyverse)
source("source.R")

# Check if hindcast data exists
if (!file.exists(here("data/summary/all_hindcasts.rds"))) {
  stop("all_hindcasts.rds not found. Please regenerate this file from the model analysis pipeline.")
}

hindcast_in <- readRDS(here("data/summary/all_hindcasts.rds"))

not_na <- hindcast_in[!is.na(hindcast_in$truth),]
top_plots <- names(tail(sort(table(not_na$plotID)), 8))
hindcast_in <- hindcast_in[hindcast_in$plotID %in% top_plots,]

# Ecto at all Harv sites
all_cov_hindcast = hindcast_in %>% filter(model_name=="all_covariates" & site_prediction != "New time x site (random effect)")

# Create plot if data is available
if (nrow(all_cov_hindcast) > 0) {
  p1 <- ggplot(all_cov_hindcast %>% filter(siteID=="CPER" & taxon_name=="ectomycorrhizal")) +
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
  cat("Ectomycorrhizal plot created successfully\n")
} else {
  cat("No data available for ectomycorrhizal plotting\n")
}

library(ggforce)

# Check if we have data for functional group plots
if (exists("fg_top_plots") && nrow(fg_top_plots) > 0) {
  for(i in 1:2){
    p2 <- ggplot(fg_top_plots) +
      geom_line(aes(x = dates, y = med), show.legend = FALSE) +
      geom_ribbon(aes(x = dates, ymin = lo, ymax = hi, fill = fcast_period), alpha=0.6, na.rm = FALSE) +
      labs(title = "Example hindcasts (calibration: 2013-2016)") +
      theme_bw()+
      scale_fill_brewer(palette = "Paired") +
      theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
            legend.position = "bottom",legend.title = element_text(NULL),
            plot.margin = unit(c(.2, .2, 2, .2), "cm")) + 
      ylab(NULL) +
      geom_point(aes(x = dates, y = as.numeric(truth))) + 
      xlab(NULL) + 
      labs(fill='') +
      facet_grid_paginate(plotID~taxon, ncol = 1, nrow = 8, page = i)
    
    print(p2)
  }
} else {
  cat("fg_top_plots data not available for functional group plotting\n")
}

# Check if relative abundances data exists
if (file.exists(here("data/clean/ecto_rel_abundances.rds"))) {
  ecto_rel <- readRDS(here("data/clean/ecto_rel_abundances.rds"))
  ecto_rel <- cbind(ecto_rel, parseNEONsampleIDs(rownames(ecto_rel)))
  
  # Check if fg_hindcast data exists for merging
  if (exists("fg_hindcast") && nrow(fg_hindcast) > 0) {
    ecto_rel_merged <- left_join(ecto_rel, fg_hindcast[,c("truth","dateID","plotID","lo","med","hi","fcast_period")])
    
    p3 <- ggplot(ecto_rel_merged) +
      geom_ribbon(aes(x = dates, ymin = lo, ymax = hi, fill = fcast_period), alpha=0.6, na.rm = FALSE) +
      labs(title = "Raw data (soil cores)") +
      theme_bw()+
      scale_fill_brewer(palette = "Paired") +
      theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
            legend.position = "bottom",legend.title = element_text(NULL),
            plot.margin = unit(c(.2, .2, 2, .2), "cm")) + 
      ylab(NULL) +
      geom_jitter(aes(x = dates, y = as.numeric(ectomycorrhizal)), alpha = .5, size = 4) + 
      xlab(NULL) + 
      labs(fill='') +
      facet_grid(rows=vars(plotID), scales="free")
    
    print(p3)
    cat("Raw data plot created successfully\n")
  } else {
    cat("fg_hindcast data not available for merging\n")
  }
} else {
  cat("ecto_rel_abundances.rds not found. Data may need to be regenerated.\n")
}

# Ecto at all Harv sites
if (exists("fg_hindcast") && nrow(fg_hindcast) > 0) {
  p4 <- ggplot(fg_hindcast[fg_hindcast$siteID=="HARV" & fg_hindcast$taxon_name=="ectomycorrhizal",]) +
    geom_ribbon(aes(x = dates, ymin = lo, ymax = hi, fill = fcast_period), alpha=0.6, na.rm = FALSE) +
    labs(title = "Raw data (soil cores)") +
    theme_bw()+
    scale_x_date(breaks = scales::date_breaks("4 month"),date_labels = "%b") +
    scale_fill_brewer(palette = "Paired") +
    theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
          legend.position = "bottom",legend.title = element_text(NULL),
          plot.margin = unit(c(.2, .2, 2, .2), "cm"),
          axis.text.x = element_text(angle = -20)) + 
    ylab(NULL) +
    geom_jitter(aes(x = dates, y = as.numeric(truth)), alpha = .5, size = 4) + 
    xlab(NULL) + 
    labs(fill='') +
    facet_grid_paginate(ncol = 1, plotID~., nrow=5, page = 1)
  
  print(p4)
  cat("HARV ectomycorrhizal plot created successfully\n")
} else {
  cat("fg_hindcast data not available for HARV plotting\n")
}








library(lubridate)

plots_keep <- c("HARV_001","HARV_004","HARV_013","HARV_016")
plots_keep <- c("STER_006","STER_005","CPER_002","CPER_003","CPER_001","CPER_004")
read_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_ITS_2021.rds")
ecto_rel <- read_in$ectomycorrhizal
ecto_rel$doy <- lubridate::yday(ecto_rel$dates)
ecto_rel$year <- lubridate::year(ecto_rel$dates)
ggplot(ecto_rel[ecto_rel$plotID %in% plots_keep,]) +
	#	facet_grid(rows=vars(example), drop=T, scales="free") +
	#geom_line(aes(x = dates, y = ectomycorrhizal), show.legend = F) +
	labs(title = "ectomycorrhizal abundances") +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	geom_point(aes(x = dates, y = as.numeric(ectomycorrhizal))) + xlab(NULL) + labs(fill='') +
	facet_grid_paginate(#cols=vars(taxon), rows=vars(plotID),
		~plotID, scales="free",
		ncol = 1, nrow = 8, page = 2)



# Do ecto abundances peak in winter for specific plots??
ggplot(ecto_rel[ecto_rel$plotID %in% plots_keep,]) +
	#	facet_grid(rows=vars(example), drop=T, scales="free") +
	#geom_line(aes(x = dates, y = ectomycorrhizal), show.legend = F) +
	labs(title = "ectomycorrhizal abundances") +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	geom_point(aes(x = doy, y = as.numeric(ectomycorrhizal),
								 color = as.factor(year)), size = 3) + xlab(NULL) + labs(fill='') +
	geom_smooth(aes(x = doy, y = as.numeric(ectomycorrhizal))) + xlab(NULL) + labs(fill='') +
	facet_grid(#cols=vars(taxon),
		rows=vars(plotID),
		#~plotID,
		scales="free")





ggplot(ecto_rel[ecto_rel$siteID %in% c("HARV","STER","CPER","DSNY","OSBS"),]) +
	labs(title = "ectomycorrhizal abundances") +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	geom_point(aes(x = doy, y = as.numeric(ectomycorrhizal),
								 color = as.factor(year)), size = 3) + xlab(NULL) + labs(fill='') +
	geom_smooth(aes(x = doy, y = as.numeric(ectomycorrhizal))) + xlab(NULL) + labs(fill='') + facet_grid(year~siteID)
	# facet_grid(#cols=vars(taxon),
	# 	rows=vars(plotID),
	# 	#~plotID,
	# 	scales="free")



cal <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_ITS_2021.rds")
val <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/val_groupAbundances_ITS_2021.rds")

ecto_rel <- rbind(cal$ectomycorrhizal,val$ectomycorrhizal)
ecto_rel$doy <- lubridate::yday(ecto_rel$dates)
ecto_rel$year <- lubridate::year(ecto_rel$dates)
ggplot(ecto_rel[ecto_rel$siteID %in% c("HARV","STER","CPER","DSNY","OSBS"),]) +
	labs(title = "ectomycorrhizal abundances") +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	geom_point(aes(x = doy, y = as.numeric(ectomycorrhizal),
								 color = as.factor(year)), size = 3) + xlab(NULL) + labs(fill='') +
	geom_smooth(aes(x = doy, y = as.numeric(ectomycorrhizal))) + xlab(NULL) + labs(fill='')



ecto_hindcast <- fg_hindcast %>% filter(taxon=="ectomycorrhizal")

ggplot(ecto_hindcast[ecto_hindcast$plotID %in% plots_keep ,]) +
	#	facet_grid(rows=vars(example), drop=T, scales="free") +
	geom_line(aes(x = dates, y = med), show.legend = F) +
	geom_ribbon(aes(x = dates, ymin = lo, ymax = hi, fill = fcast_period), alpha=0.6,  na.rm = F) +
	labs(title = "Example hindcasts (calibration: 2013-2016)") +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='') +
	facet_grid(#cols=vars(taxon), rows=vars(plotID),
		plotID~taxon)
