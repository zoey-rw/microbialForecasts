# Look at relationships between Pinaceae, soil fungi, and soil bacteria at NEON sites

# Load R packages, install if necessary.
#devtools::install_github("admahood/neondiversity") 
library(neondiversity) 
library(tidyverse)
library(zoo)
library(corrplot) 


usda_list <- read_csv("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/reference_data/NEON_pla_taxonomy.csv")

##### 1. Plant data #####
# Download data from NEON sites (the ones I have microbial data for)
# plant_info <- download_plant_div(sites = "all") 
# saveRDS(plant_info, "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEON_plant_data.rds")
plant_info <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEON_plant_data.rds")

# Calculate a bunch of metrics for Poaceae. Documentation for this is here: https://github.com/admahood/neondiversity
# Takes a couple minutes
plant_stats <- get_diversity_info(plant_info, families = "Poaceae", spp = "Poa pratensis L.")

# Let's subset:
plant_subset <- plant_stats %>% 
  dplyr::select(plotID, year, rc_Poaceae, nspp_total, rel_cover_exotic) %>% 
  mutate(siteID = substr(plotID, 1, 4)) %>% 
  mutate(year_date = as.Date(paste0(year, "-01-1"), "%Y-%m-%d"))



#### Add in % cover of mycorrhizal plants ----
# Table S2 from FungalRoot
fungalroot_genus <- read.csv("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/reference_data/fungalroot_genera.csv") 
fungalroot_genus <- fungalroot_genus %>% mutate(myco_type = recode(Mycorrhizal.type, 
																											 "EcM-AM" = "EcM")) %>% select(genus=Genus, myco_type)
## GET ECTOMYCORRHIZAL FAMILY LIST ##
# File is from Fungal Root preprint (2019)
old_fungalroot <- readxl::read_xlsx("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/reference_data/media-3.xlsx")
# This file is from GBIF - to add family info
fungalroot_spec <- read.csv("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/reference_data/fungalroot_species.csv", sep="\t")
fungalroot <- old_fungalroot %>% select(species, myco_type = `mycorrhiza type`)
fungalroot$genus <- word(fungalroot$species)
fungalroot <- left_join(fungalroot, fungalroot_spec[,c("genus","family")], by="genus")
fungalroot <- fungalroot %>% mutate(myco_type = recode(myco_type, "EcM; no others" = "EcM",
																											 "EcM; others not addressed" = "EcM",
																											 "EcM,AM" = "EcM",
																											 "ErM,EcM" = "EcM"))
EM_families <- fungalroot %>% group_by(family) %>% count(myco_type) %>% group_by(family) %>% mutate(per= prop.table(n) * 100) %>% filter(myco_type == "EcM" & per > 60) %>% select(family)


fungalroot_distinct <- fungalroot %>% distinct(.keep_all=T)

# Calculate % cover using neondiversity package
full_on_cover <- get_longform_cover(plant_info)

# for name_cleaner() to work, source the raw code:
source("https://raw.githubusercontent.com/admahood/neondiveRsity/master/R/diversity_data_prep.R")
full_on_cover <- name_cleaner(full_on_cover)
full_on_cover$genus <- word(full_on_cover$scientificName)

# Add genus when missing
full_on_cover$usda_genus <- usda_list[match(full_on_cover$taxonID, usda_list$taxonID),]$genus
full_on_cover$genus <- ifelse(is.na(full_on_cover$genus), full_on_cover$usda_genus, full_on_cover$genus)


# Add genus- myco associations
cover_with_myco <- left_join(full_on_cover, fungalroot_genus, by=c("genus"))
# Add family- myco associations
cover_with_myco[which(cover_with_myco$family %in% EM_families$family),]$myco_type <- "EcM"

myco <- cover_with_myco %>%
	dplyr::filter(!is.na(myco_type)) %>%
	dplyr::group_by(site, plotID, year) %>%
	#dplyr::group_by(site, plotID) %>%
	dplyr::mutate(total_cover = sum(cover))%>%
	dplyr::ungroup() %>%
	dplyr::group_by(site, plotID, subplotID,year, myco_type) %>%
	#dplyr::group_by(site, plotID, subplotID, myco_type) %>%
	dplyr::summarise(cover = sum(cover),
									 total_cover = first(total_cover)) %>%
	dplyr::ungroup() %>%
	dplyr::mutate(rc_EM = cover/total_cover) %>%
	dplyr::ungroup()

myco_relEM <- myco %>% filter(myco_type %in% c("EcM"))





max_myco <- cover_with_myco %>%
	dplyr::filter(!is.na(myco_type)) %>%
	dplyr::group_by(site, plotID, year) %>%
	#dplyr::group_by(site, plotID) %>%
	dplyr::mutate(total_cover = sum(cover))%>%
	dplyr::ungroup() %>%
	dplyr::group_by(site, plotID, subplotID,year, myco_type) %>%
	#dplyr::group_by(site, plotID, subplotID, myco_type) %>%
	dplyr::summarise(cover = sum(cover),
									 total_cover = first(total_cover)) %>%
	dplyr::ungroup() %>%
	dplyr::mutate(rc_EM = cover/total_cover) %>%
	dplyr::group_by(site, plotID, myco_type) %>%
	filter(myco_type == "EcM") %>% 
	filter(rc_EM == max(rc_EM, na.rm=T)) %>% 
	dplyr::ungroup()
myco_relEM <- myco %>% filter(myco_type %in% c("EcM"))


template <- full_on_cover %>%
	dplyr::select(site, plotID, year) %>% distinct(.keep_all=T)
final_table <- template %>%
	dplyr::left_join(myco_relEM, by = c("site", "plotID", "year")) %>% select(-c(subplotID, myco_type, total_cover))
final_table$rc_EM <- ifelse(is.na(final_table$rc_EM), 0, final_table$rc_EM)

plant_subset_save <- plant_subset
plant_subset <- left_join(plant_subset_save, final_table)


microbe_plots <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/microbe_plot_list.rds")
# Get all plotID x year combinations
poss_site_dates <- plant_subset %>% tidyr::expand(nesting(siteID, year_date))
poss_site_plots <- plant_subset %>% tidyr::expand(nesting(siteID, plotID))
poss_site_plots <- merge(poss_site_plots, microbe_plots, all=T)
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


# let's also use site-level means to fill in plots unique to our microbial dataset (only applies for a handful)
# This code is super clunky and dumb but w/e
# microbes1 <- rbind(readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_16S_2021.rds")[[1]],
# 									 readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/val_groupAbundances_16S_2021.rds")[[1]])
# microbes2 <- rbind(readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_ITS_2021.rds")[[1]],
# 									 readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/val_groupAbundances_ITS_2021.rds")[[1]])
# microbes3 <- do.call(rbind.data.frame, readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_16S.rds")) %>% select(plotID, siteID)
# microbes4 <- do.call(rbind.data.frame, readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_ITS.rds")) %>% select(plotID, siteID)
# microbe_plots <- unique(rbind(microbes3, microbes4))
# microbe_plots1 <- unique(data.frame(plotID = substr(rownames(microbes1), 1, 8),
# 														siteID = substr(rownames(microbes1), 1, 4)))
# microbe_plots2 <- unique(data.frame(plotID = substr(rownames(microbes2), 1, 8),
# 																		siteID = substr(rownames(microbes2), 1, 4)))
# microbe_plots <- unique(rbind(microbe_plots, microbe_plots1))
# microbe_plots <- unique(rbind(microbe_plots, microbe_plots2))
# saveRDS(microbe_plots, "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/microbe_plot_list.rds")


#plant_fill_site <- merge(plant_subset_fill, microbe_plots, all=T)
plant_fill_site <- plant_subset_fill %>% group_by(siteID) %>% mutate(site_mean_nspp_total = mean(nspp_total, na.rm=T),
																											site_mean_rel_cover_exotic = mean(rel_cover_exotic, na.rm=T),
																											site_mean_rc_Poaceae = mean(rc_Poaceae, na.rm=T),
																											site_mean_rc_EM = mean(rc_EM, na.rm=T)) %>% 
	mutate(nspp_total = ifelse(is.na(nspp_total), site_mean_nspp_total, nspp_total),
				 rel_cover_exotic = ifelse(is.na(rel_cover_exotic), site_mean_rel_cover_exotic, rel_cover_exotic),
				 rc_Poaceae = ifelse(is.na(rc_Poaceae), site_mean_rc_Poaceae, rc_Poaceae),
				 rc_EM = ifelse(is.na(rc_EM), site_mean_rc_EM, rc_EM)) %>% ungroup()
plant_fill_site[which(is.na(plant_fill_site$year_date)),]$year_date <- as.Date("2013-06-01")
				 #year_date = ifelse(is.na(year_date), as.Date("2013-06-01", ), year_date)) %>% ungroup()

# Save original values
plant_fill_site$nspp_total_out <- plant_fill_site$nspp_total
plant_fill_site$rc_Poaceae_out <- plant_fill_site$rc_Poaceae
plant_fill_site$rc_exotic_out <- plant_fill_site$rel_cover_exotic
plant_fill_site$rc_EM_out <- plant_fill_site$rc_EM

# Mean center and scale data
plant_fill_site$nspp_total <- scale(plant_fill_site$nspp_total_out)
plant_fill_site$rc_Poaceae <- scale(plant_fill_site$rc_Poaceae_out)
plant_fill_site$rc_exotic <- scale(plant_fill_site$rc_exotic_out)
plant_fill_site$rc_EM <- scale(plant_fill_site$rc_EM_out)


plant_fill_site[is.na(plant_fill_site$gap_filled),]$gap_filled <- T
saveRDS(plant_fill_site, "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/annual_plant_data.rds")





# VISUALIIIIZE
dsny_fill <- plant_fill_site[which(plant_fill_site$siteID == "DSNY"),]

p1 <- ggplot(dsny_fill) + geom_point(aes(x = year_date, y = nspp_total_out, shape = gap_filled)) + 
	facet_wrap(~plotID) + 
	labs(title = "Plant species richness, gap-filled with subsequent/previous observations", 
			 subtitle = "Example from one NEON site (DSNY)",
			 x = NULL, y = "Plant species richness", shape = "Gap-filled?") + 
	#scale_color_manual(values=c("#F8766D", "#00BA38"))
	scale_shape_manual(values = c(19, 4)) + theme (legend.position = c(0.7, 0.13), text = element_text(size = 16)) 
p1
p2 <- ggplot(dsny_fill) + geom_point(aes(x = year_date, y = rc_Poaceae_out, shape = gap_filled)) + 
	facet_wrap(~plotID) + 
	labs(title = "Plant species richness, gap-filled with subsequent/previous observations", 
			 subtitle = "Example from one NEON site (DSNY)",
			 x = "year", y = "Relative cover of Poacae (grasses(", shape = "Gap-filled?") + 
	#scale_color_manual(values=c("#F8766D", "#00BA38"))
	scale_shape_manual(values = c(19, 4)) + theme (legend.position = c(0.7, 0.04), text = element_text(size = 16)) 
p2
p3 <- ggplot(dsny_fill) + geom_point(aes(x = year_date, y = rc_exotic_out, shape = gap_filled)) + 
	facet_wrap(~plotID) + 
	labs(title = "Plant species richness, gap-filled with subsequent/previous observations", 
			 subtitle = "Example from one NEON site (DSNY)",
			 x = "year", y = "Relative cover of invasive species", shape = "Gap-filled?") + 
	#scale_color_manual(values=c("#F8766D", "#00BA38"))
	scale_shape_manual(values = c(19, 4)) + theme (legend.position = c(0.7, 0.04), text = element_text(size = 16)) 
p3
p4 <- ggplot(dsny_fill) + geom_point(aes(x = year_date, y = rc_EM_out, shape = gap_filled)) + 
	facet_wrap(~plotID) + 
	labs(title = "Plant species richness, gap-filled with subsequent/previous observations", 
			 subtitle = "Example from one NEON site (DSNY)",
			 x = "year", y = "Relative cover of invasive species", shape=NULL) + 
	#scale_color_manual(values=c("#F8766D", "#00BA38"))
	scale_shape_manual(values = c(19, 4), labels =c("Observed","Gap-filled")) + theme (legend.position = c(0.7, 0.04), text = element_text(size = 16)) 
p4

#ggsave(p, filename = "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/figures/plant_diversity_gapfill.png", device = "png", width = 12, height = 12, units = "in")




merged <- merge(plant_subset, basal_relEM)
fungal_abun$year <- substr(fungal_abun$dateID, 1, 4)
merged <- merge(merged, fungal_abun)
myco_merged <- merged %>% mutate(rel_cover = rc_EM)
dsny_fill <- myco_merged[which(myco_merged$siteID == "DSNY"),]



a <- ggscatter(myco_merged, x = "rel_cover", y = "relEM",
							 add = "reg.line")  + 
	xlab("Relative EM plant cover") + ylab("Relative EM tree basal area") + 
	stat_cor(label.y = 1.13) +
	stat_regline_equation(
		aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), label.y =1.05
	)
b <- ggscatter(myco_merged, x = "rel_cover", y = "ectomycorrhizal",
							 add = "reg.line")  + 
	xlab("Relative EM plant cover") + ylab("Ectomycorrhizal fungi") + 
	stat_cor(label.y = 1.1) +
	
	stat_regline_equation(
		aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), label.y =1)
c <- ggscatter(myco_merged, x = "relEM", y = "ectomycorrhizal",
							 add = "reg.line")  + 
	xlab("Relative EM tree basal area") + ylab("Ectomycorrhizal fungi") + 
	stat_cor(label.y = 1.1) +
	stat_regline_equation(
		aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), label.y =1)
ggpubr::ggarrange(a, b, c, nrow = 1, labels = c("A","B","C"))

