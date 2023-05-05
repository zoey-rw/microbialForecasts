# downloaded from Appears portal using NEON lat/long from site info CSV
source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")

modis_in <- read.csv("/projectnb/dietzelab/zrwerbin/microbialForecasts/data/clean/neon-allsites-MCD12Q2-006-results.csv")

modis_long <- modis_in %>% select(1:28) %>% 
	pivot_longer(cols = 6:28) %>% 
	filter(value != 32767) # remove "fill" values 

modis_long$Date = as.Date(modis_long$Date)
modis_long$year = lubridate::year(modis_long$Date)
modis_long$transition_date = as.Date(modis_long$value,
																		origin = "1970-01-01", tz="UTC")


modis_long$transition_date2 = as.Date(
	as.POSIXct(as.numeric(modis_long$transition_date)*24*60*60, origin = "1970-01-01", 
						 tz="UTC"))
modis_long$doy = lubridate::yday(modis_long$transition_date2)

modis_greenup = modis_long %>% 
	filter(name %in% c("MCD12Q2_006_MidGreenup_0","MCD12Q2_006_MidGreendown_0","MCD12Q2_006_Dormancy_0","MCD12Q2_006_Senescence_0","MCD12Q2_006_Greenup_0","MCD12Q2_006_Peak_0")) %>% 
	mutate(pheno_category = recode(name, "MCD12Q2_006_Dormancy_0" = "Dormancy",
																 "MCD12Q2_006_MidGreenup_0" = "MidGreenup",
																 "MCD12Q2_006_Peak_0" = "Peak",
																 "MCD12Q2_006_MidGreendown_0" = "MidGreendown",
																 "MCD12Q2_006_Senescence_0" = "Senescence",
																 "MCD12Q2_006_Greenup_0" = "Greenup"))

modis_greenup$transition_date2 <- ymd(modis_greenup$transition_date2)

#PROBLEM?

#modis_greenup[modis_greenup$ID=="GUAN",] %>% View
# https://judithtoms.wordpress.com/guanica/guanica-dry-forest/
# guanica phenology just makes no sense

pheno_categories = modis_greenup %>% 
	select(-c(Latitude, Longitude, Date, MODIS_Tile, name, transition_date, doy, value)) %>% 
	group_by(ID, year) %>% 
	pivot_wider(names_from = pheno_category, 
																		 values_from = transition_date2) 



# create averages per pixel in case a year's data is missing
site_pheno_avg = pheno_categories %>% 
	group_by(ID) %>% mutate_at(3:8, yday) %>% 
	mutate_at(3:8, mean) %>% 
	select(-year) %>% 
	distinct(.keep_all = T) %>% ungroup()

site_pheno_avg2 = site_pheno_avg %>% 
	mutate_at(2:7, function(x) as.Date(x, origin = "2014-01-01"))

completed = pheno_categories %>% ungroup() %>% 
	complete(year, ID) 
missing = completed %>% dplyr::anti_join(pheno_categories) %>% select(year, ID)
missing_replace = site_pheno_avg %>% filter(ID %in% missing$ID)
missing_replace2 = missing %>% left_join(missing_replace) 
missing_replace2 = missing_replace2 %>% 
	mutate_at(3:8, function(x) as.Date(x, origin = paste0(missing_replace2$`year`, "-01-01")))
missing_replace2$site_mean = TRUE
pheno_categories$site_mean = FALSE

pheno_categories <- rbind(pheno_categories, missing_replace2)

# Get greenup of subsequent year - used to create "dormancy" interval

pheno_categories = pheno_categories %>%
	group_by(ID) %>%
	mutate(subsequent_greenup = lead(Greenup),
				 previous_dormancy = lag(Dormancy))

#coalesced = coalesce(completed$Dormancy, site_pheno_avg$Dormancy)

# all(pheno_categories$Greenup < pheno_categories$Peak)
# all(pheno_categories$Peak < pheno_categories$Senescence)
# all(pheno_categories$Senescence < pheno_categories$Dormancy)
# all(pheno_categories$Dormancy < pheno_categories$Greenup)
pheno_categories = pheno_categories %>% 
	mutate(start_of_year = ymd(paste0(year, "-01-01")),
				 end_of_year = ymd(paste0(year, "-12-01")),
		greenup_interval = interval(Greenup, Peak),
				 peak_interval = interval(Peak, Senescence),
				 greendown_interval = interval(Senescence, Dormancy),
				 #dormancy_interval1 = interval(previous_dormancy, Greenup),
				 dormancy_interval = interval(Dormancy, subsequent_greenup))

pheno_categories_long = pheno_categories %>% select(-c(3:13)) %>%
pivot_longer(cols = 3:6)


pheno_categories_long2 = pheno_categories %>% select(c(1:11)) %>%
pivot_longer(cols = c(3:8,10:11)) %>%
mutate(arbitrary_greenness = recode(name, 
																		"start_of_year" = 0,
																		#"Dormancy" = 0, 
																		"Greenup" = .15, 
																		"subsequent_greenup" = .15,
																		"MidGreenup"= .5, 
																		"Peak" = 1, 
																		"Senescence" = .85,
																		"MidGreendown"= .5, 
																		"Dormancy" = 0.15, 
																		"previous_dormancy" = 0.15, 
																		"end_of_year" = 0))


pheno_categories_long2$phenophase = factor(pheno_categories_long2$name, 
																					 ordered = T, 
																					 levels = c("previous_dormancy","Greenup",
																					 					 "MidGreenup","Peak","MidGreendown",
																					 					 "Senescence","Dormancy","subsequent_greenup"))


saveRDS(list(pheno_categories, pheno_categories_long, pheno_categories_long2), 
				here("data/clean/modis_greenup.rds"))

library(ggrepel)



pheno_categories_long2$month = month(pheno_categories_long2$value)
labels = pheno_categories_long2 %>% ungroup %>%  
	select(phenophase, arbitrary_greenness, year, month, ID) %>% 
	filter(year==2018) %>% 
	select(-year) %>% 
	distinct() %>% filter(!phenophase %in% c("start_of_year","end_of_year","previous_dormancy","subsequent_greenup"))

ggplot(pheno_categories_long2 %>% filter(ID == "HARV")) +
	geom_point(aes(x = month, y = arbitrary_greenness, color = year), 
						 position=position_jitter(width = .1, height=.01), alpha=.5, size=5) + 
	theme_minimal() + 
	geom_label_repel(data = labels %>% filter(ID == "HARV"), 
									 aes(label = phenophase, y = arbitrary_greenness, x = month))


ggplot(pheno_categories_long2 %>% filter(ID == "CPER")) +
	geom_point(aes(x = month, y = arbitrary_greenness, color = year), 
						 position=position_jitter(width = .1, height=.01), alpha=.5, size=5) + 
	geom_smooth(aes(x = month, y = arbitrary_greenness), se=F) + 
	theme_minimal() + 
	geom_label_repel(data = labels %>% filter(ID == "CPER"), 
									 aes(label = phenophase, y = arbitrary_greenness, x = month), inherit.aes = F)
