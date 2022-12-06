# Creates the "site_effect_predictors.rds" file used to analyze model results


source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
pacman::p_load(stringr)


## Get site climate data from NEON
fieldsites_raw <- read.csv("https://www.neonscience.org/sites/default/files/NEON_Field_Site_Metadata_20220412.csv")
fieldsites <- fieldsites_raw %>% filter(!grepl("Aquatic", field_site_type)) %>%
	select(siteID = field_site_id,
				 MAT = field_mean_annual_temperature_C,
				 MAP = field_mean_annual_precipitation_mm,
				 nlcd = field_dominant_nlcd_classes,
				 #ecoregion = field_domain_id)
				 ecoregion = field_domain_id)

fieldsites$evergreen   <- ifelse(grepl("Evergreen", fieldsites$nlcd), "Evergreen", "")
fieldsites$deciduous   <- ifelse(grepl("Deciduous", fieldsites$nlcd), "Deciduous", "")
fieldsites$wetlands 	 <- ifelse(grepl("Wetlands", fieldsites$nlcd), "Wetlands", "")
fieldsites$herbaceous <- ifelse(grepl("Herbaceous", fieldsites$nlcd), "Herbaceous", "")
fieldsites$shrub      <- ifelse(grepl("Grassland|Shrub", fieldsites$nlcd), "Shrub", "")
fieldsites$crops      <- ifelse(grepl("Crops", fieldsites$nlcd), "Crops", "")
fieldsites$nlcd_clean <- apply(fieldsites[,c(6:11) ], 1, paste0, collapse = "")

# fieldsites$evergreen   <- ifelse(grepl("Evergreen", fieldsites$nlcd), 1, 0)
# fieldsites$deciduous   <- ifelse(grepl("Deciduous", fieldsites$nlcd), 1, 0)
# fieldsites$wetlands 	 <- ifelse(grepl("Wetlands", fieldsites$nlcd), 1, 0)
# fieldsites$herbaceous <- ifelse(grepl("Crops", fieldsites$nlcd), 1, 0)
# fieldsites$shrub      <- ifelse(grepl("Grassland|Shrub", fieldsites$nlcd), 1, 0)
# fieldsites$crops      <- ifelse(grepl("Crops", fieldsites$nlcd), 1, 0)
#
# classes <- expand(fieldsites, nesting(evergreen,
# 																			deciduous,
# 																			wetlands,
# 																			herbaceous,
# 																			shrub,
# 																			crops))
# classes1 <- classes %>%
# 	mutate(across(1:6, function(x)
# 		factor(x, c(0, 1),  c("", cur_column()))))
# classes1 <- lapply(classes1[,c(1:6)], as.character) %>% as.data.frame()
# classes1$cat  <- apply(classes1[,c(1:6) ], 1, paste0, collapse = "")
#
#
# fieldsites$nlcd_clean <- ifelse(grepl("Crops", fieldsites$nlcd), "Crops", fieldsites$nlcd)
# fieldsites$nlcd_clean <- ifelse(grepl("Evergreen Forest", fieldsites$nlcd), "Evergreen Forest", fieldsites$nlcd_clean)
# fieldsites$nlcd_clean <- ifelse(grepl("Deciduous|Mixed", fieldsites$nlcd), "Deciduous/Mixed Forest", fieldsites$nlcd_clean)
# fieldsites$nlcd_clean <- ifelse(grepl("Grassland|Scrub", fieldsites$nlcd), "Grassland/Scrub", fieldsites$nlcd_clean)
# fieldsites$nlcd_clean <- ifelse(grepl("Wetlands", fieldsites$nlcd), "Wetlands", fieldsites$nlcd_clean)


### Read in bulk density data
bd_raw <- readRDS(here("data/raw/bulkDensity_allsites.rds"))
bd <- bd_raw %>% filter(bulkDensTopDepth < 25) %>%
	group_by(siteID) %>%
	summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) #%>%
bd$bulkDens <- ifelse(is.na(bd$bulkDensThirdBar), bd$bulkDensFieldMoist, bd$bulkDensThirdBar)
bd <- bd %>% select(siteID, bulkDens)

### Read in micronutrient data
bgc_raw <- readRDS(here("data/raw/dp1.10047.rds"))
# remove deep horizons.
bgc <- bgc_raw %>% filter(biogeoTopDepth < 25)
most_relevant <- c("siteID", "plotID",
									 #"biogeoTopDepth", "biogeoBottomDepth",
									 #"gypsumConc",
									 #"caco3Conc",
									 "caNh4d",
									 #"caSatx",
									 "kNh4d",
									 #"kSatx",
									 "mgNh4d", #"mgSatx",
									 "naNh4d",
									 #"naSatx",
									 "cecdNh4",
									 #"alSatCecd33",
									 "alKcl",
									 "feKcl", "feOxalate",
									 #"mnKcl",
									 "mnOxalate",
									 "pOxalate", "siOxalate",
									 "estimatedOC", "sulfurTot",
									 #"pSatx",
									 #"ececCecd33",
									 "no2Satx", "no3Satx", "so4Satx",
									 "OlsenPExtractable", "MehlichIIITotP", "Bray1PExtractable")
bgc_relevant <- bgc %>% select(!!!most_relevant)
bgc_relevant$totalP <- ifelse(is.na(bgc_relevant$MehlichIIITotP),
															bgc_relevant$OlsenPExtractable,
															bgc_relevant$MehlichIIITotP)
bgc_relevant[,c("MehlichIIITotP","OlsenPExtractable","Bray1PExtractable")] <-NULL
bgc_relevant <- bgc_relevant %>% group_by(siteID) %>%
	summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

# Merge predictor dfs together
df_predictors <- fieldsites %>%
	merge(bgc_relevant, all=T) %>%
	merge(bd, all=T) %>%
	# Center and scale predictor values
	mutate(across(where(is.numeric), function(x) as.vector(scale(x))))

df_predictors <- merge(df_predictors, fieldsites_raw %>% select(siteID = field_site_id, latitude = field_latitude))

saveRDS(df_predictors, here("data/summary/site_effect_predictors.rds"))
