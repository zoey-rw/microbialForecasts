#### THIS SCRIPT WRITTEN BY STEVE GOUGHERTY

## call relevant packages
library(daymetr)
library(dplyr)

sites_locs <- read.csv("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data_construction/field-sites.csv")

output.path.monthly <- "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/daymet_monthly.rds"
output.path.weekly <- "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/daymet_weekly.rds"
output.path.daily <- "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/daymet_daily.rds"

sites_locs$lat <- unlist(lapply(strsplit(as.character(sites_locs$`Lat..Long.`), split = ",  ", fixed=T), "[[", 1))
sites_locs$lon <- unlist(lapply(strsplit(as.character(sites_locs$`Lat..Long.`), split = ",  ", fixed=T), "[[", 2))
sites_locs$MAT <- unlist(lapply(strsplit(as.character(sites_locs$Mean.Annual.Temperature), split = "/", fixed=T), "[[", 2))
sites_locs$MAT <- gsub("F$","", sites_locs$MAT)
sites_locs$MAP <- gsub(" mm$", "", sites_locs$Mean.Annual.Precipitation)

locs <- sites_locs %>% filter(!(Site.Type == "Relocatable Aquatic" | Site.Type == "Core Aquatic")) %>% dplyr::select(siteID = Site.ID, lat, lon, MAT, MAP) 


# #WRITE DATA TO CSV FILE
write.table(locs, paste0(tempdir(),"/locations.csv"),
            sep=",",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)

#download the data here. 
daymet_data <- download_daymet_batch(file_location = paste0(tempdir(),
                                                            "/locations.csv"),
                                     start=2013,
                                     end=2019,
                                     internal=FALSE,
                                     path =  tempdir(),
                                     silent=FALSE)

file.remove(paste0(tempdir(), "/locations.csv"))
temp = list.files(tempdir(), pattern="*.csv", full.names = T)
#temp = list.files("/scratch/Rtmpf4GoGH", pattern="*.csv", full.names = T)
site_names <- substr(basename(temp), 1, 4)

alldfs <-list()
for (i in 1:length(temp)){
  df <- read_daymet(temp[[i]], site = site_names[[i]])
  alldfs[[i]] <- df
}
all_daymet <- bind_rows(alldfs)
days <- format(strptime(all_daymet$yday, format="%j"), format="%m-%d")
dayyear <- paste(days, all_daymet$year, sep="-")
all_daymet$dayyear <-dayyear

all_daymet$Date <- as.Date(all_daymet$dayyear, format = "%m-%d-%Y")
all_daymet$dateID <- gsub("-","", substr(all_daymet$Date, 1, 7))

# Already summarized by day, just clean up
daymet.daily <- all_daymet %>% 
  select(siteID = site, dateID, Date, variable = measurement, value) %>% 
  filter(variable %in% c("prcp..mm.day.", "tmin..deg.c.", "tmax..deg.c.")) %>% 
  mutate(variable = recode(variable,
                           "prcp..mm.day." = "precip",
                           "tmin..deg.c." = "mintemp",
                           "tmax..deg.c." = "maxtemp")) %>% 
  pivot_wider(id_cols = c(siteID, dateID, Date), 
              names_from = variable, 
              values_from = value,
              values_fn = {mean})
# SAVE
saveRDS(daymet.daily, output.path.daily)


# Summarize by week
all_daymet$week <- lubridate::week(ymd(as.Date(all_daymet$Date)))
all_daymet$year <- lubridate::year(ymd(as.Date(all_daymet$Date)))
all_daymet$year_week <- paste0(all_daymet$year, "_", all_daymet$week)

daymet.weekly <- all_daymet  %>% 
  select(siteID = site, dateID, year_week, variable = measurement, value) %>% 
  filter(variable %in% c("prcp..mm.day.", "tmin..deg.c.", "tmax..deg.c.")) %>% 
  mutate(variable = recode(variable,
                           "prcp..mm.day." = "precip",
                           "tmin..deg.c." = "mintemp",
                           "tmax..deg.c." = "maxtemp")) %>% 
  group_by(siteID, year_week) %>% 
  # mutate(across(c(value, maxtemp, mintemp, mean.of.max.min), ~mean(.x, na.rm = TRUE), .names = "mean_{.col}")) %>% distinct(year_week, siteID, .keep_all = T) %>%  data.frame()
  pivot_wider(id_cols = c(siteID, dateID, year_week), 
              names_from = variable, 
              values_from = value,
              values_fn = {mean})
saveRDS(daymet.weekly, output.path.weekly)



# Summarize by month
daymet.monthly <- all_daymet %>% 
  select(siteID = site, dateID, variable = measurement, value) %>% 
  filter(variable %in% c("prcp..mm.day.", "tmin..deg.c.", "tmax..deg.c.") ) %>% 
  group_by(siteID, dateID, variable) %>% summarise(mean = mean(value, na.rm=T)) %>% 
  mutate(variable = recode(variable,
                           "prcp..mm.day." = "precip",
                           "tmin..deg.c." = "mintemp",
                           "tmax..deg.c." = "maxtemp")) %>% 
  pivot_wider(id_cols = c(siteID, dateID), 
              names_from = variable, 
              values_from = c(mean), 
              names_glue = "{variable}_{.value}")
# SAVE
saveRDS(daymet_out, output.path)
