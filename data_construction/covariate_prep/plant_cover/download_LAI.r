pacman::p_load(tidyverse)

get_LAI <- function(
	siteID = "NIWO",
	numeric_site_id = 1,
	start_year = 2013,
	end_year = 2020,
	out_dir = '/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/modis/',
	ncores=10,
	field_site_info = con
){
	#loading libraries
	#library(PEcAn.visualization)
	library(PEcAn.assim.sequential)
	library(nimble)
	library(lubridate)
	library(rgdal) 
	library(ncdf4) 
	library(purrr)
	library(listviewer)
	library(dplyr)
	library(furrr)
	library(tictoc)
	
	#modis code
	source("/projectnb/dietzelab/dongchen/Multi-site/download_500_sites/call_MODIS.R")
	cat(paste0("Downloading for site ", siteID, "\n"))
	
	#convert year to YEARDOY
	start_YEARDOY <- paste0(as.character(start_year), "001")
	end_YEARDOY <- paste0(as.character(end_year), "365")
	
	qry_results <- con[con$field_site_id == siteID,]
	site_info <- list(site_id=numeric_site_id, site_name=qry_results$field_site_id, 
										lat=qry_results$field_latitude,
										lon=qry_results$field_longitude)
	
	#download LAI
	lai = call_MODIS(outdir = NULL, var = "LAI", site_info = site_info, product_dates = c(start_YEARDOY, end_YEARDOY),
									 run_parallel = TRUE, ncores = ncores, product = "MOD15A2H", band = "Lai_500m",
									 package_method = "MODISTools", QC_filter = TRUE, progress = TRUE)
	lai_data <- lai
	
	sd = call_MODIS(outdir = NULL, var = "LAI", site_info = site_info, product_dates = c(start_YEARDOY, end_YEARDOY),
									run_parallel = TRUE, ncores = ncores, product = "MOD15A2H", band = "LaiStdDev_500m",
									package_method = "MODISTools", QC_filter = TRUE, progress = TRUE)
	lai_sd <- sd
	#export LAI data
	names(lai_sd) = c("modis_date", "calendar_date", "band", "tile", "site_id", "lat", "lon", "pixels", "sd", "qc")
	output = cbind(lai_data, lai_sd$sd)
	names(output) = c(names(lai_data), "sd")
	save(output, file = file.path(out_dir, paste0(siteID, ".rds")))#export all LAI data
	cat(paste0("Download complete for site ", siteID, "\n"))
	return(output)
}

results <- Download_LAI(settings, start_year = 2013, 
												end_year = 2020, 
												out_dir = '/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/modis/')


con <- read.csv("https://www.neonscience.org/sites/default/files/NEON_Field_Site_Metadata_20210928.csv")
con <- con %>% dplyr::filter(field_site_type %in% c("Relocatable Terrestrial", "Core Terrestrial", "Gradient Terrestrial")) 
for (s in 27:length(con$field_site_id)){
	siteID <- con[s,]$field_site_id
	if(siteID=="NIWO") next()
	print(paste("Downloading LAI for ", siteID))
	out <- get_LAI(numeric_site_id = s, siteID = siteID, 
								 out_dir = '/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/modis/', start_year = 2013, 
								 end_year = 2020,
								 ncores=10)
	print(paste("Completed LAI for ", siteID))
}



bart <- load('/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/modis_test/BART.rds')
load('/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/modis_test/BART.rds')
save(out, file = '/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/modis_test/BART2.rds')

# Turns out the RDS files are actually rdata, hence the approach below
library(miceadds)
myfiles = list.files("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/modis/","*.rds", full.names="TRUE")
in_list <- list()
for (i in 1:length(myfiles)){
	print(myfiles[i])
	siteID <- tools::file_path_sans_ext( basename(myfiles[i]))
	if(siteID=="LAI_allsites") next()
	load(myfiles[i])
	output$siteID <- siteID
	in_list[[i]] <- output
}
lai_dat <- do.call(rbind, in_list)

saveRDS(lai_dat, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/modis/LAI_allsites.rds")
