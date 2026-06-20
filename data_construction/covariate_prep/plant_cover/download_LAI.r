pacman::p_load(tidyverse)

get_LAI <- function(
	siteID = "NIWO",
	numeric_site_id = 1,
	start_year = 2013,
	end_year = 2020,
	out_dir = here::here("data/raw/modis/"),
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
	
	# MODIS download helper from the PEcAn multi-site download project; not included
	# in this repo (see Methods/Zenodo). Set MF_CALL_MODIS to the path of call_MODIS.R.
	call_modis_path <- Sys.getenv("MF_CALL_MODIS", "")
	source(call_modis_path)
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

con <- read.csv("https://www.neonscience.org/sites/default/files/NEON_Field_Site_Metadata_20210928.csv")
con <- con %>% dplyr::filter(field_site_type %in% c("Relocatable Terrestrial", "Core Terrestrial", "Gradient Terrestrial")) 
for (s in seq_along(con$field_site_id)){
	siteID <- con[s,]$field_site_id
	if(siteID=="NIWO") next()
	print(paste("Downloading LAI for ", siteID))
	out <- get_LAI(numeric_site_id = s, siteID = siteID, 
								 out_dir = here::here("data/raw/modis/"), start_year = 2013, 
								 end_year = 2020,
								 ncores=10)
	print(paste("Completed LAI for ", siteID))
}




# Turns out the RDS files are actually rdata, hence the approach below
library(miceadds)
myfiles = list.files(here::here("data/raw/modis/"),"*.rds", full.names="TRUE")
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

saveRDS(lai_dat, here::here("data/raw/modis/LAI_allsites.rds"))
