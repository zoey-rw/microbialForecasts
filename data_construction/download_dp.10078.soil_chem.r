# download NEON soil chemical properties (DP1.10078.001), build site/plot time-series
# clear environment, load packages and set file paths
rm(list=ls())
library(zoo)
library(dplyr)
library(neonUtilities)
library(scales)
source('/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/merge_left.r')

output.path <- "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/dp.10078.soil_chem.rds"

# download data from each site
all_dat <- loadByProduct(dpID="DP1.10078.001", 
                         site = "all", 
                         package = "expanded", 
                         check.size = FALSE)
soil_chem_all <- all_dat$sls_soilChemistry

# create new date columns, as both a Date class object and a dateID (YYYY-MM)
soil_chem_all$asDate <- as.Date(as.yearmon(soil_chem_all$collectDate))
soil_chem_all$dateID <- substr(soil_chem_all$asDate, 1, 7)
dp1.10078 <- soil_chem_all
# remove analytical replicates
dp1.10078 <- dp1.10078[-which(dp1.10078$remarks %in% c("Duplicate samples", "Replicate samples", "Replicate samples.")),]
dp1.10078 <- dp1.10078[-which(dp1.10078$analyticalRepNumber == 2),]

# fix chem data with N and C split up
dup <- dp1.10078[duplicated(dp1.10078$sampleID),]
dup <- dp1.10078[dp1.10078$sampleID %in% dup$sampleID,] # we want all duplicated rows
N.only <- dup[dup$acidTreatment == "N",]
C.only <- dup[dup$acidTreatment == "Y",]
C.only$nitrogenPercent <- NULL
N.only.merge <- N.only[,colnames(N.only) %in% c("sampleID","nitrogenPercent")]
collapsed.C.N <- merge(C.only, N.only.merge)

# merge back in with rest of C/N data
dp1.10078 <- dp1.10078[!(dp1.10078$sampleID %in% collapsed.C.N$sampleID),]
dp1.10078$CNvals_merged <- FALSE
collapsed.C.N$CNvals_merged <- TRUE # so we know for later that these values are different
collapsed.C.N <- collapsed.C.N[,colnames(dp1.10078)] # reorder column names
collapsed.C.N$CNratio <- round(collapsed.C.N$organicCPercent/collapsed.C.N$nitrogenPercent,1)
dp1.10078 <- rbind(dp1.10078, collapsed.C.N) # merge back together

soil_chem <- dp1.10078
# pull out columns of interest
soil_chem <- soil_chem[,c("siteID","plotID","plotType","dateID","asDate", "organicCPercent", "nitrogenPercent", "CNratio","sampleID","CNvals_merged")]

# save
saveRDS(soil_chem, output.path)

