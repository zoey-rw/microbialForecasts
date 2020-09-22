# Look at relationships between Pinaceae, soil fungi, and soil bacteria at NEON sites

# Load R packages, install if necessary.
#devtools::install_github("admahood/neondiversity") 
library(neondiversity) # This is for getting the plant datas
library(tidyverse)
library(corrplot) 
library(compositions) # Just to transform relative abundance data

##### 1. Plant data #####
# Download data from NEON sites (the ones I have microbial data for)
plant_info <- download_plant_div(sites = c("CPER", "HARV", "STER", "OSBS", "CPER")) 
# Calculate a bunch of metrics for Pinaceae. Documentation for this is here: https://github.com/admahood/neondiversity
plant_stats <- get_diversity_info(plant_info, families = "Pinaceae", spp = "Pinus strobus")
# Let's check out some of this plant data
# Histogram of the relative cover (RC) of exotic Pinaceae at all plots - almost all zeros
hist(fg_plant$rc_exotic_Pinaceae)
# But there is more variation in the RC of Pinaceae in general
hist(fg_plant$rc_Pinaceae)

## If we download data for all NEON sites, there might be more variation in exotic Pinaceae RC, 
## but this takes more than 2 min to download and I am impatient so I don't actually know.
# plant_info_allNEONsite <- download_plant_div(sites = "all") 
# plant_stats_allNEONsite <- get_diversity_info(plant_info_allNEONsite, families = "Pinaceae", spp = "Pinus strobus")
# hist(fg_plant$rc_exotic_Pinaceae)

##### 2. Soil fungi data #####
# Read in ITS functional group data (not finalized, I'll let you know when I have a more full dataset)
ITS <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/output_files/ITS/fg_abun.rds")
# Prep dataframes with fungal functional group relative abundances
fg <- ITS$rel.abundances
fg <- cbind.data.frame(fg, plotID = substr(rownames(fg), 1, 8))

##### 3. Combine! #####
# Merge with plant information
fg_plant <- merge(plant_stats, fg, by = c("plotID"))
# Remove non-numeric columns and all-zero columns, just so we can look at strong correlations
fg_plant_numeric <- fg_plant %>% dplyr::select_if(function(col) is.numeric(col) & any(col != 0))
# Create correlation matrix
cor_mat <- cor(fg_plant_numeric)
# Subset to strong correlations
strong <- abs(cor_mat[1,]) > .2
cor_mat_strong <- cor_mat[strong, strong]
# Visualize correlations
corrplot::corrplot(cor_mat_strong, method="circle", type = "lower")
# Hey look! Looks like % ecto fungi in the soil is strongly correlated with a whole bunch of stuff.

# # Just as a sanity check, do all the same stuff as above, but with ClR-transformed abundances,
# # which are better at removing spurious correlations between relative abundance data
# clr_fg <- compositions::clr(ITS$abundances[,1:5] + 1) # add 1 because CLR can't handle zeros.
# clr_fg <- cbind.data.frame(clr_fg, plotID = substr(rownames(clr_fg), 1, 8))
# clr_fg_plant <- merge(plant_stats, clr_fg, by = c("plotID"))
# clr_fg_numeric <- clr_fg_plant %>% dplyr::select_if(function(col) is.numeric(col) & any(col != 0))
# clr_cor_mat <- cor(clr_fg_numeric)
# clr_strong <- abs(clr_cor_mat[1,]) > .2
# clr_cor_mat_strong <- clr_cor_mat[clr_strong, clr_strong]
# corrplot::corrplot(clr_cor_mat_strong, method="circle", type = "lower")
# # Looks like the only difference is Saprotrophs disappear

# Higher # of exotic species, fewer ecto fungi in soil?
ggplot(fg_plant) + geom_point(aes(x = Ectomycorrhizal, y = nspp_exotic, color = site))
# Just %ecto explains 20% of variation in RC of Pinaceae! NICE
summary(lm(rc_Pinaceae ~ Ectomycorrhizal, data = fg_plant))
