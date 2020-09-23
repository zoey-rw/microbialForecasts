#### Prepare map data for Boston trees app #####
library(googlesheets)
library(tigris)   # census data for map
library(sf)
  
#### read/subset/save parcel data ####
# shapefiles downloaded from https://data.boston.gov/dataset/parcels-2016-data-full
shape <- read_sf(dsn = "data/raw/Parcels_2016_Data_Full/Parcels_2016_Data_Full.shp", layer = "Parcels_2016_Data_Full")
parcel <- shape[,c("OBJECTID","LAND_SF", "PTYPE")]
parcel <- parcel[parcel$PTYPE %in% c(900:903, 910:929, 908,973, 986, 965, 978, 984),]
saveRDS(parcel, "data/parcel.rds")

#### community organization data - compiled by hand ####
gs <- gs_title("Ecological Considerations")
community_orgs <- gs %>% gs_read("CommunityOrgs")
colnames(community_orgs) <- c("neighborhood","org","link")
community_orgs$community_org_link <- paste0("<a href='",community_orgs$link,"'>",
                                              community_orgs$org, "</a>")
community_orgs <- community_orgs[,c("neighborhood", "community_org_link")]
dup_n <- community_orgs$neighborhood[duplicated(community_orgs$neighborhood)]
dup <- community_orgs[community_orgs$neighborhood %in% dup_n,]
dup <- aggregate(community_org_link ~ neighborhood,  data = dup, paste, collapse="<br><br>")
orgs <- rbind(community_orgs[!community_orgs$neighborhood %in% dup_n,], dup)
 
# to download census tract - neighborhood mapping data
#http://bostonopendata-boston.opendata.arcgis.com/datasets/34f2c48b670d4b43a617b1540f20efe3_0.csv
neighborhood <- read_csv("data/raw/Climate_Ready_Boston_Social_Vulnerability.csv")
neighborhood_tracts <- neighborhood[,c("Name","GEOID10")]
colnames(neighborhood_tracts) <- c("neighborhood","geoid_2010")
community <- merge(neighborhood_tracts, orgs, all.x=T)
community[is.na(community$community_org_link),]$community_org_link <- ""
community$community_org_link <- paste("<img src = https://static.thenounproject.com/png/2374785-200.png width=20>
                                      Community organizations:<br>",
                                      community$community_org_link,
                                      "<br><a href='https://www.sfttbos.org/'> Speak for the Trees: tree giveaways, workshops, tree inventory project, and neighborhood-focused tree planting</a>")

#### read in vegetation data ####
asthma <- read_csv("data/raw/Chronic Asthma Prevalence.csv")    # asthma prevalence
hvi <- readxl::read_xlsx("data/raw/Boston_HVI_original_Mean_LST_joined.xlsx") # HVI 
current_percent <- readxl::read_xls("data/raw/Boston_CTs_Percent_Current_Veg.xls", sheet = "Sheet1") # current vegetation per tract, as percentage
potential_percent <- readxl::read_xls("data/raw/Boston_CTs_Percent_Potential_Veg.xls", sheet = "Sheet1") # potential vegetation per tract, as percentage

# download census tract data/shapes (using Tigris)
tracts <- tracts("massachusetts", county = "Suffolk", year=2010)

# rename some columns
colnames(current_percent) <- c("geoid_2010", "current_percent")
colnames(potential_percent) <- c("geoid_2010", "potential_percent")
colnames(current.potential_total) <- c("geoid_2010", "current_total", "potential_total")
hvi <- hvi[, c("geoid_2010", "HVI_cat", "MEAN")]
colnames(asthma) <- c("geoid_2010", "asthma")

# combine all data types
layers <- Reduce(function(...) 
  merge(..., all=TRUE), 
  list(hvi, asthma, current_percent, potential_percent, current.potential_total, community)) 
layers_boston <- layers[which(!is.na(layers$asthma)),] ## subset to Boston tracts

# assign "priority" regions based on data from "Priority Scenarios_v1_05012019_LJB"
priority_heat_trees <- data.frame(stringsAsFactors=FALSE,
                                  TRACT = c(25025000803, 25025060700, 25025040600, 25025981300),
                                  LST = c(105.6097255, 103.3075869, 101.1079007, 99.02967801),
                                  EXISTING = c(9.114305354, 4.297326291, 12.97190396, 1.309414209),
                                  POTENTIAL = c(11.80695961, 10.4565422, 12.01300784, 30.97790907),
                                  NEIGHBORHOOD = c("Allston", "South Boston", "Charlestown",
                                                   "South Boston - AIRPORT*")
)

priority_heat_other <- data.frame(stringsAsFactors=FALSE,
                                  TRACT = c(25025030300, 25025030400, 25025030500, 25025070300),
                                  LST = c(104.9576499, 104.9576499, 104.9576499, 98.44355332),
                                  EXISTING = c(3.683573139, 6.751290674, 6.654703345, 6.55424079),
                                  POTENTIAL = c(2.733153958, 3.065363109, 3.873866135, 2.124896863),
                                  NEIGHBORHOOD = c("Downtown", "North End", "North End", "South End")
)
priority_hvi_trees <- data.frame(stringsAsFactors=FALSE,
                                 TRACT = c(25025050700, 25025061101, 25025082100, 25025101102),
                                 HVI = c(21L, 21L, 21L, 20L),
                                 EXISTING = c(10.62964315, 15.74703301, 25.65414725, 23.12593173),
                                 POTENTIAL = c(11.21600552, 12.62879321, 14.82878533, 18.4756271),
                                 NEIGHBORHOOD = c("East Boston", "South Boston", "Roxbury", "Mattapan")
)

priority_hvi_other <- data.frame(stringsAsFactors=FALSE,
                                 TRACT = c(25025010403, 25025050600, 25025070402, 25025050200),
                                 HVI = c(22L, 21L, 21L, 21L),
                                 EXISTING = c(8.061162378, 9.82929216, 16.86169677, 12.20720111),
                                 POTENTIAL = c(3.609607398, 4.412372465, 2.305728283, 5.853202724),
                                 NEIGHBORHOOD = c("Fenway", "East Boston", "South End", "East Boston")
)

layers_boston$priority_heat_trees <- FALSE
layers_boston$priority_heat_trees[layers_boston$geoid_2010 %in% priority_heat_trees$TRACT] <- TRUE
layers_boston$priority_hvi_trees <- FALSE
layers_boston$priority_hvi_trees[layers_boston$geoid_2010 %in% priority_hvi_trees$TRACT] <- TRUE
layers_boston$priority_heat_other <- FALSE
layers_boston$priority_heat_other[layers_boston$geoid_2010 %in% priority_heat_other$TRACT] <- TRUE
layers_boston$priority_hvi_other <- FALSE
layers_boston$priority_hvi_other[layers_boston$geoid_2010 %in% priority_hvi_other$TRACT] <- TRUE


#### create map pop-ups ####

# pollen: CAsthma above 12
layers_boston$pollen <- NA
layers_boston[layers_boston$asthma > 12,]$pollen <- "<img src = https://static.thenounproject.com/png/250340-200.png width=20> Low-pollen tree recommended (high asthma prevalence)"

# heat reduction: hvi above 15, or temp above 36
layers_boston$heat_reduc <- NA
layers_boston[!is.na(layers_boston$HVI_cat) & layers_boston$HVI_cat > 15 | layers_boston$MEAN > 36,]$heat_reduc <- "<img src = https://static.thenounproject.com/png/1819461-200.png width=20> Heat-reducing tree recommended (high temperatures or heat vulnerability)"

# emerald ash borer: all of Boston
layers_boston$EAB <-'<img src = https://static.thenounproject.com/png/438210-200.png width=20> <a href="https://massnrc.org/pests/blog/">Emerald Ash Borer sighting (Boston-wide: click for updates)</a>'

layers_boston$popups <- apply(cbind(paste0("<b>", layers_boston$neighborhood,"</b>"),
                                           layers_boston$pollen, 
                                           layers_boston$heat_reduc, 
                                           layers_boston$EAB, 
                                           layers_boston$community_org_link), 
                                    1, function(x) paste(x[!is.na(x)],  collapse = "<br><br/>"))

tract_data <- geo_join(tracts, layers_boston, "GEOID10", "geoid_2010", how="inner")

saveRDS(tracts, "data/tracts.rds")
saveRDS(tract_data, "data/tract_data.rds")

