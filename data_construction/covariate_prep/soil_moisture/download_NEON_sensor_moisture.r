# download soil moisture, temperature, neonstore approach
# Have to use custom neonstore commands because of bugs in date handling of original package

library(neonstore)
library(neonUtilities)
library(data.table)
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/functions/custom_neonstore.r")
store_dir = "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/neon_store/"
token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJ6cndlcmJpbkBidS5lZHUiLCJzY29wZSI6InJhdGU6cHVibGljIiwiaXNzIjoiaHR0cHM6Ly9kYXRhLm5lb25zY2llbmNlLm9yZy8iLCJleHAiOjE3NTc4ODg2NzAsImlhdCI6MTYwMDIwODY3MCwiZW1haWwiOiJ6cndlcmJpbkBidS5lZHUifQ.8eW8vxUOiton-kQ_Xyvva0QSHD_BDd2E5IGeNKW3WHib-m7UpTnEhGFAUlAHGdsUyz-dKE1jMOAGS5A_NRYXGg"

# Regex to only keep the first 3 depths, and the 30-minute product
# Split by site because of timeout issues :(

##### SOIL MOISTURE ####
avail = Z10::dp.avail("DP1.00094.001")
all_sites <- unlist(unique(avail$site))
for (s in 1:length(all_sites)){
neon_download(product = "DP1.00094.001", file_regex = "50[123].030", site = all_sites[s],
              dir = store_dir, type = "basic",
              .token = token)
}

# Read in all data and save as own file
df <- neon_read_custom(product = "DP1.00094.001", dir = store_dir)
dt <- data.table(df)

# Using data.table cause it's faster even though I hate it.
dt[, `:=` (sensorID = paste0(siteID, "_", horizontalPosition, ".", verticalPosition),
                         month = substr(startDateTime, 1, 7),
                         depth = substr(verticalPosition, 3, 3),
           dateID = substr(as.character(startDateTime), 1, 7),
           day = as.Date(substr(as.character(startDateTime), 1, 10)),
           date_time = gsub(" ", "T", substr(startDateTime, 1, 13))
  )]
saveRDS(dt, "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEONSoilMoist_raw_allsites.rds")

