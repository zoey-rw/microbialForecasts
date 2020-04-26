
# Read in covariate data 
cov <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/EE509_model_covariates.rds")
weather <- cov[[1]]
plot.preds <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/soilChemPlot.rds")[[1]]
d <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/groupAbundances_EE509.rds")
rank.df <- d[[1]] # Subsetting bacterial phyla
rank.df$dates <- gsub("2016", "2015", rank.df$dates) # Replacing 2015 with 2016 to avoid a huge gap in dataset.

prepModelData <- function(weather, plot.preds, rank.df, j=1){

abun <- rank.df[, j, drop=F]
  
dat <- cbind.data.frame( # Pull out relevant data
  dateID = as.numeric(as.character(stringr::str_replace_all(substr(rank.df$dates, 1, 7), "-", ""))),
  siteID = (droplevels(factor(substr(rownames(abun), 1, 4)))),
  plotID = (droplevels(factor(substr(rownames(abun), 1, 8)))),
  coreID = rownames(abun),
  y = abun[,1])

dat <- dat[!dat$siteID %in% c("NA","NA.1"),]
dates <- sort(unique(dat$dateID)) # for later

# Expand data frame to include all possible plot-date combinations
dat <- tidyr::expand(dat, nesting(siteID, plotID), dateID) %>% 
  merge(dat, all = T) %>% arrange(dateID, siteID, plotID) 

# Calculate plot mean abundance values to check our model results against.
plot.truth <- dat %>% group_by(siteID, plotID, dateID) %>%  
  summarize(plot_mean = mean(y, na.rm=T)) %>% arrange(dateID, siteID, plotID)

# Merge with plot-level covariates
plot.truth <- merge(plot.truth, plot.preds, all.x=T, all.y=F) %>% arrange(dateID, siteID, plotID) 
all.covs <- merge(plot.truth, weather, all.x = T)

# scale and mean-center.
all.covs[,5:10] <- scale(all.covs[,5:10])

# reformat weather (temporal) covariates
temp.cov <- all.covs %>% dplyr::select(c(dateID,monthly_temp,plotID)) %>%
  pivot_wider(names_from = dateID, values_from = monthly_temp) %>% 
  dplyr::select(-c(plotID)) %>% as.matrix()
precip.cov <- all.covs %>% dplyr::select(c(dateID,monthly_precip,plotID)) %>%
  pivot_wider(names_from = dateID, values_from = monthly_precip) %>% 
  dplyr::select(-c(plotID)) %>% as.matrix()
rownames(temp.cov) <- unique(plot.truth$plotID)
rownames(precip.cov) <- unique(plot.truth$plotID)

# reformat plot-level (spatial) covariates 
temp <- all.covs[all.covs$dateID=="201401",]
space.cov <- temp %>% dplyr::select(c(pH,pC,cn,pN)) %>% as.matrix()
rownames(space.cov) <- unique(plot.truth$plotID)

time.cov <- array(c(temp.cov, precip.cov), 
                      dim=c(nrow(temp.cov), ncol(temp.cov), 2), 
                      dimnames = list(rownames(temp.cov), 
                                   colnames(temp.cov), 
                                   c("temp","precip")))

# Reformat y values so that each column is one date.
y <- dat %>% pivot_wider(names_from = dateID, values_from = y) %>% 
  dplyr::select(-c(siteID, plotID, coreID)) %>% as.data.frame()

# Remove empty rows
y <- y[rowSums(!is.na(y))>0,]
dat <- dat[!is.na(dat$y),]

return(list(dat, y, time.cov, space.cov, as.factor(temp$siteID), plot.truth))
}
