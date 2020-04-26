
# Read in covariate data 
# cov <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/EE509_model_covariates.rds")
# weather <- cov[[1]]
# plot.preds <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/soilChemPlot.rds")[[1]]
# d <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/groupAbundances_EE509.rds")
# rank.df <- d[[1]] # Subsetting bacterial phyla
# rank.df$dates <- gsub("2016", "2015", rank.df$dates) # Replacing 2015 with 2016 to avoid a huge gap in dataset.
# model.dat <- prepModelData(weather, plot.preds, rank.df, j=1) 

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

# Calculate site mean abundance values to check our model results against.
site.truth <- dat %>% group_by(siteID, dateID) %>%  
  summarize(site_mean = mean(y, na.rm=T)) %>% arrange(dateID, siteID)
# Merge with site-level covariates
site.truth <- merge(site.truth, weather, all.x=T, all.y=F) %>% arrange(dateID, siteID) 

# Calculate plot mean abundance values to check our model results against.
plot.truth <- dat %>% group_by(siteID, plotID, dateID) %>%  
  summarize(plot_mean = mean(y, na.rm=T)) %>% arrange(dateID, siteID, plotID)

# Merge with plot-level covariates
plot.truth <- merge(plot.truth, plot.preds, all.x=T, all.y=F)  %>% arrange(dateID, siteID, plotID) 

# create temperature covariate matrix
temp.cov <- site.truth %>% 
  dplyr::select(-c(monthly_precip,site_mean)) %>% 
  pivot_wider(names_from = dateID, values_from = monthly_temp) %>% 
  dplyr::select(-c(siteID)) %>% scale(scale=TRUE) %>% as.matrix()

# create temperature covariate matrix
precip.cov <- site.truth %>% 
  dplyr::select(-c(monthly_temp,site_mean)) %>% 
  pivot_wider(names_from = dateID, values_from = monthly_precip) %>% 
    dplyr::select(-c(siteID)) %>% scale(scale=TRUE) %>% as.matrix()

# create pH covariate matrix
# pH.plot <- plot.truth %>% dplyr::select(-c(plot_mean,pC,cn,pN)) %>% 
#   pivot_wider(names_from = dateID, values_from = pH)
# pH.cov <- as.numeric(scale(pH.plot$`201306`, scale=F))
# names(pH.cov) <- pH.plot$plotID
temp <- plot.truth[plot.truth$dateID=="201306",]
plot.covs <- temp %>% dplyr::select(c(pH,pC,cn,pN)) %>% scale(scale=TRUE) %>% as.matrix()
rownames(plot.covs) <- temp$plotID

# Reformat y values so that each column is one date.
y <- dat %>% pivot_wider(names_from = dateID, values_from = y) %>% 
  dplyr::select(-c(siteID, plotID, coreID)) %>% as.data.frame()

# Remove empty rows
y <- y[rowSums(!is.na(y))>0,]
dat <- dat[!is.na(dat$y),]

return(list(dat, y, temp.cov, precip.cov, plot.covs, as.factor(temp$siteID), plot.truth, site.truth))
}
