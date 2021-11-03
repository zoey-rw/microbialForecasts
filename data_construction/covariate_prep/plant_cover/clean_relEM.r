#assigning mycorrhizal status and plot level relative EM abundance, forest cover, and conifer presence.
#clear environment, source paths.
rm(list=ls())

# Specify output path
output.path <- "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEON_relEM_plot.level.rds"

#load tree data
dat <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/NEON_treedata_allsites.rds")

em_traits.path <- "/projectnb/talbot-lab-data/NEFI_data/ITS/pecan_gen/reference_data/ecto_genus_traits_hobbie_Jan2018.csv"
em_species.path <- "/projectnb/talbot-lab-data/NEFI_data/ITS/pecan_gen/reference_data/myc_assignments.rds"
poa_genera.path <- "/projectnb/talbot-lab-data/NEFI_data/ITS/pecan_gen/reference_data/poaceae_genera_wikipedia.rds"
em_genera.path <- "/projectnb/talbot-lab-data/NEFI_data/ITS/pecan_gen/reference_data/tedersoo_2017_genera.csv"
NEON_plantStatus_codes.path <- "/projectnb/talbot-lab-data/NEFI_data/ITS/pecan_gen/reference_data/NEON_DP1.10098.plantStatus_decode.csv"

#1. Assign mycorrhizal associations.----
#load lookup table for legacy plantStatus codes provided by Katie Jones at NEON.
p.codes <- read.csv(NEON_plantStatus_codes.path)

#load mycorrhizal data
myc.spp <- readRDS(em_species.path)
myc.gen <- read.csv(em_genera.path)
myc.spp$genus_spp <- myc.spp$Species
poa.gen <- readRDS(poa_genera.path)

#known US AM genera.
#AM genera based on searching Colin's super mycorrhizal database. 
#If there are at least 5 records and all are AM, then it gets assigned AM at the genus level.
#if a ton of AM and one AM_ECM, still count as AM.
am_genera <- c('Thuja','Fraxinus','Nyssa','Celtis','Cornus','Diospyros','Ilex','Lonicera','Magnolia','Viburnum', as.character(poa.gen$genus))
erm_genera <- c('Rhododendron','Vaccinium')

#Everything is ID'd at least to genus. Species has lots of qualifiers. Lets clear these up.
dat$species <- sub("^(\\S*\\s+\\S+).*", "\\1", dat$scientificName)
dat$genus   <- sub(" .*", ""                 , dat$scientificName)

#assign mycorrhizal status
dat <- merge(dat,myc.spp[,c('Species','MYCO_ASSO')], by.x = 'species', by.y = 'Species', all.x=T)
dat[dat$genus %in% myc.gen$genus,]$MYCO_ASSO <- 'ECM'
dat[dat$genus %in% am_genera,    ]$MYCO_ASSO <- 'AM'
dat[dat$genus %in% erm_genera,   ]$MYCO_ASSO <- 'ERM'

#subset to trees that have a stem diameter measurement.
dat <- dat[grep('tree',dat$growthForm),]
dat <- dat[!(is.na(dat$stemDiameter)),]

#how much of the basal area is assigned a mycorrhizal association? 96%. We good.
dat$basal <- pi * (dat$stemDiameter/2)^2
metric <- round((1 - (sum(dat[is.na(dat$MYCO_ASSO),]$basal) / sum(dat$basal)))*100, 2)
cat(paste0(metric,'% of trees assigned a mycorrhizal association.'))

#We now need to account for dead trees, insect damaged trees.
#deal with legacy codes in data using key.
for(i in 1:nrow(dat)){
  if(dat$plantStatus[i] %in% p.codes$lovElementName){
    dat$plantStatus[i] <- as.character(p.codes[p.codes$lovElementName == dat$plantStatus[i],]$lovElementCode)
  }
}


#2. aggregate to plot scale - this is where uncertainty comes in.----
output.list <- list()
n.reps <- 1000
for(i in 1:n.reps){
  sim.dat <- dat
  sim.dat$stemDiameter <- rnorm(nrow(sim.dat),log(sim.dat$stemDiameter), 0.0316)
  sim.dat$stemDiameter <- exp(sim.dat$stemDiameter)
  sim.dat$basal <- pi * (sim.dat$stemDiameter/2)^2
  
  #Assign live, live_ECM and dead basal area.
  sim.dat$basal_live <- ifelse(grepl('Live', sim.dat$plantStatus) == T, sim.dat$basal, 0)
  sim.dat$basal_dead <- ifelse(grepl('Dead', sim.dat$plantStatus) == T | grepl('dead', sim.dat$plantStatus) == T, sim.dat$basal, 0)
  sim.dat$basal_ECM  <- ifelse(sim.dat$MYCO_ASSO == 'ECM', sim.dat$basal_live, 0)
  
  #aggregate.
  plot.level            <- aggregate(basal_live ~ plotID, data = sim.dat, FUN = sum, na.rm=T, na.action = na.pass)
  plot.level$basal_ECM  <- aggregate(basal_ECM  ~ plotID, data = sim.dat, FUN = sum, na.rm=T, na.action = na.pass)[,2]
  plot.level$basal_dead <- aggregate(basal_dead ~ plotID, data = sim.dat, FUN = sum, na.rm=T, na.action = na.pass)[,2]
  plot.level <- plot.level[!(plot.level$basal_live == 0),]
  plot.level$relEM <- plot.level$basal_ECM / plot.level$basal_live
  plot.level$live_fraction <- plot.level$basal_live / (plot.level$basal_live + plot.level$basal_dead)
  
  #save to list
  output.list[[i]] <- plot.level
}
for(i in 1:length(output.list)){
  if(i == 1){
    basal_live <- output.list[[i]]$basal_live
    basal_ECM  <- output.list[[i]]$basal_ECM
    basal_dead <- output.list[[i]]$basal_dead
    relEM      <- output.list[[i]]$relEM
  }
  if(i > 1){
    basal_live <- cbind(basal_live, output.list[[i]]$basal_live)
    basal_ECM  <- cbind(basal_ECM , output.list[[i]]$basal_ECM )
    basal_dead <- cbind(basal_dead, output.list[[i]]$basal_dead)
    relEM      <- cbind(relEM     , output.list[[i]]$relEM     )
  }
}
desc.list <- list(basal_live,basal_ECM,basal_dead,relEM)
plot.mean <- list()
plot.sd   <- list()
for(i in 1:length(desc.list)){
  plot.mean[[i]] <- rowMeans(desc.list[[i]])
  plot.sd  [[i]] <- matrixStats::rowSds(desc.list[[i]])
}
plot.mean <- data.frame(plot.mean)
plot.sd   <- data.frame(plot.sd  )
colnames(plot.mean) <- c('basal_live','basal_ECM','basal_dead','relEM')
colnames(plot.sd  ) <- c('basal_live_sd','basal_ECM_sd','basal_dead_sd','relEM_sd')
plotID <- output.list[[1]]$plotID
siteID <- substring(plotID,1,4)
plot.level <- data.frame(siteID,plotID,plot.mean,plot.sd)

#3. save plot-level EM relative abundance.----
saveRDS(plot.level, output.path)
