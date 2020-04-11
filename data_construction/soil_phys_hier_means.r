# prep NEON validation covariates from 2016.
source('/usr3/graduate/zrwerbin/NEFI_microbe/NEFI_functions/hierarch_core.means_JAGS.r')
source('/usr3/graduate/zrwerbin/NEFI_microbe/NEFI_functions/pC_uncertainty_neon.r')
source('/usr3/graduate/zrwerbin/NEFI_microbe/NEFI_functions/cn_uncertainty_neon.r')

# read in soil physical data

phys <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/dp.10086.soil_phys.rds")

core.level <- phys
core.level <- core.level[core.level$sampleTiming=="peakGreenness",]
core.level <- core.level[!(duplicated(core.level)),]

core.level$pC_sd <- pC_uncertainty_neon(core.level$pC)
core.level$cn_sd <- cn_uncertainty_neon(core.level$cn)

#get higher level observations and uncertainties (plot, site and global level.)----
#pC
pC.ag <- hierarch_core.means_JAGS(core.level$pC,core_plot = core.level$plotID)
pC.plot <- pC.ag$plot.table[,c('plotID','Mean','SD')]
pC.site <- pC.ag$site.table[,c('siteID','Mean','SD')]
pC.glob <- pC.ag$glob.table[,c('Mean','SD')]
colnames(pC.plot)[2:3] <- c('pC','pC_sd')
colnames(pC.site)[2:3] <- c('pC','pC_sd')
#C:N
cn.ag <- hierarch_core.means_JAGS(core.level$cn,core_plot = core.level$plotID)
cn.plot <- cn.ag$plot.table[,c('plotID','Mean','SD')]
cn.site <- cn.ag$site.table[,c('siteID','Mean','SD')]
cn.glob <- cn.ag$glob.table[,c('Mean','SD')]
colnames(cn.plot)[2:3] <- c('cn','cn_sd')
colnames(cn.site)[2:3] <- c('cn','cn_sd')
#pH - KCl (default)
pH.ag <- hierarch_core.means_JAGS(core.level$pH,core_plot = core.level$plotID)
pH.plot <- pH.ag$plot.table[,c('plotID','Mean','SD')]
pH.site <- pH.ag$site.table[,c('siteID','Mean','SD')]
pH.glob <- pH.ag$glob.table[,c('Mean','SD')]
colnames(pH.plot)[2:3] <- c('pH','pH_sd')
colnames(pH.site)[2:3] <- c('pH','pH_sd')
