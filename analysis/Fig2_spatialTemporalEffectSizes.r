library(runjags)
library(coda)
library(forestplot)
library(gridExtra)

#Summarize performance

cov.mods.fun <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/phylum_funJAGS.rds")
cov.mods.bac <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/phylum_bacJAGS.rds")
# Placeholders:
nocov.mods.fun <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/phylum_funJAGSspatiotemporal.rds")
nocov.mods.bac <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/phylum_bacJAGSspatiotemporal.rds")

####### 1. FUNGI - with covs ######
mcmc.sum <- list()
for (j in 1:length(cov.mods.fun)){
  mod <- cov.mods.fun[[j]]
  group.mcmc <- as.data.frame(mod$JAGS$mcmc[[1]])
  group.mcmc$group <- names(cov.mods.fun)[j]
  mcmc.sum[[j]] <- group.mcmc
}
mcmc.sum.all <- do.call(rbind, mcmc.sum)

# Convert variance to standard deviation, combine into data.frame
plot_var <- sqrt(1/mcmc.sum.all[,"plot_var"])
plot_var <- cbind.data.frame("plot_var", plot_var)
colnames(plot_var) <- c("param","value")
site_var <- sqrt(1/mcmc.sum.all[,"site_var"])
site_var <- cbind.data.frame("site_var", site_var)
colnames(site_var) <- c("param","value")
time_var <- sqrt(1/mcmc.sum.all[,"time_var"])
time_var <- cbind.data.frame("time_var", time_var)
colnames(time_var) <- c("param","value")
df.covs.fungi <- rbind(site_var,plot_var, time_var)


####### 2. FUNGI - no covs ######

mcmc.sum <- list()
for (j in 1:length(nocov.mods.fun)){
  mod <- nocov.mods.fun[[j]]
  group.mcmc <- as.data.frame(mod$JAGS$mcmc[[1]])
  group.mcmc$group <- names(nocov.mods.fun)[j]
  group.mcmc <- group.mcmc[sample(1:1500, 1000, replace = F),]
  mcmc.sum[[j]] <- group.mcmc
}
mcmc.sum.all <- do.call(rbind, mcmc.sum)
# Convert variance to standard deviation, combine into data.frame
plot_var <- sqrt(1/mcmc.sum.all[,"plot_var"])
plot_var <- cbind.data.frame("plot_var", plot_var)
colnames(plot_var) <- c("param","value")
site_var <- sqrt(1/mcmc.sum.all[,"site_var"])
site_var <- cbind.data.frame("site_var", site_var)
colnames(site_var) <- c("param","value")
time_var <- sqrt(1/mcmc.sum.all[,"time_var"])
time_var <- cbind.data.frame("time_var", time_var)
colnames(time_var) <- c("param","value")
df.nocovs.fungi <- rbind(site_var,plot_var, time_var)



####### 3. BACTERIA - with covs ######
mcmc.sum <- list()
for (j in 1:length(cov.mods.bac)){
  mod <- cov.mods.bac[[j]]
  group.mcmc <- as.data.frame(mod$JAGS$mcmc[[1]])
  group.mcmc$group <- names(cov.mods.bac)[j]
  mcmc.sum[[j]] <- group.mcmc
}
mcmc.sum.all <- do.call(rbind, mcmc.sum)

# Convert variance to standard deviation, combine into data.frame
plot_var <- sqrt(1/mcmc.sum.all[,"plot_var"])
plot_var <- cbind.data.frame("plot_var", plot_var)
colnames(plot_var) <- c("param","value")
site_var <- sqrt(1/mcmc.sum.all[,"site_var"])
site_var <- cbind.data.frame("site_var", site_var)
colnames(site_var) <- c("param","value")
time_var <- sqrt(1/mcmc.sum.all[,"time_var"])
time_var <- cbind.data.frame("time_var", time_var)
colnames(time_var) <- c("param","value")
df.covs.bac <- rbind(site_var,plot_var, time_var)


####### 4. BACTERIA - no covs ######
mcmc.sum <- list()
for (j in 1:length(nocov.mods.bac)){
  mod <- nocov.mods.bac[[j]]
  group.mcmc <- as.data.frame(mod$JAGS$mcmc[[1]])
  group.mcmc$group <- names(nocov.mods.bac)[j]
  group.mcmc <- group.mcmc[sample(1:1500, 1000, replace = F),]
  mcmc.sum[[j]] <- group.mcmc
}
mcmc.sum.all <- do.call(rbind, mcmc.sum)

# Convert variance to standard deviation, combine into data.frame
plot_var <- sqrt(1/mcmc.sum.all[,"plot_var"])
plot_var <- cbind.data.frame("plot_var", plot_var)
colnames(plot_var) <- c("param","value")
site_var <- sqrt(1/mcmc.sum.all[,"site_var"])
site_var <- cbind.data.frame("site_var", site_var)
colnames(site_var) <- c("param","value")
time_var <- sqrt(1/mcmc.sum.all[,"time_var"])
time_var <- cbind.data.frame("time_var", time_var)
colnames(time_var) <- c("param","value")
df.nocovs.bac <- rbind(site_var,plot_var, time_var)



# PLOTS FOR ALL 4

# Create plot
withCovsFungi <- ggplot(df.covs.fungi, aes(x = value, y = param, fill=param)) +
  stat_density_ridges(quantile_lines = T, quantiles = 2) +
  theme_ridges() + xlab(NULL) + ylab(NULL) +
  theme(legend.position = "none") + 
  labs(subtitle = "c) WITH environmental predictors (pH, precipitation, temperature)") + 
  scale_x_continuous(limits=c(0,12)) +
  scale_y_discrete(labels = NULL) +
  scale_fill_brewer(palette = "Set3", name = NULL, labels = c("site_var"="Site variability","plot_var"="Plot variability","time_var"="Time variability"))

# Create plot
withCovsBac <- ggplot(df.covs.bac, aes(x = value, y = param, fill=param)) +
  stat_density_ridges(quantile_lines = T, quantiles = 2) +
  theme_ridges() + xlab(NULL) + ylab(NULL) +
  theme(legend.position = "none") + 
  labs(subtitle = "d) WITH environmental predictors (pH, precipitation, temperature)") + 
  scale_x_continuous(limits=c(0,12)) +
  scale_y_discrete(labels = NULL) +
  scale_fill_brewer(palette = "Set3", name = NULL, labels = c("site_var"="Site variability","plot_var"="Plot variability","time_var"="Time variability"))

# Create plot
noCovsFungi <- ggplot(df.nocovs.fungi, aes(x = value, y = param, fill=param)) +
  stat_density_ridges(quantile_lines = T, quantiles = 2) +
  theme_ridges() + xlab(NULL) + ylab(NULL) +
  theme(legend.position = "none") + 
  labs(title = "Fungi", subtitle = "a) No environmental predictors", y) + 
  scale_x_continuous(limits=c(0,12)) +
  scale_y_discrete(labels = NULL) +
  scale_fill_brewer(palette = "Set3", name = NULL, labels = c("site_var"="Site variability","plot_var"="Plot variability","time_var"="Time variability"))

# Create plot
noCovsBac <- ggplot(df.nocovs.bac, aes(x = value, y = param, fill=param)) +
  stat_density_ridges(quantile_lines = T, quantiles = 2) +
  theme_ridges() + xlab(NULL) + ylab(NULL) +
  labs(title = "Bacteria", subtitle = "b) No environmental predictors") + 
  scale_x_continuous(limits=c(0,12)) +
  scale_y_discrete(labels = NULL) +
  scale_fill_brewer(palette = "Set3", name = NULL, labels = c("site_var"="Site variability","plot_var"="Plot variability","time_var"="Time variability")) + theme(legend.position = c(0.6, 0.7))

grid.arrange(noCovsFungi, noCovsBac, 
              withCovsFungi, withCovsBac,
              top=textGrob("Variance across site, plot, and time", gp=gpar(fontsize=24)), 
              bottom=textGrob("Standard deviation of random effects", gp=gpar(fontsize=16)))

v <- arrangeGrob(noCovsFungi, noCovsBac, 
             withCovsFungi, withCovsBac,
             top=textGrob("Variance across site, plot, and time", gp=gpar(fontsize=24)), 
             bottom=textGrob("Standard deviation of random effects", gp=gpar(fontsize=16)))

ggsave(filename = "Figure2.png", plot = v, device="png")
             




# OLD - forest plots.

# sum.all <- do.call(rbind, sum.all) # combine all outputs.
# ci.data <-  cbind.data.frame(mean = sqrt(1/as.numeric(sum.all[,4])), 
#                              lower = sqrt(1/as.numeric(sum.all[,1])), 
#                              upper = sqrt(1/as.numeric(sum.all[,3])))
# ci.data <- rbind(c(NA,NA,NA), ci.data) # create empty top column to match text column titles.
# tabletext <- cbind(c("Group", sum.all$cat), 
#                    c("Parameter", rownames(sum.all)), 
#                    c("Taxon", sum.all$group), 
#                    c("EffectiveSampleSize", sum.all$SSeff))
# 
# forestplot(tabletext, 
#            ci.data,
#            col=fpColors(box="royalblue",line="darkblue"))# Extract covariate posterior estimates from summary table.
# 
# plotvar.dat <- ci.data[grep("plot_var",tabletext[,2]),]
# plotvar.dat <- rbind(c(NA,NA,NA),plotvar.dat, c(mean(plotvar.dat[,1]), mean(plotvar.dat[,2]), mean(plotvar.dat[,3])))
# plotvar.text <- tabletext[grep("plot_var",tabletext[,2]),]
# plotvar.text <- rbind(c("Group", NA,"Taxon","SSeff"),
#                       plotvar.text, c(NA,NA,"Summary",NA))
# 
# grid.newpage()
# pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 3)))
# pushViewport(viewport(layout.pos.col = 1))
# forestplot(plotvar.text[,c(3)], 
#            plotvar.dat,
#            is.summary = c(TRUE, rep(FALSE,20), TRUE),
#            hrzl_lines = gpar(col="444444"),
#            col=fpColors(box="royalblue",line="darkblue"),
#            new_page = FALSE,
#            title="Plot random effect variance")
# popViewport()
# 
# sitevar.dat <- ci.data[grep("site_var",tabletext[,2]),]
# sitevar.dat <- rbind(c(NA,NA,NA),sitevar.dat, c(mean(sitevar.dat[,1]), mean(sitevar.dat[,2]), mean(sitevar.dat[,3])))
# sitevar.text <- tabletext[grep("site_var",tabletext[,2]),]
# sitevar.text <- rbind(c("Group", NA,"Taxon","SSeff"),
#                       sitevar.text, c(NA,NA,"Summary",NA))
# 
# pushViewport(viewport(layout.pos.col = 2))
# forestplot(#sitevar.text[,c(3)], 
#            sitevar.dat,
#            is.summary = c(TRUE, rep(FALSE,20), TRUE),
#            hrzl_lines = gpar(col="444444"),
#            col=fpColors(box="royalblue",line="darkblue"),
#            new_page = FALSE,
#            title="Site random effect variance")
# popViewport(1)