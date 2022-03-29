# Calculate alpha diversity
library(tidyr)
library(phyloseq)
library(dplyr)
library(lubridate)

source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/functions/helperFunctions.r")

# Output from reformat_taxonomy.r
master_ps <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/phyloseq_ITS.rds")
# 
# rared <- rarefy_even_depth(master_ps, sample.size = 5000, rngseed = 1)
# 
# shannon <- estimate_richness(rared, measures = "Shannon")
# # chao <- estimate_richness(master_ps, measures = "chao1")
# # shannon$seqDepth <- sample_sums(master_ps)
# dat <- cbind(parseNEONsampleIDs(rownames(shannon)), shannon)
# saveRDS(dat, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_ITS_full.rds")

dat <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_ITS_full.rds")


# Mean center and scale data (within each site)
dat <- dat %>% group_by(siteID) %>% mutate(Shannon_orig = Shannon,
																					 Shannon_scale_site = scale(Shannon_orig, scale = T)) %>% 
	ungroup() %>% mutate(Shannon = scale(Shannon_orig, scale = T))

# This code is redundant
# cal_ps <- prune_samples(sample_data(rared)$asDate  < "2017-01-01", rared)
# val_ps <- prune_samples(sample_data(rared)$asDate  >= "2017-01-01", rared)
# 
# shannon_cal <- estimate_richness(cal_ps, measures = "Shannon")
# shannon_val <- estimate_richness(val_ps, measures = "Shannon")
# 
# cal_dat <- cbind(parseNEONsampleIDs(rownames(shannon_cal)), shannon_cal)
# val_dat <- cbind(parseNEONsampleIDs(rownames(shannon_val)), shannon_val)


cal_dat <- dat[dat$asDate < "2017-01-01",]
val_dat <- dat[dat$asDate >= "2017-01-01",]
saveRDS(list(cal = cal_dat, val = val_dat, full = dat), "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_ITS.rds")


#### Visualize

div_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_ITS.rds")
div_cal <- div_in$cal


keep_sites <- c("CPER", "DSNY", "HARV", "OSBS", "STER", "TALL", "WOOD", "BART", 
								"JERC", "ORNL", "SCBI", "UNDE", "BLAN", "CLBJ", "DEJU", "DELA"
)


# panels for each site
ggplot(data=div_cal[div_cal$siteID %in% keep_sites,],
			 aes(x = asDate,y = Shannon))+
	geom_point(aes(shape = as.factor(horizon), color = plotID), size = 3, show.legend=F) +
	geom_smooth() +
	labs(col = "Parameter", title = "Shannon diversity") + 
	xlab("Year")+ 
	ylab(NULL)+
	facet_wrap(~siteID) +
	theme_bw() + theme(
		text = element_text(size = 16),
		axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
		axis.title=element_text(size=22,face="bold")
	) #+ ylim(c(-1,10))


# grouped by month, colored by year
div_cal$month <- month(div_cal$asDate)
div_cal$year <- year(div_cal$asDate)
ggplot(data=div_cal,
			 aes(x = month,y = Shannon))+
	geom_point(aes(color = as.factor(year)), size = 3, show.legend=F) +
	geom_smooth() +
	labs(col = "Parameter", title = "Absolute effect size")


# grouped by month, panel by site
ggplot(data=div_cal,
			 aes(x = month,y = Shannon))+
	geom_boxplot(aes(group = month)) +
	# geom_point(aes(color = as.factor(year)), size = 3, show.legend=F) +
	facet_wrap(~siteID, nrow = 4) +
	labs(col = "Shannon diversity", title = "Alpha diversity by month")


# grouped by month, panel by site
ggplot(data=div_cal,
			 aes(x = seqDepth,y = Shannon))+
	geom_point() +
	# geom_point(aes(color = as.factor(year)), size = 3, show.legend=F) +
	facet_wrap(~siteID, nrow = 4) +
	labs(col = "Shannon diversity", title = "Alpha diversity by month")




