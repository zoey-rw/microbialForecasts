# Calculate alpha diversity
library(phyloseq)
library(lubridate)
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")


# Output from reformat_taxonomy.r
# master_ps <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/phyloseq_16S.rds")
# 
# rared <- rarefy_even_depth(master_ps, sample.size = 5000, rngseed = 1)
# 
# shannon <- estimate_richness(rared, measures = "Shannon")
# # chao <- estimate_richness(master_ps, measures = "chao1")
# # shannon$seqDepth <- sample_sums(master_ps)
# dat <- cbind(parseNEONsampleIDs(rownames(shannon)), shannon)
# saveRDS(dat, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_16S_full.rds")


dat <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_16S_full.rds")



# Mean center and scale data (within each site), for only non-legacy (recent) data
recent <- dat[dat$asDate >= "2015-11-01",]
recent <- recent %>% group_by(siteID) %>% mutate(Shannon_orig = Shannon,
																								 Shannon_scale_site = scale(Shannon_orig, scale = T)) %>% 
	ungroup() %>% mutate(Shannon = scale(Shannon_orig, scale = T))


# Mean center and scale data (within each site), for all data
dat <- dat %>% group_by(siteID) %>% mutate(Shannon_orig = Shannon,
																					 Shannon_scale_site = scale(Shannon_orig, scale = T)) %>% 
	ungroup() %>% mutate(Shannon = scale(Shannon_orig, scale = T))
cal_dat <- dat[dat$asDate < "2017-01-01",]
val_dat <- dat[dat$asDate >= "2017-01-01",]


saveRDS(list(cal = cal_dat, val = val_dat, full = dat, recent= recent), "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_16S.rds")



div_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_16S.rds")
div_cal <- div_in$cal


keep_sites <- c("CPER", "DSNY", "HARV", "OSBS", "STER", "TALL", "WOOD", "BART", 
								"JERC", "ORNL", "SCBI", "UNDE", "BLAN", "CLBJ", "DEJU", "DELA"
)

# panels for each site, with global scaling
ggplot(data=div_cal[div_cal$siteID %in% keep_sites,],
			 aes(x = asDate,y = Shannon_scale_global))+
	geom_point(aes(shape = as.factor(horizon), color = plotID), size = 3, show.legend=F) +
	geom_smooth() +
	labs(col = "Parameter", title = "Shannon diversity") + 
	xlab("Date")+ 
	ylab(NULL)+
	facet_wrap(~siteID, scales = "free") +
	theme_bw() + theme(
		text = element_text(size = 16),
		axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
		axis.title=element_text(size=22,face="bold")
	) +ylim(c(-5,5))

# panels for each site
ggplot(data=div_cal[div_cal$siteID %in% keep_sites,],
			 aes(x = asDate,y = Shannon))+
	geom_point(aes(shape = as.factor(horizon), color = plotID), size = 3, show.legend=F) +
	geom_smooth() +
	labs(col = "Parameter", title = "Shannon diversity") + 
	xlab("Taxon")+ 
	ylab(NULL)+
	facet_wrap(~siteID, scales = "free") +
	theme_bw() + theme(
		text = element_text(size = 16),
		axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
		axis.title=element_text(size=22,face="bold")
	) +ylim(c(-5,5))

# all sites together on one plot
ggplot(data=div_cal[div_cal$siteID %in% keep_sites,],
			 aes(x = asDate,y = Shannon))+
	geom_point(aes(shape = as.factor(horizon), color = siteID), size = 3, show.legend=F) +
	geom_smooth() +
	labs(col = "Parameter", title = "Absolute effect size") + 
	xlab("Taxon")+ 
	ylab(NULL)+
	theme_bw() + theme(
		text = element_text(size = 16),
		axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
		axis.title=element_text(size=22,face="bold")
	) 

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

