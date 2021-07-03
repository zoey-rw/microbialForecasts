# Calculate alpha diversity
library(tidyr)
library(phyloseq)
library(dplyr)
library(lubridate)

source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/functions/helperFunctions.r")

# Combine legacy and recent phyloseq objects
recent_ps <- readRDS("/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/data/NEON_ITS_phyloseq_subset.rds")
legacy_ps <- readRDS("/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/data/NEON_ITS_phyloseq_legacy.rds")
new_sample_dat <- parseNEONsampleIDs(as.character(sample_data(recent_ps)$dnaSampleID))
rownames(new_sample_dat) <- rownames(sample_data(recent_ps))
sample_data(recent_ps) <- new_sample_dat
master_ps <- merge_phyloseq(legacy_ps, recent_ps)
colnames(tax_table(master_ps)) <- tolower(colnames(tax_table(master_ps)))

shannon <- estimate_richness(master_ps, measures = "Shannon")

dat <- cbind(parseNEONsampleIDs(rownames(shannon)), shannon)
saveRDS(dat, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_ITS_full.rds")


cal_ps <- prune_samples(sample_data(master_ps)$asDate  < "2017-01-01", master_ps)
val_ps <- prune_samples(sample_data(master_ps)$asDate  >= "2017-01-01", master_ps)

shannon_cal <- estimate_richness(cal_ps, measures = "Shannon")
shannon_val <- estimate_richness(val_ps, measures = "Shannon")

cal_dat <- cbind(parseNEONsampleIDs(rownames(shannon_cal)), shannon_cal)
val_dat <- cbind(parseNEONsampleIDs(rownames(shannon_val)), shannon_val)

saveRDS(list(cal = cal_dat, val = val_dat), "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/alpha_div_ITS.rds")



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
	) + ylim(c(-1,10))


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

