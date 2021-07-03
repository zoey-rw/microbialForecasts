# Script for preparing bacterial abundance data for phylogeny models
# Input: phyloseq objects. Combines with functional group data and outputs abundances.
# Output is genus-level abundances for 30 genera

setwd("/projectnb2/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/")

library(tidyr)
library(phyloseq)
library(padr)
library(ggplot2)
library(dplyr)
source("/projectnb2/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/binTaxGroups.r")
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/functions/helperFunctions.r")

#### 16S #####

# Combine legacy and recent phyloseq objects
recent_ps <- readRDS("/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/data/NEON_16S_phyloseq_subset.rds")
legacy_ps <- readRDS("/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/data/NEON_16S_phyloseq_legacy.rds")
new_sample_dat <- parseNEONsampleIDs(as.character(sample_data(recent_ps)$dnaSampleID))
rownames(new_sample_dat) <- rownames(sample_data(recent_ps))
sample_data(recent_ps) <- new_sample_dat
master_ps <- merge_phyloseq(legacy_ps, recent_ps)
colnames(tax_table(master_ps)) <- tolower(colnames(tax_table(master_ps)))

master_ps.filt = filter_taxa(master_ps, function(x) sum(x) > 20, TRUE)


# Make relative
ps.rel <- transform_sample_counts(master_ps, function(x) x/sum(x))

ps.filt = filter_taxa(ps.rel, function(x) sum(x) > .005, TRUE)


out <- get_tax_level_abun(master_ps, 
													tax_rank_list = "genus", 
													min_seq_depth = 5000)

# prev <- out$genus$prevalence
# prev_top <- prev[order(prev$prevalence, decreasing = T),]



# Now go through ranks to get top abundances
tax_rank <- "genus"

n.taxa <- 30

	rank_abun <- out[[tax_rank]]$rel.abundances
	prev_top <- out[[tax_rank]]$prevalence[order(out[[tax_rank]]$prevalence$prevalence, decreasing = T),]
	most_abundant_taxa <- prev_top[,1]
	most_abundant_taxa <- most_abundant_taxa[!most_abundant_taxa=="other"][1:n.taxa]
	most_abundant_taxa <- gsub("\\-|\\(|\\)", "\\.", most_abundant_taxa)
	out_top10 <- rank_abun[,colnames(rank_abun) %in% most_abundant_taxa, drop=F]
	seqDepth <- rowSums(rank_abun)
	out_top10$other <- 1-rowSums(out_top10)
	# Remove samples with a sum of above one (not sure why they exist)
	out_top10 <- out_top10[which(!rowSums(out_top10) > 1),]
	
	
	ps.rank.filt <- prune_samples(sample_names(master_ps) %in% rownames(out_top10), master_ps)
	rank.df <- cbind(sample_data(ps.rank.filt)[,c("siteID","plotID","dateID","sampleID","dates","plot_date")], out_top10)
	
	# organize by date
	rank.df$dates <- as.Date(as.character(rank.df$dates), "%Y%m%d")
	rank.df <- rank.df[order(rank.df$dates),]
	# subset by date (till 2016).
	val <- rank.df[rank.df$dates >= "2017-01-01",]
	cal <- rank.df[rank.df$dates < "2017-01-01",]
	
	dim(cal)
	dim(val)

saveRDS(cal, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_phylo_30tax.rds")
saveRDS(val, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/val_groupAbundances_phylo_30tax.rds")
