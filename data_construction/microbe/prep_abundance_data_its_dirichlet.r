# Script for preparing fungal abundance data for models
# Input: phyloseq objects. Combines with functional group data and outputs abundances.
# Output is a list of taxonomic rank/functional groups
setwd("/projectnb2/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/")

library(tidyr)
library(phyloseq)
library(padr)
library(ggplot2)
library(dplyr)
library(data.table)
source("/projectnb2/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/binTaxGroups.r")
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/functions/helperFunctions.r")
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/functions/assign_funfun_new.r")


#### 16S #####

# Combine legacy and recent phyloseq objects
recent_ps <- readRDS("/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/data/NEON_ITS_phyloseq_subset.rds")
legacy_ps <- readRDS("/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/data/NEON_ITS_phyloseq_legacy.rds")
new_sample_dat <- parseNEONsampleIDs(as.character(sample_data(recent_ps)$dnaSampleID))
rownames(new_sample_dat) <- rownames(sample_data(recent_ps))
sample_data(recent_ps) <- new_sample_dat
master_ps <- merge_phyloseq(legacy_ps, recent_ps)
colnames(tax_table(master_ps)) <- tolower(colnames(tax_table(master_ps)))

# Assign functional groups 
tax_df = as.data.frame(as(tax_table(master_ps), "matrix"))
tax.fun <- assign_fungal_guilds(tax_table = tax_df, n.cores = 16)
#tax_table(master_ps) <- as.matrix(new_tax)
na <- tax.fun[which(is.na(tax.fun$guild)),]

#Build phylo-functional group taxonomy tables.----
fg <- data.table(tax.fun)
groups <- c('Arbuscular','Animal Pathogen','Plant Pathogen','Saprotroph','Wood Saprotroph','Ectomycorrhizal')
fg[,groups] <- 0
for (g in 1:length(groups)) {
	group <- groups[[g]]
	fg <- fg %>% 
		dplyr::mutate(!!group := dplyr::case_when(grepl(!!group, guild) ~ !!group,
																							TRUE ~ "other")) %>% as.data.frame()
}
fg <- fg[,colnames(fg) %in% c('phylum','class','order','family','genus','species', groups)]

# replace white spaces with underscores.
fg <- apply(fg, 2, function(x) gsub('\\s+', '_',x))
colnames(fg) <- lapply(colnames(fg), function(x) gsub('\\s+', '_',x))

# reassign to phyloseq object.
rownames(fg) <- taxa_names(master_ps)
tax_table(master_ps) <- fg



out <- get_tax_level_abun(master_ps, 
													tax_rank_list = colnames(tax_table(master_ps)), 
													min_seq_depth = 5000)

# Now go through ranks to get top abundances
n.taxa <- 10
cal.out.bac <- list()
val.out.bac <- list()

for (tax_rank in names(out)){
	
	rank_abun <- out[[tax_rank]]$rel.abundances
	
	prev_top <- out[[tax_rank]]$prevalence[order(out[[tax_rank]]$prevalence$prevalence, decreasing = T),]
	
	if (tax_rank=="phylum") prev_top <- prev_top[prev_top$prevalence > .3,]
	most_abundant_taxa <- prev_top[1:n.taxa, 1]
	
	most_abundant_taxa <- most_abundant_taxa[!most_abundant_taxa=="other"]
	
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
	
	cal.out.bac[[tax_rank]]	<- cal
	val.out.bac[[tax_rank]]	<- val
}
names(cal.out.bac)[1:6] <- paste0(names(cal.out.bac)[1:6], "_fun")
names(val.out.bac)[1:6] <- paste0(names(val.out.bac)[1:6], "_fun")

saveRDS(cal.out.bac, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_ITS_10tax.rds")
saveRDS(val.out.bac, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/val_groupAbundances_ITS_10tax.rds")


#saveRDS(cal.out.bac, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_ITS.rds")
