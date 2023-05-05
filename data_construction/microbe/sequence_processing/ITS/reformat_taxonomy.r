# Converting ITS outputs into formats suitable for modeling.
source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")

library(phyloseq)
library(dplyr)
library(tidyr)
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/functions/assign_fungal_guilds.r")
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/functions/helperFunctions.r")
source("/projectnb2/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/binTaxGroups.r")

# RECENT DATA #
# Load combined sequence table and taxonomic table
seqtab_orig <- readRDS("/projectnb/microbiome/zrwerbin/NEON_amplicon/ITS/NEON_ITS_ASV_subset.rds")
taxa_orig <- read.table("/projectnb/microbiome/zrwerbin/NEON_amplicon/ITS/NEON_ITS_taxonomy/final_tax_table.csv")

# Prep for phyloseq
new_tax <- do.call(rbind, (lapply(taxa_orig$V2, parse_taxonomy_qiime)))
rownames(new_tax) <- taxa_orig$V1
sample_dat <- parseNEONsampleIDs(rownames(seqtab_orig))

# Assign functional groups
tax <- as.data.frame(new_tax)
tax.fun <- assign_fungal_guilds(tax_table = tax, n.cores = 8)
#na <- tax.fun[which(is.na(tax.fun$guild_collapse)),]
ps_recent <- phyloseq(otu_table(seqtab_orig, taxa_are_rows = F), tax_table(as.matrix(tax.fun)), sample_data(sample_dat))


# LEGACY DATA #
# Load combined sequence table and taxonomic table
seqtab_legacy <- readRDS("/projectnb/microbiome/zrwerbin/NEON_amplicon/ITS/NEON_ITS_ASV_collapsed_legacy.rds")
taxa_legacy <- read.table("/projectnb/microbiome/zrwerbin/NEON_amplicon/ITS/NEON_ITS_taxonomy_legacy/final_tax_table.csv")
# Prep for phyloseq
new_tax_legacy <- do.call(rbind, (lapply(taxa_legacy$V2, parse_taxonomy_qiime)))
rownames(new_tax_legacy) <- taxa_legacy$V1
sample_dat_legacy <- parseNEONsampleIDs(rownames(seqtab_legacy))
# Assign functional groups
tax_legacy <- as.data.frame(new_tax_legacy)
tax.fun_legacy <- assign_fungal_guilds(tax_table = tax_legacy, n.cores = 8)
#na <- tax.fun[which(is.na(tax.fun$guild_collapse)),]
ps_legacy <- phyloseq(otu_table(seqtab_legacy, taxa_are_rows = F), tax_table(as.matrix(tax.fun_legacy)), sample_data(sample_dat_legacy))



# Combine legacy and recent!
ps_its <- merge_phyloseq(ps_recent, ps_legacy)
saveRDS(ps_its, "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/phyloseq_ITS.rds")



ps_its <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/phyloseq_ITS.rds")

## Get abundances
out <- get_tax_level_abun(ps_its,
													tax_rank_list = colnames(tax_table(ps_its))[1:5],
													min_seq_depth = 3000)

master_ps <- ps_its

# Now go through ranks to get top abundances
n.taxa <- 20
cal.out.bac <- list()
val.out.bac <- list()
out.fun = list()
ranks <- names(out)
ranks <- ranks[ranks != "kingdom"]
for (tax_rank in ranks){

	rank_abun <- out[[tax_rank]]$rel.abundances

	prev_top <- out[[tax_rank]]$prevalence[order(out[[tax_rank]]$prevalence$prevalence, decreasing = T),]

	if (tax_rank %in% c("phylum","genus")) prev_top <- prev_top[prev_top$prevalence > .2,]
	most_abundant_taxa <- prev_top[, 1]
	most_abundant_taxa <- most_abundant_taxa[!grepl("unassigned|other|genus|class|family|order|^fungi",most_abundant_taxa)][1:n.taxa]


	out_top10 <- rank_abun[,colnames(rank_abun) %in% most_abundant_taxa, drop=F]
	seqDepth <- rowSums(rank_abun)
	out_top10$other <- 1-rowSums(out_top10)
	# Remove samples with a sum of above one (not sure why they exist)
	out_top10 <- out_top10[which(!rowSums(out_top10) > 1),]


	# ps.rank.filt <- prune_samples(sample_names(master_ps) %in% rownames(out_top10), master_ps)
	# rank.df <- cbind(sample_data(ps.rank.filt)[,c("siteID","plotID","dateID","sampleID","dates","plot_date")], out_top10)
	sample_dat <- parseNEONsampleIDs(rownames(out_top10)) %>% select(c("siteID","plotID","dateID","sampleID","dates","plot_date"))
	rank.df <- cbind(sample_dat, out_top10)

	# organize by date
	rank.df$dates <- as.Date(as.character(rank.df$dates), "%Y%m%d")
	rank.df <- rank.df[order(rank.df$dates),]
	out.fun[[tax_rank]] <- rank.df
}


names(out.fun)[1:5] <- paste0(names(out.fun)[1:5], "_fun")


saveRDS(out.fun_save, here("data/clean/groupAbundances_ITS_2023.rds"))




# Append previous fg abundances since they didn't change
old_cal = readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_ITS_2021.rds")
old_val = readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/val_groupAbundances_ITS_2021.rds")

fg_list = list()
for (i in 6:length(old_cal)){
	fg_list[[i]] <- rbind(old_cal[[i]], old_val[[i]])
}
fg_list = fg_list[6:length(old_cal)]
names(fg_list) <- names(old_cal)[6:length(old_cal)]
out.fun_save <- c(out.fun, fg_list)

saveRDS(out.fun_save, here("data/clean/groupAbundances_ITS_2023.rds"))



