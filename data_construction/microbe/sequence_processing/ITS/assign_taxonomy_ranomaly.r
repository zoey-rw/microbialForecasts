## Assigning taxonomy to NEON ITS dataset using the UNITE and UTOPIA databases.

library(ranomaly)
library(phyloseq)
library(dplyr)

source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/functions/helperFunctions.r")

#### Prep full ASV table (from Clara Qin's Google drive) ####
library(dada2)

orig.seqtab <- readRDS("/projectnb/microbiome/zrwerbin/NEON_amplicon/ITS/NEON_ITS_ASV_joinedruns.rds")


# Remove low-prevalence taxa and low-read count samples.
seqtab_subset <- orig.seqtab[which(rowSums(orig.seqtab) > 3000),]
present <- colSums(seqtab_subset > 0)
keep <- which(present > 2)
seqtab_subset <- seqtab_subset[,keep]
dim(orig.seqtab)
dim(seqtab_subset)
# Collapse ASVs from different runs
seqtab_collapse <- collapseNoMismatch(seqtab_subset)
saveRDS(seqtab_collapse, "/projectnb/microbiome/zrwerbin/NEON_amplicon/ITS/NEON_ITS_ASV_collapsed.rds")

# ### keep <- which(present > 5)
# seqtab_subset <- orig.seqtab[,keep]
# print(dim(seqtab_subset))
# saveRDS(seqtab_subset, "/projectnb/microbiome/zrwerbin/NEON_amplicon/ITS/NEON_ITS_ASV_subset.rds")
#### ####
#### Read in smaller dataset ####
# seqtab_subset <- readRDS("/projectnb/microbiome/zrwerbin/NEON_amplicon/ITS/NEON_ITS_ASV_subset.rds")
out.file <- "/projectnb/microbiome/zrwerbin/NEON_amplicon/ITS/NEON_ITS_taxonomy"
dada_res <- list()
dada_res$seqtab.export <- seqtab_collapse
dada_res$seqtab.nochim <- seqtab_collapse


# # # Testing with Harvard site.
# out.file <- "/projectnb/microbiome/zrwerbin/NEON_amplicon/ITS/NEON_ITS_taxonomy_HARV"
# seqtab_HARV <- seqtab_subset[grepl("HARV",rownames(seqtab_subset)),]
# # taxa to remove
# present <- colSums(seqtab_HARV > 0)
# keep <- which(present > 5)
# seqtab_HARV_subset <- seqtab_HARV[,keep]
# print(dim(seqtab_HARV_subset))
# dada_res$seqtab.export <- seqtab_HARV_subset
# dada_res$seqtab.nochim <- seqtab_HARV_subset


tax1 <- "/projectnb/microbiome/ref_db/taxonomy_db/unite_fungi_v8.2_idtaxa.rdata"
tax2 <- "/projectnb/microbiome/ref_db/taxonomy_db/utopia90_201908.rdata"
tax.table = assign_taxo_fun(dada_res = dada_res, id_db = c(tax1,tax2), 
														verbose = 3, output = out.file)

## Combine with metadata into phyloseq
sample_dat <- read.csv("/projectnb/microbiome/zrwerbin/NEON_amplicon/ITS/mmg_metadata_ITS_QCd_20210331.csv")
sample_dat_subset <- sample_dat %>% 
	select(dnaSampleID, geneticSampleID, sampleTotalReadNumber, sampleFilteredReadNumber) %>% 
	filter(dnaSampleID %in% rownames(seqtab_subset)) %>% 
	distinct(dnaSampleID)
rownames(sample_dat_subset) <- sample_dat_subset$dnaSampleID

#new_sample_dat <- parseNEONsampleIDs(rownames(sample_dat_subset))

ps = phyloseq::phyloseq(otu_table = otu_table(dada_res$seqtab.nochim, taxa_are_rows = F), 
												tax.table = tax.table, 
												sample_data = sample_dat_subset)


new_sample_dat <- parseNEONsampleIDs(rownames(otu_table(ps)))
rownames(new_sample_dat) <- new_sample_dat$sampleID
sample_data(ps) <- new_sample_dat

saveRDS(ps, "/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/data/outputs/anomaly_tax/ps_ITS_HARV.rds")

# out <- get_tax_level_abun(ps, 
# 													tax_rank_list = "Genus", 
# 													min_seq_depth = 5000)
# 
# rel_abun <- out$Genus$rel.abundances
# rel_abun <- cbind.data.frame(rel_abun, parseNEONsampleIDs(rownames(rel_abun)))
# rel_abun$year <- substr(rel_abun$dateID, 1, 4)
# ggplot(rel_abun) + geom_jitter(aes(y = other, x = reorder(siteID,other), color = year), width=.2) + ggtitle("Relative abundance of unidentified genera") + theme(axis.text.x = element_text(angle = 270,vjust = 0.5, hjust=1))
# 
# 
# ps <- harv_ps
# 
# 
# top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:200]
# prune.dat_top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
# prune.dat_top20 <- prune_taxa(top20, prune.dat_top20)
# 
# dat.aglo = tax_glom(prune.dat_top20, taxrank = "Genus")
# dat.dataframe = psmelt(dat.aglo)
# ggplot(dat.dataframe, aes(x=siteID, y=Abundance, fill=Genus)) + geom_bar(stat="identity", position="fill") + facet_grid(~siteID, scale="free")
# 
# 
# N <- 50
# barplot(sort(taxa_sums(dat.aglo), TRUE)[1:N]/nsamples(dat.aglo), las=2)
# 
# barplot(dat.dataframe$Genus ~ dat.dataframe$Abundance)
# 
# means <- dat.dataframe %>% group_by(Genus) %>% summarize(mean = mean(Abundance, na.rm=T))
# ggplot(dat.dataframe) + geom_point(aes(y = Abundance, x = Genus)) + coord_flip()
# 
# library(FUNGuildR)
# funguild_assign()
