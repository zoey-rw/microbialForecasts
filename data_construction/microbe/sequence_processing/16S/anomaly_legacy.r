library(ranomaly)
library(dada2)

source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/functions/helperFunctions.r")

# legacy.seqtab <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEON_16S_ITS_data_construction/data/raw_seqs/16S/raw_from_neon/seq_tables/legacyOtuTable_16S.rds")
# legacy.seqtab.df <- as.data.frame(legacy.seqtab)
# legacy.seqtab.df$siteID <- substr(rownames(legacy.seqtab.df), 1, 4)
# 
# each_site_otu <- list()
# for(site in (unique(legacy.seqtab.df$siteID))) {
# 	print(site)
# 	site_dat <- legacy.seqtab.df %>% dplyr::filter(siteID==!!site) 
# 	if(nrow(site_dat)==0) next()
# 	site_dat$siteID <- NULL
# 	site_otu <- site_dat %>% as.matrix()
# 	print(dim(site_otu))
# 	site_otu_collapse <- collapseNoMismatch(site_otu, verbose=T, identicalOnly = T)
# 	print(dim(site_otu_collapse))
# 	# taxa to remove
# 	present <- colSums(site_otu_collapse > 0, na.rm=T)
# 	keep <- which(present > 3)
# 	site_otu_subset <- site_otu_collapse[,keep] %>% as.matrix()
# 	print(dim(site_otu_subset))
# 	
# 	each_site_otu[[site]] <- site_otu_subset
# }
# 
# seqtab_joined <- mergeSequenceTables(tables = each_site_otu, tryRC = T)
# dim(seqtab_joined)
# 



# sampledata_full <- read.csv("/projectnb/talbot-lab-data/zrwerbin/NEON_soil_microbe_processing/data/mmg_soilMetadata_16S_2020-09-08.csv")
# # Create table to link internalLabID with sampleID
# link_ID <- sampledata_full %>% select(dnaSampleID, geneticSampleID, internalLabID) %>%	distinct() 
# link_ID$geneticSampleID <- gsub("-DNA[12]", "", link_ID$dnaSampleID)
# # Fix rownames on sequence table
# seqtab.nochim <- seqtab_joined
# seq.internalLabID <-  gsub("_filt.fastq.gz", "", rownames(seqtab.nochim))
# 
# # Sample names are wrong for bacteria; this creates a table to link deprecatedVialID and geneticSampleID
# map <- readRDS("/projectnb/talbot-lab-data/zrwerbin/NEFI_microbe/pre-release-map_16S.rds")[,c("geneticSampleID","sampleID")]
# map <- map[!duplicated(map),,drop=F]
# new_genIDs <- map[match(seq.internalLabID, map$sampleID),]$geneticSampleID
# new_dnaIDs <- as.character(link_ID[match(new_genIDs, link_ID$geneticSampleID),]$dnaSampleID)
# new_rownames <- ifelse(is.na(new_dnaIDs), new_genIDs, new_dnaIDs) # switch to dnaSampleIDs *if possible*
# keep <- which(!is.na(new_rownames) & !duplicated(new_rownames))
# seqtab <- seqtab.nochim[keep,]
# rownames(seqtab) <- new_rownames[keep]



# saveRDS(seqtab, "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/MCC_otu_16S_legacy.rds")


seqtab_joined <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/MCC_otu_16S_legacy.rds")

out.file <- "/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/data/outputs/legacy_anomaly_tax_16S"

dada_res <- list()
dada_res$seqtab.export <- seqtab_joined
dada_res$seqtab.nochim <- seqtab_joined
dada_res$otu.table <- t(seqtab_joined)

tax1 <- "/projectnb/microbiome/ref_db/taxonomy_db/SILVA_SSU_r138_2019.RData"
tax2 <- "/projectnb/microbiome/ref_db/taxonomy_db/gtbd120_idtaxa.rdata"

tax.table = assign_taxo_fun(dada_res = dada_res, id_db = c(tax1,tax2), 
														verbose = 3, output = out.file)
# 
# tax_table <- as.matrix(cbind.data.frame(apply(X = tax.table, 2, function(x) gsub("[kpcofgs]\\_\\_","", x))))
# 
# 
# tree = generate_tree_fun(dada_res)
# 
# 
# metadata <- parseNEONsampleIDs(rownames(dada_res$seqtab.export))
# rownames(metadata) <- metadata$sampleID
# metadata$sample.id <- colnames(dada_res$otu.table)
# write.csv(metadata, "/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/data/outputs/anomaly_tax/metadata_HARV.csv")
# 
# ps = generate_phyloseq_fun(dada_res = dada_res, 
# 													 tax.table = tax_table, 
# 													 tree = tree,metadata = "/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/data/outputs/anomaly_tax/metadata_HARV.csv")
# 
# phylo_tree <- phy_tree(tree)
# taxa_names(phylo_tree) <- taxa_names(tax_table(tax_table))
# library(phyloseq)
# ps = phyloseq::phyloseq(otu_table = otu_table(dada_res$seqtab.nochim, taxa_are_rows = F), 
# 												tax_table(tax_table),
# 												sample_data(metadata), 
# 												phylo_tree)
# 
# saveRDS(ps, "/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/data/outputs/anomaly_tax/ps_16S_HARV.rds")
# 
# library(ggtree)
# library(scales)
# p <- ggtree(phylo_tree)
# ggtree(ps) + geom_text2(aes(subset=!isTip, label=label), hjust=-.2, size=2) +
# 	geom_tiplab(aes(label=Genus, color = Phylum), hjust=-.3) +
# 	#geom_point(aes(x=x+hjust, color=Phylum, size=Abundance),na.rm=TRUE) +
# 	scale_size_continuous(trans=log_trans(5)) +
# 	theme(legend.position="right")
# 
# bars_fun(data = ps, top = 20, Ord1 = "plot_date", Fact1 = "plot_date", rank="Genus", relative = TRUE)
# 
# heatmap_plot = heatmap_fun(data = ps, column1 = "plotID", top = 20, output = "./plot_heatmap/", rank = "Genus")
# 
# out <- get_tax_level_abun(ps, 
# 													tax_rank_list = "Genus", 
# 													min_seq_depth = 5000)
# 
# rel_abun <- out$Genus$rel.abundances
# rel_abun <- cbind.data.frame(rel_abun, parseNEONsampleIDs(rownames(rel_abun)))
# rel_abun$year <- substr(rel_abun$dateID, 1, 4)
# ggplot(rel_abun) + geom_jitter(aes(y = other, x = reorder(siteID,other), color = year), width=.2, height=0) + ggtitle("Relative abundance of unidentified genera") + theme(axis.text.x = element_text(angle = 270,vjust = 0.5, hjust=1))
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
