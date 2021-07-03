library(ranomaly)
library(easycsv)
library(tidyverse)

# all_mcc <- fread_folder(directory = "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/mcc",
# 						 extension = "csv", showProgress=T, combine="data.frame")
# 
# mcc_16S <- all_mcc %>% dplyr::filter(kingdom=="Bacteria") 
# mcc_16S$siteID <- substr(mcc_16S$dnaSampleID, 1, 4)
# 
# library(dada2)
# each_site_otu <- list()
# for(site in (unique(mcc_16S$siteID))) {
# 	print(site)
# 	site_dat <- mcc_16S %>% dplyr::filter(siteID==!!site) 
# 	site_otu <- site_dat %>% 
# 		pivot_wider(id_cols = dnaSampleID, names_from = taxonSequence, 
# 								values_from = individualCount, values_fn = mean) %>% 
# 		column_to_rownames("dnaSampleID") %>% as.matrix()
# 	print(dim(site_otu))
# 	site_otu_collapse <- collapseNoMismatch(site_otu, verbose=T, identicalOnly = T)
# 	print(dim(site_otu_collapse))
# 	# taxa to remove
# 	present <- colSums(site_otu_collapse > 0, na.rm=T)
# 	keep <- which(present > 5)
# 	site_otu_subset <- site_otu_collapse[,keep]
# 	print(dim(site_otu_subset))
# 	site_otu_subset[is.na(site_otu_subset)] <- 0
# 	each_site_otu[[site]] <- site_otu_subset
# 	
# 	
# }
# 
# seqtab_joined <- mergeSequenceTables(tables = each_site_otu, tryRC = T)
# saveRDS(seqtab_joined, "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/MCC_otu_16S.rds")


seqtab_joined <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/MCC_otu_16S.rds")
dada_res <- list()
dada_res$seqtab.export <- seqtab_joined
dada_res$seqtab.nochim <- seqtab_joined
dada_res$otu.table <- t(seqtab_joined)

tax1 <- "/projectnb/microbiome/ref_db/taxonomy_db/SILVA_SSU_r138_2019.RData"
tax2 <- "/projectnb/microbiome/ref_db/taxonomy_db/gtbd120_idtaxa.rdata"
out.file <- "/projectnb/microbiome/zrwerbin/NEON_amplicon/16S/NEON_16S_taxonomy"

tax.table = assign_taxo_fun(dada_res = dada_res, id_db = c(tax1,tax2), 
														verbose = 3, output = out.file)


saveRDS(each_site_otu, "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/MCC_otu_16S.rds")
