# visualizing & comparing Kraken classification success with different databases

library(data.table)
library(tidyverse)

sample_names <- c("HARV_002-O-20190716-COMP-DNA1","HARV_004-O-31-8-20131122-gen","HARV_005-O-20190716-COMP-DNA1","HARV_010-O-20160803-COMP-DNA1")
sample_names <- c("HARV_002-O-20190716-COMP-DNA1")

all_files <- list.files("/projectnb/microbiome/dgolden/Struo2/kraken2_testing")
sample_names <- gsub("DB_original_|DB_updated_|.kreport", "", all_files) %>% unique()

#sample_names <- c("HARV_002-O-20190716-COMP-DNA1","HARV_004-O-31-8-20131122-gen")

original_db_results <- list()
updated_db_results <- list()
updated_db_classified <- list()
updated_db_unclassified <- list()
original_db_classified <- list()
original_db_unclassified <- list()

# Loop through all sample files and add to lists
for (s in sample_names) {
	message("Reading in sample: ", s)
	path1 <- paste0("/projectnb2/microbiome/dgolden/Struo2/kraken2_testing/DB_original_", s, ".kreport")
	path2 <- paste0("/projectnb2/microbiome/dgolden/Struo2/kraken2_testing/DB_updated_", s, ".kreport")
	original_db_results[[s]] <- fread(path1, col.names = c("classified_status","seq_id","tax_id","seq_length","kmer_mapping"), nThread = 8, colClasses = c("character","character","character","character","character"))
	updated_db_results[[s]] <- fread(path2, col.names = c("classified_status","seq_id","tax_id","seq_length","kmer_mapping"), nThread = 8, colClasses = c("character","character","character","character","character"))
	updated_db_classified[[s]] <- dplyr::filter(updated_db_results[[s]], classified_status=="C")
	updated_db_unclassified[[s]] <- dplyr::filter(updated_db_results[[s]], classified_status=="U")
	original_db_classified[[s]] <- dplyr::filter(original_db_results[[s]], classified_status=="C")
	original_db_unclassified[[s]] <- dplyr::filter(original_db_results[[s]], classified_status=="U")

}

original_taxids <- unique(original_db_results[[s]]$tax_id)
updated_taxids <- unique(updated_db_results[[s]]$tax_id)

only_update <- setdiff(updated_taxids, original_taxids) %>% unique()
only_original <- setdiff(original_taxids, updated_taxids) %>% unique()

options(digits = 5)


# calculate pct classified per sample
pct_classified_list <- list()
for (s in sample_names){
	original_pct <- table(original_db_results[[s]]$classified_status) %>%
		as.data.frame() %>%
		mutate(percent = Freq/sum(Freq),
					 database_version = "original")
	updated_pct <- table(updated_db_results[[s]]$classified_status) %>%
		as.data.frame() %>%
		mutate(percent = Freq/sum(Freq),
					 database_version = "updated")
	pct_classified_list[[s]] <- rbind(original_pct, updated_pct) %>% mutate(sample = !!s)
}
pct_classified <- rbindlist(pct_classified_list) %>%
	pivot_wider(id_cols = c("sample"), names_from = c("Var1", 'database_version'), values_from = "percent") %>%
	mutate(improvement = C_updated - C_original)


# Improvement in classification: values very close to 0
ggplot(pct_classified) +
	geom_point(aes(x = sample, y = improvement)) +
	geom_hline(yintercept = 0, linetype=2) +
	ylab("% of additional reads classified") +
	ggtitle("Improvement in read classification with JGI genomes") +
	theme_bw() +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05)) +
	scale_y_sqrt()

# Actual points pretty much overlap
ggplot(pct_classified) +
	geom_point(aes(x = sample, y = C_original), alpha = .5) +
	geom_point(aes(x = sample, y = C_updated), color = "red", alpha = .5) +
	ylab("Reads classified") +
	ggtitle("Read classificiation with JGI genomes") +
	theme_bw() +
	theme(axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05)) #+
	#scale_y_log10()

##### LINK TAXONOMY #####

# Fix missing taxonomy info

# prelim_map
taxonomy_mapping <- fread("/projectnb/microbiome/dgolden/Struo2/custom_dbs/kraken2/taxonomy/prelim_map.txt",
													nThread = 8, sep = "\t")

# Splitting the column makes it freeze (due to size?)
#taxonomy_mapping1 <- taxonomy_mapping %>% separate(col="V2", into = c(NA, "tax_id1", "tax_id2"), sep = "|")

tax_map = stringi::stri_split_fixed(taxonomy_mapping$V2, "|", n = 3)
tax1 <- lapply(tax_map, "[[", 2) %>% unlist()
tax2 <- lapply(tax_map, "[[", 3)  %>% unlist()
tax_map_df =  cbind.data.frame(tax1, tax2, taxonomy_mapping$V3)

# Check for overlap
length(intersect(tax_map_df$tax1, original_db_results[[1]]$tax_id))
length(intersect(tax_map_df$tax1, updated_db_results[[1]]$tax_id))
tax_map_df[grepl("609216830", tax_map_df$`taxonomy_mapping$V3`),]
#nothing :(


#### Another attempt but with seqid2taxid
seqid2taxid <- fread("/projectnb/microbiome/dgolden/Struo2/custom_dbs/kraken2/seqid2taxid.map", nThread = 8, sep = "\t")
seqid2taxid = stringi::stri_split_fixed(seqid2taxid$V1, "|", n = 3)

tax1 <- lapply(seqid2taxid, "[[", 2) %>% unlist()
tax2 <- lapply(seqid2taxid, "[[", 3)  %>% unlist()
tax_map_df =  cbind.data.frame(tax1, tax2)

# Check for overlap
# This works! Over 50K intersecting
length(intersect(tax_map_df$tax1, original_db_results[[1]]$tax_id))
length(intersect(tax_map_df$tax1, updated_db_results[[1]]$tax_id))
tax_map_df[grepl("609216830", tax_map_df$tax1),]


tax_map_init[grepl("1000232", tax_map_init$ncbi_species_taxid)]

updated_db_results[[1]]$taxonomy <- tax_map_df[match(updated_db_results[[1]]$tax_id, tax_map_df$tax1),]$tax2
original_db_results[[1]]$taxonomy <- tax_map_df[match(original_db_results[[1]]$tax_id, tax_map_df$tax1),]$tax2

newly_classed <- updated_db_results[[1]] %>% filter(tax_id %in% only_update)
original_unique_classed <- original_db_results[[1]] %>% filter(tax_id %in% only_original)
updated_db_results[[1]]$taxonomy


tax_map_init <- read_tsv("/projectnb/microbiome/dgolden/Struo2/custom_dbs/JGI_downloads_data_distinct_name_acc_premade_filter.tsv")

tax_map_init2 <- read_tsv("/projectnb/microbiome/dgolden/Struo2/custom_dbs/JGI_downloads_data_1_200.tsv")
tax_map_init$ncbi_species_taxid
length(intersect(tax_map_init$ncbi_species_taxid, updated_db_results[[1]]$tax_id))
length(intersect(tax_map_init$ncbi_taxid_nonspecies, updated_db_results[[1]]$tax_id))
length(intersect(tax_map_init$`ORGANISM NCBI TAX ID`, updated_db_results[[1]]$tax_id))
609216830

tax_map_init[grepl("609216", tax_map_init$ncbi_species_taxid),]
tax_map_init[grepl("609216", tax_map_init$ncbi_species_taxid),]



tax_IDs_kreport <- updated_db_results[[1]]$tax_id
tax_IDs <- tax_map_init$ncbi_species_taxid

tax_IDs_slice <- substr(tax_IDs, start=1, stop=7)
tax_IDs_kreport <- as.character(tax_IDs_kreport)
tax_IDs_kreport_deduplicate <- tax_IDs_kreport[!duplicated(tax_IDs_kreport)]
tax_IDs_kreport_slice_A <- substr(tax_IDs_kreport_deduplicate, start=1, stop=7)
shared_tax_IDs_A <- tax_IDs_kreport_slice_A[tax_IDs_kreport_slice_A %in% tax_IDs_slice]
tax_IDs_kreport_slice_B <- substr(tax_IDs_kreport_deduplicate, start=2, stop=8)
shared_tax_IDs_B <- tax_IDs_kreport_slice_B[tax_IDs_kreport_slice_B %in% tax_IDs_slice]
tax_IDs_kreport_slice_C <- substr(tax_IDs_kreport_deduplicate, start=3, stop=9)
shared_tax_IDs_C <- tax_IDs_kreport_slice_C[tax_IDs_kreport_slice_C %in% tax_IDs_slice]
tax_IDs_kreport_slice_D <- substr(tax_IDs_kreport_deduplicate, start=4, stop=10)
shared_tax_IDs_D <- tax_IDs_kreport_slice_D[tax_IDs_kreport_slice_D %in% tax_IDs_slice]
tax_IDs_kreport_slice_E <- substr(tax_IDs_kreport_deduplicate, start=5, stop=11)
shared_tax_IDs_E <- tax_IDs_kreport_slice_E[tax_IDs_kreport_slice_E %in% tax_IDs_slice]
shared_tax_IDs_AB <- append(shared_tax_IDs_A, shared_tax_IDs_B)
shared_tax_IDs_CD <- append(shared_tax_IDs_C, shared_tax_IDs_D)
shared_tax_IDs_ABCD <- append(shared_tax_IDs_AB, shared_tax_IDs_CD)
shared_tax_IDs <- append(shared_tax_IDs_ABCD, shared_tax_IDs_E)
