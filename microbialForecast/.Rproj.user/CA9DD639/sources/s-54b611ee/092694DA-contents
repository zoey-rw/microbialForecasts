#remotes::install_github("fbreitwieser/pavian")
library(pavian)
library(tidyverse)
library(data.table)
source("/projectnb/talbot-lab-data/zrwerbin/soil_genome_db/custom_pavian.r")


sample_names <- c("HARV_002-O-20190716-COMP-DNA1")

sample_names <- "HARV_037-O-20190715-COMP-DNA1"
all_files <- list.files("/projectnb/microbiome/dgolden/Struo2/kraken2_testing2")
sample_names <- gsub(pattern = "DB_original_output_|DB_updated_output_|DB_original_report_|DB_updated_report_|.kreport|_bracken_species",
										 replace = "",
										 all_files) %>% unique()

#sample_names <- c("HARV_002-O-20190716-COMP-DNA1","HARV_004-O-31-8-20131122-gen")

original_db_results <- list()
updated_db_results <- list()
updated_db_classified <- list()
updated_db_unclassified <- list()
original_db_classified <- list()
original_db_unclassified <- list()
original_reports <- list()
updated_reports <- list()
summary_original <- list()
summary_updated <- list()
# Loop through all sample files and add to lists
for (s in sample_names) {
	message("Reading in sample: ", s)
	path1 <- paste0("/projectnb2/microbiome/dgolden/Struo2/kraken2_testing2/DB_original_output_", s, ".kreport")
	path2 <- paste0("/projectnb2/microbiome/dgolden/Struo2/kraken2_testing2/DB_updated_output_", s, ".kreport")
	path3 <- paste0("/projectnb2/microbiome/dgolden/Struo2/kraken2_testing2/DB_original_report_", s, ".kreport")
	path4 <- paste0("/projectnb2/microbiome/dgolden/Struo2/kraken2_testing2/DB_updated_report_", s, ".kreport")
	original_db_results[[s]] <- fread(path1, col.names = c("classified_status","seq_id","tax_id","seq_length","kmer_mapping"), nThread = 8, colClasses = c("character","character","character","character","character"))
	updated_db_results[[s]] <- fread(path2, col.names = c("classified_status","seq_id","tax_id","seq_length","kmer_mapping"), nThread = 8, colClasses = c("character","character","character","character","character"))
	updated_db_classified[[s]] <- dplyr::filter(updated_db_results[[s]], classified_status=="C")
	updated_db_unclassified[[s]] <- dplyr::filter(updated_db_results[[s]], classified_status=="U")
	original_db_classified[[s]] <- dplyr::filter(original_db_results[[s]], classified_status=="C")
	original_db_unclassified[[s]] <- dplyr::filter(original_db_results[[s]], classified_status=="U")

	original_reports[[s]] <- read_report3(path3) # edited pavian function, sourced at top of script
	updated_reports[[s]] <- read_report3(path4) # edited pavian function, sourced at top of script


	summary_original[[s]] <- summarize_report_custom(original_reports[[s]]) %>% mutate(database = "GTDB")
	summary_updated[[s]] <- summarize_report_custom(updated_reports[[s]]) %>% mutate(database = "GTDB_JGI")

}


seqid2taxid_updated <- fread("/projectnb/microbiome/dgolden/Struo2/custom_dbs/GTDB_release207_updated_fulltax/kraken2/seqid2taxid.map",
nThread = 8, sep = "\t")
#### Another attempt but with seqid2taxid
seqid2taxid_updated <- fread("/projectnb/microbiome/dgolden/Struo2/custom_dbs/GTDB_release207_updated_fulltax/kraken2/seqid2taxid.map",
														 nThread = 8, sep = "\t")
#seqid2taxid <- fread("/projectnb/microbiome/dgolden/Struo2/custom_dbs/kraken2/seqid2taxid.map", nThread = 8, sep = "\t")
seqid2taxid_updated = stringi::stri_split_fixed(seqid2taxid_updated$V1, "|", n = 3)

tax1 <- lapply(seqid2taxid_updated, "[[", 2) %>% unlist()
tax2 <- lapply(seqid2taxid_updated, "[[", 3)  %>% unlist()
tax_map_df =  cbind.data.frame(tax1, tax2)

# Check for overlap
# This works! Almost 55K intersecting, which leaves 4k unclassified
length(intersect(tax_map_df$tax1, updated_db_results[[1]]$tax_id))


# This should be added to the loop once it is working properly

tax_map_df[match("2440106587", tax_map_df$tax1),]$tax2

updated_db_results[[1]]$taxonomy <- tax_map_df[match(updated_db_results[[1]]$tax_id, tax_map_df$tax1),]$tax2

report_mapping <- updated_reports[[1]] %>% select(taxRank, taxID, name, taxLineage) %>% distinct()

# this didn't work too well because the taxon names are funky for subspecies
# in the original pavian::read_report2 function had the "collapse.taxRanks" command that i took out because it was v slow,
# i guess we need it though.

# > tail(updated_reports[[s]])
# percentage cladeReads taxonReads taxRank      taxID                   name depth             taxLineage
# 126178          0          1          0       G   38762556             g_JAAOXW01     3             g_JAAOXW01
# 126179          0          1          0       S 2330724092 s_JAAOXW01 sp018221005     4 s_JAAOXW01 sp018221005
# 126180          0          1          1      S1 2577974019           s1_018221005     5           s1_018221005
# 126181          0          1          0       G 3828307224             g_JAAOXX01     3             g_JAAOXX01
# 126182          0          1          0       S 2030629558 s_JAAOXX01 sp018220935     4 s_JAAOXX01 sp018220935
# 126183          0          1          1      S1 1270471748           s1_018220935     5           s1_018220935

# Using the report file to assign taxonomy to the tax_map
tax_map_df$taxonomy <- report_mapping[match(as.numeric(tax_map_df$tax1), report_mapping$taxID),]$taxLineage

# Using the report file to assign taxonomy to the full-output
updated_db_results[[1]]$taxonomy <- report_mapping[match(updated_db_results[[1]]$tax_id, as.character(report_mapping$taxID)),]$taxLineage






#### This part goes outside of loop
seqid2taxid_original <- read_tsv("/projectnb2/microbiome/dgolden/Struo2/custom_dbs/GTDB_release207/kraken2/seqid2taxid.map")

#seqid2taxid <- fread("/projectnb/microbiome/dgolden/Struo2/custom_dbs/kraken2/seqid2taxid.map", nThread = 8, sep = "\t")
seqid2taxid_original_repaired = stringi::stri_split_fixed(seqid2taxid_original$`X2 X3`, "|", n = 3)
seqid2taxid_original_repaired = do.call(rbind, seqid2taxid_original_repaired)
seqid2taxid_original_repaired = seqid2taxid_original_repaired %>% as.data.frame()
colnames(seqid2taxid_original_repaired) = c("V1", "tax2","tax1")
seqid2taxid_original_repaired = seqid2taxid_original_repaired %>% separate(col = "tax2", into=c(NA, tax2), sep = "\t")
tax_map_df_original =  cbind.data.frame(tax1 = seqid2taxid_original_repaired$tax1, tax2 = seqid2taxid_original_repaired$tax2)

tax_map_df_original = stringi::stri_split_fixed(seqid2taxid_original$`X2 X3`, "|", n = 3)


# This would go inside a loop to fix each file
report_mapping_original <- original_reports[[1]] %>% select(taxRank, taxID, name, taxLineage) %>% distinct()
# Using the report file to assign taxonomy to the tax_map
tax_map_df_original$taxonomy <- report_mapping_original[match(as.numeric(tax_map_df_original$tax1), report_mapping_original$taxID),]$taxLineage
# Using the report file to assign taxonomy to the full-output
original_db_results[[1]]$taxonomy <- report_mapping[match(original_db_results[[1]]$tax_id, as.character(report_mapping$taxID)),]$taxLineage

# Testing the --use-names flag: WORKS TO ADD THE TAXONOMY!!

new_original_path <- "/projectnb2/microbiome/dgolden/Struo2/kraken2_testing/names_DB_original_output_HARV_037-O-20190715-COMP-DNA1.kreport"
new_original_outputs = fread(new_original_path, col.names = c("classified_status","seq_id","tax_id","seq_length","kmer_mapping"), nThread = 8, colClasses = c("character","character","character","character","character"))

path2="/projectnb2/microbiome/dgolden/Struo2/kraken2_testing2/DB_updated_output_HARV_037-O-20190715-COMP-DNA1.kreport"
updated_outputs <- fread(path2, col.names = c("classified_status","seq_id","tax_id","seq_length","kmer_mapping"), nThread = 8, colClasses = c("character","character","character","character","character"))
updated_outputs$taxonomy <- report_mapping[match(updated_outputs$tax_id, as.character(report_mapping$taxID)),]$taxLineage

colnames(updated_outputs)[c(1,3:6)] <- paste0(colnames(updated_outputs)[c(1,3:6)], "_updated")

head(updated_outputs)
merged_full_outputs <- merge(updated_outputs, new_original_outputs)

# GOAL: figure out if any JGI organism names are in the updated Kraken output files - doesn't necessarily have to be from this df
differing_class_status <- merged_full_outputs %>% filter(classified_status != classified_status_updated)
table(differing_class_status$taxonomy_updated) %>% sort %>% tail(20)

unique_classes <- differing_class_status$taxonomy_updated %>% unique()
output_taxonomic_classes <- lapply(unique_classes, function(x) {
	str_split(x, "_" )

	})

head(sample_tsv)
