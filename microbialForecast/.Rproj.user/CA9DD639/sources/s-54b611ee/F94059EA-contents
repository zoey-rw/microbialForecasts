# Add GTDB taxonomy to existing JGI dataset

# Metadata downloaded from IMG using the following query:
# MER-FS Metagenome: Query: (Sequencing Assembly Annotation -- High Quality [ Yes ]) AND (Environmental Classification -- Ecosystem Type [ Soil ])
# (Sequencing Assembly Annotation -- High Quality [ Yes ]): 84637 count(s).
#	(Environmental Classification -- Ecosystem Type [ Soil ]): 13124 count(s).
# Final Combination: 3738 count(s).
# Then added GTDB-Tk column, selected "all", and exported to Excel

hq_img_metadata <- data.table::fread("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/highquality_IMG_soil.txt")

hq_img_simple <- hq_img_metadata %>% select(taxon_oid, `IMG Genome ID`, `GTDB-tk Lineage`, `GOLD Analysis Project ID`, `GOLD Sequencing Project ID`, `NCBI Taxon ID`, `NCBI Assembly Accession`) %>% unique()

### Testing: using the genome_info TSV created with create_GOLD_download_tsv.r ###
#length(intersect(genome_info$NCBI_accession, hq_img_metadata$`NCBI Assembly Accession`)) #2879
#length(intersect(genome_info$`AP GOLD ID`, hq_img_metadata$`GOLD Analysis Project ID`)) #3738
#
#genome_info$GTDB_lineage1 <- hq_img_simple[match(genome_info$`AP GOLD ID`, hq_img_metadata$`GOLD Analysis #Project ID`),]$`GTDB-tk Lineage`
#genome_info$GTDB_lineage2 <- hq_img_simple[match(genome_info$NCBI_accession, hq_img_metadata$`NCBI Assembly #Accession`),]$`GTDB-tk Lineage`

# GTDB Taxonomy added for 2407 out of 4304 taxa
ncbi_ids_dan <- data.table::fread("/projectnb/microbiome/dgolden/Struo2/ncbi_ids.txt", header = F)
ncbi_ids_dan$GTDB_lineage1 <- hq_img_simple[match(ncbi_ids_dan$V1, hq_img_metadata$`NCBI Taxon ID`),]$`GTDB-tk Lineage`
table(is.na(ncbi_ids_dan$GTDB_lineage1))

# OKAY NOW IMPLEMENT!!!

# Not sure which TSV is best to use
JGI_tsv_dan <- data.table::fread("/projectnb/microbiome/dgolden/Struo2/custom_dbs/JGI_downloads_gtdb_ids_prokaryotes.tsv")

# Match based on NCBI Taxon ID: Missing 1183 taxa
JGI_tsv_dan$GTDB_lineage1 <- hq_img_simple[match(JGI_tsv_dan$ncbi_species_taxid,
																								 hq_img_metadata$`NCBI Taxon ID`),]$`GTDB-tk Lineage`
table(is.na(JGI_tsv_dan$GTDB_lineage1))

# Match based on AP GOLD ID: Missing 1455 taxa
JGI_tsv_dan$GTDB_lineage2 <- hq_img_simple[match(JGI_tsv_dan$`AP GOLD ID`,
																								 hq_img_metadata$`GOLD Analysis Project ID`),]$`GTDB-tk Lineage`
table(is.na(JGI_tsv_dan$GTDB_lineage2))

# Inspect rows with no taxonomy - which approach is better? Idk
missing_gtdb1 <- JGI_tsv_dan %>% filter(is.na(GTDB_lineage1))
table(missing_gtdb1$`ORGANISM NCBI PHYLUM`)

missing_gtdb2 <- JGI_tsv_dan %>% filter(is.na(GTDB_lineage2))
table(missing_gtdb2$`ORGANISM NCBI PHYLUM`)

# Save to use as Struo2 input TSV



# All environmental-associated metadata from JGI
# 56k rows, but 33k have nothing or 0 in the NCBI Taxon ID column
all_JGI_metadata <- data.table::fread("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/all_JGI_env_metadata.txt")
all_JGI_metadata_simple <- all_JGI_metadata %>%
	select(taxon_oid, `IMG Genome ID`, `GTDB-tk Lineage`, `GOLD Analysis Project ID`,
				 `GOLD Sequencing Project ID`, `NCBI Taxon ID`, `NCBI Assembly Accession`) %>%
	unique() %>% filter(`NCBI Taxon ID` != 0)

JGI_tsv_dan$GTDB_lineage3 <- all_JGI_metadata_simple[match(JGI_tsv_dan$`ORGANISM NCBI TAX ID`,
																													 all_JGI_metadata_simple$`NCBI Taxon ID`),]$`GTDB-tk Lineage`
table(is.na(JGI_tsv_dan$GTDB_lineage3))
#FALSE  TRUE
#2778   812
