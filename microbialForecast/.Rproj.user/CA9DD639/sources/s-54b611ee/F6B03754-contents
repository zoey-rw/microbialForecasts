# For GOLD data from 20 Jun, 2022
# Downloaded from https://gold.jgi.doe.gov/downloads

library(readxl)
library(tidyverse)

gold_data_readme <- readxl::read_xlsx("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/goldData.xlsx", sheet = 1)
gold_data_study <- readxl::read_xlsx("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/goldData.xlsx", sheet = 2)
gold_data_biosample <- readxl::read_xlsx("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/goldData.xlsx",  sheet = 3)
gold_data_organism <- readxl::read_xlsx("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/goldData.xlsx",  sheet = 4)
gold_data_sequencing <- readxl::read_xlsx("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/goldData.xlsx",  sheet = 5)
gold_data_analysis <- readxl::read_xlsx("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/goldData.xlsx",  sheet = 6)

ecosystem_paths <- readxl::read_xlsx("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/GOLDs5levelEcosystemClassificationPaths.xlsx")

soil_paths <- ecosystem_paths %>% filter(grepl("Soil", `ECOSYSTEM TYPE`))
soil_path_ids <- soil_paths %>% select(`ECOSYSTEM PATH ID`) %>% unlist()

# Filter organisms to organisms found in soil
soil_organisms <- gold_data_organism %>% filter(`ORGANISM ECOSYSTEM PATH ID` %in% soil_path_ids)
# 16213 organisms
soil_organism_ids <- soil_organisms %>% select(`ORGANISM GOLD ID`) %>% unique() %>% unlist()

# Just testing the removal of fungi
no_fungi <- soil_organisms %>% filter(`ORGANISM NCBI SUPERKINGDOM` != "Eukaryota")

# Filter to genomes (analysis projects) from soil organisms
soil_sequencing <- gold_data_sequencing %>% filter(`ORGANISM GOLD ID` %in% soil_organism_ids)
soil_sequencing_ids <- soil_sequencing %>% select(`PROJECT GOLD ID`) %>% unlist()

# Filter to genomes (sequencing metadata) from soil organisms
# This file does not have all the completeness columns we want :(
soil_ap <- gold_data_analysis %>% filter(`AP PROJECT GOLD IDS` %in% soil_sequencing_ids)

genome_info <- merge(soil_ap, soil_sequencing, by.x = "AP PROJECT GOLD IDS", by.y = "PROJECT GOLD ID")
# Most of these columns can be removed, either before or after merging


# Next: add accession info as a column

accession_info <- genome_info$`AP GENBANK`

accession_info_test <- "[{\"genbankId\":\"ARFV00000000\",\"assemblyAccession\":\"GCA_000374445.1\"}]"
get_accession <- function(accession_info_test) {
	accession_info <- jsonlite::parse_json(accession_info_test)
	if (length(accession_info)==0) {
		return(NA)
	}
	out <- accession_info[[1]]
	print(out$assemblyAccession)
	to_return <- ifelse(is.null(out$assemblyAccession), NA, out$assemblyAccession)
	return(to_return)
}
get_accession(accession_info_test)

# Approach 1
genome_info$NCBI_accession <- lapply(genome_info$`AP GENBANK`, get_accession)

# Approach 2
genome_info$NCBI_accession <- NA
for (i in 1:nrow(genome_info)){
	print(i)
	genome_info$NCBI_accession[i] <- get_accession(genome_info$`AP GENBANK`[i])
}

# Check why ~1000 are NA
accession_NAs <- genome_info[is.na(genome_info$NCBI_accession),]

# read in taxmap file from GTDB207 release of gtdb-taxdump
# https://github.com/shenwei356/gtdb-taxdump/releases
taxmap = data.table::fread("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/gtdb-taxdump/R207/taxid.map", data.table=F, col.names = c("NCBI_accession_ID", "GTDB_ID"))
genome_info$custom_GTDB_ID = taxmap[match(genome_info$NCBI_accession, taxmap$NCBI_accession_ID),]$GTDB_ID

unique_NCBI_accession_list = unique(genome_info$NCBI_accession)

unique_NCBI_accession_list <- ifelse(!grepl(".1$", unique_NCBI_accession_list), paste0(unique_NCBI_accession_list, ".1"), unique_NCBI_accession_list)

NCBI_GTDB_key = cbind.data.frame(unique_NCBI_accession_list,
																 GTDB_ID = taxmap[match(unique_NCBI_accession_list, taxmap$NCBI_accession_ID),]$GTDB_ID)


taxmap[grepl("GCA_000092925", taxmap$NCBI_accession_ID),]



taxmap_2 = data.table::fread("/projectnb/microbiome/dgolden/Struo2/custom_dbs/GTDB_release207/kraken2/seqid2taxid.map", data.table=F, col.names = c("NCBI_accession_ID", "GTDB_ID"))





taxdump_changelog = data.table::fread("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/taxid-changelog.csv.gz", data.table=F)#, col.names = c("NCBI_accession_ID", "GTDB_ID"))

taxdump_changelog_not_deleted <- taxdump_changelog %>% filter(change != "DELETE")


JGI_input_TSV = data.table::fread("/projectnb/microbiome/dgolden/Struo2/custom_dbs/JGI_input_unique_species.tsv", data.table=F)

JGI_input_TSV_simple = JGI_input_TSV %>% select(partial_accession, full_accession, ncbi_taxonomy,
																								`ORGANISM NCBI SPECIES`, `ORGANISM NCBI TAX ID`,
																								unique_ncbi_organism_name, fasta_file_path) #%>% distinct(.keep_all = T)
write.csv(JGI_input_TSV_simple, "/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/JGI_TSV_simple.csv")


taxmap3 = data.table::fread("/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/new_JGI_taxdump/taxid.map", data.table=F, header = F)
taxmap3 <- taxmap3 %>% separate(col="V8", into = c("fasta_file_path", "GTDB_taxID"), sep = "\t")

