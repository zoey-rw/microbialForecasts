library(pavian)
library(tidyverse)


fp = 
  "/projectnb2/microbiome/dgolden/Struo2/7_19_test_report.kreport"

# Read in reports
report_orig <- read_report2(fp)

report_mapping <- report_orig %>% select(taxRank, taxID, name, taxLineage) %>% distinct()

#### Another attempt but with seqid2taxid
seqid2taxid <- fread("/projectnb/microbiome/dgolden/Struo2/custom_dbs/kraken2/seqid2taxid.map", nThread = 8, sep = "\t")
seqid2taxid = stringi::stri_split_fixed(seqid2taxid$V1, "|", n = 3)
tax1 <- lapply(seqid2taxid, "[[", 2) %>% unlist()
tax2 <- lapply(seqid2taxid, "[[", 3)  %>% unlist()
tax_map_df =  cbind.data.frame(tax1, tax2)

# Using the report file to assign taxonomy to the tax_map
tax_map_df$taxonomy <- report_mapping[match(as.numeric(tax_map_df$tax1), report_mapping$taxID),]$taxLineage

# Using the report file to assign taxonomy to the full-output
updated_db_results[[1]]$taxonomy <- report_mapping[match(updated_db_results[[1]]$tax_id, as.character(report_mapping$taxID)),]$taxLineage

# Input TSV 
tax_map_init <- read_tsv("/projectnb/microbiome/dgolden/Struo2/custom_dbs/JGI_downloads_data_distinct_name_acc_premade_filter.tsv")

