pacman::p_load(tidyverse, minval, sybil,readxl)
library(stringr)
source("/projectnb2/talbot-lab-data/zrwerbin/interactions/source.R")

#all_models_orig <- read.csv("/projectnb2/talbot-lab-data/zrwerbin/interactions/models/n_cycle_models.csv")
all_models_orig <- readxl::read_xls("/projectnb2/talbot-lab-data/metabolic_models/curated_models/N_cycle/orig/AOB+NOB8-2A.xls", .name_repair = "universal")

all_mets_orig <- readxl::read_xls("/projectnb2/talbot-lab-data/metabolic_models/curated_models/N_cycle/orig/AOB+NOB8-2A.xls", .name_repair = "universal", sheet = 2)

all_models <- all_models_orig %>% filter(!grepl("ECO.", Rxn.description)) %>% rename(Equation = Formula)
eco <- all_models_orig %>% filter(grepl("ECO", Rxn.description)) %>% rename(Equation = Formula)

repaired_models <- all_models %>% mutate(Equation = gsub("<==>", "<=>", Equation),
																				 Equation = gsub("-->", "=>", Equation),
																				 Equation = gsub("<--", "<=", Equation))
repaired_models <- repaired_models %>%
	separate(Rxn.description, into = c("Organism","Description"), sep="\\. ") %>%
	filter(!is.na(Description)) %>% filter(Equation != "")
repaired_models$Equation <- str_replace_all(repaired_models$Equation, "\\[[abcdifgh]\\]", "\\[e\\]")
repaired_models$Equation <- str_replace_all(repaired_models$Equation, "\\[[jklmnopq]\\]", "\\[p\\]")
repaired_models$Equation <- str_replace_all(repaired_models$Equation, "\\[[rstuvwxy]\\]", "\\[c\\]")
repaired_models$Equation <- str_replace_all(repaired_models$Equation, "[[:space:]]\\.5|^\\.5", "0\\.5")



eco <- eco %>% mutate(Equation = gsub("<==>", "<=>", Equation),
																				 Equation = gsub("-->", "=>", Equation),
																				 Equation = gsub("<--", "<=", Equation))
eco <- eco %>%
	separate(Rxn.description, into = c("Organism","Description"), sep="\\. ") %>%
	filter(!is.na(Description)) %>% filter(Equation != "")
eco$Equation <- str_replace_all(eco$Equation, "\\[[abcdifgh]\\]", "\\[e\\]")
eco$Equation <- str_replace_all(eco$Equation, "\\[[jklmnopq]\\]", "\\[p\\]")
eco$Equation <- str_replace_all(eco$Equation, "\\[[rstuvwxy]\\]", "\\[c\\]")
eco$Equation <- str_replace_all(eco$Equation, "[[:space:]]\\.5|^\\.5", "0\\.5")
eco <- eco %>% select(-c(Rxn.name, Description)) %>% distinct(.keep_all = T)

model_list <- repaired_models %>% group_by(Organism) %>% group_split()
names(model_list) <- c("EX","NDE","NET","NEU","NHA","NMU","NOC","NSP","NWI")

#idk why rep wasn't working here
model_list$EX <- as.data.frame(model_list$EX)
env_list <- list(model_list$EX, model_list$EX, model_list$EX, model_list$EX,
							model_list$EX, model_list$EX, model_list$EX, model_list$EX)
names(env_list) <- c("NDE","NET","NEU","NHA","NMU","NOC","NSP","NWI")


for (i in 2:9){
	to_write <- model_list[[i]] %>% select("ID" = "Rxn.name",
	                                       "SPECIES" = "Organism",
	                                       "DESCR" = "Description",
	                                       "REACTION" = "Equation",
	                                       "GPR" = "Genes",
	                                       "REF" = "References",
	                                       "LOWER.BOUND" = "LB",
	                                       "UPPER.BOUND" = "UB",
	                                       "OBJECTIVE" = "Objective") %>% as.data.frame()
	spec <- unique(to_write$SPECIES)[[1]]
	print(spec)

	env = env_list[[spec]] %>% select("ID" = "Rxn.name",
																							 "SPECIES" = "Organism",
																							 "DESCR" = "Description",
																							 "REACTION" = "Equation",
																							 "GPR" = "Genes",
																							 "REF" = "References",
																							 "LOWER.BOUND" = "LB",
																							 "UPPER.BOUND" = "UB",
																							 "OBJECTIVE" = "Objective") %>% as.data.frame()
	env$DESCR = paste0(env$DESCR,"_exchange")
	env$OBJECTIVE=0
	to_write = rbind(to_write, env)
	#colnames(to_write) <- c("ID","SPECIES","DESCR","REACTION","GPR","REF")
	# to_write <- to_write %>% mutate(LOWER.BOUND = NA,
	# 																UPPER.BOUND = 1000,
	# 																OBJECTIVE = NA) %>% as.data.frame()
	to_write$GPR <- gsub("α|ε", "", to_write$GPR)
	to_write$GPR <- gsub(",[[:blank:]][[:alnum:]]"," or ",to_write$GPR)
	to_write$GPR <- gsub(",[[:blank:]]$","",to_write$GPR)

	to_write$GPR <- sapply(to_write$GPR, janitor::make_clean_names, use_make_names = F, case= "none")
	to_write$GPR <- gsub("or", "", to_write$GPR)
	to_write$GPR <- gsub("[[:blank:]][[:blank:]]", " ", to_write$GPR)
	to_write$GPR <- gsub("[[:blank:]][[:alnum:]]*", "", to_write$GPR)


	out_fp = paste0("/projectnb2/talbot-lab-data/metabolic_models/curated_models/N_cycle/",spec,"_v2.xml")

	fixed <- repair_SBML_model(modelData = to_write,
														 modelID = spec,
														 boundary = "e")

	# fixed[grepl("M_cyt550_c", "15991" fixed)]
	# fixed <- gsub("M_cyt550e_c","MNXM731975", fixed)

	writeLines(text = fixed, con = out_fp, sep = "\n")
}


metanetx_chem <- data.table::fread("/projectnb2/talbot-lab-data/metabolic_models/scripts/metabolic_model_curation/reference_data/chem_xref.tsv", fill=T, skip=351, nThread = 10)


kegg_id = "C00126"
all_mets_orig$metanet_id <- NA
for (i in 1:nrow(all_mets_orig)){
  kegg_id = all_mets_orig[i,]$Metabolite.KeggID %>% str_trim()
  if (is.na(kegg_id)) next()
  print(kegg_id)
  metanetx_kegg = metanetx_chem[grepl("kegg", metanetx_chem$`#source`),]
  metanet_info <- metanetx_kegg[grepl(kegg_id, metanetx_kegg$`#source`),]
  metanet_id <- unique(metanet_info$ID)
  print(metanet_id)
  all_mets_orig[i,]$metanet_id = metanet_id
}
write.csv(all_mets_orig, "/projectnb2/talbot-lab-data/metabolic_models/curated_models/N_cycle/metanet_key.csv")

ChEBI[grepl("ferrocytochrome",ChEBI$FORMULA)]

library("BiGGR")


fp <- '/projectnb2/talbot-lab-data/metabolic_models/curated_models/N_cycle/metanet_cobra_sbml/NDE5.COBRA-sbml3.xml'
out_fp <- '/projectnb2/talbot-lab-data/metabolic_models/curated_models/N_cycle/metanet_cobra_sbml/NDE_comps_combined.xml'

nde = readLines(fp)
to_write$GPR <- gsub("MNXC2", "MNXC19", to_write$GPR)


writeLines(text = fixed, con = out_fp, sep = "\n")

