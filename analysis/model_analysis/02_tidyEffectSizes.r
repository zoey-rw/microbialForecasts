# Combine effect size estimates (beta covariates) from all model
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
pacman::p_load(stringr, forestplot, gridExtra) 

# Read in summaries
sum.all.fg <- readRDS("./data/summary/fg_summaries.rds")
sum.fg <- sum.all.fg$full_uncertainty$summary_df %>% filter(grepl("beta", rowname)) 

sum.all.tax <- readRDS("./data/summary/taxa_summaries.rds")
sum.tax <- sum.all.tax$summary_df %>% filter(taxon != "other" & grepl("beta", rowname))

sum.div.all <- readRDS("./data/summary/div_summaries.rds")
sum.div <- sum.div.all$summary_df %>% filter(grepl("beta", rowname)) %>% 
	mutate(rank = paste0("diversity_", group))

df <- data.table::rbindlist(list(sum.tax, 
																 sum.fg, 
																 sum.div), fill=TRUE)

df$pretty_group <- ifelse(df$group=="16S", "Bacteria", "Fungi")

df$pretty_name <- recode(df$rank, !!!pretty_rank_names)

df$only_rank <- sapply(str_split(df$rank, "_",  n = 2), `[`, 1)
df$only_rank <- ordered(tolower(df$only_rank), levels = c("genus",
																								 "family",
																								 "order", 
																								 "class",
																								 "phylum", "functional", "diversity"))
df$rank <- ordered(df$rank, levels = c("genus_bac","genus_fun",
																			 "family_bac","family_fun",
																			 "order_bac", "order_fun",
																			 "class_bac", "class_fun",
																			 "phylum_bac","phylum_fun",
																			 "functional_group", "diversity"))

df$beta <- ordered(df$beta, levels = c("sin", "cos", 
																			 "Ectomycorrhizal trees",
																			 "LAI",
																			 "pC", 
																			 "pH",
																			 "Temperature", 
																			 "Moisture","rho"))
levels(df$beta)[levels(df$beta)=="Ectomycorrhizal trees"] <- "Ectomycorrhizal\ntrees"


saveRDS(df, "./data/summary/all_fcast_effects.rds")
