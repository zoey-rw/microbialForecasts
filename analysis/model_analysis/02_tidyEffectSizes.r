# Combine effect size estimates (beta covariates) from all model
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
pacman::p_load(stringr, forestplot, gridExtra) 

# Read in summaries
sum.all.fg <- readRDS("./data/summary/fg_summaries.rds")
sum.fg <- sum.all.fg$full_uncertainty$summary_df %>% mutate(rank = "functional_group")

sum.all.tax <- readRDS("./data/summary/taxa_summaries.rds")
sum.tax <- sum.all.tax$summary_df %>% filter(taxon != "other")

sum.div.all <- readRDS("./data/summary/div_summaries.rds")
sum.div <- sum.div.all$summary_df %>% 
	mutate(rank = paste0("diversity_", group))

df <- data.table::rbindlist(list(sum.tax, 
																 sum.fg, 
																 sum.div), fill=TRUE)


# Add prettier data values 
df$pretty_group <- ifelse(df$group=="16S", "Bacteria", "Fungi")
df$pretty_name <- recode(df$rank, !!!pretty_rank_names)
df$only_rank <- sapply(str_split(df$rank, "_",  n = 2), `[`, 1)
df$only_rank <- ordered(df$only_rank, levels = c("genus",
																																			 "family",
																																			 "order", 
																																			 "class",
																																			 "phylum", "functional", "diversity"))
df$rank <- ordered(df$rank, levels = c("genus_bac","genus_fun",
																														 "family_bac","family_fun",
																														 "order_bac", "order_fun",
																														 "class_bac", "class_fun",
																														 "phylum_bac","phylum_fun",
																														 "functional_group", "diversity_16S", "diversity_ITS"))
df$pretty_name <- ordered(df$pretty_name, levels = c("Genus",
																																					 "Family",
																																					 "Order", 
																																					 "Class",
																																					 "Phylum", "Functional group", "Diversity"))


# Filter by effect type

# Linear model beta (covariate) effects
beta_effects <- df %>% filter(grepl("beta", rowname)) 
beta_effects$beta <- ordered(beta_effects$beta, levels = c("sin", "cos", 
																			 "Ectomycorrhizal trees",
																			 "LAI",
																			 "pC", 
																			 "pH",
																			 "Temperature", 
																			 "Moisture","rho"))
levels(beta_effects$beta)[levels(beta_effects$beta)=="Ectomycorrhizal trees"] <- "Ectomycorrhizal\ntrees"
saveRDS(beta_effects, "./data/summary/all_fcast_effects.rds")

# Site effects
site_effects <- df %>% filter(grepl("site", rowname)) 
saveRDS(site_effects, "./data/summary/site_effects.rds")
