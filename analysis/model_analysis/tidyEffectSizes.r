# Combine effect size estimates (beta covariates) from all model
library(stringr)
library(forestplot)
library(ggplot2)
library(gridExtra)
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/helperFunctions.r")

sum.all.fg <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/fg_summaries.rds")
sum.fg <- sum.all.fg$full_uncertainty$summary_df %>% filter(grepl("beta", rowname)) 
sum.all.fg_refit <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/refit_fg_summaries.rds")
sum.fg_refit <- sum.all.fg_refit$full_uncertainty$summary_df %>% filter(grepl("beta", rowname)) 

sum.all.tax <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/allrank_summaries.rds")
sum.tax <- sum.all.tax$summary_df %>% filter(taxon != "other" & grepl("beta", rowname)) 
sum.all.tax_refit <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/refit_allrank_summaries.rds")

# TODO: put this into summary code
sum.tax_refit <- sum.all.tax_refit$summary_df %>% filter(taxon != "other" & grepl("beta", rowname)) 

sum.div.all <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/div_summaries.rds")
sum.div <- sum.div.all$summary_df %>% filter(uncert == "full_uncertainty" & grepl("beta", rowname))
sum.div.all_refit <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/refit_div_summaries.rds")
sum.div_refit <- sum.div.all_refit$summary_df %>% filter(uncert == "full_uncertainty" & grepl("beta", rowname))

df <- data.table::rbindlist(list(sum.tax, 
																 sum.fg, 
																 sum.div, 
																 sum.fg_refit, 
																 sum.tax_refit,
																 sum.div_refit), fill=TRUE)

df$pretty_group <- ifelse(df$group=="16S", "Bacteria", "Fungi")

df$pretty_name <- recode(df$rank,
												 "genus_bac" = "Genus",
												 "family_bac" = "Family",
												 "order_bac" = "Order", 
												 "class_bac" = "Class", 
												 "phylum_bac" = "Phylum",
												 "genus_fun" = "Genus",
												 "family_fun" = "Family",
												 "order_fun" = "Order", 
												 "class_fun" = "Class", 
												 "phylum_fun" = "Phylum",
												 "functional_group" = "Functional group")

df$only_rank <- sapply(str_split(df$rank, "_",  n = 2), `[`, 1)
df$only_rank <- ordered(df$only_rank, levels = c("genus",
																								 "family",
																								 "order", 
																								 "class",
																								 "phylum", "functional"))
df$rank <- ordered(df$rank, levels = c("genus_bac","genus_fun",
																			 "family_bac","family_fun",
																			 "order_bac", "order_fun",
																			 "class_bac", "class_fun",
																			 "phylum_bac","phylum_fun",
																			 "functional_group"))

df$beta <- ordered(df$beta, levels = c("Ectomycorrhizal trees",
																			 "Plant species richness",
																			 "pC", 
																			 "pH",
																			 "Temperature", "Moisture"))
levels(df$beta)[levels(df$beta)=="Ectomycorrhizal trees"] <- "Ectomycorrhizal\ntrees"
levels(df$beta)[levels(df$beta)=="Plant species richness"] <- "Plant species\nrichness"


saveRDS(df, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/all_fcast_effects.rds")