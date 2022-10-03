# Combine effect size estimates (beta covariates) from all workflows
# This script must be run before other analysis scripts

source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
pacman::p_load(stringr, forestplot, gridExtra)

# Read in summaries and combine into fewer dfs for parameter effects

# Functional groups
sum.all.fg <- readRDS(here("data", paste0("summary/beta_fg_summaries_20151101_20180101.rds")))
sum.fg <- sum.all.fg$summary_df %>% mutate(rank = "functional_group")

sum.all.fg_refit <- readRDS(here("data", paste0("summary/beta_fg_summaries_20151101_20200101.rds")))
sum.fg_refit <- sum.all.fg_refit$summary_df %>% mutate(rank = "functional_group")

# Taxonomic groups
sum.all.tax <- readRDS(here("data/summary/single_taxon_summaries_201511_201801.rds"))
sum.tax <- sum.all.tax$summary_df %>% mutate(rank_taxon = rank, rank = rank_name)

sum.all.tax_refit <- readRDS(here("data/summary/single_taxon_summaries_201511_202001.rds"))
sum.tax_refit <- sum.all.tax_refit$summary_df  %>% mutate(rank_taxon = rank, rank = rank_name)

# Diversity models
sum.div.all <- readRDS(here("data", "/summary/div_summaries_effects.rds"))
sum.div <- sum.div.all %>%
	filter(time_period %in% c("2015-11_2018-01","2015-11_2020-01"))
sum.div$group <- sapply(str_split(sum.div$taxon, "_",  n = 2), `[`, 2)
sum.div$rank <- paste0(sum.div$rank, "_", sum.div$group)


df <- data.table::rbindlist(list(sum.tax,
																 sum.tax_refit,
																 sum.fg,
																 sum.fg_refit,
																 sum.div), fill=TRUE) %>%
	mutate(pretty_group = ifelse(group == "16S", "Bacteria", "Fungi"))


# Add prettier data values
df$pretty_name <- recode(df$rank, !!!microbialForecast:::pretty_rank_names)
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


# For saving: filter by effect type

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
saveRDS(beta_effects, here("data", "summary/all_fcast_effects.rds"))

# Site effects
site_effects <- df %>% filter(grepl("site", rowname))
saveRDS(site_effects, here("data", "summary/site_effects.rds"))
