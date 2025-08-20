# Combine & separately save model parameter and effect size estimates (beta covariates) from all models

source("source.R")
pacman::p_load(stringr, forestplot, gridExtra)

# Read in summaries and combine into fewer dfs for parameter effects

# Functional groups

dirichlet_summaries <- readRDS(here("data/summary/dirichlet_regression_summaries.rds"))

sum.all <- dirichlet_summaries$summary_df %>% filter(time_period == "2015-11_2018-01") %>% 
	mutate(tax_rank = rank,
				 time_period = recode(time_period, !!!microbialForecast:::date_recode))
df <- sum.all %>%
	mutate(pretty_group = ifelse(group %in% c("16S","bac"), "Bacteria", "Fungi"))

# Add prettier data values
df$pretty_name <- recode(df$rank_only, !!!microbialForecast:::pretty_rank_names) %>%
	ordered(levels = c("Genus","Family","Order","Class","Phylum","Functional group","Diversity"))

df$only_rank <- sapply(str_split(df$rank_only, "_",  n = 2), `[`, 1) %>%
	ordered(levels = c("genus","family","order","class","phylum","functional","diversity"))

# df$tax_rank <- ordered(df$tax_rank, levels = c("genus_bac","genus_fun",
# 																														 "family_bac","family_fun",
# 																														 "order_bac", "order_fun",
# 																														 "class_bac", "class_fun",
# 																														 "phylum_bac","phylum_fun",
# 																														 "functional_group", "diversity_16S", "diversity_ITS"))


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
saveRDS(beta_effects, here("data", "summary/dirichlet_predictor_effects.rds"))
