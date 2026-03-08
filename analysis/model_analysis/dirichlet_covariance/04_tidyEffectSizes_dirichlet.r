# Combine & separately save model parameter and effect size estimates (beta covariates) from all models

source("source.R")
pacman::p_load(stringr, forestplot, gridExtra)

# Read in summaries and combine into fewer dfs for parameter effects

# Functional groups

dirichlet_summaries <- readRDS(here("data/summary/dirichlet_regression_summaries.rds"))

sum.all <- dirichlet_summaries$summary_df %>%
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

# 1. Linear model beta (covariate) effects
# Our updated summarize script cleanly labels covariates in the 'beta' column
beta_effects <- df %>%
	filter(!is.na(beta) & beta != "UNKNOWN")

beta_effects$beta <- ordered(beta_effects$beta,
                             levels = c("sin", "cos",
                                        "Ectomycorrhizal trees",
                                        "LAI", "pC", "pH",
                                        "Temperature", "Moisture"))

levels(beta_effects$beta)[levels(beta_effects$beta)=="Ectomycorrhizal trees"] <- "Ectomycorrhizal\ntrees"

# 2. Temporal persistence (rho) effects
# Isolate the AR(1) parameter into its own dataframe for separate analysis
rho_effects <- df %>%
	# Safely check parameter or rowname depending on how extract_summary_row formatted it
	filter(if("parameter" %in% names(.)) grepl("rho", parameter) else grepl("rho", rowname))

# Save the tidy outputs
saveRDS(beta_effects, here("data", "summary/dirichlet_predictor_effects.rds"))
saveRDS(rho_effects, here("data", "summary/dirichlet_rho_effects.rds"))
