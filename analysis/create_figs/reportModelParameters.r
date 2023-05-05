# Summarize all the estimated model parameters for paper supplement.

source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")
library(kableExtra)


sum.in <- readRDS(here("data", paste0("summary/logit_beta_regression_summaries.rds")))
sum.all <- sum.in$summary_df  %>% mutate(tax_rank = rank,
																				 time_period = recode(time_period, !!!microbialForecast:::date_recode))
df <- sum.all %>%
	mutate(pretty_group = ifelse(group %in% c("16S","bac"), "Bacteria", "Fungi"))

# Add prettier data values
df$pretty_name <- recode(df$rank_only, !!!microbialForecast:::pretty_rank_names) %>%
	ordered(levels = c("Genus","Family","Order","Class","Phylum","Functional group","Diversity"))

df$only_rank <- sapply(str_split(df$rank_only, "_",  n = 2), `[`, 1) %>%
	ordered(levels = c("genus","family","order","class","phylum","functional","diversity"))

df = df %>% mutate(parameter_name = ifelse(is.na(beta), rowname, beta))
df_wide = df  %>% filter(time_period ==  "2015-11_2018-01" & parameter_name != "site_effect") %>%
	pivot_wider(id_cols = c("model_id","taxon","model_name","fcast_type","pretty_group","only_rank"),
							names_from = parameter_name,
							values_from = "Mean", values_fn = list) %>% as.data.frame() %>%
	rename("Covariates included" = "model_name",
				 "Site effect variation" = "sig",
				 "Kingdom" = "pretty_group",
				 "Microbial group type" = "fcast_type",
				 "Group" = "taxon",
				 "Rank" = "only_rank")


kable(df_wide, "html") %>%
	kable_styling(bootstrap_options = c("striped", "hover")) %>%
	cat(., file = here("figures", "model_params.html"))


# Check number of plots per site - usually 10, some have fewer or more
aggregate(plotID ~ siteID, hindcast_data, function(x) length(unique(x))) 
