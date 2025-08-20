# Visualize effect size estimates (beta covariates) from all model
source("source.R")
pacman::p_load(stringr, forestplot, gridExtra, ggpubr) 

sum.all <- readRDS(here("data/summary/all_fcast_effects.rds"))
df_refit <- sum.all %>% filter((time_period == "2013-06_2017-01" & 
																	fcast_type == "Taxonomic") |
															 	(time_period == "2015-11_2018-01")) %>% 
	filter(model_name == "all_covariates")



tukey_list <- list()
group_out <- list()
beta_names <- c(#"sin", "cos", 
	"Ectomycorrhizal\ntrees", "LAI", "pC", 
	"pH", "Temperature", "Moisture")

for (group in c("16S","ITS")){
	group_df <- sum.all %>% filter(group == !!group)
	
for (b in beta_names) {
	beta_df <- group_df %>% filter(beta == !!b)
	out <- tukey(beta_df$pretty_name, beta_df$effSize, y.offset = .1) 
	tukey_list[[b]] <- out %>% mutate(beta = !!b)
}
	group_out[[group]] <- data.table::rbindlist(tukey_list) %>% mutate(group = !!group)
}
tukey_list_tax <- data.table::rbindlist(group_out)
colnames(tukey_list_tax)[1] <- "pretty_name"
tukey_list_tax$pretty_group <- ifelse(tukey_list_tax$group == "16S", "Bacteria", "Fungi")
tukey_list_tax$pretty_name <- ordered(tukey_list_tax$pretty_name, levels = c("Genus",
																										 "Family",
																										 "Order", 
																										 "Class",
																										 "Phylum", "Functional group", "Diversity"))


ggplot(df_refit %>% filter(!beta %in% c("sin","cos"))) +
	geom_jitter(aes(x = pretty_name,y = effSize,
									fill = pretty_group), 
							shape=21, color="black", width=.1, height = 0, size=4) + 
	labs(title = "Absolute effect size") + 
	xlab("Rank")+ 
	ylab(NULL) +
	facet_grid(rows = vars(beta), cols = vars(pretty_group), drop = T,
						 scales = "free", 
						 space = "free_x") + #,strip.position="bottom",nrow=2) +
	theme_bw() + theme(#axis.ticks.x=element_blank(),
		text = element_text(size = 16),
		axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
			angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=22,face="bold"),
		strip.text.y = element_text(size=12,face="bold")
	) + 
	geom_smooth(aes(x = as.numeric(pretty_name), y = effSize), show.legend = F) +
	geom_text(data = tukey_list_tax, 
						aes(x = as.numeric(pretty_name), y = tot + .2, 
								label = Letters_Tukey), show.legend = F, color = 1, size =5) 
