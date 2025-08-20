# Within-site heterogeneity (core_sd)
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")


sum.all.fg <- readRDS("./data/summary/fg_summaries.rds")
sum.fg <- sum.all.fg$full_uncertainty$summary_df %>% filter(!grepl("beta", rowname)) 

sum.all.tax <- readRDS("./data/summary/taxa_summaries.rds")
sum.tax <- sum.all.tax$summary_df %>% filter(!grepl("beta", rowname))

sum.div.all <- readRDS("./data/summary/div_summaries.rds")
sum.div <- sum.div.all$summary_df %>% filter(!grepl("beta", rowname))
df <- data.table::rbindlist(list(sum.tax, 
																 sum.fg, 
																 sum.div), fill=TRUE)

df_refit <- df %>% filter(time_period=="refit") 

core_sd <- df_refit %>% filter(grepl("core", rowname))


ggplot(core_sd,
			 aes(x = fg_cat,y = Mean)) +
	geom_jitter(aes(color = group), size = 5, height = 0, width=.1, alpha = .5,
		shape = 16, show.legend = F) + #ylim(c(0,.9)) +
	labs(col = "Parameter", title = "Core_sd size") + 
	xlab("Parameter")+ 
	ylab(NULL) +
	facet_grid(rows = vars(model_name),  drop = T,
						 scales = "free_x", space = "free_x") + #,strip.position="bottom",nrow=2) +
	theme_bw() + theme(#axis.ticks.x=element_blank(),
		text = element_text(size = 16),
		axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
			angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=22,face="bold")
	)  + geom_text(aes(label = outlier), na.rm = TRUE, hjust = -0.3) 

core_sd <- core_sd %>% 
	mutate(outlier = ifelse(core_sd > .15, taxon, as.numeric(NA))) 
