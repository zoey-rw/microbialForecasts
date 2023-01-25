source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
pacman::p_load(stringr, forestplot, gridExtra, ggpubr, rstatix)
sum.all <- readRDS("./data/summary/all_fcast_effects.rds")

#pass_filter <- readRDS(here("data/summary/tax_filter_pass.rds"))

sum.all_params <- sum.all %>% filter(beta %in% c("sin","cos") & !grepl("other", taxon))
seas_vals <- sum.all_params %>% pivot_wider(id_cols = c("taxon","model_name","fcast_type","time_period",
																												"pretty_name","pretty_group","rank","only_rank"),
																						names_from = beta,
																						values_from = "Mean")
# Couldn't figure out how to vectorize.
out <- list()
for (i in 1:nrow(seas_vals)) {
#	out[[i]] <- sin_cos_to_seasonality(seas_vals$sin[[i]], seas_vals$cos[[i]])
	out[[i]] <- sin_cos_to_seasonality(exp(seas_vals$sin[[i]]), exp(seas_vals$cos[[i]]))

}
out <- rbindlist(out)
seas_vals <- cbind.data.frame(seas_vals, out) %>%
	mutate(beta = "residual\nseasonal\namplitude", effSize= amplitude) %>%
	select(-c(sin, cos, max, amplitude_orig,amplitude)) %>%
	filter(fcast_type != "Diversity")


df_cal <- sum.all %>%
	filter(beta %in% beta_names & #(workaround since the model names aren't all saving)
				 	#model_name == "all_covariates" &
				 	!grepl("other", taxon) &
				 	#time_period == "2015-11_2018-01" &
				 	!time_period %in% c("2015-11_2020-01"))

df_fg <- df_cal %>% filter(fcast_type == "Functional")

df_tax <- df_cal %>% filter(fcast_type == "Taxonomic")
#df_tax_pass <- df_tax %>% merge(pass_filter, all.y=T)

df_cal_fg_tax <- rbind(df_tax, df_fg) %>% filter(!beta %in% c("sin","cos")) %>%
	rbind(seas_vals %>% filter(model_name=="all_covariates"), fill=T)



# Subset to converged only
keep = readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/converged_taxa_list.rds")
unconverged <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/unconverged_taxa_list.rds")
df_cal_fg_tax_converged = df_cal_fg_tax %>%
	mutate(taxon_model_rank = paste(taxon, model_name, rank)) #%>%
	#filter(taxon_model_rank %in% keep$taxon_model_rank) %>%
	#filter(!taxon_model_rank %in% unconverged$taxon_model_rank)
seas_vals_converged = seas_vals %>%
	mutate(taxon_model_rank = paste(taxon, model_name, rank)) #%>%
	#filter(!taxon_model_rank %in% unconverged$taxon_model_rank)

plotting_df = df_cal_fg_tax_converged %>%
	filter(!beta %in% c("sin","cos")) %>%
	filter(fcast_type == "Taxonomic") %>%
	mutate(exp_Mean = exp(Mean),
				 beta = "residual\nseasonal\namplitude", effSize= amplitude)

tukey_list <- list()
beta_names <- c("Ectomycorrhizal\ntrees", "LAI", "pC",
	"pH", "Temperature", "Moisture","residual\nseasonal\namplitude")
for(b in beta_names) {
	x = "only_rank"
	y = "effSize"
	df <- plotting_df[which(plotting_df$beta == b),]
	new.df <- cbind.data.frame(x = df$only_rank, y = df$effSize)
	abs_max <- max(new.df[,"y"], na.rm = T)
	maxs <- new.df %>% group_by(x) %>%
		summarise(tot=max(y, na.rm=T)+ 0.2 * abs_max)
	Tukey_test <- aov(y ~ x, data=new.df) %>%
		agricolae::HSD.test("x", group=TRUE) %>%
		.$groups %>%
		as_tibble(rownames="x") %>%
		rename("Letters_Tukey"="groups") %>%
		dplyr::select(-y) %>%
		left_join(maxs, by="x") %>%
		rename("only_rank"="x")
	Tukey_test$beta <- b
	tukey_list[[b]] <- Tukey_test
}

tukey_list_tax_fg <- data.table::rbindlist(tukey_list)



###### # Plot with Tukey, effect sizes across ranks ----
ranks_beta_plot <- ggplot(plotting_df) +
	geom_jitter(aes(x = only_rank,y = effSize,
									color = pretty_group),
							width=.3, height = 0, size=4, alpha = .5, show.legend = F) +
	labs(title = "Absolute effect size") +
	theme_minimal(base_size = 18) +
	xlab("Rank")+
	ylab(NULL) +
	facet_grid(rows = vars(beta), cols = vars(pretty_group), drop = T,
						 scales = "free", space = "free_x") +
	theme(axis.text.x=element_text(
		angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=22,face="bold"),
		strip.text.y = element_text(size=12,face="bold"))


###### # Plot with Tukey, effect sizes across ranks ----
ranks_beta_plot <- ggplot(plotting_df, aes(x = only_rank,y = effSize,
																					 color = pretty_group)) +
	geom_jitter(
							width=.3, height = 0, size=4, alpha = .5, show.legend = F) +
	geom_violin(draw_quantiles = c(.5), show.legend = F) +
	labs(title = "Absolute effect size") +
	theme_minimal(base_size = 18) +
	xlab("Rank")+
	ylab(NULL) +
	facet_grid(cols = vars(beta),
						 rows = vars(pretty_group), drop = T,
						 scales = "free", space = "free") +
	theme(axis.text.x=element_text(
		angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=22,face="bold"),
		strip.text.y = element_text(size=12,face="bold"))


ranks_beta_plot <- ranks_beta_plot  +
	stat_compare_means(data = plotting_df, aes(label = paste0("p = ", ..p.format..)),
										 method = "anova", size=5, label.y.npc = .5) +
	geom_text(data = tukey_list_tax_fg,
																							 aes(x = only_rank, y = tot,
																							 		label = Letters_Tukey), show.legend = F, color = 1, size =4)
ranks_beta_plot






ggplot(plotting_df,
			 aes(x = only_rank,y = effSize,
			 																	 color = pretty_group)) +
	geom_jitter(width=.3, height = 0, size=4, alpha = .5, show.legend = F) +
	labs(title = "Absolute effect size") +
	theme_minimal(base_size = 18) +
	facet_grid(rows = vars(beta),
						 cols = vars(pretty_group), drop = T,
						 scales = "free", space = "free_x") +
	stat_compare_means(aes(label = paste0("p = ", ..p.format..))) +
	geom_text(data = tukey_list_tax_fg,
						aes(x = only_rank, y = tot,
								label = Letters_Tukey), show.legend = F, color = 2, size =6)

anova_fit <- plotting_df %>%
	group_by(pretty_group, beta) %>%
	anova_test(effSize ~ only_rank) %>%
	add_significance()

#### Run Tukey ###
tukey <- plotting_df %>%
	group_by(pretty_group, beta) %>%
	tukey_hsd(effSize ~ only_rank) %>%
	add_significance() %>%
	add_xy_position()



ggplot(plotting_df, aes(x = only_rank,y = effSize,
												color = pretty_group)) +
	geom_jitter(width=.3, height = 0, size=4, alpha = .5, show.legend = F) +
	labs(title = "Absolute effect size") +
	theme_minimal(base_size = 18) +
	facet_grid(rows = vars(beta),
						 cols = vars(pretty_group), drop = T,
						 scales = "free", space = "free_x") +
	stat_pvalue_manual(tukey,
										 hide.ns = T)+
	labs(subtitle = get_test_label(anova_fit,
																 detailed = TRUE),
			 caption = get_pwc_label(tukey))



## RESIDUAL SEASONALITY ONLY

ggplot(plotting_df %>% filter(beta=="residual\nseasonal\namplitude"),
			 aes(x = only_rank,y = effSize,
			 		color = pretty_group)) +
	geom_jitter(width=.3, height = 0, size=4, alpha = .5, show.legend = F)  +
	geom_violin(draw_quantiles = c(.5), show.legend = F) +
	labs(title = "Absolute effect size") +
	theme_minimal(base_size = 18) +
	facet_grid(rows = vars(pretty_group), drop = T,
						 scales = "free", space = "free_x") +
	stat_compare_means(aes(label = paste0("p = ", ..p.format..)))
