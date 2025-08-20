source("source.R")
pacman::p_load(stringr, forestplot, gridExtra, ggpubr, rstatix)
#remotes::install_github('rpkgs/gg.layers')


sum.all <- readRDS(here("data", "summary/predictor_effects.rds"))
seasonal_amplitude_in = readRDS(here("data/summary/seasonal_amplitude.rds"))
cycl_vals_scores = seasonal_amplitude_in[[6]]
cycl_vals_y = seasonal_amplitude_in[[1]]
#pass_filter <- readRDS(here("data/summary/tax_filter_pass.rds"))


converged <- readRDS(here("data/summary/weak_converged_taxa_list.rds"))
converged_strict <- readRDS(here("data/summary/converged_taxa_list.rds"))

beta_names <- c("Ectomycorrhizal\ntrees", "LAI", "pC",
								"pH", "Temperature", "Moisture"#,"residual\nseasonal\namplitude"
								)

df_cal <- sum.all %>%
	filter(beta %in% beta_names & #(workaround since the model names aren't all saving)
				 	model_name %in% c("env_cov","env_cycl") &
				 	!grepl("other", taxon) &
				 	time_period == "2015-11_2018-01" &
				 	!time_period %in% c("2015-11_2020-01")) %>% mutate(rank_only=only_rank)
df_cal$beta <- droplevels(df_cal$beta)

df_cal_fg_tax <- df_cal
# df_cal_fg_tax <-  rbindlist(list(df_cal,cycl_vals_scores), fill=T) %>%
# 	filter(model_id %in% converged)

df_cyclical <- merge(cycl_vals_y, df_cal) %>%
	filter(model_id %in% converged)

plotting_df = df_cal_fg_tax %>%
	filter(!beta %in% c("sin","cos") & model_name %in% c("env_cov")) %>%
	filter(time_period == "2015-11_2018-01")
	# mutate(exp_Mean = exp(Mean),
	# 			 beta = "residual\nseasonal\namplitude", effSize= amplitude)

tukey_list <- list()


for (pretty_group in c("Bacteria","Fungi")){
	#group_df <- sum.all %>% filter(group == !!group)
	tukey_group_list = list()
	group_df = plotting_df %>% filter(pretty_group==!!pretty_group)
for(b in beta_names) {
	x = "rank_only"
	y = "effSize"
	df <- group_df[which(group_df$beta == b),]
	
	new.df <- cbind.data.frame(x = df$rank_only, y = df$effSize)
	new.df <- new.df[!is.na(new.df$y), ]
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
		rename("rank_only"="x")
	Tukey_test$beta <- b
	Tukey_test$pretty_group <- pretty_group
	tukey_group_list[[b]] <- Tukey_test
}
	tukey_list[[pretty_group]] = data.table::rbindlist(tukey_group_list)
}

tukey_list_tax_fg <- data.table::rbindlist(tukey_list)

# Initialize test result variables


###### # Plot with Tukey, effect sizes across ranks ----
ranks_beta_plot <- ggplot(plotting_df, aes(x = rank_only,y = effSize,
																					 color = pretty_group)) +
	geom_jitter(aes(shape=as.factor(significant)),
							width=.1, height = 0, size=3, alpha = .2, show.legend = F) +
	geom_violin(draw_quantiles = c(.5), show.legend = F) +
	labs(title = "Absolute effect size") +
	theme_minimal(base_size = 18) +
	xlab("Rank")+
	ylab(NULL) +
	facet_grid(cols = vars(beta),
						 rows = vars(pretty_group),
						 drop = T,
						 scales = "free", space = "free") +
	theme(axis.text.x=element_text(
		angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=22),
		strip.text.y = element_text(size=12)) +
	scale_shape_manual(values = c(21, 16), name = NULL,																																						 labels = c("Not significant","Significant")) + guides(color="none")


ranks_beta_plot <- ranks_beta_plot  +
	stat_compare_means(data = plotting_df, aes(label = paste0("p = ", ..p.format..)),
										 method = "anova", size=5, label.y.npc = .5) +
	geom_text(data = tukey_list_tax_fg,
																							 aes(x = rank_only, y = tot, color =pretty_group,
																							 		label = Letters_Tukey), show.legend = F, color = 1, size =4)
ranks_beta_plot






ggplot(plotting_df,
			 aes(x = rank_only,y = effSize,
			 																	 color = pretty_group)) +
	geom_jitter(width=.3, height = 0, size=4, alpha = .5, show.legend = F) +
	labs(title = "Absolute effect size") +
	theme_minimal(base_size = 18) +
	facet_grid(rows = vars(beta),
						 cols = vars(pretty_group), drop = T,
						 scales = "free", space = "free_x") +
	stat_compare_means(aes(label = paste0("p = ", ..p.format..))) +
	geom_text(data = tukey_list_tax_fg,
						aes(x = rank_only, y = tot,
								label = Letters_Tukey), show.legend = F, color = 2, size =6)

anova_fit <- plotting_df %>%
	group_by(pretty_group, beta) %>%
	anova_test(effSize ~ rank_only) %>%
	add_significance()

#### Run Tukey ###
tukey <- plotting_df %>%
	group_by(pretty_group, beta) %>%
	tukey_hsd(effSize ~ rank_only) %>%
	add_significance() %>%
	add_xy_position()

ggplot(plotting_df, aes(x = rank_only,y = effSize,
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



ggplot(plotting_df %>% filter(beta=="residual\nseasonal\namplitude"),
			 aes(x = only_rank,y = effSize,
			 		color = pretty_group)) +
	geom_jitter(width=.3, height = 0, size=4, alpha = .5, show.legend = F)  +
	geom_violin(draw_quantiles = c(.5), show.legend = F) +
	labs(title = "Absolute effect size") +
	theme_minimal(base_size = 18) +
	facet_grid(rows = vars(pretty_group),
						 cols = vars(model_name),
						 drop = T,
						 scales = "free", space = "free_x") +
	stat_compare_means(aes(label = paste0("p = ", ..p.format..)))

ggplot(df_cyclical %>%
			 	filter(fcast_type != "Diversity" & time_period == "2015-11_2018-01"),
			 aes(x = only_rank,y = effSize,
												color = pretty_group)) +

	geom_violin(draw_quantiles = c(.5), show.legend = F) +
	geom_jitter(
		width=.3, height = 0, size=4, alpha = .5, show.legend = F) +
	labs(title = "Seasonal amplitude") +
	theme_minimal(base_size = 18) +
	xlab("Rank")+
	ylab(NULL) +
	facet_grid(
						 rows = vars(pretty_group), drop = T,
						 scales = "free", space = "free") +
	theme(axis.text.x=element_text(
		angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=22,face="bold"),
		strip.text.y = element_text(size=12,face="bold"))
