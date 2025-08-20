# Visualize effect size estimates (beta covariates) from all model
source("source.R")
pacman::p_load(stringr, forestplot, gridExtra, ggpubr)

sum.all <- readRDS(here("data/summary/predictor_effects.rds"))

#pass_filter <- readRDS(here("data/summary/tax_filter_pass.rds"))


converged <- readRDS(here("data/summary/weak_converged_taxa_list.rds"))
converged_strict <- readRDS(here("data/summary/converged_taxa_list.rds"))

seas_vals =readRDS(here("data/summary/seasonal_amplitude.rds"))

all_seas_vals = seas_vals[[5]] %>% filter(time_period == "2015-11_2018-01")

df_cal_fg_tax <- sum.all %>% filter(#(time_period == "2013-06_2017-01" & fcast_type == "Taxonomic") |
															 	(time_period == "2015-11_2018-01")) %>%
	filter(model_name == "env_cycl" & model_id %in% converged)  %>%
	filter(!beta %in% c("sin","cos"))
df_refit_fg_tax <- df_cal_fg_tax

# df_fg <- df_cal %>% filter(fcast_type == "Functional")
#
# df_tax <- df_cal %>% filter(fcast_type == "Taxonomic")
#df_tax_pass <- df_tax %>% merge(pass_filter, all.y=T)

# df_cal_fg_tax <- rbind(df_tax, df_fg) %>% filter(!beta %in% c("sin","cos")) %>%
# 	rbind(all_seas_vals %>% filter(model_name=="env_cycl"), fill=T)



####### Absolute effect sizes per rank ------
a <- ggplot(df_cal_fg_tax,
						aes(x = beta,y = effSize)) +
	geom_jitter(aes(shape = as.factor(significant),
		color = beta), size = 5, height = 0, width=.4, alpha = .3,
		shape = 16, show.legend = F) + ylim(c(0,.9)) +
	labs(col = "Parameter", title = "Absolute effect size") +
	xlab("Parameter")+
	ylab(NULL) +
	facet_grid(rows = vars(fcast_type), cols = vars(pretty_group), drop = F,
						 scales = "free_x", space = "free_x") + #,strip.position="bottom",nrow=2) +
	theme_bw() + theme(#axis.ticks.x=element_blank(),
		text = element_text(size = 16),
		axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
			angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=22,face="bold")
		#strip.text.y = element_text(size=24,hjust=0,vjust = 1,angle=180,face="bold")
	) + scale_shape_manual(values = c(21, 16), name = NULL,
										 labels = c("Not significant","Significant"))

###### Relative effect sizes per rank -----
b <- ggplot(df_cal_fg_tax,
						aes(x = beta,y = Mean)) +
	geom_jitter(aes(#shape = as.factor(significant),
		color = beta), size = 5, height = 0, width=.4, alpha = .3,
		shape = 16, show.legend = F) + ylim(c(-.9,.9)) +
	labs(col = "Parameter", title = "Effect size (refit)") +
	xlab("Parameter")+
	ylab(NULL) +
	geom_hline(yintercept = 0, linetype=2) +
	facet_grid(rows = vars(fcast_type), cols = vars(pretty_group), drop = T,
						 scales = "free_x", space = "free_x") + #,strip.position="bottom",nrow=2) +
	theme_bw() + theme(#axis.ticks.x=element_blank(),
		text = element_text(size = 16),
		axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
			angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=22,face="bold")
		#strip.text.y = element_text(size=24,hjust=0,vjust = 1,angle=180,face="bold")
	) + scale_shape_manual(values = c(21, 16), name = NULL,
												 labels = c("Not significant","Significant"))
#####
ggarrange(a,b)

#
# keep = readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/converged_taxa_list.rds")
# keep = keep[keep$median_gbr < 3,]
#
# unconverged <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/unconverged_taxa_list.rds")



###### # Plot with Tukey, fungi vs bacteria effect sizes -----
b_vs_f_fcast_type_plot <- ggplot(data=df_cal_fg_tax,
																 aes(x = pretty_group,
																 		color = pretty_group, y = effSize)) +

	geom_violin(draw_quantiles = c(.5), show.legend = F, color = 1) +
	geom_jitter(aes(#shape = as.factor(significant),
		color = pretty_group), size = 5, height = 0, width=.2, alpha = .2,
		shape = 16, show.legend = F) +
	#labs(col = "Parameter", title = "Effect sizes of environmental predictors ") +
	xlab("Kingdom")+
	ylab("Absolute effect size") +
	facet_grid(rows = vars(fcast_type),
						 cols = vars(beta), #as.table = T,
						 drop = T,
						 scales = "free") +
	theme_minimal() +
	theme(text = element_text(size = 22),
				axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
					angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=20)
	) + stat_compare_means(method="t.test", hide.ns = T, show.legend = F)  +
	scale_y_log10()# +coord_flip()
b_vs_f_fcast_type_plot  + stat_compare_means(
	aes(label = after_stat(p.signif)),
	method = "t.test",
	label.y.npc = .8, hide.ns = T, size=10, show.legend = F)

tukey_list <- list()
beta_names <- c(#"sin", "cos",
	"Ectomycorrhizal\ntrees", "LAI", "pC",
	"pH", "Temperature", "Moisture")
for(b in beta_names) {
	x = "pretty_group"
	y = "effSize"
	df <- df_tax[which(df_tax$beta == b),]
	new.df <- cbind.data.frame(x = df$pretty_group, y = df$effSize)
	abs_max <- max(new.df[,"y"], na.rm = T)
	maxs <- new.df %>% group_by(x) %>%
		summarise(tot=max(y, na.rm=T)+ 0.3 * abs_max)
	Tukey_test <- aov(y ~ x, data=new.df) %>%
		agricolae::HSD.test("x", group=TRUE) %>%
		.$groups %>%
		as_tibble(rownames="x") %>%
		rename("Letters_Tukey"="groups") %>%
		dplyr::select(-y) %>%
		left_join(maxs, by="x") %>%
		rename("pretty_group"="x")
	Tukey_test$beta <- b
	tukey_list[[b]] <- Tukey_test
}

tukey_list_tax <- data.table::rbindlist(tukey_list)
tukey_list_tax$fcast_type <- "Taxonomic"

tukey_list <- list()
for(b in beta_names) {
	x = "pretty_group"
	y = "effSize"
	df <- df_fg[which(df_fg$beta == b),]
	new.df <- cbind.data.frame(x = df$pretty_group, y = df$effSize)
	abs_max <- max(new.df[,"y"], na.rm = T)
	maxs <- new.df %>% group_by(x) %>%
		summarise(tot=max(y, na.rm=T)+ 0.3 * abs_max)
	Tukey_test <- aov(y ~ x, data=new.df) %>%
		agricolae::HSD.test("x", group=TRUE) %>%
		.$groups %>%
		as_tibble(rownames="x") %>%
		rename("Letters_Tukey"="groups") %>%
		dplyr::select(-y) %>%
		left_join(maxs, by="x") %>%
		rename("pretty_group"="x")
	Tukey_test$beta <- b
	tukey_list[[b]] <- Tukey_test
}
tukey_list_fg <- data.table::rbindlist(tukey_list)
tukey_list_fg$fcast_type <- "Functional"
tukey <- rbind(tukey_list_fg, tukey_list_tax)

tukey$beta <- order_betas(tukey$beta)

b_vs_f_fcast_type_plot <- b_vs_f_fcast_type_plot + geom_text(data = tukey, aes(x = pretty_group, y = tot+.1,
																										 label = Letters_Tukey), show.legend = F, color = 1, size = 5)
######
b_vs_f_fcast_type_plot # + ggtitle("Full dataset effect sizes")

###### # Plot with Tukey, effect sizes across ranks ----
ranks_beta_plot <- ggplot(df_cal_fg_tax %>% filter(!beta %in% c("sin","cos")),
													aes(x = only_rank,y = effSize,
															color = pretty_group)) +
	geom_jitter(aes(shape=as.factor(significant)),

							width=.3, height = 0, size=2, alpha = .5) +
	labs(title = "Absolute effect size") +
	geom_violin(draw_quantiles = c(.5), show.legend = F) +
	xlab("Rank")+
	ylab(NULL) +
	facet_grid(cols = vars(beta), rows = vars(pretty_group), drop = T,
						 scales = "free", space = "free_x") + #,strip.position="bottom",nrow=2) +
	theme_minimal() + theme(#axis.ticks.x=element_blank(),
		text = element_text(size = 16),
		axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
			angle = 320, vjust=1, hjust = -0.05)
	)  +
	scale_shape_manual(values = c(21, 16), name = NULL,																																						 labels = c("Not significant","Significant")) + guides(color="none")


tukey_list <- list()
beta_names <- c(#"sin", "cos",
	"Ectomycorrhizal\ntrees", "LAI", "pC",
	"pH", "Temperature", "Moisture")
for(b in beta_names) {
	x = "only_rank"
	y = "effSize"
	df <- df_refit_fg_tax[which(df_refit_fg_tax$beta == b),] #%>% filter(significant==1)
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

ranks_beta_plot <- ranks_beta_plot + geom_text(data = tukey_list_tax_fg,
														aes(x = only_rank, y = tot,
																label = Letters_Tukey), show.legend = F, color = 1, size =5)
######
ranks_beta_plot


tukey_eff_rank = df_refit_fg_tax %>% filter(!is.na(rank_only)) %>%
	group_by(rank_only, beta, .drop=T) %>%
	summarize(tukey(x = rank_only, y = effSize))


# View actual taxa and effects
df_refit_fg <- df_refit_fg_tax %>% filter(fcast_type == "Functional")
ggplot(data=df_refit_fg %>% filter(!beta %in% c("sin","cos")),
			 aes(x = fg_cat,y = Mean)) +
	geom_jitter(aes(color = beta, shape = as.factor(significant)), size = 5, height = 0, width=.1, alpha = .8) +
	labs(col = "Parameter", title = "Effect size (refit)") +
	xlab("Taxon")+
	ylab(NULL) +
	geom_hline(yintercept = 0, linetype=2) +
	facet_grid(#rows = vars(only_rank),
		cols = vars(pretty_group), drop = T,
		scales = "free", space = "free_x") + #,strip.position="bottom",nrow=2) +
	theme_bw() + theme(#axis.ticks.x=element_blank(),
		text = element_text(size = 16),
		axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
		#angle = 320, vjust=1, hjust = -0.05
		#),
		axis.title=element_text(size=22,face="bold")
		#strip.text.y = element_text(size=24,hjust=0,vjust = 1,angle=180,face="bold")
	) + scale_shape_manual(values = c(21, 16), name = NULL,
												 labels = c("Not significant","Significant")) +
	# ggrepel::geom_text_repel(data =df[df$beta_num==1, ], aes(label=taxon, x = taxon),
	# 												 nudge_y = -0.6, #direction="y",
	# 												 min.segment.length = Inf) +
	geom_hline(yintercept = 0, linetype=2) + ylim(c(-1,1))


##### Cyclical parameter effects (cycl_only) -----
df_cycl <- sum.all %>% filter(model_name == "cycl_only")
f_vs_b_cycl <- ggplot(seas_vals_converged,
						aes(x = pretty_group,y = effSize)) +
	geom_violin(draw_quantiles = c(.5)) +
	geom_jitter(aes(#shape = as.factor(significant),
		color = pretty_group), size = 5, height = 0, width=.1, alpha = .1,
		shape = 16, show.legend = F) + ylim(c(0,.9)) +
	xlab("Domain")+
	ylab("Seasonal amplitude") +
	facet_grid(rows = vars(fcast_type),
						 #cols = vars(group),
						 #drop = F,
						 scales = "free", space = "free") + #,strip.position="bottom",nrow=2) +
	theme_minimal(base_size = 18) + theme(#axis.ticks.x=element_blank(),
		axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
			angle = 320, vjust=1, hjust = -0.05))  +
	stat_compare_means(
										 method = "t.test", size=5, label.y.npc = .75, label.x.npc = .4) + scale_y_log10()
f_vs_b_cycl




# Do functional groups show stronger effects than taxa?
ggplot(df_refit_fg_tax %>% filter(!beta %in% c("sin","cos")), aes(x = fcast_type,y = effSize,
																																	color = pretty_group)) +
	geom_violin(draw_quantiles = c(.5)) +
	geom_jitter(
							# shape=21,
							# color="black",
							width=.3, height = 0, size=4, alpha = .2) +
	labs(title = "Absolute effect size") +
	xlab("Rank")+
	ylab(NULL) +
	facet_grid(rows = vars(beta), cols = vars(pretty_group), drop = T,
						 scales = "free", space = "free_x") + #,strip.position="bottom",nrow=2) +
	theme_bw() + theme(#axis.ticks.x=element_blank(),
		text = element_text(size = 16),
		axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
			angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=22,face="bold"),
		strip.text.y = element_text(size=12,face="bold")
	)
