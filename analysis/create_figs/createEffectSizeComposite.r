# Visualize effect size estimates (beta covariates) from all model
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
pacman::p_load(stringr, forestplot, gridExtra, ggpubr) 

sum.all <- readRDS("./data/summary/all_fcast_effects.rds")
df_refit <- sum.all %>% filter(time_period == "refit" & model_name == "all_covariates")
df_refit_fg_tax <- df_refit %>% filter(fcast_type != "Diversity")

####### Absolute effect sizes per rank ------
a <- ggplot(df_refit,
						aes(x = beta,y = effSize)) +
	geom_jitter(aes(#shape = as.factor(significant), 
		color = beta), size = 5, height = 0, width=.1, alpha = .5,
		shape = 16, show.legend = F) + ylim(c(0,.9)) +
	labs(col = "Parameter", title = "Absolute effect size (refit)") + 
	xlab("Parameter")+ 
	ylab(NULL) +
	facet_grid(rows = vars(fcast_type), cols = vars(group), drop = F,
						 scales = "free_x", space = "free_x") + #,strip.position="bottom",nrow=2) +
	theme_bw() + theme(#axis.ticks.x=element_blank(),
		text = element_text(size = 16),
		axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
			angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=22,face="bold")
		#strip.text.y = element_text(size=24,hjust=0,vjust = 1,angle=180,face="bold")
	) #+ scale_shape_manual(values = c(21, 16), name = NULL, 
#										 labels = c("Not significant","Significant")) 

###### Relative effect sizes per rank -----
b <- ggplot(df_refit,
						aes(x = beta,y = Mean)) +
	geom_jitter(aes(#shape = as.factor(significant), 
		color = beta), size = 5, height = 0, width=.1, alpha = .5,
		shape = 16, show.legend = F) + ylim(c(-.9,.9)) +
	labs(col = "Parameter", title = "Effect size (refit)") + 
	xlab("Parameter")+ 
	ylab(NULL) +
	geom_hline(yintercept = 0, linetype=2) +
	facet_grid(rows = vars(fcast_type), cols = vars(group), drop = T,
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


###### # Plot with Tukey, fungi vs bacteria effect sizes -----
b_vs_f_fcast_type_plot <- ggplot(data=df_refit %>% filter(!beta %in% c("sin","cos")),
																 aes(x = pretty_group, 
																 		color = pretty_group, y = effSize)) +
	geom_jitter(aes(#shape = as.factor(significant), 
		color = pretty_group), size = 5, height = 0, width=.2, alpha = .5,
		shape = 16, show.legend = F) +
	geom_boxplot(color = 1) +
	labs(col = "Parameter", title = "Parameter effects across forecast types") + 
	xlab("Domain")+ 
	ylab("Absolute effect size") +
	facet_grid(rows = vars(fcast_type),
						 cols = vars(beta), as.table = T,
						 drop = T,
						 scales = "free") + #,strip.position="bottom",nrow=2) +
	theme_bw() + theme(#axis.ticks.x=element_blank(),
		text = element_text(size = 22),
		axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
			angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=22,face="bold")
		#strip.text.y = element_text(size=24,hjust=0,vjust = 1,angle=180,face="bold")
	) 

tukey_list <- list()
beta_names <- c(#"sin", "cos", 
	"Ectomycorrhizal\ntrees", "LAI", "pC", 
	"pH", "Temperature", "Moisture")
for(b in beta_names) {
	x = "pretty_group"
	y = "effSize"
	df <- df_refit[which(df_refit$beta == b),]
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
	df <- df_refit[which(df_refit$beta == b),]
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
tukey_list_fg$fcast_type <- "Functional group"
tukey <- rbind(tukey_list_fg, tukey_list_tax)

b_vs_f_fcast_type_plot <- b_vs_f_fcast_type_plot + geom_text(data = tukey, aes(x = pretty_group, y = tot+.1, 
																										 label = Letters_Tukey), show.legend = F, color = 1) 
######
b_vs_f_fcast_type_plot

###### # Plot with Tukey, effect sizes across ranks ----
ranks_beta_plot <- ggplot(df_refit_fg_tax %>% filter(!beta %in% c("sin","cos"))) +
	geom_jitter(aes(x = only_rank,y = effSize,
									fill = pretty_group), 
							shape=21, color="black", width=.1, height = 0, size=4) + 
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
	) + geom_smooth(aes(x = as.numeric(only_rank), y = effSize), show.legend = F) + scale_fill_manual(values = c("grey30","grey90"))

tukey_list <- list()
beta_names <- c(#"sin", "cos", 
	"Ectomycorrhizal\ntrees", "LAI", "pC", 
	"pH", "Temperature", "Moisture")
for(b in beta_names) {
	x = "only_rank"
	y = "effSize"
	df <- df_refit_fg_tax[which(df_refit_fg_tax$beta == b),]
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
																label = Letters_Tukey), show.legend = F, color = 2, size =6) 
######
ranks_beta_plot


# View actual taxa and effects
df_refit_fg <- df_refit_fg_tax %>% filter(fcast_type == "Functional group")
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
df_cycl <- sum.all %>% filter(time_period == "refit" & model_name == "cycl_only")
f_vs_b_cycl <- ggplot(df_cycl,
						aes(x = beta,y = effSize)) +
	geom_jitter(aes(#shape = as.factor(significant), 
		color = beta), size = 5, height = 0, width=.1, alpha = .5,
		shape = 16, show.legend = F) + ylim(c(0,.9)) +
	labs(col = "Parameter", title = "Absolute effect size (refit)") + 
	xlab("Parameter")+ 
	ylab(NULL) +
	facet_grid(rows = vars(fcast_type), cols = vars(group), drop = F,
						 scales = "free_x", space = "free_x") + #,strip.position="bottom",nrow=2) +
	theme_bw() + theme(#axis.ticks.x=element_blank(),
		text = element_text(size = 16),
		axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
			angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=22,face="bold"))
