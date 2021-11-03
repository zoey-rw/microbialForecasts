# Visualize effect size estimates (beta covariates) from all model

library(stringr)
library(forestplot)
library(ggplot2)
library(gridExtra)
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/helperFunctions.r")


sum.all <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/all_fcast_effects.rds")


library(ggpubr)
a <- ggplot(data=sum.all[sum.all$fcast_period=="calibration",],
			 aes(x = beta,y = effSize)) +
		geom_jitter(aes(#shape = as.factor(significant), 
										color = beta), size = 5, height = 0, width=.1, alpha = .5,
								shape = 16, show.legend = F) + ylim(c(0,.9)) +
	labs(col = "Parameter", title = "Absolute effect size (calibration)") + 
	xlab("Parameter")+ 
	ylab(NULL) +
	facet_grid(rows = vars(fcast_type), cols = vars(group), drop = T,
						 scales = "free_x", space = "free_x") + #,strip.position="bottom",nrow=2) +
	theme_bw() + theme(#axis.ticks.x=element_blank(),
		text = element_text(size = 16),
		axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
			angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=22,face="bold")
		#strip.text.y = element_text(size=24,hjust=0,vjust = 1,angle=180,face="bold")
	) 

b <- ggplot(data=sum.all[sum.all$fcast_period=="refit",],
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

ggarrange(a,b)






c <- ggplot(data=sum.all[sum.all$fcast_period=="calibration",],
						aes(x = beta,y = Mean)) +
	geom_jitter(aes(#shape = as.factor(significant), 
		color = beta), size = 5, height = 0, width=.1, alpha = .5,
		shape = 16, show.legend = F) + ylim(c(-.9,.9)) +
	labs(col = "Parameter", title = "Effect size (calibration)") + 
	xlab("Taxon")+ 
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

d <- ggplot(data=sum.all[sum.all$fcast_period=="refit",],
						aes(x = beta,y = Mean)) +
	geom_jitter(aes(#shape = as.factor(significant), 
		color = beta), size = 5, height = 0, width=.1, alpha = .5,
		shape = 16, show.legend = F) + ylim(c(-.9,.9)) +
	labs(col = "Parameter", title = "Effect size (refit)") + 
	xlab("Taxon")+ 
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

ggarrange(c,d)

ggplot(data=sum.all,
			 aes(x = beta,y = Mean)) +
	geom_jitter(aes(color = fcast_period), size = 5, height = 0, width=.1, alpha = .5,
		shape = 16, show.legend = F) +
	labs(col = "Parameter", title = "Effect size (refit)") + 
	xlab("Taxon")+ 
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












# 
# 
# ## Violin plots of parameter estimates
# div_ITS_samps <- sum.div.all$samples$full_uncertainty_ITS
# div_ITS_samps_long <- do.call(rbind.data.frame, 
# 															div_ITS_samps[,grep("beta|rho", colnames(div_ITS_samps[[1]]))]) %>% 
# 	rownames_to_column("rowname") %>% pivot_longer(2:8) %>% 
# 	mutate(beta_num = gsub("beta\\[|\\]", "", name)) %>% 
# 	mutate(beta = recode(beta_num,
# 											 "1" = "Temperature",
# 											 "2" = "Moisture",
# 											 "3" = "pH",
# 											 "4" = "pC",
# 											 "5" = "Plant species richness",
# 											 "6" = "Ectomycorrhizal trees",
# 											 "rho" = "Autocorrelation",
# 											 .missing = "Autocorrelation"),
# 				 fcast_period = "calibration")
# 
# div_ITS_samps_refit <- sum.div.all_refit$samples$full_uncertainty_ITS
# div_ITS_samps_refit_long <- do.call(rbind.data.frame, div_ITS_samps_refit[,grep("beta|rho", colnames(div_ITS_samps_refit[[1]]))]) %>% 
# 	rownames_to_column("rowname") %>% pivot_longer(2:8) %>% 
# 	mutate(beta_num = gsub("beta\\[|\\]", "", name)) %>% 
# 	mutate(beta = recode(beta_num,
# 											 "1" = "Temperature",
# 											 "2" = "Moisture",
# 											 "3" = "pH",
# 											 "4" = "pC",
# 											 "5" = "Plant species richness",
# 											 "6" = "Ectomycorrhizal trees",
# 											 "rho" = "Autocorrelation",
# 											 .missing = "Autocorrelation"),
# 				 fcast_period = "refit")
# 
# ggplot(data=samps_long,
# 			 aes(x = reorder(beta, value),y = value)) +
# 	geom_violin(aes(fill = beta), trim=FALSE, show.legend = F) + 
# 	ylab("Effect size") + 
# 	xlab("Parameter")+ ggtitle("Drivers of fungal evenness (calibration: 2013-2016)") + 
# 	theme_minimal(base_size = 16) + geom_hline(aes(yintercept=0), linetype=2) + 
# 	theme(
# 		axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
# 			angle = 320, vjust=1, hjust = -0.05),
# 		axis.title=element_text(size=18,face="bold")
# 	) 


# Refit parameter estimates, row = fcast type, col = parameter
df_refit <- df[df$fcast_period=="refit",]

b_vs_f_fcast_type_plot <- ggplot(data=df_refit,
			 aes(x = pretty_group, 
			 		color = pretty_group, y = effSize)) +
	geom_jitter(aes(#shape = as.factor(significant), 
		color = pretty_group), size = 5, height = 0, width=.2, alpha = .5,
		shape = 16, show.legend = F) +
	geom_boxplot(color = 1) +
	labs(col = "Parameter", title = "Parameter effects across forecast types") + 
	xlab("Kingdom")+ 
	ylab("Absolute effect size") +
	facet_grid(rows = vars(fcast_type),
		cols = vars(beta), 
		drop = T,
						 scales = "free") + #,strip.position="bottom",nrow=2) +
	theme_bw() + theme(#axis.ticks.x=element_blank(),
		text = element_text(size = 22),
		axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
			angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=22,face="bold")
		#strip.text.y = element_text(size=24,hjust=0,vjust = 1,angle=180,face="bold")
	) 



sum_fg_tax <- df[df$fcast_period=="refit" & df$fcast_type != "Diversity",]
sum_fg_tax$pretty_group <- ifelse(sum_fg_tax$group=="16S", "Bacteria", "Fungi")

# Compare bac and fun tax & functional ranks for each beta
ggplot(data=sum_fg_tax[grepl("beta", sum_fg_tax$rowname),],
			 aes(x = as.numeric(only_rank),y = effSize, color = pretty_group)) +
	geom_jitter(aes(#shape = as.factor(significant), 
		x = only_rank, color = pretty_group), size = 5, height = 0, width=.1, alpha = .5,
		shape = 16, show.legend = F) +
	labs(col = "", title = "Absolute effect size") + 
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
	) + geom_smooth(show.legend = F)



df_refit

df <- sum.all[sum.all$fcast_period=="refit" & sum.all$fcast_type == "Taxonomic",]

tukey_list <- list()
for(b in 1:6) {
	x = "pretty_group"
	y = "effSize"
	beta_name <-  as.character(unique(sum.all[sum.all$beta_num == b,]$beta))
 df <- sum.all[sum.all$fcast_period=="refit" & sum.all$fcast_type == "Taxonomic" & sum.all$beta_num == b,]
new.df <- cbind.data.frame(x = df$pretty_group, y = df$effSize)
abs_max <- max(new.df[,"y"], na.rm = T)
maxs <- new.df %>%
	group_by(x) %>%
	summarise(tot=max(y, na.rm=T)+ 0.3 * abs_max)
Tukey_test <- aov(y ~ x, data=new.df) %>%
	agricolae::HSD.test("x", group=TRUE) %>%
	.$groups %>%
	as_tibble(rownames="x") %>%
	rename("Letters_Tukey"="groups") %>% 
	dplyr::select(-y) %>%
	left_join(maxs, by="x") %>% 
	rename("pretty_group"="x")
Tukey_test$beta <- beta_name
tukey_list[[b]] <- Tukey_test
}

tukey_list_tax <- data.table::rbindlist(tukey_list)
tukey_list_tax$fcast_type <- "Taxonomic"

tukey_list <- list()
for(b in 1:6) {
	x = "pretty_group"
	y = "effSize"
	beta_name <-  as.character(unique(sum.all[sum.all$beta_num == b,]$beta))
	df <- sum.all[sum.all$fcast_period=="refit" & sum.all$fcast_type == "Functional group" & sum.all$beta_num == b,]
	new.df <- cbind.data.frame(x = df$pretty_group, y = df$effSize)
	abs_max <- max(new.df[,"y"], na.rm = T)
	maxs <- new.df %>%
		group_by(x) %>%
		summarise(tot=max(y, na.rm=T)+ 0.3 * abs_max)
	Tukey_test <- aov(y ~ x, data=new.df) %>%
		agricolae::HSD.test("x", group=TRUE) %>%
		.$groups %>%
		as_tibble(rownames="x") %>%
		rename("Letters_Tukey"="groups") %>% 
		dplyr::select(-y) %>%
		left_join(maxs, by="x") %>% 
		rename("pretty_group"="x")
	Tukey_test$beta <- beta_name
	tukey_list[[b]] <- Tukey_test
}


tukey_list_fg <- data.table::rbindlist(tukey_list)
tukey_list_fg$fcast_type <- "Functional group"

tukey <- rbind(tukey_list_fg, tukey_list_tax)

b_vs_f_fcast_type_plot + geom_text(data = tukey, aes(x = pretty_group, y = tot+.1, 
																										 label = Letters_Tukey), show.legend = F, color = 1) 
																	 	




# View actual taxa and effects

df <- sum.all[sum.all$fcast_period=="refit" & sum.all$fcast_type == "Taxonomic",]
ggplot(data=df,
			 aes(x = taxon_num,y = Mean)) +
	geom_jitter(aes(color = beta, shape = as.factor(significant)), size = 5, height = 0, width=.1, alpha = .8) +
	labs(col = "Parameter", title = "Effect size (refit)") + 
	xlab("Taxon")+ 
	ylab(NULL) +
	geom_hline(yintercept = 0, linetype=2) +
	facet_grid(rows = vars(only_rank), cols = vars(pretty_group), drop = T,
						 scales = "free", space = "free_x") + #,strip.position="bottom",nrow=2) +
	theme_bw() + theme(#axis.ticks.x=element_blank(),
		text = element_text(size = 16),
		axis.text.x=element_text(NULL),
#angle = 45, hjust = 1, vjust = 1),
			#angle = 320, vjust=1, hjust = -0.05
		#),
		axis.title=element_text(size=22,face="bold")
		#strip.text.y = element_text(size=24,hjust=0,vjust = 1,angle=180,face="bold")
	) + scale_shape_manual(values = c(21, 16), name = NULL, 
												 labels = c("Not significant","Significant")) + 
	ggrepel::geom_text_repel(data =df[df$beta_num==1, ], aes(label=taxon, x = taxon_num),
													 nudge_y = -0.6, #direction="y",
													 min.segment.length = Inf) + 
	geom_hline(yintercept = 0, linetype=2) 
