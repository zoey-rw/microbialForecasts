# Comparing functional groups vs taxonomy
source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")
#p_load(multcomp)
library(ggallin)
library(emmeans)


in_list <- readRDS(here("data/summary/fcast_horizon_input.rds"))
to_plot <- in_list[[1]]
fcast_horizon_null_site <-  in_list[[3]]
model_id_list <- unique(to_plot$model_id)

weak_converged <- readRDS(here("data/summary/weak_converged_taxa_list.rds"))
strict_converged <- readRDS(here("data/summary/converged_taxa_list.rds"))
converged = weak_converged
#converged = strict_converged
#converged =  scores_list$converged_strict_list

sum.all <- readRDS(here("data", "summary/predictor_effects.rds"))


df_refit = sum.all %>% filter(time_period=="2015-11_2020-01" & model_name=="env_cycl") %>%
	filter(model_id %in% converged)
df_cal = sum.all %>% filter(time_period=="2015-11_2018-01" & model_name=="env_cov"
														) %>%
	filter(model_id %in% converged)



df_cal <- merge(df_cal, fcast_horizon_null_site[,c("model_id","abundance")], all.x=T)

df_cal$adjEffectSize = df_cal$effSize/df_cal$abundance

plant_fg_groups = c("chitin_complex","cellulolytic",#"assim_nitrate_reduction",  
										"oligotroph", "copiotroph", "herbicide_stress","cellobiose_complex","cellulose_complex",
										"chitinolytic",
										"ectomycorrhizal", "lignolytic","saprotroph","lichenized",
										"plant_pathogen", "endophyte",
										"phylum_bac", "genus_bac", "class_bac", 
										"family_bac", "order_bac", "order_fun", "class_fun", "phylum_fun", 
										"family_fun", "genus_fun")
cal_for_stats = df_cal #%>% 
	#filter(significant == 1)
	#filter(rank %in% plant_fg_groups)
#cal_for_stats = df_cal 


#df_cal = df_cal %>% filter(rank_only %in% c("phylum","functional"))

beta_names <- c("Ectomycorrhizal\ntrees", "LAI", "pC",
								"Temperature", "Moisture","pH")
# Test for differences in effect sizes between fcast types
beta_tukey_group = list()
for (group in c("Fungi", "Bacteria")){
	beta_tukey = list()
	for (beta in beta_names) {
		df_group=cal_for_stats %>% filter(pretty_group==!!group)
		out = df_group %>% filter(beta %in% !!beta) %>%
			#group_by(beta) %>%
			#aov(adjEffectSize~fcast_type,.) %>%
			aov(effSize~fcast_type,.) %>%
			emmeans::emmeans(object = ., pairwise ~ "fcast_type", adjust = "tukey") %>% .$emmeans %>%
			multcomp::cld(object = ., Letters = letters) %>% as.data.frame() %>% rename(Letters_Tukey = `.group`) %>%
			#rownames_to_column("beta") %>%
			mutate(pretty_group = !!group, beta = !!beta)
		beta_tukey[[beta]] = out
	}
	beta_tukey_group[[group]] = do.call(rbind, beta_tukey)
}
tukey_cal_beta = do.call(rbind.data.frame, beta_tukey_group)
tukey_cal_beta$tot = tukey_cal_beta$upper.CL + .1

# Subset to only predictors that have a difference between types
diff_letters <- tukey_cal_beta %>% group_by(beta) %>% #any(Letters_Tukey != "a")  %>%
	mutate(eq = replace(Letters_Tukey, n_distinct(Letters_Tukey)==1, '') ) #%>% filter(eq != "")
df_cal_diff = merge(cal_for_stats, diff_letters, all.x = T) %>% filter(eq != "")


tukey_cal_beta2 = cal_for_stats %>% filter(!beta %in% c("rho","sin","cos")) %>%
	group_by(pretty_group, beta) %>%
	summarize(tukey(x = fcast_type, y = effSize)) %>%
	rename(fcast_type = x)

tukey_beta_diff = merge(tukey_cal_beta2,
												diff_letters[,c("fcast_type","pretty_group","beta","eq")], all.x = T) %>% filter(eq != "")

###### # Plot with Tukey, effect sizes across fcast types ----
c <- ggplot(df_cal_diff %>%
							filter(!beta %in% c("sin","cos","rho")),
						aes(x = fcast_type,y = effSize,color = pretty_group)) +
	xlab(NULL)+
	ylab("Absolute effect size") +
	facet_grid(rows = vars(beta),
						 cols = vars(pretty_group),
						 drop = T,
						 scales = "free", space = "free_x") +
	theme_minimal() +
	theme(text = element_text(size = 16),
				axis.text.x=element_text(
					angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=18),
				strip.text.y = element_text(size=16)) +
	geom_violin(draw_quantiles = c(.5), show.legend = F) +
	geom_jitter(width=.2, height = 0, size=4, alpha = .1, show.legend = F)

fcast_type_beta_plot <- c +
	geom_text(data = tukey_beta_diff,
						
						aes(x = fcast_type, y = tot, label = Letters_Tukey), show.legend = F, color = 1, size =6)
######
fcast_type_beta_plot




library(rstatix)
model_stat_pvalue <- cal_for_stats %>% 
	group_by(pretty_group, beta) %>% 
	compare_means(formula=effSize ~ fcast_type, data=., group.by = c("pretty_group","beta")) %>% 
	add_y_position(data=., x = fcast_type, y = effSize, test="wilcoxon", group=c("pretty_group","beta"))

model_stat_pvalue <- cal_for_stats %>% 
	group_by(pretty_group, beta) %>% 
	rstatix::tukey_hsd(effSize ~ fcast_type) %>%
	#filter(p.adj < 0.05) %>% 
	rstatix::add_y_position() #%>% 
	# mutate(y.position = seq(min(y.position), max(y.position),
	# 												length.out = n()))
ggplot(cal_for_stats %>%
			 	filter(!beta %in% c("sin","cos","rho")),
			 aes(x = fcast_type,y = effSize,color = pretty_group)) +
	xlab(NULL)+
	ylab("Absolute effect size") +
	facet_grid(cols = vars(beta),
						 rows = vars(pretty_group),
						 drop = T,
						 scales = "free", space = "free_x") +
	theme_minimal(base_size = 22) +
	theme(text = element_text(size = 20),
				axis.text.x=element_text(
					angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=18),
				strip.text.y = element_text(size=16)) +
	geom_boxplot(show.legend = F, outlier.shape = NA) +
#	geom_violin(draw_quantiles = c(.5), show.legend = F) +
	geom_jitter(width=.1, height = 0, size=3, alpha = .5, show.legend = F) +
	stat_pvalue_manual(model_stat_pvalue, 
										 size=9,
										 bracket.nudge.y = -0.1, 
										 hide.ns = T)












model_stat_pvalue2 <- cal_for_stats %>% 
	filter(model_id %in% strict_converged) %>% 
	group_by(beta) %>% 
	rstatix::tukey_hsd(effSize ~ fcast_type) %>%
	#filter(p.adj < 0.05) %>% 
	rstatix::add_y_position() #%>% 
# mutate(y.position = seq(min(y.position), max(y.position),
# 												length.out = n()))
ggplot(cal_for_stats %>% 
			 	filter(model_id %in% strict_converged)  %>%
			 	filter(rank %in% plant_fg_groups) %>% 
			 #	filter(rank_only %in% c("phylum","functional"))  %>%
			 	filter(!beta %in% c("sin","cos","rho")),
			 aes(x = fcast_type,y = effSize)) +
	xlab(NULL)+
	ylab("Absolute effect size") +
	facet_grid(cols = vars(beta),
						 drop = T,
						 scales = "free", space = "free_x") +
	theme_minimal(base_size = 22) +
	theme(text = element_text(size = 20),
				axis.text.x=element_text(
					angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=18),
				strip.text.y = element_text(size=16)) +
	geom_boxplot(show.legend = F, outlier.shape = NA) +
	#	geom_violin(draw_quantiles = c(.5), show.legend = F) +
	geom_jitter(width=.1, height = 0, size=3, alpha = .5, show.legend = F) +
	stat_pvalue_manual(model_stat_pvalue2, 
										 size=9,
										 bracket.nudge.y = -0.1, 
										 hide.ns = T)


## phenology


seas_in <- readRDS("/projectnb/dietzelab/zrwerbin/microbialForecasts/data/summary/seasonal_amplitude.rds")

df_cal = seas_in[[6]] %>% filter(time_period=="2015-11_2018-01"# & model_name=="env_cycl"
) %>%
	filter(model_id %in% strict_converged)

df_cal$only_rank <- sapply(str_split(df_cal$rank_only, "_",  n = 2), `[`, 1) %>%
	ordered(levels = c("genus","family","order","class","phylum","functional","diversity"))

pheno_beta = df_cal %>% mutate(any_sig = ifelse(significant_sin|significant_cos, 1, 0)) %>% 
	filter(model_id %in% converged) %>% 
	group_by(model_name, pretty_group, fcast_type, only_rank, any_sig) %>% 
	summarise(n = n()) %>%
	mutate(freq = n / sum(n)) %>% 
	filter(any_sig==1) %>% 
	rename(proportion_sig = freq) %>% mutate(beta = "seasonal effect") 

non_pheno_beta = sum.all %>% 
	filter(time_period=="2015-11_2018-01") %>%
	filter(!beta %in% c("sin","cos")) %>% 
	filter(model_id %in% converged) %>% 
	group_by(model_name, pretty_group, fcast_type, only_rank, beta, significant) %>% 
	summarise(n = n()) %>%
	mutate(freq = n / sum(n)) %>% 
	#filter(significant==1) %>% 
	rename(proportion_sig = freq)
non_pheno_beta[non_pheno_beta$beta=="Moisture" & 
							 	non_pheno_beta$model_name=="env_cycl" & 
							 	non_pheno_beta$fcast_type=="Functional" & 
							 	non_pheno_beta$pretty_group=="Fungi",]$significant = 1
non_pheno_beta[non_pheno_beta$beta=="Moisture" & 
							 	non_pheno_beta$model_name=="env_cycl" & 
							 	non_pheno_beta$fcast_type=="Functional" & 
							 	non_pheno_beta$pretty_group=="Fungi",]$proportion_sig = 0
non_pheno_beta[non_pheno_beta$beta=="LAI" & 
							 	non_pheno_beta$model_name=="env_cycl" & 
							 	non_pheno_beta$fcast_type=="Functional" & 
							 	non_pheno_beta$pretty_group=="Fungi",]$significant = 1
non_pheno_beta[non_pheno_beta$beta=="LAI" & 
							 	non_pheno_beta$model_name=="env_cycl" & 
							 	non_pheno_beta$fcast_type=="Functional" & 
							 	non_pheno_beta$pretty_group=="Fungi",]$proportion_sig = 0

to_plot = rbindlist(list(pheno_beta, non_pheno_beta), fill = T)

ggplot(to_plot %>% 
			 	filter(only_rank %in% c("functional") &
			 				 	beta %in% c(#"seasonal effect",
			 				 							"Temperature","LAI","Moisture")),
			 aes(x = only_rank,y = proportion_sig, color=model_name)) +
	xlab(NULL)+
	ylab("Absolute effect size") +
	facet_grid(cols = vars(pretty_group),
						 rows=vars(beta),
						 drop = T,
						 scales = "free", space = "free_x") +
	theme_minimal(base_size = 22) +
	theme(text = element_text(size = 20),
				axis.text.x=element_text(
					angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=18),
				strip.text.y = element_text(size=16)) +
	geom_boxplot(show.legend = F, outlier.shape = NA) +
	geom_jitter(width=.1, height = 0, size=3, alpha = .5)#, show.legend = F)









pheno_df = df_cal %>% mutate(significant = ifelse(significant_sin|significant_cos, 1, 0),
														 beta = "seasonal effect",
														 effSize=amplitude)
effSizes = rbindlist(list(pheno_df,sum.all %>% 
															filter(time_period=="2015-11_2018-01") %>%
															filter(!beta %in% c("sin","cos"))), fill=T) %>% 
	filter(model_id %in% converged) 
ggplot(effSizes %>% 
			 	filter(#model_id %in% strict_converged & 
			 		rank_only %in% c("functional","phylum") &
			 				 		beta %in% c("seasonal effect","Temperature","LAI","Moisture")),
			 aes(x = rank_only,y = effSize, color=model_name)) +
	xlab(NULL)+
	ylab("Absolute effect size") +
	facet_grid(cols = vars(pretty_group),
						 rows=vars(beta),
						 drop = T,
						 scales = "free", space = "free_x") +
	theme_minimal(base_size = 22) +
	theme(text = element_text(size = 20),
				axis.text.x=element_text(
					angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=18),
				strip.text.y = element_text(size=16)) +
	geom_boxplot(show.legend = F, outlier.shape = NA, position = position_dodge(width = 1)) +
	geom_point(aes(group=model_name, shape=as.factor(significant)), #show.legend = F, 
						 size=3,
						 alpha=.5, 
						 #position=position_dodge(width = .9)) + #, 
						 #position=position_jitter(height=0, width=.1))   + 
	 position=position_jitterdodge(jitter.height = 0, jitter.width=0.01, 
	 															dodge.width = .9))   + 
	scale_shape_manual(values = c(21, 16), name = NULL,																																						 labels = c("Not significant","Significant")) 

