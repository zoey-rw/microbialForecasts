library(ggstatsplot)

source("source.R")
library(ggpubr)
library(scoringRules)

clr_summaries <- readRDS(here("data/summary/clr_regression_summaries.rds"))
beta_summaries <- readRDS(here("data/summary/logit_beta_regression_summaries.rds"))

beta_mean_abun = beta_summaries$plot_est %>% filter(!is.na(truth)) %>% 
	group_by(model_id) %>% summarize(beta_mean_abun = mean(truth))

clr_mean_abun = clr_summaries$plot_est %>% filter(!is.na(truth)) %>% 
	group_by(model_id) %>% summarize(clr_mean_abun = mean(truth))

clr_keep_list = readRDS(here("data/summary/clr_converged_taxa_list.rds"))
clr_predictor_effects <- readRDS(here("data/summary/clr_predictor_effects.rds")) %>% 
	filter(time_period == "2015-11_2018-01")  %>% 
	mutate(model_type="clr") %>% 
	select(beta, taxon, Mean, rank, model_name, model_id, model_type) %>% 
	filter(model_name=="env_cycl")
clr_predictor_effects$converged = ifelse(clr_predictor_effects$model_id %in% clr_keep_list, 1, 0)

beta_keep_list = readRDS(here("data/summary/converged_taxa_list.rds"))
beta_predictor_effects <- readRDS(here("data/summary/predictor_effects.rds"))  %>% 
	filter(time_period == "2015-11_2018-01") %>% 
	mutate(model_type="beta_reg") %>% 
	select(beta, taxon, Mean, rank, model_name, model_id, model_type) %>% 
	filter(model_name=="env_cycl")
beta_predictor_effects$converged = ifelse(beta_predictor_effects$model_id %in% beta_keep_list, 1, 0)


dirichlet_predictor_effects <- readRDS(here("data/summary/dirichlet_predictor_effects.rds")) %>% 
	mutate(model_type="dirichlet") %>% mutate(model_name = ifelse(model_name== "all_covariates","env_cycl",model_name),
																						time_per = gsub("-","",time_period),
																						time_per = ifelse(time_per== "201511_201801","20151101_20180101",time_per),
																					
																						model_id = paste0(model_name,"_",taxon,"_",time_per)) %>% 
	select(beta, taxon, Mean, rank, model_name, model_id, model_type)
#beta_predictor_effects$converged = ifelse(beta_predictor_effects$model_id %in% beta_keep_list, 1, 0)


# all_models <- merge(clr_predictor_effects, beta_predictor_effects, by=c("beta", "taxon", "Mean", "rank", "model_name", "model_id"))

all_models_long <- data.table::rbindlist(list(clr_predictor_effects, beta_predictor_effects, dirichlet_predictor_effects), 
																				 use.names = T, fill = T)  #%>% 
#	filter(rank %in% c("phylum_fun","genus_fun","phylum_bac","genus_bac"))
# 
# filter(rank %in% c("phylum_fun","genus_fun","phylum_bac","class_bac","class_fun","order_bac", "order_fun","family_bac","family_fun"))

all_models <- all_models_long %>% 
	#select(beta, taxon, Mean, rank, model_name, group, model_id, model_type) %>% 
	group_by(beta) %>% 
	pivot_wider(names_from="model_type",
							values_from = c("Mean", "converged")) 
	#pivot_wider(names_from="beta",values_from = "Mean")

all_models <- merge(all_models, beta_mean_abun, all.x=T)
all_models <- merge(all_models, clr_mean_abun, all.x=T)
all_models$mean_abun = ifelse(is.na(all_models$beta_mean_abun), NA, all_models$beta_mean_abun)


all_models$both_converged = ifelse(all_models$converged_beta_reg==1 & 
																	 	all_models$converged_clr==1, 1, 0)
for_cor = all_models  %>% 
	filter(both_converged==1) %>% 
	pivot_wider(names_from = c("beta"), values_from = c(Mean_clr, Mean_beta_reg)) 
	
	
	ggcorrmat(
		data = for_cor,
		cor.vars = c(Mean_clr_cos:Mean_beta_reg_LAI)
	)



fig_param_estimates = ggplot(all_models %>% 	filter(both_converged==1),
														 aes(x = Mean_clr, y = Mean_beta_reg)) + 
	geom_point(data=all_models, aes(x = Mean_clr, y = Mean_beta_reg, alpha = both_converged), show.legend = F) + 
	#geom_point(aes(color=converged)) + 
	facet_wrap(~beta, scales="free") + 
	geom_smooth(method="lm", se=F, show.legend = F) + 
	stat_cor(method="spearman", 
					 aes(label = ..r.label..)) + 
	theme_classic(base_size = 14) + 
	xlab("CLR parameter estimate") + 
	ylab("Beta regression parameter estimate") + theme(strip.text = element_text( margin = margin( b = 0, t = 0) ) )
fig_param_estimates

convergence=all_models %>% filter(!is.na(converged_clr)) %>%  
	select(model_id, mean_abun, converged_clr, converged_beta_reg) %>% distinct()

convergence_sum_betareg=all_models %>% 
	filter(!is.na(converged_clr)) %>%  
	select(model_id, mean_abun, converged_beta_reg) %>% 
	distinct() %>% group_by(converged_beta_reg) %>% 
	tally() %>% 
	mutate(model_type="Beta regression (rel. abundance)",
				 converged = n/187)

convergence_sum_clr=all_models %>% 
	filter(!is.na(converged_clr)) %>%  
	select(model_id, mean_abun, converged_clr) %>% 
	distinct() %>% 
	group_by(converged_clr) %>% 
	tally() %>% 
	mutate(model_type="CLR transformation", converged = n/187)

overall_convergence_sum = rbindlist(list(convergence_sum_betareg, 
																				 convergence_sum_clr), fill = T) %>% 
	filter(converged_beta_reg == 1 | converged_clr == 1)

fig_overall_convergence = ggplot(overall_convergence_sum, 
			 aes(y = converged, x = model_type, 
			 		color = model_type)) + 
	geom_point(show.legend = F) + 
	ylab("Percent of full-covariate models converged") + 
	theme_classic(base_size = 14)  + ylim(c(0,1)) + xlab("Model type")





convergence$rare = ifelse(convergence$mean_abun < .03, 1, 0)
convergence$binned_abun = round(convergence$mean_abun, 2)

	
	convergence_binned_sum_clr1=convergence %>%  
	select(model_id, binned_abun, converged_clr) %>% 
	distinct()   %>% 
	group_by(binned_abun, converged_clr) %>% 
	add_tally(name="group count") %>% 
	ungroup %>% 
	group_by(binned_abun) %>% 
	add_tally(name="total") 

to_add_clr = convergence_binned_sum_clr1[convergence_binned_sum_clr1$converged_clr==0 & convergence_binned_sum_clr1$total==1,] %>% 
	mutate(pct_converged = 0) %>% 
	select(pct_converged, binned_abun, total) %>% distinct() %>% 
	mutate(model_type="CLR transformation")

	
	convergence_binned_sum_clr2 = convergence_binned_sum_clr1 %>% filter(converged_clr == 1) %>% 
	mutate(pct_converged = `group count`/total) %>% ungroup %>% 
	select(pct_converged, binned_abun, total) %>% distinct() %>% 
	mutate(model_type="CLR transformation")

convergence_binned_sum_beta_reg1=convergence %>%  
	select(model_id, binned_abun, converged_beta_reg) %>% 
	distinct()   %>% 
	group_by(binned_abun, converged_beta_reg) %>% 
	add_tally(name="group count") %>% ungroup %>% 
	group_by(binned_abun) %>% 
	add_tally(name="total")


to_add_beta_reg = convergence_binned_sum_beta_reg1[convergence_binned_sum_beta_reg1$converged_beta_reg==0 & convergence_binned_sum_beta_reg1$total==1,] %>% 
	mutate(pct_converged = 0) %>% 
	select(pct_converged, binned_abun, total) %>% distinct() %>% 
	mutate(model_type="Beta regression (rel. abundance)")
	
convergence_binned_sum_beta_reg2 = convergence_binned_sum_beta_reg1 %>% 
	filter(converged_beta_reg == 1) %>% 
	mutate(pct_converged = `group count`/total) %>% ungroup %>% 
	select(pct_converged, binned_abun, total) %>% distinct() %>% 
	mutate(model_type="Beta regression (rel. abundance)")

percent_convergence = rbind(convergence_binned_sum_clr, convergence_binned_sum_beta_reg, to_add_beta_reg, to_add_clr)
fig_rarity_convergence = ggplot(percent_convergence %>% filter(total > 1), 
			 aes(y = pct_converged, x = binned_abun, color=model_type)) + 
	#geom_point(alpha=.5) +
	geom_jitter(height = 0, width=0, aes(size=total),
							alpha=.5) + 
	geom_smooth(method = "glm", se=F, 
							method.args = list(family = "binomial")) + 
	theme_classic(base_size = 14) + theme(legend.position="top", #, #c(.8,.8),
																			legend.direction = "vertical") + 
	#scale_x_sqrt() + 
	xlab("Mean relative abundance of microbial group") + 
	ylab("Percent of full-covariate models converged") + 
	labs(color='Model type', size="Count") 
fig_rarity_convergence

library(ggpubr)

BC = ggarrange(fig_rarity_convergence, fig_overall_convergence, widths = c(3,1))
ggarrange(fig_param_estimates, BC, nrow = 2)


ggarrange(fig_param_estimates, fig_rarity_convergence, nrow = 1)



ggplot(convergence, aes(y = converged_beta_reg, x = mean_abun)) + 
	geom_jitter(height = 0.01, width=.011)


betareg_glm = glm(formula = converged_beta_reg ~ mean_abun, family = binomial, data = convergence)
clr_glm = glm(formula = converged_clr ~ mean_abun, family = binomial, data = convergence)



#:::::::::::::::::::::::::::::::::::::::::::::::::
ggboxplot(convergence, x = "converged_clr", y = "mean_abun",
							 add = "jitter") + stat_compare_means()
ggboxplot(convergence, x = "converged_beta_reg", y = "mean_abun",
					add = "jitter") + stat_compare_means()



ggplot(all_models %>% filter(converged==1), aes(x = clr, y = beta_reg)) + geom_point(aes(color=converged)) + facet_wrap(~beta) + geom_smooth(method="lm", se=F)


correlation_data <- all_models %>%
	
	group_by(beta) %>%
	
	summarise(beta_clr_correlation = cor(Mean_beta_reg, Mean_clr, use="complete.obs"),
						beta_dirichlet_correlation = cor(Mean_beta_reg, Mean_dirichlet, use="complete.obs"),
						clr_dirichlet_correlation = cor(Mean_clr, Mean_dirichlet, use="complete.obs"))
library(corrplot)
for_corrplot = as.matrix(correlation_data$beta_clr_correlation)
rownames(for_corrplot) <- correlation_data$beta
corrplot(for_corrplot, cl.pos="n")
library(RColorBrewer)
scalebluered <- colorRampPalette(brewer.pal(8, "RdBu"))(8)
colorlegend(xlim=c(-1,0), ylim=c(7,8), 
						scalebluered, 
						c(seq(-1,1,.25)), align="l", vertical=TRUE, addlabels=TRUE)

