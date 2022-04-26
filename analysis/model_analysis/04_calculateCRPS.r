# Calculate and visualize CRPS for diversity forecasts
pacman::p_load(scoringRules, reshape2, parallel, lubridate, data.table) 
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

# Combine forecast types
div_hindcast_data <- readRDS("./data/summary/hindcast_div.rds") %>% 
	mutate(fcast_type = "Diversity", rank = paste0("diversity_", group),
				 taxon = rank)
fg_hindcast_data <- readRDS("./data/summary/hindcast_fg.rds") %>% 
	mutate(fcast_type = "Functional group")
tax_hindcast_data <- readRDS("./data/summary/hindcast_tax_test.rds") %>% 
	mutate(fcast_type = "Taxonomic") # replace
hindcast_data <- rbindlist(list(div_hindcast_data, 
																fg_hindcast_data, 
																tax_hindcast_data), fill = T)

hindcast_data$only_rank <- sapply(str_split(hindcast_data$rank, "_",  n = 2), `[`, 1)
hindcast_data$only_rank <- ordered(hindcast_data$only_rank, levels = c("genus",
																								 "family",
																								 "order", 
																								 "class",
																								 "phylum", "functional", "diversity"))
hindcast_data$rank <- ordered(hindcast_data$rank, levels = c("genus_bac","genus_fun",
																			 "family_bac","family_fun",
																			 "order_bac", "order_fun",
																			 "class_bac", "class_fun",
																			 "phylum_bac","phylum_fun",
																			 "functional_group", "diversity_16S", "diversity_ITS"))


# Add prettier data values 
hindcast_data$newsite <- ifelse(hindcast_data$new_site, "New site", "Observed site")
hindcast_data$pretty_group <- ifelse(hindcast_data$group=="16S", "Bacteria", "Fungi")
hindcast_data$pretty_name <- recode(hindcast_data$rank, !!!pretty_rank_names)
hindcast_data$pretty_name <- ordered(hindcast_data$pretty_name, levels = c("Genus",
																																			 "Family",
																																			 "Order", 
																																			 "Class",
																																			 "Phylum", "Functional group", "Diversity"))

saveRDS(hindcast_data, "./data/summary/all_hindcasts.rds")



# Subset to hindcast observations
hindcast_data <- hindcast_data %>% filter(!is.na(hindcast_data$truth) & dates > "2016-12-31")
# Score using continuous ranked probability score with estimated mean and SD
scored_hindcasts <- hindcast_data %>% 
	mutate(truth = as.numeric(truth)) %>% 
	mutate(crps = crps_norm(truth, mean, sd))


# Compare means by fcast type/group/newsite
scored_hindcasts_mean <- scored_hindcasts %>% 
	group_by(fcast_type, pretty_group, model_name, newsite) %>% 
	dplyr::summarize(crps_mean = mean(crps, na.rm=T))
# Compare means by plot
scored_hindcasts_plot <- scored_hindcasts %>% 
	group_by(fcast_type, pretty_group, siteID, plotID, model_name, pretty_name, newsite) %>% 
	summarize(crps_mean = mean(crps, na.rm=T))
# Compare means by siteID
scored_hindcasts_site <- scored_hindcasts %>% 
	group_by(fcast_type, pretty_group, siteID, model_name, pretty_name, newsite) %>% 
	summarize(crps_mean = mean(crps, na.rm=T))
# Compare means by only-rank
scored_hindcasts_rank <- scored_hindcasts %>% 
	group_by(fcast_type, pretty_group, model_name, newsite, pretty_name) %>% 
	dplyr::summarize(crps_mean = mean(crps, na.rm=T))
# Compare means by specific group/rank
scored_hindcasts_taxon <- scored_hindcasts %>% 
	group_by(fcast_type, pretty_group, model_name, newsite, pretty_name, taxon) %>% 
	dplyr::summarize(crps_mean = mean(crps, na.rm=T))


# Skill score from Scavia 2021
skill_score_rank <- scored_hindcasts_rank %>% 
	pivot_wider(id_cols = c("fcast_type","pretty_group","model_name","pretty_name"), 
							values_from = "crps_mean", names_from = "newsite") %>% 
	mutate(skill_score = (1 - (`New site`/`Observed site`)))
skill_score_taxon <- scored_hindcasts_taxon %>% 
	pivot_wider(id_cols = c("fcast_type","pretty_group","model_name","pretty_name","taxon"), 
							values_from = "crps_mean", names_from = "newsite") %>% 
	mutate(skill_score = (1 - (`New site`/`Observed site`)))


to_save <- list("scored_hindcasts" = scored_hindcasts, 
								"scored_hindcasts_mean" = scored_hindcasts_mean, 
								"scored_hindcasts_plot" = scored_hindcasts_plot, 
								"scored_hindcasts_site" = scored_hindcasts_site, 
								"scored_hindcasts_taxon" = scored_hindcasts_taxon,
								"skill_score_rank" = skill_score_rank,
								"skill_score_taxon" = skill_score_taxon)
# Save output before visualizing
saveRDS(to_save, "./data/summary/CRPS_hindcasts.rds")



# View raw CRPS scores by group and model_name: 
# Fungal forecasts are better at both observed and new sites
ggplot(scored_hindcasts) + 
	coord_trans(y = "log10") +
	geom_jitter(aes(x = model_name, y = crps, color = group), #width=.2, #height = 0, 
							alpha = .1, size=4, 
							position=position_jitterdodge(dodge.width = 1)) + 
	geom_violin(aes(x = model_name, y = crps, color = group), 
							draw_quantiles = c(0.5), show.legend=F) + 
	facet_wrap(~newsite) +  
	ylab("Continuous ranked probability score (CRPS)") + xlab(NULL) + 
	theme_minimal(base_size=18)

# Compare mean skill scores
# skill_score$skill_score <- ifelse(is.na(skill_score$skill_score), 0, skill_score$skill_score)
ggplot(skill_score %>% filter(model_name == "all_covariates"), 
			 aes(x = only_rank, y = skill_score, 
			 		color = pretty_group, shape = fcast_type)) + 
#	geom_violin(aes(fill = uncert), draw_quantiles = c(0.5), show.legend=F) + 
	#facet_grid(~model_name, drop = T) +  
	geom_jitter(aes(x = only_rank, y = skill_score), width=.1, 
							height = 0, alpha = .8, size=4) + 
	ylab("Skill score (% change in CRPS)") + xlab(NULL) + 
	theme_minimal(base_size=18) + ggtitle("Change in predictability at new sites") + 
	theme(text = element_text(size = 20),
		axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=22)) 


# Compare plot-mean CRPS across ranks
ggplot(scored_hindcasts_plot  %>% filter(model_name == "all_covariates")) + 
	geom_point(aes(x = only_rank, y = crps_mean, color = pretty_group), 
						 alpha = .1, size=4,
						 position=position_jitterdodge(dodge.width = 1)) + 
	geom_violin(aes(x = only_rank, y = crps_mean, color = pretty_group), 
							draw_quantiles = c(0.5), show.legend=F) + 
	coord_trans(y = "log10") +
	facet_grid(cols=vars(fcast_type), drop = T, scales = "free", space = "free", as.table = T) +  
	ylab("Continuous ranked probability score (CRPS)") + xlab(NULL) + 
	theme_minimal(base_size=20) + ggtitle("Predictability at observed sites") + 
	theme(text = element_text(size = 22),
				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=24), legend.position = c(.9,1.1)) + 
	guides(color=guide_legend(title=NULL))



# Do T-tests

# F vs. B difference between observed site forecasts
# plot means
oldsites_allcov <- scored_hindcasts_plot %>% filter(newsite=="Observed site" & model_name == "all_covariates")
oldsites_bac <- oldsites_allcov[oldsites_allcov$group=="16S",]
oldsites_fun <- oldsites_allcov[oldsites_allcov$group=="ITS",]
t.test(oldsites_bac$crps_mean, oldsites_fun$crps_mean)
# all values
oldsites_allcov <- scored_hindcasts %>% filter(newsite=="Observed site" & model_name == "all_covariates")
oldsites_bac <- oldsites_allcov[oldsites_allcov$group=="16S",]
oldsites_fun <- oldsites_allcov[oldsites_allcov$group=="ITS",]
t.test(oldsites_bac$crps, oldsites_fun$crps)



# F vs. B difference between new site forecasts
# plot means
newsites_allcov_plot <- scored_hindcasts_plot %>% filter(newsite=="New site" & model_name == "all_covariates")
newsites_bac <- newsites_allcov_plot[newsites_allcov_plot$group=="16S",]
newsites_fun <- newsites_allcov_plot[newsites_allcov_plot$group=="ITS",]
t.test(newsites_bac$crps_mean, newsites_fun$crps_mean)
# site means
newsites_allcov_site <- scored_hindcasts_site %>% filter(newsite=="New site" & model_name == "all_covariates")
newsites_bac <- newsites_allcov_site[newsites_allcov_site$group=="16S",]
newsites_fun <- newsites_allcov_site[newsites_allcov_site$group=="ITS",]
t.test(newsites_bac$crps_mean, newsites_fun$crps_mean)
# all values
newsites_allcov <- scored_hindcasts %>% filter(newsite=="New site" & model_name == "all_covariates")
newsites_bac <- newsites_allcov %>% filter(group=="16S")
newsites_fun <- newsites_allcov %>% filter(group=="ITS")
t.test(newsites_bac$crps, newsites_fun$crps)



ggplot(skill_score, aes(x = pretty_group, y = skill_score, color = pretty_group)) + 
	ggtitle("Skill score: forecasting to new sites") +
	facet_wrap(~model_name) +  
	geom_point(size=4, show.legend = F) + 
	#ylim(c( -.25, 0.)) +
	ylab("Skill score (decrease in mean CRPS)") + xlab(NULL) + 
	theme_minimal(base_size=18)






its_bart_fcast <- hindcast_data_oldsites[hindcast_data_oldsites$plotID %in% c("CPER_001","CPER_002","CPER_004"),] %>% 
	filter(scenario == "full_uncertainty_ITS")

plot_summaries <- summaries$plot_est %>% filter(time_period == "calibration" & scenario == "full_uncertainty_ITS")
its_bart_est <- plot_summaries %>% filter(plotID %in% c("CPER_001","CPER_002","CPER_004"))  %>% arrange(siteID,plotID,dateID) #,"BART_002","BART_004"),]



ggplot() +
	facet_grid(rows=vars(plotID), 
						 cols=vars(model_name), 
						 drop=T, scales="free", space="free") +
	geom_line(data = its_bart_fcast, aes(x = dates, y = `50%`), show.legend = F, na.rm = T) +
	geom_ribbon(data = its_bart_fcast, aes(x = dates, ymin = `2.5%`, ymax = `97.5%`), 
							alpha=0.7,  na.rm = T, show.legend = F)  + 
	geom_point(data = its_bart_fcast, aes(x = dates, y = as.numeric(truth)))  +
	scale_x_date(date_breaks = "6 month", 
							 labels=date_format("%b-%Y"),
							 limits = as.Date(c('2014-01-01','2021-01-01'))) + #ylim(c(-1.1,1.1)) + 
	geom_line(data = its_bart_est, aes(x = dates, y = `50%`), show.legend = F, color = 2, na.rm = T) +
	geom_ribbon(data = its_bart_est, aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill = 2,
							alpha=0.7,  na.rm = T, show.legend = F)  
