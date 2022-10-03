source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
library(ggpubr)
library(scoringRules)

dirichlet_summaries <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/dirichlet_taxa_summaries.rds")
dirichlet_refit <- dirichlet_summaries$plot_est %>% filter(time_period == "2015-11_2020-01")

beta_summaries <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/single_taxon_summaries_201511_202001.rds")
)
beta_refit <- beta_summaries$plot_est %>% filter(time_period == "2015-11_2020-01")
beta_refit_not_na <- beta_refit %>% filter(!is.na(truth))
beta_refit_long <- beta_refit %>% pivot_longer()



dirichlet_refit <- dirichlet_refit %>% rename("2.5%_dirichlet"= "2.5%" ,
																											"25%_dirichlet"="25%"  ,
																											"50%_dirichlet"="50%"  ,
																											"75%_dirichlet"="75%"  ,
																											"97.5%_dirichlet"="97.5%",
																											#"lo_dirichlet"="lo"   ,
																											#"med_dirichlet"="med"  ,
																											#"hi_dirichlet"="hi"   ,
																											"mean_dirichlet"= "Mean" ,
																											"sd_dirichlet"="SD",
																											"truth_dirichlet" = "truth") %>%
	select(-c(date_num, species_num, plot_num, site_num, timepoint, rank)) %>% filter(taxon != "other")


beta_refit <- beta_refit %>% select(-c(date_num, species_num, plot_num, site_num, timepoint, rank)) %>% filter(taxon != "other")

all_hindcasts <- full_join(dirichlet_refit, beta_refit, all=T)

library(dplyr)
all_models <- left_join(dirichlet_refit, beta_refit,
					#by = c("cyl"),
					suffix = c("_dirichlet", "_beta"))
all_models <- inner_join(dirichlet_refit, beta_refit,
												#by = c("cyl"),
												suffix = c("_dirichlet", "_beta"))

all_models$discrepancy = all_models$truth_dirichlet - all_models$truth
all_models[which(all_models$discrepancy > .1),]
head(all_models)

beta_basid <- beta_refit %>% filter(taxon=="basidiomycota" &
															 	group=="ITS" & model_name == "all_covariates")

dirichlet_basid <- dirichlet_refit %>% filter(taxon=="basidiomycota" &
																			group=="ITS" & model_name == "all_covariates")

ggplot() +
	geom_line(data = dirichlet_basid %>% filter(siteID == "HARV"), aes(x = dates, y = `50%_dirichlet`), show.legend = F) +
	geom_ribbon(data = dirichlet_basid %>% filter(siteID == "HARV"), aes(x = dates, ymin = `2.5%_dirichlet`, ymax = `97.5%_dirichlet`),fill="red", alpha=0.6) +
	geom_line(data = beta_basid %>% filter(siteID == "HARV"), aes(x = dates, y = `50%`), show.legend = F) +
	geom_ribbon(data = beta_basid %>% filter(siteID == "HARV"), aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill="blue", alpha=0.6) +
	geom_point(data = dirichlet_basid %>% filter(siteID == "HARV"), aes(x = dates, y = as.numeric(truth_dirichlet))) + xlab(NULL) + #labs(fill='') +
	theme_bw()+
	facet_grid(rows=vars(plotID),
						 #cols = vars(calibration_label),
						 drop=T, scales="free")

# Calculate CRPS scores
models_scored <- all_models %>%
	filter(!is.na(mean_dirichlet) & !is.na(Mean)) %>%
	mutate(truth_dirichlet = as.numeric(truth_dirichlet),
				 truth = as.numeric(truth),
				 crps_dirichlet = crps_norm(truth_dirichlet, mean_dirichlet, sd_dirichlet),
				 crps_beta = crps_norm(truth, Mean, SD))

scored_rank <- models_scored %>%
	group_by(rank_name, taxon) %>%
	summarise(mean_CRPS_dirichlet = mean(crps_dirichlet, na.rm=T),
						mean_CRPS_single = mean(crps_beta, na.rm=T)) %>%
	pivot_longer(cols=3:4)

ggplot(scored_rank) + geom_point(aes(x = taxon, y = value, color = name)) +
	facet_wrap(~rank_name, drop = T, scales = "free")


# View actual models & hindcasts for a single site
a <- ggplot(all_models %>% filter(taxon=="basidiomycota" &
	group=="ITS" & model_name == "all_covariates" & siteID == "HARV")) +
	facet_grid(rows=vars(plotID),
						 #cols = vars(calibration_label),
						 drop=T, scales="free") +
	geom_line(aes(x = dates, y = `50%`), show.legend = F) +
	geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill="red", alpha=0.6) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='') +
	xlim(c(as.Date("2013-06-01"), as.Date("2020-01-01"))) + ggtitle("Fungal diversity hindcasts at NEON site: CPER")


# View actual models & hindcasts for a single site
b <- ggplot(all_models %>% filter(taxon=="basidiomycota" &
																		group=="ITS" & model_name == "all_covariates" & siteID == "HARV")) +
	facet_grid(rows=vars(plotID),
						 #cols = vars(calibration_label),
						 drop=T, scales="free") +
	geom_line(aes(x = dates, y = `50%`), show.legend = F) +
	geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill="red", alpha=0.6) +
	theme_bw()+
	scale_fill_brewer(palette = "Paired") +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='') +
	xlim(c(as.Date("2013-06-01"), as.Date("2020-01-01"))) + ggtitle("Fungal diversity hindcasts at NEON site: CPER")
ggarrange(a,b)
