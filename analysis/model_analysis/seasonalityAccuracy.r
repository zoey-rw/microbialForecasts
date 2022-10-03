# from 05_explainSiteEffects.r script
site_descr <- readRDS(here("data/summary/site_effect_predictors.rds"))
site_descr$latitude <- round(site_descr$latitude, 1)

# from 04_seasonalEffSize.r script
seasonality_in <- readRDS(here("data/summary/seasonalAmplitude.rds"))
all_cov_vals <- seasonality_in[[1]]
cycl_vals <- seasonality_in[[2]]

season_site_df_cycl <- merge(site_descr, seasonality_in[[2]])
season_site_df_all_cov <- merge(site_descr, seasonality_in[[1]])


predictability_scores <- hindcast_filter
cycl_vals_scores <- merge(cycl_vals %>% select(-model_name), predictability_scores)
all_cov_vals_scores <- merge(all_cov_vals %>% select(-model_name), predictability_scores)


sum.all.tax_refit <- readRDS(here("data/summary/single_taxon_summaries_201511_202001.rds"))
tax_plot_est <- sum.all.tax_refit[[2]]
tax_plot_est$month <- lubridate::month(tax_plot_est$dates)
tax_plot_est$month_label <- lubridate::month(tax_plot_est$dates, label = T)
tax_refit_site <- merge(tax_plot_est, site_descr)

sum.all.tax_cal <- readRDS(here("data/summary/single_taxon_summaries_201511_201801.rds"))
tax_plot_est <- sum.all.tax_cal[[2]]
tax_plot_est$month <- lubridate::month(tax_plot_est$dates)
tax_plot_est$month_label <- lubridate::month(tax_plot_est$dates, label = T)
tax_cal_site <- merge(tax_plot_est, site_descr)


# Merge seasonality with hindcast data
# All cov model
phylum_all_cov <-  tax_refit_site %>% filter(model_name=="all_covariates") %>%
	filter(rank_only == "phylum") %>%
	distinct()
phylum_cycl_only <-  tax_cal_site  %>% filter(model_name=="cycl_only") %>%
	filter(rank_only == "phylum") %>%
	distinct()
genus_cycl_only <-  tax_cal_site  %>% filter(model_name=="cycl_only") %>%
	filter(rank_only == "genus") %>%
	distinct()


# Check for natural break in seasonality: looks like .12
ggplot(cycl_vals,
					aes(x=amplitude)) +
	theme_bw(base_size = 18) +
	geom_histogram() +
	geom_vline(xintercept=.12, linetype =2, color = 2)

high_seasonality_tax <- cycl_vals %>% filter(amplitude > .12)


pass_filter = readRDS(here("data/summary/tax_filter_pass.rds"))

tax_refit_filter = tax_refit_site %>% merge(pass_filter, all.y=T)


mean_abundance <- tax_cal_site %>% filter(model_name=="all_covariates") %>%
	select(-model_name) %>%
	filter(!is.na(truth)) %>%
	group_by(group, taxon, rank_only) %>% summarize(mean_abundance = mean(truth),
																									median_abundance = median(truth))

high_seasonality_tax <-  high_seasonality_tax %>%
	merge(mean_abundance)

mean_abundance_cycl <- mean_abundance %>%
	merge(cycl_vals) %>%
	mutate(species = taxon) %>% select(-model_name)

mean_abundance_all_cov <- mean_abundance %>%
	merge(all_cov_vals) %>%
	mutate(species = taxon)


ggplot(mean_abundance_cycl %>% filter(only_rank %in% c("phylum","genus")),
			 aes(x = mean_abundance, y = amplitude,
			 		#group = ecoregion,
			 		color = only_rank)) +
	geom_point() +
	geom_smooth(method="lm",
							size=1, alpha = .2, na.rm = T, se = F) +
	scale_x_log10() +
	ylab("Seasonal amplitude") +
	xlab("Average abundance") +
	#facet_grid(rows = vars(taxon)) +
	facet_wrap(~pretty_group, scales="free") +
	theme_bw(base_size = 18) +
	stat_cor() + labs(color="Taxonomic rank")


#tax_refit_site_abun <- merge(tax_refit_site, mean_abundance)
ggplot(genus_cycl_only %>% filter(taxon %in% high_seasonality_tax$taxon),
			 aes(x = month, y = `50%`,

			 		group = taxon,
			 		color = pretty_group), alpha=.2) +
	geom_smooth(size=1, alpha = .2, na.rm = T, se = F, span=.99) +
	#scale_x_date(breaks = "2 month") +
	xlab("Month") +
	ylab("Abundances (modeled)") +
	#facet_grid(pretty_group~taxon, scales="free") +
	ggtitle("Seasonal genera are more abundant for bacteria than fungi",
					"Genera with strong seasonality peak in spring or late fall") +
	theme_bw(base_size = 18) + labs(color = "Domain")


#scores_list = readRDS(here("data", paste0("summary/scoring_metrics_cv.rds")))

# hindcast_wide_site from fig1 script
site_scores_descr <- merge(site_descr, hindcast_wide_site)
site_scores_descr <- merge(site_scores_descr, cycl_vals %>% select(-model_name))

site_effects_all <- readRDS(here("data/summary/site_effects.rds"))
site_effects_refit_only <- site_effects_all %>% filter(time_period == "2015-11_2020-01" &
																											 	model_name == "all_covariates") %>% mutate(site_effect = Mean)

site_scores_descr <- merge(site_scores_descr, site_effects_refit_only)


ggplot(site_scores_descr,
			 aes(x = latitude, y = CRPS_truncated)) +
	geom_point(aes(color = pretty_name), position = position_jitter(height=0, width=.1)) +
	geom_smooth(method="lm",
		size=1, alpha = .2) +
	facet_wrap(~pretty_name) +
	xlab("Site latitude") +
	ylab("Site hindcast accuracy") +
	theme_bw(base_size = 18) +labs(colour="Rank") +
	ggtitle("Hindcasts improve with proximity to equator") +
	stat_cor(aes(color = amplitude), label.x = 30)


ggscatter(site_scores_descr, #%>%
						#filter(pretty_name == "Functional group"),
					x = "latitude", y = "CRPS_truncated",
					add = "reg.line",                         # Add regression line
					color = "pretty_name" , fullrange = F,           # Color by groups "cyl"
)+
	stat_cor(aes(color = amplitude), label.x = 20) +
	labs(colour="Rank") +
	ggtitle("Hindcasts improve with proximity to equator") +
	xlab("Site latitude") +
	ylab("Site hindcast accuracy")




summary(lm(CRPS_truncated ~ latitude*amplitude, site_scores_descr))

ggplot(site_scores_descr,
			 aes(x = amplitude, y = CRPS_truncated, color = pretty_name)) +
	geom_point(show.legend = F) +
	geom_smooth(method="lm",
							size=1, alpha = .2, show.legend = F) +
	facet_wrap(~pretty_name) +
	xlab("Seasonal amplitude of taxon") +
	ylab("Site CRPS (hindcast accuracy)") +
	stat_cor(color=1) +
	theme_bw(base_size = 18) +
	ggtitle("Seasonality is associated with worse forecasts?")


ggplot(site_scores_descr,
			 aes(x = amplitude, y = latitude)) +
	geom_jitter(aes(color = CRPS_truncated, size = CRPS_truncated), alpha=.3) +
	xlab("Seasonal amplitude of taxon") +
	ylab("Latitude") +
	theme_bw(base_size = 18) +
	ggtitle("Seasonality & proximity from equator \nare associated with worse forecasts") +
	labs(color = "Hindcast accuracy", size = "Hindcast accuracy")




# Site effect ~ latitude
ggplot(site_scores_descr,
			 aes(x = site_effect, y = latitude)) +
	geom_jitter(aes(color = pretty_name), alpha=.3)  +
	facet_wrap(~pretty_name)
