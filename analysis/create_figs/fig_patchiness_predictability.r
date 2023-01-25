# View predictability by beta effects
library(cowplot)

source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
library(ggrepel)

scores_list <- readRDS(here("data/summary/scoring_metrics_cv.rds"))
cv_metric_scaled <- scores_list$cv_metric_scaled

unconverged = scores_list$unconverged_list %>% filter(median_gbr > 3)

scores_list$unconverged_list
tax_scores <- scores_list$scoring_metrics_cv %>%
	mutate(taxon_model_rank = paste(taxon, model_name, rank)) %>%
	filter(!taxon_model_rank %in% unconverged$taxon_model_rank) %>%
	filter(model_name == "all_covariates" &
				 	metric %in% c("CRPS_truncated","RSQ","RSQ.1") &
				 	site_prediction == "New time (observed site)" &
				 	fcast_type != "Diversity")  %>%
	filter(metric == "RSQ")


ggplot(tax_scores, aes(x = score, y = per_site_cv)) +
	geom_point(aes( color = pretty_group), alpha=.5) +
	geom_smooth(span=1, position=position_jitter(height = .01, width = .01), method="loess") +
	#	facet_grid(~metric) +
	stat_cor(
		aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
		p.digits = 1, label.x.npc = .65, label.y.npc = .65) +
	#geom_rug(aes( color = pretty_group), alpha=.05, sides="t", show.legend = F) +
	geom_rug(aes( color = pretty_group), alpha=.05, sides="r") +
	theme_minimal(base_size = 18) +
	ylab("Patchiness (variation within site)") +
	xlab("Predictability (RSQ at new times)") +
	labs(color="Domain") + coord_flip()

