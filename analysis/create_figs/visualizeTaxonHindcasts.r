
# View forecasts for individual taxonomic groups
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

pacman::p_load(Rfast, moments, ggpubr)

hindcasts_raw <- readRDS(here("data/summary/all_hindcasts.rds"))

hindcasts <- hindcasts_raw %>% mutate(taxon = species)
scoring_metrics_cv <- readRDS(here("./data/summary/scoring_metrics_cv.rds"))
phy_scores <- scoring_metrics_cv$scoring_metrics %>% filter(model_name == "all_covariates" &
																													pretty_name %in% c("Phylum","Functional group"))
genus_scores <- scoring_metrics_cv$scoring_metrics %>% filter(model_name == "all_covariates" &
																															pretty_name %in% c("Genus"))

fg_scores <- scoring_metrics_cv$scoring_metrics %>% filter(model_name == "all_covariates" &
																															pretty_name %in% c("Functional group"))


calibration_only = hindcasts %>% filter(fcast_period=="calibration") %>% filter(timepoint > plot_start_date)
asco_calibration <- calibration_only %>% filter(taxon == "ascomycota") %>% filter(!is.na(truth) & !is.na(mean))
asco_calibration_wide <- asco_calibration %>% filter(!siteID %in% c("DEJU","LENO")) %>%
	select(siteID, plotID, dateID, model_name, taxon, mean, truth) %>%
	pivot_wider(names_from = model_name, values_from = c(mean))


# For the calibration, the first observed date per plot is always wonky due to model structure, so we leave it out of our scoring metrics


taxon_scores = phy_scores %>%
	group_by(taxon, model_name, site_prediction) %>%
	summarize(mean_RSQ = mean(RSQ)) %>%
	pivot_wider(names_from = site_prediction, values_from ="mean_RSQ")

worse_at_new_sites = c("acidobacteriota",#"candidatus.solibacter",
											 "penicillium")
fine_at_new_sites = c("ascomycota",#"n_fixation",
											#"oligotroph",
											"copiotroph")



plot_model(hindcasts, taxon = "acidobacteriota", siteID = "BART")

new_plots = table(hindcasts[hindcasts$new_site==T & !is.na(hindcasts$truth) & hindcasts$pretty_group !="Fungi",]$plotID) %>% sort %>% head(30)
obs_plots = table(hindcasts[hindcasts$fcast_period == "calibration" & !is.na(hindcasts$truth) & hindcasts$pretty_group !="Fungi",]$plotID) %>% sort %>% head(30)
new_plots_fun = table(hindcasts[hindcasts$new_site==T & !is.na(hindcasts$truth) & hindcasts$pretty_group=="Fungi",]$plotID) %>% sort %>% head(30)
obs_plots_fun = table(hindcasts[hindcasts$fcast_period == "calibration" & hindcasts$pretty_group=="Fungi" & !is.na(hindcasts$truth),]$plotID) %>% sort %>% head(30)

intersect(names(obs_plots),names(obs_plots_fun))
intersect(names(new_plots),names(new_plots_fun))
plots_to_viz = c(#"CLBJ_006","DELA_038",
	"HARV_034","CPER_001","DSNY_001",
								 #"DSNY_044","TALL_001","YELL_001","MLBS_001",
								 "BONA_004",#"MLBS_002",
	"SOAP_031")

input_df <- hindcasts %>% mutate(fcast_period_alpha = ifelse(fcast_period == "calibration", .8, .4)) %>%
																 filter(model_name =="all_covariates" &
																 			 	taxon %in% c(worse_at_new_sites,fine_at_new_sites) &
																 			 	#siteID %in% c("HARV","DSNY","ORNL"))
																 			 	plotID %in% plots_to_viz) %>%
	filter(!grepl("random",site_prediction))


acido <- hindcasts %>% mutate(fcast_period_alpha = ifelse(fcast_period == "calibration", .8, .4)) %>%
	filter(model_name =="all_covariates" &
				 	taxon %in% c("acidobacteriota")) %>%
	filter(!grepl("random",site_prediction))



penicillium <- hindcasts %>% mutate(fcast_period_alpha = ifelse(fcast_period == "calibration", .8, .4)) %>%
	filter(model_name =="all_covariates" &
				 	taxon %in% c("penicillium")) %>%
	filter(!grepl("random",site_prediction)) %>% group_by(siteID) %>% summarize(avg = mean(truth, na.rm=T)) %>% arrange(avg)


input_df <- hindcasts %>% mutate(fcast_period_alpha = ifelse(fcast_period == "calibration", .8, .4)) %>%
	filter(model_name =="all_covariates" &
				 	taxon %in% fine_at_new_sites &
				 	#siteID %in% c("HARV","DSNY","ORNL"))
				 	plotID %in% c("HARV_001","YELL_001")) %>%
	filter(!grepl("random",site_prediction))

ggplot(input_df, aes(x = dates, group=plotID))  +
	ylab("Abundance") +
	theme_bw() +
	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
				legend.position = "bottom",legend.title = element_text(NULL),
				plot.margin = unit(c(.2, .2, 2, .2), "cm")) +
	facet_grid(newsite~taxon, scales="free", drop=T) +
	geom_ribbon(aes(ymin = lo, ymax = hi, fill = newsite), alpha=0.3) +
	geom_ribbon(aes(ymin = lo_25, ymax = hi_75, fill = newsite), alpha=0.5) +
	geom_line(aes(y = med, color = newsite), show.legend = F) +
	geom_jitter(aes(y = as.numeric(truth)), width = .2, height = 0) +
	xlab(NULL)

# Maybe bad forecasts are just unconverged models????
single_tax_summaries <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/single_taxon_summaries_201511_201801.rds")

all_covariates <- hindcasts_raw %>% filter(model_name == "all_covariates" & fcast_type != "Diversity")

p1 <- plot_model(all_covariates, taxon = "acidobacteriota", siteID = "BART")




plot_model(all_out, taxon = "actinobacteriota", siteID = "CPER", site_plots = "facet")
all_rank_list <- list()
for (rank in names(keep_list)){
	print(rank)
	rank_taxa <- keep_list[[rank]]$taxon.name
	rank_plot_list <- list()
	for (i in rank_taxa){
		rank_plot_list[[i]] <- plot_model(all_out, taxon = i, siteID = "BART", model_type = "all") +
			scale_y_sqrt()
	}
	all_rank_list[[rank]] <- rank_plot_list
}

ggpubr::ggarrange(plotlist = all_rank_list[["order_fun"]][1:10], common.legend = T)
ggpubr::ggarrange(plotlist = rank_plot_list[1:10], common.legend = T)

input_df <- all_out %>% filter(species == i)

input_df$fcast_period <- ifelse(input_df$dates > "2018-01-01", "hindcast", "calibration")
not_na_hindcast <- input_df %>% filter(fcast_period == "hindcast")# & !is.na(truth))




pdf("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/figures/hindcasts_good_example.pdf", height = 6)
for(i in 1:length(top_plots)){
}
	dev.off()
