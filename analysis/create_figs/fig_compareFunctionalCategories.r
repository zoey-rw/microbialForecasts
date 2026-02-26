# Compare functional group forecasts by category/kingdom

pacman::p_load(scoringRules, reshape2, parallel, lubridate, data.table, ggforce, ggrepel)
source("source.R")

# Increase memory limit to handle large data files
tryCatch({
  mem.maxVSize(Inf)
  cat("Memory limit increased to unlimited\n")
}, error = function(e) {
  cat("Note: Could not increase memory limit:", e$message, "\n")
})

# Read in hindcasts - try Parquet first for memory efficiency
parquet_file <- here("data/summary/parquet/all_hindcasts_plsr2.parquet")
rds_file <- here("data/summary/all_hindcasts_plsr2.rds")

if (file.exists(parquet_file)) {
  # Try arrow first, then nanoparquet, then fallback to RDS
  if (requireNamespace("arrow", quietly = TRUE)) {
    cat("Using Parquet file with arrow for memory efficiency...\n")
    hindcast_data <- arrow::read_parquet(parquet_file)
  } else if (requireNamespace("nanoparquet", quietly = TRUE)) {
    cat("Using Parquet file with nanoparquet for memory efficiency...\n")
    hindcast_data <- nanoparquet::read_parquet(parquet_file)
  } else {
    cat("Parquet file exists but neither arrow nor nanoparquet available, using RDS file...\n")
    if (file.exists(rds_file)) {
      hindcast_data <- readRDS(rds_file)
    } else {
      stop("RDS file not found!")
    }
  }
} else if (file.exists(rds_file)) {
  cat("Parquet file not found, using RDS file...\n")
  hindcast_data <- readRDS(rds_file)
} else {
  stop("Neither Parquet nor RDS hindcast files found!")
}

# Read in hindcast scores
scores_list = readRDS(here("data", paste0("summary/scoring_metrics_plsr2.rds")))

# Check available columns first
cat("Available columns in scoring_metrics_long:\n")
print(colnames(scores_list$scoring_metrics_long))

# Subset to functional groups - only select existing columns
available_cols <- colnames(scores_list$scoring_metrics_long)
cols_to_select <- c("fcast_type", "pretty_group", "model_name", "pretty_name", "taxon")
cols_to_select <- cols_to_select[cols_to_select %in% available_cols]

fcast_info_simple <- scores_list$scoring_metrics_long %>% ungroup %>%
	select(all_of(cols_to_select)) %>%
	distinct()

# Get all functional taxa from scoring metrics (40+ functional groups)
# Define functional taxa based on what's available in the data
functional_taxa <- c("cellulose_complex", "acetogen_anaerobic", "assim_nitrate_reduction", 
                    "assim_nitrite_reduction", "benomyl_antibiotic", "cellobiose_complex",
                    "chitin_complex", "chitinolytic", "copiotroph", "dissim_nitrate_reduction",
                    "dissim_nitrite_reduction", "erythromycin_antibiotic", "gentamycin_antibiotic",
                    "glucose_simple", "glycerol_simple", "heat_stress", "herbicide_stress",
                    "lignolytic", "n_fixation", "oligotroph", "animal_pathogen", "lichenized",
                    "streptomycin_antibiotic", "sucrose_complex", "talaromyces")

# Check if taxon column exists, if not use rank_name or species
if("taxon" %in% colnames(scores_list$scoring_metrics_long)) {
  taxon_col <- "taxon"
} else if("rank_name" %in% colnames(scores_list$scoring_metrics_long)) {
  taxon_col <- "rank_name"
} else if("species" %in% colnames(scores_list$scoring_metrics_long)) {
  taxon_col <- "species"
} else {
  stop("No suitable taxon column found")
}

cat("Using column:", taxon_col, "for functional taxa filtering\n")

# Get all models that have functional taxa
all_models_with_fg <- scores_list$scoring_metrics_long %>% 
  filter(!!sym(taxon_col) %in% functional_taxa) %>% 
  pull(model_id) %>% 
  unique()

converged_models <- all_models_with_fg
cat("Using all models with functional taxa:", length(converged_models), "\n")

# Get functional groups - these are taxa with functional names

fg_rsq <- scores_list$scoring_metrics_long %>%
	filter(model_id %in% converged_models) %>%
	filter(metric %in% c("RSQ.1","RSQ","RMSE.norm")) %>%
	mutate(score = ifelse(score < 0, 0, score)) %>%
	filter(!!sym(taxon_col) %in% functional_taxa) %>%
	distinct()  %>%
	merge(fcast_info_simple, all.x=T, all.y=F)

cat("Functional group data loaded:", nrow(fg_rsq), "rows\n")
cat("Unique functional groups:", length(unique(fg_rsq[[taxon_col]])), "\n")

# Check if we have enough data for analysis
if (nrow(fg_rsq) == 0) {
  stop("No functional group data available for analysis")
}

if (length(unique(fg_rsq[[taxon_col]])) < 3) {
  stop("Insufficient functional groups for statistical analysis. Need at least 3 groups, found: ", length(unique(fg_rsq[[taxon_col]])))
}

fg_rsq$fg_category <- microbialForecast:::assign_fg_categories(fg_rsq[[taxon_col]]) %>% make.names
fg_rsq$fg_source <- assign_fg_sources(fg_rsq[[taxon_col]]) #%>% make.names



pretty_names <- list("cellulolytic" = "Cellulose degraders",
										 "assim_nitrite_reduction" = "Assimilatory nitrite reducers",
										 "dissim_nitrite_reduction" = "Dissimilatory nitrite reducers",
										 "assim_nitrate_reduction" = "Assimilatory nitrate reducers",
										 "n_fixation" = "Nitrogen fixers",
										 "dissim_nitrate_reduction" = "Dissimilatory nitrate reducers",
										 "nitrification" = "Nitrifiers",
										 "denitrification" = "Denitrifiers",
										 "chitinolytic" = "Chitin degraders",
										 "lignolytic" = "Lignin degraders",
										 "methanotroph" = "Methanotrophs",
										 "copiotroph" = "Copiotrophs",
										 "oligotroph" = "Oligotrophs",
										 "benomyl_antibiotic" = "Benomyl-resistant",
										 "glucose_simple" = "Glucose-enriched",
										 "pyruvate_simple" = "Pyruvate-enriched",
										 "streptomycin_antibiotic" = "Streptomycin-resistant",
										 "sucrose_complex"  = "Sucrose-enriched",
										 "acetogen_anaerobic" = "Acetogen anaerobic",
										 "chloramphenicol_antibiotic"  = "Chloramphenicol-resistant",
										 "erythromycin_antibiotic"  = "Erythromycin-resistant",
										 "gentamycin_antibiotic"  = "Gentamycin-resistant",
										 "glycerol_simple"   = "Glycerol-enriched",
										 "acetate_simple"  = "Acetate-enriched",
										 "acidic_stress"   = "Acidic stress-tolerant",
										 "cellobiose_complex"   = "Cellobiose-enriched",
										 "cellulose_complex"   = "Cellulose-enriched",
										 "chitin_complex"   = "Chitin-enriched",
										 "galactose_simple"   = "Galactose-enriched",
										 "xylose_simple"   = "Xylose-enriched",
										 "salt_stress" = "Salt stress-tolerant",
										 "herbicide_stress" = "Herbicide stress-tolerant",
										 "osmotic_stress" = "Osmotic stress-tolerant",
										 "heat_stress" = "Heat stress-tolerant",
										 "light_stress" = "Light stress-tolerant",
										 "arbuscular" = "Arbuscular mycorrhizae",
										 "endophyte" = "Endophyte",
										 "litter_saprotroph" = "Litter saprotrophs",
										 "lichenized" = "Lichenized fungi",
										 "animal_pathogen" = "Animal pathogens",
										 "plant_pathogen" = "Plant pathogens",
										 "saprotroph" = "Saprotrophs",
										 "wood_saprotroph" = "Wood saprotrophs",
										 "ectomycorrhizal" = "Ectomycorrhizae"
)
fg_rsq$pretty_fg_names <- recode(fg_rsq[[taxon_col]], !!!pretty_names)# %>% make.names
#fg_rsq$pretty_fg_names <- recode(fg_rsq$taxon, !!!microbialForecast:::pretty_names)


# Statistical tests - use all available data
stat_pvalue_fg_source <- fg_rsq %>%
	filter(metric %in% "RMSE.norm" &
				 	site_prediction == "New time (observed site)") %>%
	rstatix::tukey_hsd(score ~ fg_source) %>%
	#filter(p.adj < 0.05) %>%
	rstatix::add_y_position(step.increase = .02) #%>%
#	mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))

stat_pvalue_fg_category <- fg_rsq %>%
	filter(metric %in% "RMSE.norm" &
				 	site_prediction == "New time (observed site)") %>%
	rstatix::tukey_hsd(score ~ fg_category) %>%
	#filter(p.adj < 0.05) %>%
	rstatix::add_y_position(step.increase = .02) #%>%
	#mutate(y.position = seq(min(y.position), max(y.position),length.out = n()))

pos <- position_jitter(width = 0.1, height=0, seed = 1)

# Looking at assignment method
fg_source = ggplot(fg_rsq %>%
			 				 	filter(metric %in% "RMSE.norm" &
			 				 				 	site_prediction == "New time (observed site)"),
			 aes(x = fg_source, y = as.numeric(score)))  +
	#geom_boxplot(alpha=.5) +
	geom_jitter(aes(color=pretty_group),
							size=4, #width = .3, height=0,
							alpha=.5,
							position=pos) +
	xlab(NULL) + labs(fill='') +
	theme_classic(base_size = 22) +
	theme(axis.text.x=element_text(size=20,
			angle = 320, vjust=1, hjust = -0.05)	) +
	ylab("Forecast error (nRMSE)") + labs(color=NULL) +
	geom_text_repel(aes(label=pretty_fg_names), size=6,
									 max.overlaps = 10,
									 position=pos) +
	ggpubr::	stat_pvalue_manual(stat_pvalue_fg_source,
														 label = "p.adj.signif",
														 #bracket.nudge.y = -1.1,
size=7, hide.ns = T) +
	scale_y_continuous(trans = "log10")


png(here("figures","functional_group_source.png"), width = 800, height=1000)
print(fg_source)
dev.off()

# Looking at functional group category
fg_category_data <- fg_rsq %>%  filter(metric %in% "RMSE.norm" &
														site_prediction == "New time (observed site)")
if (nrow(fg_category_data) > 0 && length(unique(fg_category_data$fg_category)) > 0) {
  # After filtering to single metric and site_prediction, we don't need faceting
  # Check if we have both metric and site_prediction columns with multiple values
  n_metrics <- if ("metric" %in% names(fg_category_data)) length(unique(fg_category_data$metric)) else 0
  n_site_preds <- if ("site_prediction" %in% names(fg_category_data)) length(unique(fg_category_data$site_prediction)) else 0
  
  p_fg_category <- ggplot(fg_category_data,
  			 aes(x = fg_category, y = as.numeric(score)))  +
  	geom_boxplot() +
  	geom_jitter( size=3, width = .1, height=0, alpha=.1) +
  	xlab(NULL) + labs(fill='') +
  	stat_compare_means() + theme_bw()  + theme_bw(base_size = 18) +
  	theme(
  		axis.text.x=element_text(
  			angle = 320, vjust=1, hjust = -0.05)	) + ylab("RSQ") +
  	#geom_text(aes(label=pretty_fg_names), hjust=0, vjust=0) +
  	ggpubr::	stat_pvalue_manual(stat_pvalue_fg_category,
  														 label = "p.adj.signif", #, bracket.nudge.y = -.4,
  														 size=7, hide.ns = T)
  
  # Only add faceting if both columns exist and have multiple values
  # Since we filtered to single metric and site_prediction, we skip faceting here
  # (faceting would fail with only one value in each)
  ggsave(here("figures", "fg_category_rsq.png"), p_fg_category, width = 8, height = 6, dpi = 200)
} else {
  cat("No data available for fg_category plot\n")
}


# Not restricted by site-pred type or metric
# Check if we have data and required columns before plotting
if (nrow(fg_rsq) > 0) {
  has_metric <- "metric" %in% names(fg_rsq) && length(unique(fg_rsq$metric)) > 0
  has_site_pred <- "site_prediction" %in% names(fg_rsq) && length(unique(fg_rsq$site_prediction)) > 0
  
  n_metrics <- if (has_metric) length(unique(fg_rsq$metric)) else 0
  n_site_preds <- if (has_site_pred) length(unique(fg_rsq$site_prediction)) else 0
  
  p1 <- ggplot(fg_rsq, aes(x = fg_category, y = as.numeric(score)))  +
  	geom_boxplot() +
  	geom_jitter(aes(color=pretty_group), size=3, width = .1, height=0, alpha=.1) +
  	xlab(NULL) + labs(fill='') +
  	stat_compare_means() + theme_bw()
  
  # Only add faceting if we have multiple values in both columns
  if (n_metrics > 1 && n_site_preds > 1) {
    p1 <- p1 + facet_grid(metric~site_prediction)
  } else if (n_metrics > 1) {
    p1 <- p1 + facet_wrap(~metric)
  } else if (n_site_preds > 1) {
    p1 <- p1 + facet_wrap(~site_prediction)
  }
  # Skip printing if no faceting variables (would cause error)
  if (n_metrics > 1 || n_site_preds > 1) {
    print(p1)
  }

  p2 <- ggplot(fg_rsq, aes(x = fg_source, y = as.numeric(score)))  +
  	geom_boxplot() +
  	geom_jitter(aes(color=pretty_group), size=3, width = .1, height=0, alpha=.1) +
  	xlab(NULL) + labs(fill='') +
  	stat_compare_means() + theme_bw()
  
  # Only add faceting if we have multiple values in both columns
  if (n_metrics > 1 && n_site_preds > 1) {
    p2 <- p2 + facet_grid(metric~site_prediction)
  } else if (n_metrics > 1) {
    p2 <- p2 + facet_wrap(~metric)
  } else if (n_site_preds > 1) {
    p2 <- p2 + facet_wrap(~site_prediction)
  }
  # Skip printing if no faceting variables (would cause error)
  if (n_metrics > 1 || n_site_preds > 1) {
    print(p2)
  }
} else {
  cat("No data available for unrestricted functional group plots\n")
}



# Use current hindcast data structure
both_hindcast_periods <- hindcast_data

fungi_allcov <- both_hindcast_periods  %>%
	filter(taxon %in% c("animal_pathogen", "ectomycorrhizal", "endophyte", "lichenized",
												"plant_pathogen", "saprotroph"),
				 model_name == "all_covariates")

fungi_cycl <- both_hindcast_periods  %>%
	filter(taxon %in% c("animal_pathogen", "ectomycorrhizal", "endophyte", "lichenized",
											"plant_pathogen", "saprotroph"))


# Get best plots for examples (most calibration AND validation points)
not_na_hindcast <- both_hindcast_periods %>% filter(fcast_period == "hindcast" & !is.na(truth))
not_na_calibration <- both_hindcast_periods %>% filter(fcast_period == "calibration" & !is.na(truth))
plot_hindcast_counts <- sort(table(not_na_hindcast$plotID))
plot_calibration_counts <- sort(table(not_na_calibration$plotID))
top_hindcast_plots <- names(tail(plot_hindcast_counts, 30))
top_calibration_plots <- names(tail(plot_calibration_counts, 30))
top_plots <- intersect(top_hindcast_plots, top_calibration_plots)
top_plots



# Check if data exists and faceting variables have values
fungi_harv <- fungi_allcov %>% filter(plotID=="HARV_013")
if (nrow(fungi_harv) > 0) {
  p_fungi_harv <- ggplot(fungi_harv) +
  	geom_line(aes(x = dates, y = med), show.legend = F, linetype=2) +
  	geom_line(aes(x = dates, y = `50%`), show.legend = F) +
  	geom_ribbon(aes(x = dates, ymin = lo, ymax = hi), alpha=0.6, fill="blue") +
  	geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill="red", alpha=0.6) +
  	theme_bw()+
  	scale_fill_brewer(palette = "Paired") +
  	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
  			legend.position = "bottom",legend.title = element_text(NULL),
  			plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
  	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='')
  
  has_taxon <- "taxon" %in% names(fungi_harv) && length(unique(fungi_harv$taxon)) > 0
  has_time_period <- "time_period" %in% names(fungi_harv) && length(unique(fungi_harv$time_period)) > 0
  
  if (has_taxon && has_time_period && length(unique(fungi_harv$taxon)) > 1 && length(unique(fungi_harv$time_period)) > 1) {
    p_fungi_harv <- p_fungi_harv + facet_grid(rows=vars(taxon), cols = vars(time_period), drop=T, scales="free")
  } else if (has_taxon && length(unique(fungi_harv$taxon)) > 1) {
    p_fungi_harv <- p_fungi_harv + facet_wrap(~taxon)
  } else if (has_time_period && length(unique(fungi_harv$time_period)) > 1) {
    p_fungi_harv <- p_fungi_harv + facet_wrap(~time_period)
  }
  ggsave(here("figures", "hindcasts_fungal_harv.png"), p_fungi_harv, width = 10, height = 8, dpi = 200)
} else {
  cat("No data available for HARV_013 plot\n")
}




top_plots
fungi_allcov_top_plots <- fungi_allcov %>%
	filter(plotID %in% top_plots & time_period == "2015-11_2018-01")

# Check if we have data and required columns before plotting
if (nrow(fungi_allcov_top_plots) > 0 && length(top_plots) > 0) {
  has_taxon <- "taxon" %in% names(fungi_allcov_top_plots) && length(unique(fungi_allcov_top_plots$taxon)) > 0
  has_plotID <- "plotID" %in% names(fungi_allcov_top_plots) && length(unique(fungi_allcov_top_plots$plotID)) > 0
  
  if (has_taxon && has_plotID && length(unique(fungi_allcov_top_plots$taxon)) > 1 && length(unique(fungi_allcov_top_plots$plotID)) > 1) {
    png(here("figures","hindcasts_fungal_functional.png"), height = 10, width = 12, units = "in", res = 200)
    for(i in 1:length(top_plots)){
    	p <- ggplot(fungi_allcov_top_plots) +
    					geom_line(aes(x = dates, y = med), show.legend = F, linetype=2) +
    					geom_line(aes(x = dates, y = `50%`), show.legend = F) +
    					geom_ribbon(aes(x = dates, ymin = lo, ymax = hi), alpha=0.6, fill="blue") +
    					geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill="red", alpha=0.6) +
    					theme_bw()+
    					scale_fill_brewer(palette = "Paired") +
    					theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
    								legend.position = "bottom",legend.title = element_text(NULL),
    								plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
    					geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='')
    	
    	# Only add faceting if we have multiple values
    	if (length(unique(fungi_allcov_top_plots$taxon)) > 1 && length(unique(fungi_allcov_top_plots$plotID)) > 1) {
    	  p <- p + facet_grid_paginate(
    						taxon~plotID,
    						drop=T, scales="free",
    						ncol = 1, nrow = 6, page = i)
    	} else if (length(unique(fungi_allcov_top_plots$taxon)) > 1) {
    	  p <- p + facet_wrap(~taxon)
    	} else if (length(unique(fungi_allcov_top_plots$plotID)) > 1) {
    	  p <- p + facet_wrap(~plotID)
    	}
    	print(p)
    }
    dev.off()
  } else {
    cat("Insufficient data for faceting in fungal functional plots\n")
  }
} else {
  cat("No data available for fungal functional plots\n")
}



library(scoringRules)
both_hindcast_periods$newsite <- ifelse(both_hindcast_periods$new_site, "New site", "Observed site")

scored_hindcasts <- both_hindcast_periods %>%
	filter(!is.na(truth) & fcast_period == "hindcast") %>%
	mutate(truth = as.numeric(truth)) %>% mutate(crps = crps_norm(truth, mean, sd))


# Check which grouping columns exist
group_cols <- c("fcast_type", "category", "group", "time_period", "model_name", "newsite", "taxon")
available_group_cols <- group_cols[group_cols %in% names(scored_hindcasts)]
scored_hindcasts_mean <- scored_hindcasts %>% 
  group_by(across(any_of(available_group_cols))) %>%
	dplyr::summarize(crps_mean = mean(crps, na.rm=T), .groups = "drop")

# Check if model_name exists and has values before faceting
scored_hindcasts_obs <- scored_hindcasts_mean %>% filter(newsite=="Observed site")
if (nrow(scored_hindcasts_obs) > 0 && "model_name" %in% names(scored_hindcasts_obs) && length(unique(scored_hindcasts_obs$model_name)) > 0) {
  # Use available grouping column (pretty_group, group, or taxon) instead of category
  color_col <- if ("pretty_group" %in% names(scored_hindcasts_obs)) "pretty_group"
               else if ("group" %in% names(scored_hindcasts_obs)) "group"
               else if ("taxon" %in% names(scored_hindcasts_obs)) "taxon"
               else NULL
  
  p_scored <- ggplot(scored_hindcasts_obs)
  if (!is.null(color_col)) {
    p_scored <- p_scored + geom_violin(aes(x = time_period, y = crps_mean), draw_quantiles = c(.5)) +
      geom_jitter(aes_string(x = "time_period", y = "crps_mean", color = color_col))
  } else {
    p_scored <- p_scored + geom_violin(aes(x = time_period, y = crps_mean), draw_quantiles = c(.5)) +
      geom_jitter(aes(x = time_period, y = crps_mean))
  }
  
  if (length(unique(scored_hindcasts_obs$model_name)) > 1) {
    p_scored <- p_scored + facet_grid(~model_name)
  }
  ggsave(here("figures", "scored_hindcasts_by_period.png"), p_scored, width = 10, height = 6, dpi = 200)
} else {
  cat("No data available for scored hindcasts plot\n")
}


# Skip this plot if category column doesn't exist
scored_filtered <- scored_hindcasts_mean %>% filter(model_name=="all_covariates" &
																						time_period == "2015-11_2018-01" &
																						newsite=="Observed site")
if (nrow(scored_filtered) > 0 && "category" %in% names(scored_filtered)) {
  p_cat <- ggplot(scored_filtered) +
  	geom_violin(aes(x = category, y = crps_mean), draw_quantiles = c(.5)) +
  	geom_jitter(aes(x = category, y = crps_mean, color = category))
  ggsave(here("figures", "scored_by_category.png"), p_cat, width = 8, height = 5, dpi = 200)
} else {
  cat("Skipping category plot - category column not available\n")
}

# Only 16S groups since ITS is own category
scored_16s <- scored_hindcasts_mean  %>%  filter(model_name=="all_covariates" &
                                            time_period == "2015-11_2018-01" &
                                            newsite=="Observed site" &
                                            group == "16S")
if (nrow(scored_16s) > 0 && "category" %in% names(scored_16s)) {
  p_16s <- ggplot(scored_16s) +
    geom_violin(aes(x = reorder(category, crps_mean), y = crps_mean), draw_quantiles = c(.5)) +
    geom_jitter(aes(x = reorder(category, crps_mean), y = crps_mean, color = category), size=3, alpha = .8, show.legend = F) +
    theme_bw(base_size = 18) +xlab("Functional category") +ylab("CRPS score") +
    ggtitle("Predictability of bacterial functional categories")
  ggsave(here("figures", "scored_16s_category.png"), p_16s, width = 10, height = 6, dpi = 200)
} else {
  cat("Skipping 16S category plot - category column not available\n")
}


# Skill score from Scavia 2021
# Only calculate if category column exists
if ("category" %in% names(scored_hindcasts_mean)) {
  id_cols <- c("fcast_type","category","group","model_name","time_period","taxon")
  id_cols <- id_cols[id_cols %in% names(scored_hindcasts_mean)]
  
  skill_score <- scored_hindcasts_mean %>%
  	pivot_wider(id_cols = all_of(id_cols), values_from = "crps_mean", names_from = "newsite") %>%
  	mutate(skill_score = (1 - (`New site`/`Observed site`)))

  skill_score_filtered <- skill_score %>% filter(model_name=="all_covariates" & time_period == "2015-11_2018-01")
  if (nrow(skill_score_filtered) > 0) {
    p_skill <- ggplot(skill_score_filtered) +
    	geom_violin(aes(x = category, y = skill_score), draw_quantiles = c(.5)) +
    	geom_jitter(aes(x = category, y = skill_score, color = category))
    ggsave(here("figures", "skill_score_by_category.png"), p_skill, width = 8, height = 5, dpi = 200)
  }
} else {
  cat("Skipping skill score calculation - category column not available\n")
}






# Check if data exists and model_name has values before faceting
hindcast_oligo <- hindcast_data %>% filter(plotID=="BART_002" & taxon == "oligotroph")
if (nrow(hindcast_oligo) > 0) {
  p_oligo <- ggplot(hindcast_oligo) +
  	geom_line(aes(x = dates, y = med), show.legend = F, linetype=2) +
  	geom_line(aes(x = dates, y = `50%`), show.legend = F) +
  	geom_ribbon(aes(x = dates, ymin = lo, ymax = hi), alpha=0.6, fill="blue") +
  	geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill="red", alpha=0.6) +
  	theme_bw()+
  	scale_fill_brewer(palette = "Paired") +
  	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
  			legend.position = "bottom",legend.title = element_text(NULL),
  			plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
  	geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='')
  
  if ("model_name" %in% names(hindcast_oligo) && length(unique(hindcast_oligo$model_name)) > 1) {
    p_oligo <- p_oligo + facet_grid(cols = vars(model_name), drop=T, scales="free")
  }
  ggsave(here("figures", "hindcast_oligo.png"), p_oligo, width = 10, height = 6, dpi = 200)
} else {
  cat("No data available for oligotroph plot\n")
}

