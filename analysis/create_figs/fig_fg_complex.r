# read in seasonal values
library(lubridate)
source("source.R")
source("microbialForecast/R/load_plot_estimates.r")

scores_list = readRDS(here("data/summary/scoring_metrics_plsr2.rds"))
converged = scores_list$converged_list
#converged_strict = scores_list$converged_strict_list

sum.in <- readRDS(here("data", paste0("summary/logit_beta_fixed_priors_summaries.rds")))

# Load plot estimates using new memory-efficient functions
cat("Loading plot estimates for functional group analysis...\n")
plot_estimates_env_cycl <- load_plot_estimates_for_phenology("env_cycl")
plot_estimates_cycl_only <- load_plot_estimates_for_phenology("cycl_only")

# Combine plot estimates from both model types
plot_estimates <- rbind(plot_estimates_env_cycl, plot_estimates_cycl_only, fill = TRUE)

# Filter to only converged models
plot_estimates <- plot_estimates %>% filter(model_id %in% converged)
plot_estimates$month = lubridate::month(plot_estimates$dates)
cycl_only_est = plot_estimates %>% filter(grepl("cycl_only",model_name))
cycl_only_est$fg_category <- microbialForecast:::assign_fg_categories(cycl_only_est$taxon)


env_cycl_est = plot_estimates %>% filter(grepl("env_cycl",model_name))
env_cycl_est$fg_category <- microbialForecast:::assign_fg_categories(env_cycl_est$taxon)

seas_in = readRDS(here("data/summary/seasonal_amplitude.rds"))


seas_vals_long <- seas_in[[1]] %>% filter(model_id %in% converged)
seas_vals_wide <- seas_in[[6]] %>% filter(model_id %in% converged)

# Extract groups with "max" dates in winter
seas_vals_long$max_month = month(seas_vals_long$max_y_date)
seas_vals_short <- seas_vals_long %>% select(-c(dates,y_cycl)) %>% distinct(.keep_all = T)

cycl_only_vals = seas_vals_short %>%
	filter(model_name == "cycl_only")

cat("cycl_only_vals rows:", nrow(cycl_only_vals), "\n")
if(nrow(cycl_only_vals) > 0) {
  cat("cycl_only_vals fcast_type values:", paste(unique(cycl_only_vals$fcast_type), collapse = ", "), "\n")
  cat("cycl_only_vals columns:", paste(colnames(cycl_only_vals), collapse = ", "), "\n")
}

winter_groups = cycl_only_vals %>% filter(
	time_period=="2015-11_2018-01") %>%
	filter(max_month %in% c(1,2,12))
unique(winter_groups$taxon)


fg_seasonal = cycl_only_vals %>% filter(fcast_type %in% c("Functional", "functional"))
cat("fg_seasonal rows after filter:", nrow(fg_seasonal), "\n")
if(nrow(fg_seasonal) > 0) {
  fg_seasonal$fg_category <- microbialForecast:::assign_fg_categories(fg_seasonal$taxon)
  cat("fg_seasonal unique fg_category:", paste(unique(fg_seasonal$fg_category), collapse = ", "), "\n")
} else {
  # Try without filtering by fcast_type to see what we have
  cat("Trying without fcast_type filter...\n")
  if("fcast_type" %in% colnames(cycl_only_vals)) {
    cat("All fcast_type values in cycl_only_vals:", paste(unique(cycl_only_vals$fcast_type), collapse = ", "), "\n")
  } else {
    cat("fcast_type column not found in cycl_only_vals\n")
  }
}

if (nrow(fg_seasonal) > 0 && length(unique(fg_seasonal$fg_category)) > 0) {
  p_fg_seasonal <- ggplot(fg_seasonal,
  			 aes(x=max_y_date, y=fg_category)) +
  	geom_jitter(aes(colour = pretty_group), width = 0.1, height=.1, alpha=.5, show.legend = F) +
  	theme_bw(base_size = 20) +
  	ggtitle(paste0("Peak month of seasonal trend")) +
  	xlab(NULL) + ylab(NULL)
  print(p_fg_seasonal)
  tryCatch({
    png(here("figures","fg_complex_seasonal_peak.png"), width = 1200, height = 800)
    print(p_fg_seasonal)
    dev.off()
  }, error = function(e) {
    cat("Error saving fg_complex_seasonal_peak:", e$message, "\n")
  })
} else {
  cat("No data available for fg_seasonal plot\n")
}
#facet_grid(rows=vars(taxon_name)) +



cycl_only_est$month_date = as.Date(paste0(cycl_only_est$month, "-01-2016"), format = "%m-%d-%Y")

# "pretty_names" is from globalVariables.r but isn't exporting...
cycl_only_est$pretty_fg_names <- recode(cycl_only_est$taxon, !!!microbialForecast:::pretty_names)

simple_sugars = cycl_only_est %>% filter(fg_category=="Simple substrates" & fcast_type=="Functional")
complex_sugars = cycl_only_est %>% filter(fg_category=="Complex substrates" & fcast_type=="Functional")

simple_complex = cycl_only_est %>% filter(fg_category %in% c("Simple substrates","Complex substrates") & fcast_type=="Functional") %>%
	filter(!taxon %in% c("chitinolytic","cellulolytic","lignolytic"))



copio_oligo_est = cycl_only_est  %>%
	filter(grepl("otroph", taxon) &  time_period=="2015-11_2018-01")



pheno_categories_in <- readRDS(here("data/clean/modis_greenup.rds"))
bart_harv_pheno = pheno_categories_in[[1]] %>% filter(ID %in% c("BART","HARV") & year == "2016") %>% ungroup

northern_sites <- c("ABBY", "BART", "BONA",  "DEJU", "HARV",
										"HEAL", "NIWO", "ONAQ", "RMNP", "STEI", "TOOL",
										"TREE", "UNDE", "WREF", "YELL")

northern_sites <- c("BART", "BONA", "HARV",
										"HEAL", "STEI", "TOOL",
										"TREE", "UNDE", "WREF", "YELL")

northern_sites <- c("HARV",
										"BART")

plant_associated = cycl_only_est  %>%
	filter(time_period=="2015-11_2018-01")  %>%
	filter(taxon %in% c("plant_pathogen","oligotroph","heat_stress","copiotroph","saprotroph"))

plant_associated_filtered <- plant_associated %>% filter(siteID %in% northern_sites) %>% filter(taxon != "heat_stress")
if (nrow(plant_associated_filtered) > 0 && length(unique(plant_associated_filtered$pretty_fg_names)) > 0) {
  ggplot(plant_associated_filtered,
  			 aes(x = month_date, color=siteID)) +
  
  	geom_smooth(aes(y = `Mean`), method="loess", span=1, se=F) +
  	theme_classic()+
  	scale_fill_brewer(palette = "Paired") +
  	theme(text = element_text(size = 18), panel.spacing = unit(.2, "cm"),
  				legend.position = "bottom",legend.title = element_text(NULL),
  				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
  	xlab(NULL) + labs(fill='') +
  	facet_wrap(~pretty_fg_names, ncol=1, scales="free") +
	annotate(geom = 'rect', xmin=as.Date("2016-01-01"), xmax=as.Date("2016-05-02"), ymin=-Inf, ymax=Inf, alpha=.2, fill='lightgray', ) +
	annotate(geom = 'rect', xmin=as.Date("2016-05-02"), xmax=as.Date("2016-06-26"), ymin=-Inf, ymax=Inf, alpha=.2, fill='lightgreen') +
	annotate(geom = 'rect', xmin=as.Date("2016-06-26"), xmax=as.Date("2016-08-14"), ymin=-Inf, ymax=Inf, alpha=.2, fill='green') +
	annotate(geom = 'rect', xmin=as.Date("2016-08-14"), xmax=as.Date("2016-12-31"), ymin=-Inf, ymax=Inf, alpha=.2, fill='lightgray') +
	scale_x_date(date_labels = "%B") + labs(color = "NEON site")
} else {
  cat("No data available for plant_associated plot\n")
}

simple_complex_filtered <- simple_complex %>% filter(siteID %in% northern_sites) %>% filter(taxon != "heat_stress")
if (nrow(simple_complex_filtered) > 0 && length(unique(simple_complex_filtered$fg_category)) > 0 && 
    length(unique(simple_complex_filtered$pretty_fg_names)) > 0) {
  ggplot(simple_complex_filtered,
  			 aes(x = month_date, color=siteID)) +
  
  	geom_smooth(aes(y = `Mean`), method="loess", span=1, se=F) +
  	theme_classic()+
  	scale_fill_brewer(palette = "Paired") +
  	theme(text = element_text(size = 18), panel.spacing = unit(.2, "cm"),
  				legend.position = "bottom",legend.title = element_text(NULL),
  				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
  	xlab(NULL) + labs(fill='') +
  	facet_wrap(fg_category~pretty_fg_names,  scales="free") +
	annotate(geom = 'rect', xmin=as.Date("2016-01-01"), xmax=as.Date("2016-05-02"), ymin=-Inf, ymax=Inf, alpha=.2, fill='lightgray', ) +
	annotate(geom = 'rect', xmin=as.Date("2016-05-02"), xmax=as.Date("2016-06-26"), ymin=-Inf, ymax=Inf, alpha=.2, fill='lightgreen') +
	annotate(geom = 'rect', xmin=as.Date("2016-06-26"), xmax=as.Date("2016-08-14"), ymin=-Inf, ymax=Inf, alpha=.2, fill='green') +
	annotate(geom = 'rect', xmin=as.Date("2016-08-14"), xmax=as.Date("2016-12-31"), ymin=-Inf, ymax=Inf, alpha=.2, fill='lightgray') +
	scale_x_date(date_labels = "%B") + labs(color = "NEON site")
} else {
  cat("No data available for simple_complex plot\n")
}

simple_sugars_filtered <- simple_sugars %>% filter(siteID %in% northern_sites)
if (nrow(simple_sugars_filtered) > 0 && "taxon" %in% names(simple_sugars_filtered) && 
    length(unique(simple_sugars_filtered$taxon)) > 0) {
  ggplot(simple_sugars_filtered,
  			 aes(fill=species, x = as.numeric(month))) +
  	#geom_point(aes(y = `Mean`), show.legend = F, color="red",
  	#					 alpha=.1, position=position_jitter(height=0)) +
  	#geom_point(aes(y = as.numeric(truth)), alpha = .3, position=position_jitter(height=0)) +
  	geom_smooth(aes(y = `Mean`), method="loess", span=1, se=F) +
  	theme_bw()+
  	scale_fill_brewer(palette = "Paired") +
  	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
  				legend.position = "bottom",legend.title = element_text(NULL),
  				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
  	xlab(NULL) + labs(fill='') +
  	facet_wrap(~taxon, scales="free")
} else {
  cat("No data available for simple_sugars plot\n")
}

complex_sugars_filtered <- complex_sugars %>% filter(siteID %in% northern_sites)
if (nrow(complex_sugars_filtered) > 0 && "taxon" %in% names(complex_sugars_filtered) && 
    length(unique(complex_sugars_filtered$taxon)) > 0) {
  ggplot(complex_sugars_filtered,
  			 aes(fill=species, x = as.numeric(month))) +
  	geom_point(aes(y = `Mean`), show.legend = F, color="red",
  						 alpha=.1, position=position_jitter(height=0)) +
  	geom_point(aes(y = as.numeric(truth)), alpha = .5, position=position_jitter(height=0)) +
  	geom_smooth(aes(y = `Mean`), method="loess", span=1, se=F) +
  	theme_bw()+
  	scale_fill_brewer(palette = "Paired") +
  	theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
  				legend.position = "bottom",legend.title = element_text(NULL),
  				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
  	xlab(NULL) + labs(fill='') +
  	facet_wrap(~taxon,nrow=2, scales="free")
} else {
  cat("No data available for complex_sugars plot\n")
}


copio_oligo_env_cycl_est = env_cycl_est  %>%
	filter(grepl("otroph", taxon) | taxon %in% c("glucose_simple",
																							 "cellulose_complex","acetate_simple") &  time_period=="2015-11_2018-01")

copio_oligo_filtered <- copio_oligo_env_cycl_est %>% filter(siteID %in% northern_sites)
if (nrow(copio_oligo_filtered) > 0 && length(unique(copio_oligo_filtered$fg_category)) > 0 && 
    length(unique(copio_oligo_filtered$taxon)) > 0) {
  # Use Mean instead of 50% if 50% doesn't exist
  y_col_copio <- if ("50%" %in% names(copio_oligo_filtered)) "50%" else if ("Mean" %in% names(copio_oligo_filtered)) "Mean" else "Mean"
  ggplot(copio_oligo_filtered,
  			 aes(x = as.numeric(month))) +
  	geom_point(aes_string(y = y_col_copio, color = "siteID"), show.legend = F,
  						 alpha=.1, position=position_jitter(height=0)) +
  	geom_point(aes(y = as.numeric(truth), color = siteID), alpha = .5,
  						 position=position_jitter(height=0)) +
  	geom_smooth(aes_string(y = y_col_copio), se = F, color=1) +
  	theme_bw( base_size = 22) +
  	scale_fill_brewer(palette = "Paired") +
  	theme(text = element_text(size = 18), panel.spacing = unit(.2, "cm"),
  				legend.position = "bottom",legend.title = element_text(NULL),
  				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
  	xlab(NULL) + labs(fill='') +
  	facet_wrap(fg_category~taxon, scales="free")
} else {
  cat("No data available for copio_oligo_env_cycl plot\n")
}


all_decomposers <- cycl_only_est %>% filter(fcast_type=="Functional") %>%
	filter(fg_category %in% c("Complex substrates", "Simple substrates") | taxon %in% c("saprotroph","copiotroph","oligotroph"))

all_decomposers_filtered <- all_decomposers %>% filter(siteID %in% northern_sites)
if (nrow(all_decomposers_filtered) > 0 && length(unique(all_decomposers_filtered$fg_category)) > 0 && 
    length(unique(all_decomposers_filtered$taxon)) > 0) {
  # Use Mean instead of 50% if 50% doesn't exist
  y_col_decomp <- if ("50%" %in% names(all_decomposers_filtered)) "50%" else if ("Mean" %in% names(all_decomposers_filtered)) "Mean" else "Mean"
  ggplot(all_decomposers_filtered,
  			 aes(x = as.numeric(month))) +
  	geom_point(aes_string(y = y_col_decomp, color = "siteID"), show.legend = F,
  						 alpha=.1, position=position_jitter(height=0)) +
  	geom_point(aes(y = as.numeric(truth), color = siteID), alpha = .5,
  						 position=position_jitter(height=0)) +
  	geom_smooth(aes_string(y = y_col_decomp), se = F, color=1) +
  	theme_bw( base_size = 22) +
  	scale_fill_brewer(palette = "Paired") +
  	theme(text = element_text(size = 18), panel.spacing = unit(.2, "cm"),
  				legend.position = "bottom",legend.title = element_text(NULL),
  				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
  	xlab(NULL) + labs(fill='') +
  	facet_wrap(fg_category~taxon, scales="free")
} else {
  cat("No data available for all_decomposers plot\n")
}


all_decomposers_env <- env_cycl_est %>% filter(fcast_type=="Functional") %>%
	filter(fg_category %in% c("Complex substrates", "Simple substrates") | taxon %in% c("saprotroph","copiotroph","oligotroph"))

all_decomposers_env_filtered <- all_decomposers_env %>% filter(siteID %in% northern_sites)
if (nrow(all_decomposers_env_filtered) > 0 && length(unique(all_decomposers_env_filtered$fg_category)) > 0 && 
    length(unique(all_decomposers_env_filtered$taxon)) > 0) {
  # Use Mean instead of 50% if 50% doesn't exist
  y_col_decomp_env <- if ("50%" %in% names(all_decomposers_env_filtered)) "50%" else if ("Mean" %in% names(all_decomposers_env_filtered)) "Mean" else "Mean"
  ggplot(all_decomposers_env_filtered,
  			 aes(x = as.numeric(month))) +
  	geom_point(aes_string(y = y_col_decomp_env, color = "siteID"), show.legend = F,
  						 alpha=.1, position=position_jitter(height=0)) +
  	geom_point(aes(y = as.numeric(truth), color = siteID), alpha = .5,
  						 position=position_jitter(height=0)) +
  	geom_smooth(aes_string(y = y_col_decomp_env), se = F, color=1) +
  	theme_bw( base_size = 22) +
  	scale_fill_brewer(palette = "Paired") +
  	theme(text = element_text(size = 18), panel.spacing = unit(.2, "cm"),
  				legend.position = "bottom",legend.title = element_text(NULL),
  				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
  	xlab(NULL) + labs(fill='') +
  	facet_wrap(fg_category~taxon, scales="free") + ggtitle("Full models, estimates plotted by month ")
} else {
  cat("No data available for all_decomposers_env plot\n")
}


all_decomposers_env_harv_bart <- all_decomposers_env %>% filter(siteID %in% c("HARV","BART"))
if (nrow(all_decomposers_env_harv_bart) > 0 && length(unique(all_decomposers_env_harv_bart$fg_category)) > 0 && 
    length(unique(all_decomposers_env_harv_bart$taxon)) > 0) {
  # Use Mean instead of 50% if 50% doesn't exist
  y_col_harv_bart <- if ("50%" %in% names(all_decomposers_env_harv_bart)) "50%" else if ("Mean" %in% names(all_decomposers_env_harv_bart)) "Mean" else "Mean"
  ggplot(all_decomposers_env_harv_bart,
  			 aes(x = as.numeric(month))) +
  	geom_point(aes_string(y = y_col_harv_bart, color = "siteID"), show.legend = F,
  						 alpha=.1, position=position_jitter(height=0)) +
  	geom_point(aes(y = as.numeric(truth), color = siteID), alpha = .5,
  						 position=position_jitter(height=0)) +
  	geom_smooth(aes_string(y = y_col_harv_bart), se = F, color=1) +
  	theme_bw( base_size = 22) +
  	scale_fill_brewer(palette = "Paired") +
  	theme(text = element_text(size = 18), panel.spacing = unit(.2, "cm"),
  				legend.position = "bottom",legend.title = element_text(NULL),
  				plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
  	xlab(NULL) + labs(fill='') +
  	facet_wrap(fg_category~taxon, scales="free") + ggtitle("Full models, estimates plotted by month ")
} else {
  cat("No data available for all_decomposers_env HARV/BART plot\n")
}
