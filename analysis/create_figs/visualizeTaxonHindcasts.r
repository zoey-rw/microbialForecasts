
# View forecasts for individual taxonomic groups
source("source.R")

pacman::p_load(Rfast, moments, ggpubr)

# Check if hindcast data exists
if (!file.exists(here("data/summary/all_hindcasts.rds"))) {
  stop("all_hindcasts.rds not found. Please regenerate this file from the model analysis pipeline.")
}

hindcasts_raw <- readRDS(here("data/summary/all_hindcasts.rds"))
hindcasts <- hindcasts_raw %>% mutate(taxon = species)

# Check if scoring metrics exist
if (!file.exists(here("data/summary/scoring_metrics_cv.rds"))) {
  stop("scoring_metrics_cv.rds not found. Please regenerate this file from the model analysis pipeline.")
}

scoring_metrics_cv <- readRDS(here("data/summary/scoring_metrics_cv.rds"))
phy_scores <- scoring_metrics_cv$scoring_metrics %>% 
  filter(model_name == "all_covariates" &
         pretty_name %in% c("Phylum","Functional group"))
genus_scores <- scoring_metrics_cv$scoring_metrics %>% 
  filter(model_name == "all_covariates" &
         pretty_name %in% c("Genus"))
fg_scores <- scoring_metrics_cv$scoring_metrics %>% 
  filter(model_name == "all_covariates" &
         pretty_name %in% c("Functional group"))

calibration_only = hindcasts %>% filter(fcast_period=="calibration") %>% filter(timepoint > plot_start_date)
asco_calibration <- calibration_only %>% filter(taxon == "ascomycota") %>% filter(!is.na(truth) & !is.na(mean))
asco_calibration_wide <- asco_calibration %>% filter(!siteID %in% c("DEJU","LENO")) %>%
  select(siteID, plotID, dateID, model_name, taxon, mean, truth) %>%
  pivot_wider(names_from = model_name, values_from = c(mean))

# For the calibration, the first observed date per plot is always wonky due to model structure, so we leave it out of our scoring metrics

if (nrow(phy_scores) > 0) {
  taxon_scores = phy_scores %>%
    group_by(taxon, model_name, site_prediction) %>%
    summarize(mean_RSQ = mean(RSQ)) %>%
    pivot_wider(names_from = site_prediction, values_from ="mean_RSQ")
} else {
  cat("No phy_scores data available for taxon_scores analysis.\n")
  taxon_scores <- data.frame()
}

worse_at_new_sites = c("acidobacteriota", "penicillium")
fine_at_new_sites = c("ascomycota", "copiotroph")

# Check if plot_model function exists
if (exists("plot_model")) {
  plot_model(hindcasts, taxon = "acidobacteriota", siteID = "BART")
} else {
  cat("plot_model function not available.\n")
}

# Get plot information
if (nrow(hindcasts) > 0) {
  not_na <- hindcasts[!is.na(hindcasts$truth),]
  top_plots <- names(tail(sort(table(not_na$plotID)), 8))
  hindcast_in <- hindcasts[hindcasts$plotID %in% top_plots,]
  
  new_plots = table(hindcasts[hindcasts$new_site==TRUE & !is.na(hindcasts$truth) & hindcasts$pretty_group !="Fungi",]$plotID) %>% sort %>% head(30)
  obs_plots = table(hindcasts[hindcasts$fcast_period == "calibration" & !is.na(hindcasts$truth) & hindcasts$pretty_group !="Fungi",]$plotID) %>% sort %>% head(30)
  new_plots_fun = table(hindcasts[hindcasts$new_site==TRUE & !is.na(hindcasts$truth) & hindcasts$pretty_group=="Fungi",]$plotID) %>% sort %>% head(30)
  obs_plots_fun = table(hindcasts[hindcasts$fcast_period == "calibration" & hindcasts$pretty_group=="Fungi" & !is.na(hindcasts$truth),]$plotID) %>% sort %>% head(30)
  
  intersect(names(obs_plots),names(obs_plots_fun))
  intersect(names(new_plots),names(new_plots_fun))
  plots_to_viz = c("HARV_034","CPER_001","DSNY_001", "BONA_004", "SOAP_031")
  
  input_df <- hindcasts %>% mutate(fcast_period_alpha = ifelse(fcast_period == "calibration", .8, .4)) %>%
    filter(model_name =="all_covariates" &
           taxon %in% c(worse_at_new_sites,fine_at_new_sites) &
           plotID %in% plots_to_viz) %>%
    filter(!grepl("random",site_prediction))
  
  acido <- hindcasts %>% mutate(fcast_period_alpha = ifelse(fcast_period == "calibration", .8, .4)) %>%
    filter(model_name =="all_covariates" &
           taxon %in% c("acidobacteriota")) %>%
    filter(!grepl("random",site_prediction))
  
  penicillium <- hindcasts %>% mutate(fcast_period_alpha = ifelse(fcast_period == "calibration", .8, .4)) %>%
    filter(model_name =="all_covariates" &
           taxon %in% c("penicillium")) %>%
    filter(!grepl("random",site_prediction)) %>% 
    group_by(siteID) %>% 
    summarize(avg = mean(truth, na.rm=TRUE)) %>% 
    arrange(avg)
  
  input_df <- hindcasts %>% mutate(fcast_period_alpha = ifelse(fcast_period == "calibration", .8, .4)) %>%
    filter(model_name =="all_covariates" &
           taxon %in% fine_at_new_sites &
           plotID %in% c("HARV_001","YELL_001")) %>%
    filter(!grepl("random",site_prediction))
  
  # Create plot if data is available
  if (nrow(input_df) > 0) {
    ggplot(input_df, aes(x = dates, group=plotID))  +
      ylab("Abundance") +
      theme_bw() +
      theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
            legend.position = "bottom",legend.title = element_text(NULL),
            plot.margin = unit(c(.2, .2, 2, .2), "cm")) +
      facet_grid(newsite~taxon, scales="free", drop=TRUE) +
      geom_ribbon(aes(ymin = lo, ymax = hi, fill = newsite), alpha=0.3) +
      geom_ribbon(aes(ymin = lo_25, ymax = hi_75, fill = newsite), alpha=0.5) +
      geom_line(aes(y = med, color = newsite), show.legend = FALSE) +
      geom_jitter(aes(y = as.numeric(truth)), width = .2, height = 0) +
      xlab(NULL)
  } else {
    cat("No data available for plotting.\n")
  }
} else {
  cat("No hindcast data available for analysis.\n")
}

# Check if single taxon summaries exist
if (file.exists(here("data/summary/single_taxon_summaries_201511_201801.rds"))) {
  single_tax_summaries <- readRDS(here("data/summary/single_taxon_summaries_201511_201801.rds"))
  cat("Single taxon summaries loaded successfully.\n")
} else {
  cat("single_taxon_summaries_201511_201801.rds not found.\n")
}
