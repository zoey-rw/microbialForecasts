
# Visualize phenological parameters from all models
source("source.R")
pacman::p_load(ggrepel, ggforce, gridExtra, ggpubr)


#######
# convert sin and cos to maximum and minimum of amplitude
# from Stolwijk 1999

# Agrees with results from Stolwijk 1999
# getMaxMin(sin = 0.009822, cos = 0.06929, max_only = F)

##### Create season data frame for plotting #####
mo <- 1:12
mo_sin <- sin((2*pi*mo)/12)
mo_cos <- cos((2*pi*mo)/12)
seasons <- cbind.data.frame(mo, mo_sin, mo_cos)
seasons$season <- "winter"
seasons$season[3:5] <- "spring"
seasons$season[6:8] <- "summer"
seasons$season[9:11] <- "fall"

if(exists("model_name") && model_name == "cycl_only") {
  cov_key <- microbialForecast:::cycl_only_key
}

#####

# Check if forecast effects data exists
if (!file.exists(here("data/summary/all_fcast_effects.rds"))) {
  stop("all_fcast_effects.rds not found. Please regenerate this file from the model analysis pipeline.")
}

df_orig <- readRDS(here("data/summary/all_fcast_effects.rds"))

# Filter for cyclical effects
df <- df_orig %>%
  filter(beta %in% c("cos", "sin") &
         time_period == "2015-11_2018-01") %>%
  group_by(taxon, time_period, model_name) %>%
  mutate(cyc_significant = ifelse(any(significant==1), 1, 0))

# Get cyclical values
vals <- df %>% filter(cyc_significant == 1) %>%
  pivot_wider(id_cols = c("taxon","model_name","fcast_type","pretty_group"),
              values_from = "Mean", names_from = "beta")

# Calculate amplitude and peak timing
if (nrow(vals) > 0) {
  vals$amplitude <- sqrt(vals$sin^2 + vals$cos^2)
  vals$max <- getMaxMin(vals$sin, vals$cos)
  vals$max <- unlist(vals$max)
  
  # Create predictability quadrants plot
  p1 <- ggplot(vals, aes(x = amplitude, y = max)) +
    geom_point(aes(color = pretty_group)) +
    geom_hline(yintercept = 6, linetype = "dashed") +
    geom_vline(xintercept = median(vals$amplitude, na.rm = TRUE), linetype = "dashed") +
    theme_bw() +
    labs(title = "Phylum Predictability Quadrants",
         x = "Seasonal Amplitude",
         y = "Peak Month",
         color = "Group") +
    scale_y_continuous(breaks = 1:12, labels = month.abb) +
    theme(text = element_text(size = 14))
  
  print(p1)
  cat("Phylum predictability quadrants plot created successfully\n")
} else {
  stop("No significant cyclical data available. Please check data generation.")
}

# Check if functional group summaries exist
if (file.exists(here("data/summary/fg_summaries.rds"))) {
  fg_in <- readRDS(here("data/summary/fg_summaries.rds"))
  fg_effects <- fg_in$summary_df %>%
    filter(beta %in% c("cos", "sin") &
           taxon != "other" &
           model_name == "cycl_only") %>%
    group_by(group, fg_cat, rank, taxon) %>%
    mutate(cyc_significant = ifelse(any(significant==1), 1, 0),
           rank_only = "functional")
  
  cat("Functional group phenology data loaded successfully\n")
} else {
  cat("fg_summaries.rds not found. Data may need to be regenerated.\n")
}

# Check if convergence summaries exist
if (file.exists(here("data/summary/taxon_convergence_summaries_201511_201801.rds"))) {
  convergence_in <- readRDS(here("data/summary/taxon_convergence_summaries_201511_201801.rds"))
  cat("Taxon convergence summaries loaded successfully\n")
} else {
  cat("taxon_convergence_summaries_201511_201801.rds not found. Data may need to be regenerated.\n")
}
