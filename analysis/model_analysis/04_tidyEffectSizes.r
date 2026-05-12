# Combine & separately save model parameter and effect size estimates (beta covariates) from all models

source("../../source.R")
pacman::p_load(stringr, forestplot, gridExtra, dplyr)

# Source helper functions if package not available
if (!requireNamespace("microbialForecast", quietly = TRUE)) {
  if (file.exists(here("microbialForecast/R/statsFunctions.r"))) {
    source(here("microbialForecast/R/statsFunctions.r"))
  }
  if (file.exists(here("microbialForecast/R/helperFunctions.r"))) {
    source(here("microbialForecast/R/helperFunctions.r"))
  }
}

# Read in summaries and combine into fewer dfs for parameter effects

# Use pre-combined summaries from 03b_combine_summaries.r (much faster than reading individual files)
combined_summary_file <- here("data", "summary/logit_beta_fixed_priors_summaries.rds")

if (file.exists(combined_summary_file)) {
  cat("Loading pre-combined summaries from", basename(combined_summary_file), "...\n")
  sum.in <- readRDS(combined_summary_file)
  
  if (!is.list(sum.in) || !"summary_df" %in% names(sum.in)) {
    stop("Combined summary file has incorrect structure. Please regenerate using 03b_combine_summaries.r")
  }
  
  cat("Loaded", nrow(sum.in$summary_df), "rows from combined summary file\n")
  if (!is.null(sum.in$plot_est) && nrow(sum.in$plot_est) > 0) {
    cat("Plot estimates available:", nrow(sum.in$plot_est), "rows\n")
  }
} else {
  cat("Pre-combined summary file not found. Please run 03b_combine_summaries.r first.\n")
  cat("Expected file:", combined_summary_file, "\n")
  stop("Combined summary file not found. Run 03b_combine_summaries.r to create it.")
}

# Plot estimates are already in sum.in from the combined file
if (is.null(sum.in$plot_est) || nrow(sum.in$plot_est) == 0) {
  cat("No plot estimates in combined file. Proceeding without plot estimates.\n")
  sum.in$plot_est <- data.frame()
} else {
  cat("Plot estimates available from combined file:", nrow(sum.in$plot_est), "rows\n")
}

sum.all <- sum.in$summary_df  %>% filter(model_name != "all_covariates") %>% 
	mutate(tax_rank = rank,
																				 time_period = recode(time_period, !!!microbialForecast:::date_recode))

# Check if we have cloglog driver uncertainty models
has_cloglog_driver_uncertainty <- any(grepl("20130601_20180101_with_legacy_covariate_beta_regression", sum.all$model_id))

# Handle cloglog driver uncertainty models
if (has_cloglog_driver_uncertainty) {
  # Add a flag to identify driver uncertainty models
  # Model IDs have pattern: model_name_species_20130601_20180101_with_legacy_covariate_beta_regression
  sum.all <- sum.all %>%
    mutate(
      is_driver_uncertainty = grepl("20130601_20180101_with_legacy_covariate_beta_regression", model_id),
      # Ensure model_name is correctly set (should already be cycl_only, env_cov, or env_cycl)
           model_name = case_when(
        grepl("^cycl_only_", model_id) ~ "cycl_only",
        grepl("^env_cov_", model_id) ~ "env_cov",
        grepl("^env_cycl_", model_id) ~ "env_cycl",
        TRUE ~ model_name  # Keep original if pattern doesn't match
      )
    )
  cat("Detected cloglog driver uncertainty models\n")
  cat("  Driver uncertainty models:", sum(sum.all$is_driver_uncertainty), "\n")
  cat("  Model types:", paste(unique(sum.all$model_name[sum.all$is_driver_uncertainty]), collapse=", "), "\n")
} else {
  # Add flag for regular models
  sum.all <- sum.all %>%
    mutate(is_driver_uncertainty = FALSE)
  cat("Processing regular (non-driver uncertainty) models\n")
}
df <- sum.all %>%
	mutate(pretty_group = ifelse(group %in% c("16S","bac"), "Bacteria", "Fungi"))

# Add prettier data values
# Define pretty_rank_names locally if package not available
if (!exists("pretty_rank_names") || is.null(pretty_rank_names)) {
  if (requireNamespace("microbialForecast", quietly = TRUE)) {
    pretty_rank_names <- microbialForecast:::pretty_rank_names
  } else {
    # Fallback: define a simple pretty_rank_names
    pretty_rank_names <- c(
      "genus" = "Genus",
      "family" = "Family", 
      "order" = "Order",
      "class" = "Class",
      "phylum" = "Phylum",
      "functional" = "Functional group"
    )
  }
}

df$pretty_name <- recode(df$rank_only, !!!pretty_rank_names) %>%
	ordered(levels = c("Genus","Family","Order","Class","Phylum","Functional group","Diversity"))

df$only_rank <- sapply(str_split(df$rank_only, "_",  n = 2), `[`, 1) %>%
	ordered(levels = c("genus","family","order","class","phylum","functional","diversity"))

# df$tax_rank <- ordered(df$tax_rank, levels = c("genus_bac","genus_fun",
# 																														 "family_bac","family_fun",
# 																														 "order_bac", "order_fun",
# 																														 "class_bac", "class_fun",
# 																														 "phylum_bac","phylum_fun",
# 																														 "functional_group", "diversity_16S", "diversity_ITS"))


# For saving: filter by effect type

# Linear model beta (covariate) effects
beta_effects <- df %>% filter(grepl("beta", rowname))
# Replace "Ectomycorrhizal trees" with "Ectomycorrhizal\ntrees" for plotting
beta_effects$beta[beta_effects$beta == "Ectomycorrhizal trees"] <- "Ectomycorrhizal\ntrees"
# Make beta column an ordered factor for proper plotting order
beta_effects$beta <- ordered(beta_effects$beta, levels = c("Ectomycorrhizal\ntrees", "LAI", "pC", "pH", "Temperature", "Moisture", "sin", "cos"))
# Make rank_only column an ordered factor for proper plotting order
beta_effects$rank_only <- ordered(beta_effects$rank_only, levels = c("genus", "family", "order", "class", "phylum", "functional"))
saveRDS(beta_effects, here("data", "summary/predictor_effects.rds"))

# Site effects
site_effects_raw <- df %>% filter(grepl("site", rowname))

cat("Site effects before aggregation: ", nrow(site_effects_raw), " rows\n")

# Group by essential columns only - site effects should be unique per model/taxon/site/time_period
# Aggregate numeric columns, use [1] for other columns to avoid nested structure issues
site_effects <- site_effects_raw %>%
  group_by(model_id, model_name, taxon, siteID, time_period) %>%
  summarize(
    Mean = mean(Mean, na.rm = TRUE),
    SD = mean(SD, na.rm = TRUE),
    Median = mean(Median, na.rm = TRUE),
    # Use [1] instead of first() to avoid nested structure issues
    significant = if("significant" %in% names(site_effects_raw)) significant[1] else NA,
    effSize = if("effSize" %in% names(site_effects_raw)) effSize[1] else NA,
    site_num = if("site_num" %in% names(site_effects_raw)) site_num[1] else NA,
    param = if("param" %in% names(site_effects_raw)) param[1] else NA,
    rowname = if("rowname" %in% names(site_effects_raw)) rowname[1] else NA,
    beta_num = if("beta_num" %in% names(site_effects_raw)) beta_num[1] else NA,
    beta = if("beta" %in% names(site_effects_raw)) beta[1] else NA,
    group = if("group" %in% names(site_effects_raw)) group[1] else NA,
    fcast_type = if("fcast_type" %in% names(site_effects_raw)) fcast_type[1] else NA,
    pretty_group = if("pretty_group" %in% names(site_effects_raw)) pretty_group[1] else NA,
    rank = if("rank" %in% names(site_effects_raw)) rank[1] else NA,
    rank_only = if("rank_only" %in% names(site_effects_raw)) rank_only[1] else NA,
    pretty_name = if("pretty_name" %in% names(site_effects_raw)) pretty_name[1] else NA,
    only_rank = if("only_rank" %in% names(site_effects_raw)) only_rank[1] else NA,
    is_driver_uncertainty = if("is_driver_uncertainty" %in% names(site_effects_raw)) is_driver_uncertainty[1] else NA,
    tax_rank = if("tax_rank" %in% names(site_effects_raw)) tax_rank[1] else NA,
    .groups = "drop"
  )

cat("Site effects after aggregation: ", nrow(site_effects), " rows\n")

saveRDS(site_effects, here("data", "summary/site_effects.rds"))

# Linear model beta (covariate) effects
rho_effects <- df %>% filter(grepl("rho", rowname) | grepl("precision", rowname))
saveRDS(rho_effects, here("data", "summary/rho_core_sd_effects.rds"))

# Seasonality (cos/sin) effects
seas_params <- df %>% filter(beta %in% c("sin","cos") & !grepl("other", taxon))
if (nrow(seas_params) > 0) {
  seas_vals <- seas_params %>% pivot_wider(id_cols = c("model_id","taxon","model_name","fcast_type","time_period",
																												#"pretty_name","
																										 "pretty_group","rank","rank_only"),
																						names_from = beta,
																						values_from = c("Mean","significant"),
																						values_fn = mean) %>% rename(sin="Mean_sin", cos = "Mean_cos")
} else {
  cat("No seasonal parameters found - creating empty seasonal data\n")
  seas_vals <- data.frame()  # Empty dataframe for seasonal values
}

if (nrow(seas_vals) > 0) {
  # Convert to amplitude and max
  # Handle list columns by taking the first value and ensure single values
  out <- list()
  for (i in 1:nrow(seas_vals)) {
    sin_val <- if(is.list(seas_vals$sin[[i]])) seas_vals$sin[[i]][[1]] else seas_vals$sin[[i]]
    cos_val <- if(is.list(seas_vals$cos[[i]])) seas_vals$cos[[i]][[1]] else seas_vals$cos[[i]]
    
    # Ensure single values (not vectors)
    if(length(sin_val) > 1) sin_val <- sin_val[1]
    if(length(cos_val) > 1) cos_val <- cos_val[1]
    
    # Handle NA values
    if(is.na(sin_val)) sin_val <- 0
    if(is.na(cos_val)) cos_val <- 0
    
    out[[i]] <- sin_cos_to_seasonality(sin_val, cos_val)
  }
  out <- rbindlist(out)
  seas_vals <- cbind.data.frame(seas_vals, out)

  # Primary calibration period in current production data is 2013-06_2018-01 (recoded form).
  # "env_cycl" is the model carrying BOTH sin/cos AND env covariates, so its sin/cos are
  # the "residual seasonal trend" after env covariates are accounted for.
  # "env_cov" has NO sin/cos terms and is intentionally excluded here.
  primary_period <- "2013-06_2018-01"
  refit_period   <- "2015-11_2020-01"

  cycl_vals    <- seas_vals %>% filter(time_period == primary_period & model_name == "cycl_only")
  all_cov_vals <- seas_vals %>% filter(time_period == primary_period & model_name == "env_cycl")

  cycl_vals_refit    <- seas_vals %>% filter(time_period == refit_period & model_name == "cycl_only")
  all_cov_vals_refit <- seas_vals %>% filter(time_period == refit_period & model_name == "env_cycl")

  input_dateID = c("201401","201402","201403","201404","201405","201406","201407","201408","201409","201410","201411","201412")
  dates = fixDate(input_dateID)
  input_date_df = data.frame(x = lubridate::month(dates),
                             dates = dates)
  out_seas_vals =list()

  seas_vals_only = seas_vals %>% filter(grepl("cycl",model_name))
  for (row in 1:nrow(seas_vals_only)){
    alpha = seas_vals_only[row,]$sin
    beta = seas_vals_only[row,]$cos
    
    # Handle list columns by taking the first value
    if(is.list(alpha)) alpha <- alpha[[1]]
    if(is.list(beta)) beta <- beta[[1]]
    
    # Ensure single values
    if(length(alpha) > 1) alpha <- alpha[1]
    if(length(beta) > 1) beta <- beta[1]
    
    # Handle NA values
    if(is.na(alpha)) alpha <- 0
    if(is.na(beta)) beta <- 0
    
    df = input_date_df
    y_cycl = alpha * sin(2*pi*input_date_df$x/12) + beta * cos(2*pi*input_date_df$x/12)
    names(y_cycl) = input_date_df$dates
    out_seas_vals[[row]] <- data.frame(y_cycl) %>% t %>% as.data.frame()
  }
  seas_vals_to_plot = rbindlist(out_seas_vals, fill=T)
  date_cols <- names(seas_vals_to_plot)  # column names = ISO date strings produced by fixDate()
  seas_vals_to_plot <- cbind.data.frame(seas_vals_only, seas_vals_to_plot)
  seas_vals_long = seas_vals_to_plot %>%
    pivot_longer(cols = all_of(date_cols), names_to = "dates", values_to = "y_cycl") %>%
    mutate(dates = as.Date(dates))

  max_vals =	seas_vals_long %>% group_by(model_name, taxon, time_period, model_id) %>%
    filter(y_cycl == max(y_cycl, na.rm=T)) %>%
    mutate(max_y_date = dates) %>% select(-c(dates, y_cycl))
  seas_vals_long <- merge(seas_vals_long, max_vals, all=T)

  # Saved as a named list to make slot meanings explicit. Positional access still works,
  # so existing consumers (fig3_f_b_seasonality.r uses [[6]], etc.) continue to function.
  # Slot 2 is env_cycl (sin/cos AFTER env covariates) — i.e. the "residual seasonal trend".
  # Slot 3 is cycl_only (pure seasonal model). env_cov has no sin/cos and is not present.
  saveRDS(list(
    seas_vals_long      = seas_vals_long,
    env_cycl_vals       = all_cov_vals,
    cycl_only_vals      = cycl_vals,
    env_cycl_vals_refit = all_cov_vals_refit,
    cycl_only_vals_refit = cycl_vals_refit,
    seas_vals           = seas_vals
  ), here("data/summary/seasonal_amplitude.rds"))
} else {
  cat("No seasonal data to process - creating empty seasonal amplitude file\n")
  # Create empty seasonal amplitude data for compatibility
  saveRDS(list(data.frame(), data.frame(), data.frame(), data.frame(), data.frame(), data.frame()), 
          here("data/summary/seasonal_amplitude.rds"))
}

#
# ggplot(rho_effects %>% filter(model_name=="all_covariates"),
# 			 aes(x = pretty_group,y = abs(Mean), color=pretty_group)) +
# 	geom_boxplot() +
# 	geom_jitter( size = 5, height = 0, width=.4, alpha = .3,
# 							 shape = 16, show.legend = F) +
# 	ylab(NULL) + stat_compare_means()
#
#
#
# 	select_plots <- c("HARV_033","OSBS_026","WOOD_044","KONZ_001")
#
# 	examples = ggplot(sum.in$plot_est %>% filter(taxon %in% c("basidiomycota") & plotID %in% select_plots),
# 																							 aes(fill=species, x = dates, group=plotID)) +
# 		#rows = vars(fcast_type), drop=T, scales="free") +
# 		geom_ribbon(
# 								aes(x = dates,ymin = `2.5%`, ymax = `97.5%`), alpha=0.2) +
# 		geom_ribbon(
# 								aes(x = dates,ymin = `25%`, ymax = `75%`), alpha=.5) +
# 		geom_line(aes(x = dates, y = `50%`), alpha=0.3) +
# 		geom_point(aes(y = as.numeric(truth)), position = position_jitter()) +
# 		xlab(NULL) + labs(fill='') +
# 		scale_fill_brewer(palette = "Set2") +
# 		scale_color_brewer(palette = "Set2") +
# 		theme(panel.spacing = unit(.2, "cm"),
# 					legend.position = "bottom",legend.title = element_text(NULL),
# 					plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
# 		ggtitle("Hindcasts at 4 plots") +
# 		theme_minimal(base_size = 20) +
# 		scale_y_sqrt() +
# 		theme(legend.position = "none") + facet_grid(plotID~model_name, scales="free")
# 	examples
#
# 	seasonal_cycl = ggplot(df, aes(x=dates, y=y_cycl)) +
# 		geom_line(colour = "red") +
# 		theme_bw(base_size = 20) +
# 		ggtitle(paste0("Seasonal trend in abundance: ", taxon_name)) + ylab("Cyclic component") +
# 		#ylim(c(-.2,.2))
# 		ylim(c(-.1,.1))
#
#
# 	# Check the predicted~observed time series
# ggplot(sum.in$plot_est) +
# geom_point(aes( x = Mean, y = as.numeric(truth),
# 										color = siteID), show.legend = F, alpha=.5)  +
# 		ylab("Observed") + xlab("Predicted") +
# # 		theme_minimal(base_size = 18) +
# 		ggtitle("Overestimating plot means") +
# 		geom_abline(slope=1, intercept = 0) + facet_wrap(~taxon)

