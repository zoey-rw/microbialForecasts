# Combine & separately save model parameter and effect size estimates (beta covariates) from all truncated normal models

source("../../source.R")
pacman::p_load(stringr, forestplot, gridExtra)

# Read in summaries and combine into fewer dfs for parameter effects

# Functional groups
sum.in <- readRDS(here("data", "summary/truncated_normal_summaries.rds"))
sum.all <- sum.in$summary_df  %>% filter(model_name != "all_covariates") %>% 
	mutate(tax_rank = rank,
																				 time_period = recode(time_period, !!!microbialForecast:::date_recode))
df <- sum.all %>%
	mutate(pretty_group = ifelse(group %in% c("16S","bac"), "Bacteria", "Fungi"))

# Add prettier data values
df$pretty_name <- recode(df$rank_only, !!!microbialForecast:::pretty_rank_names) %>%
	ordered(levels = c("Genus","Family","Order","Class","Phylum","Functional group","Diversity"))

df$only_rank <- sapply(str_split(df$rank_only, "_",  n = 2), `[`, 1) %>%
	ordered(levels = c("genus","family","order","class","phylum","functional","diversity"))

# For saving: filter by effect type

# Create parameter type column for truncNorm models
df$parameter_type <- ifelse(grepl("beta", df$rowname), "beta", 
                           ifelse(grepl("core_sd|sigma", df$rowname), "variance", 
                                 ifelse(grepl("intercept|legacy_effect|rho", df$rowname), "effect", "other")))

# Linear model beta (covariate) effects
beta_effects <- df %>% filter(grepl("beta", rowname))
# For truncNorm models, we need to extract the beta parameter number
beta_effects$beta_param <- sapply(strsplit(beta_effects$rowname, "\\[|\\]"), function(x) {
  if(length(x) > 1) paste0("beta_", x[2]) else "beta"
})
saveRDS(beta_effects, here("data", "summary/truncated_normal_predictor_effects.rds"))

# Site effects
site_effects <- df %>% filter(grepl("site", rowname))
saveRDS(site_effects, here("data", "summary/truncated_normal_site_effects.rds"))

# Linear model beta (covariate) effects - truncNorm specific parameters
rho_effects <- df %>% filter(grepl("rho", rowname) | grepl("core_sd", rowname) | grepl("sigma", rowname))
saveRDS(rho_effects, here("data", "summary/truncated_normal_rho_core_sd_effects.rds"))

# Seasonality (cos/sin) effects - for truncNorm models, we need to check if there are seasonal parameters
# Since truncNorm models don't have sin/cos parameters in the same way, we'll create empty seasonal data
cat("TruncNorm models don't have sin/cos seasonal parameters - creating empty seasonal data\n")
seas_params <- data.frame()  # Empty dataframe for seasonal parameters
seas_vals <- data.frame()  # Empty dataframe for seasonal values

if (nrow(seas_vals) > 0) {
  # Convert to amplitude and max
  # Didn't vectorize this function, oops
  out <- list()
  for (i in 1:nrow(seas_vals)) {
    out[[i]] <- sin_cos_to_seasonality(seas_vals$sin[[i]], seas_vals$cos[[i]])
  }
  out <- rbindlist(out)
  seas_vals <- cbind.data.frame(seas_vals, out)

  cycl_vals <- seas_vals %>% filter(time_period == "2015-11_2018-01" &
                                    model_name == "cycl_only")
  all_cov_vals <- seas_vals %>% filter(time_period == "2015-11_2018-01" &
                                       model_name == "all_covariates|env_cov")

  cycl_vals_refit <- seas_vals %>% filter(time_period == "2015-11_2020-01" &
                                          model_name == "cycl_only")
  all_cov_vals_refit <- seas_vals %>% filter(time_period == "2015-11_2020-01" &
                                             model_name == "all_covariates|env_cov")

  input_dateID = c("201401","201402","201403","201404","201405","201406","201407","201408","201409","201410","201412")
  dates = fixDate(input_dateID)
  input_date_df = data.frame(x = lubridate::month(dates),
                             dates = dates)
  out_seas_vals =list()

  seas_vals_only = seas_vals %>% filter(grepl("cycl",model_name))
  for (row in 1:nrow(seas_vals_only)){
    alpha = seas_vals_only[row,]$sin
    beta = seas_vals_only[row,]$cos
    df = input_date_df
    y_cycl = alpha * sin(2*pi*input_date_df$x/12) + beta * cos(2*pi*input_date_df$x/12)
    names(y_cycl) = input_date_df$dates
    out_seas_vals[[row]] <- data.frame(y_cycl) %>% t %>% as.data.frame()
  }

  seas_vals_only$y_cycl = out_seas_vals
  seas_vals_only$y_cycl = lapply(seas_vals_only$y_cycl, function(x) x %>% mutate(date = rownames(x)) %>% gather(key = "plot", value = "y_cycl", -date))

  saveRDS(seas_vals_only, here("data", "summary/truncated_normal_seasonal_effects.rds"))
} else {
  cat("No seasonal data to save\n")
}

# Save the main dataframe
saveRDS(df, here("data", "summary/truncated_normal_all_effects.rds"))

cat("Truncated normal effect sizes processing completed!\n")
cat("Files saved:\n")
cat("  - truncated_normal_predictor_effects.rds\n")
cat("  - truncated_normal_site_effects.rds\n")
cat("  - truncated_normal_rho_core_sd_effects.rds\n")
if (nrow(seas_params) > 0) cat("  - truncated_normal_seasonal_effects.rds\n")
cat("  - truncated_normal_all_effects.rds\n")
