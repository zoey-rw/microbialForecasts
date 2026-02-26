# Explore drivers of site-level forecast skill
# Why do some sites have no skilled forecasts?
# Analyzes within-model variation across sites.

library(ggplot2)
library(ggpubr)
library(dplyr)
library(data.table)
library(ggrepel)
source("source.R")

# ============================================================
# 1. Compute site-level skill from hindcast scoring data
# ============================================================

fi <- readRDS(here("data/summary/fcast_horizon_input.rds"))
model_mean <- as.data.table(fi[[4]])  # per model x site x months_since_obs
null_site <- as.data.table(fi[[3]])   # null model per model x site

# Focus on env_cycl
scoring <- model_mean[model_name == "env_cycl"]
null_rsq <- null_site[model_name == "env_cycl",
                       .(null_RSQ = mean(RSQ.1, na.rm = TRUE)), by = model_id]
scoring <- merge(scoring, null_rsq, by = "model_id", all.x = TRUE)

# Classify each site as having hindcast data or not
site_has_hindcast <- scoring[months_since_obs > 0 & !is.na(RSQ.1),
                              .(has_hindcast = TRUE), by = siteID]
all_sites_in_data <- unique(scoring$siteID)

cat("Total sites in scoring data:", length(all_sites_in_data), "\n")
cat("Sites with hindcast evaluations:", nrow(site_has_hindcast), "\n")
cat("Sites without hindcast evaluations:",
    length(setdiff(all_sites_in_data, site_has_hindcast$siteID)), "\n")
cat("Sites without hindcast data:\n")
cat(paste(sort(setdiff(all_sites_in_data, site_has_hindcast$siteID)), collapse = ", "), "\n\n")

# ============================================================
# 2. Per-model, per-site horizon: max months beating null
# ============================================================

site_hz <- scoring[months_since_obs > 0 & RSQ.1 > null_RSQ & is.finite(RSQ.1),
                    .(site_horizon = max(months_since_obs, na.rm = TRUE)),
                    by = .(model_id, species, siteID, pretty_group)]

# Include model x site combos that never beat null
all_combos <- unique(scoring[siteID %in% site_has_hindcast$siteID,
                              .(model_id, species, siteID, pretty_group)])
site_hz <- merge(all_combos, site_hz,
                 by = c("model_id", "species", "siteID", "pretty_group"), all.x = TRUE)
site_hz[is.na(site_horizon), site_horizon := 0]

cat("Site-level horizon data:", nrow(site_hz), "model x site combos\n")

# ============================================================
# 3. Count calibration observations per site
# ============================================================

# Use all_truth_data for observation counts per site
truth_file <- here("data/clean/all_truth_data.rds")
if (file.exists(truth_file)) {
  truth_data <- as.data.table(readRDS(truth_file))

  # Count per site: number of timepoints and plot-level observations
  cal_counts <- truth_data[!is.na(truth),
                            .(n_obs = .N,
                              n_timepoints = uniqueN(dateID),
                              n_plots = uniqueN(plotID),
                              date_min = min(dates, na.rm = TRUE),
                              date_max = max(dates, na.rm = TRUE)),
                            by = siteID]
  cal_counts[, cal_span_months := as.numeric(difftime(date_max, date_min, units = "days")) / 30.44]

  cat("\nCalibration observation counts per site:\n")
  print(cal_counts[order(n_timepoints)])

  rm(truth_data); gc(verbose = FALSE)
} else {
  cat("all_truth_data.rds not found, trying calibration_only_processed.rds\n")
  cal_file <- here("data/summary/calibration_only_processed.rds")
  if (file.exists(cal_file)) {
    cal_data <- as.data.table(readRDS(cal_file))
    cal_counts <- cal_data[!is.na(truth) & !is.nan(truth),
                            .(n_obs = .N,
                              n_timepoints = uniqueN(dateID),
                              n_plots = uniqueN(plotID)),
                            by = siteID]
    cal_counts[, cal_span_months := NA_real_]
    rm(cal_data); gc(verbose = FALSE)
  } else {
    cat("No calibration data found. Creating placeholder.\n")
    cal_counts <- data.table(siteID = character(), n_obs = integer(),
                             n_timepoints = integer(), n_plots = integer(),
                             cal_span_months = numeric())
  }
}

# ============================================================
# 4. Load site-level environmental predictors
# ============================================================

site_env <- as.data.table(readRDS(here("data/clean/site_effect_predictors.rds")))
cat("\nSite environmental predictors:", nrow(site_env), "sites,", ncol(site_env), "variables\n")

# ============================================================
# 5. Aggregate site-level skill metrics
# ============================================================

# Per site: mean horizon, proportion skilled, mean RSQ
site_skill <- site_hz[, .(mean_horizon = mean(site_horizon, na.rm = TRUE),
                           med_horizon = median(site_horizon, na.rm = TRUE),
                           prop_skilled = mean(site_horizon > 0),
                           n_models = .N),
                        by = siteID]

# Add which sites have no hindcast data at all
no_hindcast_sites <- setdiff(all_sites_in_data, site_has_hindcast$siteID)
no_hindcast_dt <- data.table(siteID = no_hindcast_sites,
                              mean_horizon = NA_real_,
                              med_horizon = NA_real_,
                              prop_skilled = NA_real_,
                              n_models = 0L)
site_skill_all <- rbind(site_skill, no_hindcast_dt)
site_skill_all[, has_hindcast := siteID %in% site_has_hindcast$siteID]

# Merge everything
site_df <- merge(site_skill_all, cal_counts, by = "siteID", all.x = TRUE)
site_df <- merge(site_df, site_env[, .(siteID, latitude, MAT, MAP, bulkDens, estimatedOC)],
                 by = "siteID", all.x = TRUE)

cat("\nCombined site data:", nrow(site_df), "sites\n")
cat("With hindcast:", sum(site_df$has_hindcast), "\n")
cat("Without hindcast:", sum(!site_df$has_hindcast), "\n")

# ============================================================
# 6. What distinguishes sites with/without hindcast data?
# ============================================================

cat("\n=== SITES WITHOUT HINDCAST DATA ===\n")
cat("These sites have no post-calibration observations for evaluation.\n")
print(site_df[has_hindcast == FALSE,
              .(siteID, n_timepoints, n_plots, n_obs, cal_span_months, latitude)])

cat("\n=== SITES WITH HINDCAST DATA ===\n")
print(site_df[has_hindcast == TRUE,
              .(siteID, prop_skilled, mean_horizon, n_timepoints, n_plots, cal_span_months, latitude)][
                order(prop_skilled)])

# ============================================================
# Panel A: Calibration time points for sites with vs without hindcast data
# ============================================================

site_df$hindcast_status <- ifelse(site_df$has_hindcast,
                                   "Has hindcast evaluation",
                                   "No hindcast data")

panel_a <- ggplot(site_df %>% filter(!is.na(n_timepoints)),
                  aes(x = n_timepoints, y = hindcast_status)) +
  geom_point(aes(color = hindcast_status), size = 3, alpha = .6,
             position = position_jitter(height = .15, width = .3)) +
  geom_text_repel(aes(label = siteID), size = 2.5, max.overlaps = 20) +
  stat_summary(fun = median, geom = "crossbar", width = 0.4,
               color = "black", show.legend = FALSE) +
  theme_bw(base_size = 14) +
  xlab("Number of calibration time points") +
  ylab(NULL) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("Has hindcast evaluation" = "steelblue",
                                 "No hindcast data" = "grey50"))

# ============================================================
# Panel B: Calibration span (months) vs whether site has hindcast
# ============================================================

panel_b <- ggplot(site_df %>% filter(!is.na(cal_span_months)),
                  aes(x = cal_span_months, y = hindcast_status)) +
  geom_point(aes(color = hindcast_status), size = 3, alpha = .6,
             position = position_jitter(height = .15, width = .3)) +
  geom_text_repel(aes(label = siteID), size = 2.5, max.overlaps = 20) +
  stat_summary(fun = median, geom = "crossbar", width = 0.4,
               color = "black", show.legend = FALSE) +
  theme_bw(base_size = 14) +
  xlab("Calibration span (months)") +
  ylab(NULL) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("Has hindcast evaluation" = "steelblue",
                                 "No hindcast data" = "grey50"))

# ============================================================
# Panels C-F: Within hindcast-evaluable sites only
# ============================================================

site_eval <- site_df[has_hindcast == TRUE & !is.na(prop_skilled)]

# Panel C: Calibration timepoints vs proportion skilled
panel_c <- ggplot(site_eval,
                  aes(x = n_timepoints, y = prop_skilled)) +
  geom_point(size = 3, alpha = .7, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue", linewidth = 1) +
  geom_text_repel(aes(label = siteID), size = 2.5, max.overlaps = 20) +
  theme_bw(base_size = 14) +
  xlab("Number of calibration time points") +
  ylab("Proportion of taxa\nwith skilled forecasts")

tryCatch({
  ct <- cor.test(site_eval$n_timepoints, site_eval$prop_skilled)
  panel_c <- panel_c +
    labs(subtitle = paste0("r = ", round(ct$estimate, 2), ", p = ", round(ct$p.value, 3)))
}, error = function(e) NULL)

# Panel D: Latitude vs proportion skilled (for evaluable sites only)
panel_d <- ggplot(site_eval,
                  aes(x = latitude, y = prop_skilled)) +
  geom_point(size = 3, alpha = .7, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue", linewidth = 1) +
  geom_text_repel(aes(label = siteID), size = 2.5, max.overlaps = 20) +
  theme_bw(base_size = 14) +
  xlab("Latitude (°N)") +
  ylab("Proportion of taxa\nwith skilled forecasts")

tryCatch({
  ct2 <- cor.test(site_eval$latitude, site_eval$prop_skilled)
  panel_d <- panel_d +
    labs(subtitle = paste0("r = ", round(ct2$estimate, 2), ", p = ", round(ct2$p.value, 3)))
}, error = function(e) NULL)

# Panel E: MAT vs proportion skilled
panel_e <- ggplot(site_eval %>% filter(!is.na(MAT)),
                  aes(x = MAT, y = prop_skilled)) +
  geom_point(size = 3, alpha = .7, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue", linewidth = 1) +
  geom_text_repel(aes(label = siteID), size = 2.5, max.overlaps = 20) +
  theme_bw(base_size = 14) +
  xlab("Mean Annual Temperature (scaled)") +
  ylab("Proportion of taxa\nwith skilled forecasts")

tryCatch({
  ct3 <- with(site_eval %>% filter(!is.na(MAT)), cor.test(MAT, prop_skilled))
  panel_e <- panel_e +
    labs(subtitle = paste0("r = ", round(ct3$estimate, 2), ", p = ", round(ct3$p.value, 3)))
}, error = function(e) NULL)

# Panel F: Within-model site variation — distribution of per-model site SD
model_site_var <- site_hz[, .(site_sd = sd(site_horizon, na.rm = TRUE),
                               site_range = max(site_horizon) - min(site_horizon),
                               mean_hz = mean(site_horizon, na.rm = TRUE)),
                            by = .(model_id, species, pretty_group)]

panel_f <- ggplot(model_site_var %>% filter(!is.na(site_sd)),
                  aes(x = mean_hz, y = site_sd, color = pretty_group)) +
  geom_point(size = 2, alpha = .4) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  theme_bw(base_size = 14) +
  xlab("Mean forecast horizon across sites (months)") +
  ylab("SD of horizon across sites") +
  labs(color = "Kingdom") +
  theme(legend.position = "top")

# ============================================================
# Compose and save
# ============================================================

fig_top <- ggarrange(panel_a, panel_b, labels = c("A", "B"))
fig_mid <- ggarrange(panel_c, panel_d, labels = c("C", "D"))
fig_bot <- ggarrange(panel_e, panel_f, labels = c("E", "F"))
fig_all <- ggarrange(fig_top, fig_mid, fig_bot, nrow = 3)

png(here("figures", "site_skill_drivers.png"), width = 1400, height = 1500, res = 150)
print(fig_all)
dev.off()

cat("\nFigure saved to figures/site_skill_drivers.png\n")

# ============================================================
# Summary correlation table
# ============================================================

cat("\n=== CORRELATION TABLE (evaluable sites only) ===\n")
cor_vars <- c("prop_skilled", "mean_horizon", "n_timepoints", "n_plots",
              "cal_span_months", "latitude", "MAT", "MAP", "bulkDens", "estimatedOC")
available_vars <- cor_vars[cor_vars %in% names(site_eval)]
cor_data <- site_eval[, ..available_vars]
cor_data <- cor_data[complete.cases(cor_data)]

if (nrow(cor_data) > 5) {
  cor_mat <- cor(cor_data, use = "pairwise.complete.obs")
  cat("\nCorrelations with prop_skilled:\n")
  cors_with_skill <- cor_mat["prop_skilled", ]
  print(sort(cors_with_skill, decreasing = TRUE))
}
