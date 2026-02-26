source("source.R")

library(scales)
library(spgs) # for uniformity test

# run siteEffectVariogram.r to generate
variograms = readRDS(here("data/summary/site_effect_variograms.rds"))

# Use plsr2 converged list
scores_list <- readRDS(here("data/summary/scoring_metrics_plsr2.rds"))
converged_strict <- scores_list$converged_strict_list
cat("Loaded converged strict list with", length(converged_strict), "models\n")

# look at example
variograms[[2]][12]

sig_results = variograms[[1]]

# CRITICAL FIX: Fill in missing model_name values by extracting from model_id
# Some models have NA model_name but it can be extracted from model_id pattern
sig_results <- sig_results %>%
  mutate(
    model_name = ifelse(
      is.na(model_name),
      case_when(
        grepl("^env_cycl", model_id) ~ "env_cycl",
        grepl("^env_cov", model_id) ~ "env_cov",
        grepl("^cycl_only", model_id) ~ "cycl_only",
        TRUE ~ NA_character_
      ),
      model_name
    )
  )

sig_results_long = sig_results %>%
	pivot_longer(cols=5:6) %>%
	mutate(value = as.numeric(value))

# Check for overlap between converged models and variogram data
overlap <- intersect(converged_strict, sig_results$model_id)
cat("Overlap between converged models and variogram data:", length(overlap), "models\n")

if (length(overlap) > 0) {
  # Filter to converged models only
  sig_results_long <- sig_results_long %>% filter(model_id %in% converged_strict)
  cat("Filtered to", nrow(sig_results_long), "rows using converged models\n")
} else {
  cat("No overlap between converged models and variogram data.\n")
  cat("Variogram data uses different model ID format than current converged models.\n")
  cat("This suggests the variogram data is from an older analysis.\n")
  stop("Cannot proceed: no matching models between converged list and variogram data")
}

var_plot = ggplot(sig_results_long) +
	#geom_point(position = position_jitterdodge(jitter.width = .05, jitter.height = 0)) +
	geom_density(aes(x = value), show.legend = F) +
	facet_grid(rows=vars(name), scales="free") +
	geom_vline(xintercept = .05, color=2, linetype=2) +
	theme_classic(base_size = 20)+
	ylab("Variogram p-value") +
	xlab(NULL) +
	ggtitle("Site effect autocorrelation is mostly captured by PLSR") +
	#annotate("text", x = .1, y = 1.5, label = "p.value < .05 indicates \nspatial autocorrelation")  +
	geom_rug(aes(x = value), alpha=.3,length= unit(0.06, "npc"), show.legend = F) + xlab("variogram p-value")

var_plot <- tag_facet(var_plot)

# Neither set of values follows a uniform distribution
# CRITICAL FIX: Add error handling for insufficient data
site_effect_vals <- sig_results_long[sig_results_long$name=="site effect",]$value
if(length(site_effect_vals) >= 2 && sum(!is.na(site_effect_vals)) >= 2) {
  tryCatch({
    chisq.unif.test(site_effect_vals)
    # X-squared = 387.23, df = 2, a = 0, b = 1, p-value < 2.2e-16
  }, error = function(e) {
    cat("Warning: chisq.unif.test failed for site effect:", conditionMessage(e), "\n")
  })
} else {
  cat("Warning: Insufficient data for site effect uniformity test\n")
}

residual_vals <- sig_results_long[sig_results_long$name=="site effect residuals",]$value
if(length(residual_vals) >= 2 && sum(!is.na(residual_vals)) >= 2) {
  tryCatch({
    chisq.unif.test(residual_vals)
    # X-squared = 39.581, df = 17, a = 0, b = 1, p-value = 0.001482
  }, error = function(e) {
    cat("Warning: chisq.unif.test failed for residuals:", conditionMessage(e), "\n")
  })
} else {
  cat("Warning: Insufficient data for residual uniformity test\n")
}

# Non-uniform residuals are driven by env_cov model (without seasonality predictors)
# CRITICAL FIX: Add error handling for each model type
for(model_type in c("cycl_only", "env_cov", "env_cycl")) {
  model_vals <- sig_results_long %>% 
    filter(name=="site effect residuals" & model_name==model_type) %>%
    select(value) %>% unlist()
  
  if(length(model_vals) >= 2 && sum(!is.na(model_vals)) >= 2) {
    tryCatch({
      result <- chisq.unif.test(model_vals)
      cat("Model", model_type, "uniformity test completed\n")
    }, error = function(e) {
      cat("Warning: chisq.unif.test failed for", model_type, ":", conditionMessage(e), "\n")
    })
  } else {
    cat("Warning: Insufficient data for", model_type, "uniformity test\n")
  }
}
# Expected results (when data is sufficient):
# cycl_only: X-squared = 7.1852, df = 4, a = 0, b = 1, p-value = 0.1264
# env_cov: X-squared = 11.27, df = 4, a = 0, b = 1, p-value = 0.02369
# env_cycl: X-squared = 6.9259, df = 5, a = 0, b = 1, p-value = 0.2262

png(here("figures","variogram_p_val.png"), width = 800, height=1000)
print(var_plot)
dev.off()

# CRITICAL FIX: Properly handle variogram example plot (recordedplot object)
if(length(variograms) >= 2 && length(variograms[[2]]) >= 12) {
  example_plot <- variograms[[2]][[12]]
  if(inherits(example_plot, "recordedplot")) {
    png(here("figures","variogram_example.png"), width = 600, height=800)
    replayPlot(example_plot)
    dev.off()
    cat("Generated variogram_example.png\n")
  } else if(inherits(example_plot, "ggplot")) {
    png(here("figures","variogram_example.png"), width = 600, height=800)
    print(example_plot)
    dev.off()
    cat("Generated variogram_example.png\n")
  } else {
    cat("Warning: variograms[[2]][[12]] is not a recognized plot object, skipping example figure\n")
  }
} else {
  cat("Warning: variograms structure insufficient for example plot\n")
}


by_model = ggplot(sig_results_long) +
	#geom_point(position = position_jitterdodge(jitter.width = .05, jitter.height = 0)) +
	geom_density(aes(x = value), show.legend = F) +
	facet_grid(rows=vars(name), cols=vars(model_name), scales="free") +
	geom_vline(xintercept = .05, color=2, linetype=2) +
	theme_classic(base_size = 20)+
	ylab("Variogram p-value") +
	xlab(NULL) +
	ggtitle("Site effect autocorrelation is mostly captured by PLSR") +
	#annotate("text", x = .1, y = 1.5, label = "p.value < .05 indicates \nspatial autocorrelation")  +
	geom_rug(aes(x = value), alpha=.3,length= unit(0.06, "npc"), show.legend = F) + xlab("variogram p-value")

png(here("figures","variogram_p_val_by_model.png"), width = 1200, height=1000)
print(by_model)
dev.off()

# ============================================================================
# Two-panel summary figure: % significant + ECDF
# ============================================================================

library(ggpubr)

# Pretty labels
model_labs <- c(cycl_only = "Seasonal only",
                env_cov   = "Environment only",
                env_cycl  = "Environment + seasonal")
stage_labs <- c("site effect"           = "Raw site effects",
                "site effect residuals" = "After PLSR modeling")

# Ensure model_name is an ordered factor for consistent panel order
sig_results_long <- sig_results_long %>%
  mutate(
    model_name  = factor(model_name, levels = c("cycl_only", "env_cov", "env_cycl")),
    stage_label = factor(stage_labs[name], levels = stage_labs),
    significant = value < 0.05
  )

# ── Panel A: % of taxa with significant spatial autocorrelation ──────────────
pct_sig <- sig_results_long %>%
  group_by(model_name, stage_label) %>%
  summarise(
    n_total = n(),
    n_sig   = sum(significant, na.rm = TRUE),
    pct_sig = 100 * n_sig / n_total,
    .groups = "drop"
  )

pA_var <- ggplot(pct_sig,
                 aes(x = model_name, y = pct_sig, fill = stage_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_text(aes(label = paste0(round(pct_sig, 0), "%")),
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 4) +
  scale_x_discrete(labels = model_labs) +
  scale_fill_manual(values = c("Raw site effects"       = "#D55E00",
                                "After PLSR modeling"   = "#0072B2"),
                    name = NULL) +
  coord_cartesian(ylim = c(0, max(pct_sig$pct_sig) * 1.15)) +
  labs(x = NULL,
       y = "Taxa with significant\nspatial autocorrelation (%)",
       tag = "A") +
  theme_bw(base_size = 14) +
  theme(strip.background   = element_rect(fill = "grey92", color = NA),
        strip.text         = element_text(face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        legend.position    = "top",
        axis.text.x        = element_text(angle = 20, hjust = 1))

# ── Panel B: ECDF of p-values vs uniform reference ──────────────────────────
pB_var <- ggplot(sig_results_long,
                 aes(x = value, color = stage_label, linetype = stage_label)) +
  stat_ecdf(linewidth = 0.8, pad = FALSE) +
  geom_abline(slope = 1, intercept = 0, color = "grey50",
              linetype = "dotted", linewidth = 0.6) +
  facet_wrap(~model_name, labeller = labeller(model_name = model_labs)) +
  scale_color_manual(values = c("Raw site effects"     = "#D55E00",
                                "After PLSR modeling"  = "#0072B2"),
                     name = NULL) +
  scale_linetype_manual(values = c("Raw site effects"     = "solid",
                                   "After PLSR modeling"  = "solid"),
                        name = NULL) +
  annotate("text", x = 0.65, y = 0.15, label = "uniform\nreference",
           color = "grey50", size = 3.5, fontface = "italic") +
  labs(x = "Variogram p-value",
       y = "Cumulative proportion of taxa",
       tag = "B") +
  theme_bw(base_size = 14) +
  theme(strip.background   = element_rect(fill = "grey92", color = NA),
        strip.text         = element_text(face = "bold"),
        panel.grid.minor   = element_blank(),
        legend.position    = "top")

fig_variogram_summary <- ggarrange(pA_var, pB_var, nrow = 1, widths = c(1, 1.5),
                                   common.legend = TRUE, legend = "top")

ggsave(here("figures", "fig_variogram_summary.png"), fig_variogram_summary,
       width = 12, height = 5.5, dpi = 300)
cat("Saved: figures/fig_variogram_summary.png\n")

# Print summary table
cat("\n-- Spatial autocorrelation summary --\n")
print(as.data.frame(pct_sig))
