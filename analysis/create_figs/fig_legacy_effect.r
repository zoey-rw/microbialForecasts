# Compare legacy_effect parameter across ranks and domains
# This script extracts and visualizes the legacy_effect parameter from model summaries

source("source.R")
pacman::p_load(stringr, forestplot, gridExtra, ggpubr, rstatix)

# Read in the main summary data
sum.in <- readRDS(here("data", "summary/logit_beta_fixed_priors_summaries.rds"))
sum.all <- sum.in$summary_df %>% filter(model_name != "all_covariates") %>% 
  mutate(tax_rank = rank,
         time_period = recode(time_period, !!!microbialForecast:::date_recode))

# Use existing pretty_group column (already has correct Bacteria/Fungi values)
df <- sum.all

# Add prettier data values - handle legacy_effect parameters specially
df$pretty_name <- case_when(
  df$rank_only == "genus" ~ "Genus",
  df$rank_only == "family" ~ "Family", 
  df$rank_only == "order" ~ "Order",
  df$rank_only == "class" ~ "Class",
  df$rank_only == "phylum" ~ "Phylum",
  df$rank_only == "functional" ~ "Functional group",
  TRUE ~ as.character(df$rank_only)
) %>%
  ordered(levels = c("Genus","Family","Order","Class","Phylum","Functional group","Diversity"))

df$only_rank <- df$rank_only %>%
  ordered(levels = c("genus","family","order","class","phylum","functional","diversity"))

# Filter for legacy_effect parameters
legacy_effects <- df %>% 
  filter(grepl("legacy_effect", rowname)) %>%
  mutate(effSize = abs(Mean))  # Absolute effect size
  # Note: significance is already calculated in the data

# Remove extreme outlier that dominates the plot
# assim_nitrite_reduction has an extreme legacy effect of -26.19
legacy_effects <- legacy_effects %>%
  filter(!(taxon == "assim_nitrite_reduction" & abs(Mean) > 20))

cat("Removed extreme outlier (assim_nitrite_reduction) from analysis\n")
cat("Remaining legacy effect parameters:", nrow(legacy_effects), "\n")

# Check if we have any legacy effect data
if (nrow(legacy_effects) == 0) {
  cat("No legacy_effect parameters found in the data.\n")
  cat("Available rownames containing 'legacy':\n")
  print(unique(df$rowname[grepl("legacy", df$rowname, ignore.case = TRUE)]))
  stop("No legacy_effect data available for plotting")
}

cat("Found", nrow(legacy_effects), "legacy_effect parameter estimates\n")

# Create the main comparison plot
legacy_plot <- ggplot(legacy_effects, 
                      aes(x = pretty_name, y = Mean, color = pretty_group)) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.7) +
  geom_point(aes(shape = as.factor(significant)), 
             size = 4, alpha = 0.7, 
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2)) +
  geom_violin(aes(fill = pretty_group), alpha = 0.3, 
              draw_quantiles = c(0.5), show.legend = FALSE) +
  facet_grid(rows = vars(pretty_group), scales = "free_y") +
  labs(title = "Legacy Effect Parameter Estimates",
       subtitle = "Effect of legacy sampling period (2013-2015) on microbial abundance",
       x = "Taxonomic Rank",
       y = "Legacy Effect Estimate",
       color = "Domain",
       shape = "Significant") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14)) +
  scale_shape_manual(values = c(21, 16), 
                     labels = c("Not significant", "Significant")) +
  guides(color = guide_legend(override.aes = list(shape = 16)),
         shape = guide_legend(override.aes = list(color = "black")))

# Add simple statistical comparisons if we have enough data
if (nrow(legacy_effects) > 10) {
  # Test for differences between domains (simplified)
  tryCatch({
    domain_test <- t.test(Mean ~ pretty_group, data = legacy_effects)
    cat("Domain comparison p-value:", domain_test$p.value, "\n")
  }, error = function(e) {
    cat("Could not perform domain comparison:", e$message, "\n")
  })
}

# Create a summary table
legacy_summary <- legacy_effects %>%
  group_by(pretty_group, pretty_name) %>%
  summarise(
    n_models = n(),
    mean_effect = mean(Mean, na.rm = TRUE),
    sd_effect = sd(Mean, na.rm = TRUE),
    median_effect = median(Mean, na.rm = TRUE),
    prop_significant = mean(significant, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(pretty_group, pretty_name)

cat("\nLegacy Effect Summary by Rank and Domain:\n")
print(legacy_summary)

# Create a second plot showing absolute effect sizes
legacy_abs_plot <- ggplot(legacy_effects, 
                          aes(x = pretty_name, y = effSize, color = pretty_group)) +
  geom_point(aes(shape = as.factor(significant)), 
             size = 4, alpha = 0.7,
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2)) +
  geom_violin(aes(fill = pretty_group), alpha = 0.3, 
              draw_quantiles = c(0.5), show.legend = FALSE) +
  facet_grid(rows = vars(pretty_group), scales = "free_y") +
  labs(title = "Legacy Effect Magnitude",
       subtitle = "Absolute magnitude of legacy sampling period effects",
       x = "Taxonomic Rank",
       y = "Absolute Legacy Effect Size",
       color = "Domain",
       shape = "Significant") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14)) +
  scale_shape_manual(values = c(21, 16), 
                     labels = c("Not significant", "Significant")) +
  # scale_y_log10() +  # Removed log scale to avoid warnings
  guides(color = guide_legend(override.aes = list(shape = 16)),
         shape = guide_legend(override.aes = list(color = "black")))

# Create a combined plot
combined_plot <- ggarrange(legacy_plot, legacy_abs_plot, 
                          ncol = 2, nrow = 1,
                          common.legend = TRUE, legend = "bottom",
                          labels = c("A", "B"))

# Save the plots
ggsave(here("figures", "legacy_effect_comparison.png"), 
       legacy_plot, width = 10, height = 8, dpi = 300)

ggsave(here("figures", "legacy_effect_magnitude.png"), 
       legacy_abs_plot, width = 10, height = 8, dpi = 300)

ggsave(here("figures", "legacy_effect_combined.png"), 
       combined_plot, width = 16, height = 8, dpi = 300)

# Save summary data
write.csv(legacy_summary, here("figures", "legacy_effect_summary.csv"), row.names = FALSE)

cat("\nPlots saved to figures/ directory:\n")
cat("- legacy_effect_comparison.png\n")
cat("- legacy_effect_magnitude.png\n") 
cat("- legacy_effect_combined.png\n")
cat("- legacy_effect_summary.csv\n")

# Identify taxa driving patterns
cat("\n=== TAXA DRIVING PATTERNS ===\n")

# Top 10 largest absolute legacy effects
top_effects <- legacy_effects %>%
  arrange(desc(effSize)) %>%
  head(10) %>%
  select(taxon, pretty_group, pretty_name, Mean, effSize)

cat("Top 10 largest absolute legacy effects:\n")
print(top_effects)

# Summary by domain and rank
cat("\nLegacy effect summary by domain and rank (excluding outlier):\n")
legacy_summary_clean <- legacy_effects %>%
  group_by(pretty_group, pretty_name) %>%
  summarise(
    n_models = n(),
    mean_effect = mean(Mean, na.rm = TRUE),
    sd_effect = sd(Mean, na.rm = TRUE),
    median_effect = median(Mean, na.rm = TRUE),
    max_abs_effect = max(effSize, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(pretty_group, desc(max_abs_effect))

print(legacy_summary_clean)

# Print key findings
cat("\n=== KEY FINDINGS ===\n")
cat("Total legacy effect parameters analyzed:", nrow(legacy_effects), "\n")
cat("Bacteria models:", sum(legacy_effects$pretty_group == "Bacteria"), "\n")
cat("Fungi models:", sum(legacy_effects$pretty_group == "Fungi"), "\n")
cat("Mean legacy effect (Bacteria):", round(mean(legacy_effects$Mean[legacy_effects$pretty_group == "Bacteria"], na.rm = TRUE), 4), "\n")
cat("Mean legacy effect (Fungi):", round(mean(legacy_effects$Mean[legacy_effects$pretty_group == "Fungi"], na.rm = TRUE), 4), "\n")
cat("Domain comparison p-value:", round(domain_test$p.value, 4), "\n")
cat("========================\n")
