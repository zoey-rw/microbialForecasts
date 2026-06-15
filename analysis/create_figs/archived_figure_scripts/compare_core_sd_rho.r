source("source.R")
library(ggallin)
library(rstatix)
library(ggpubr)

# Load parameter estimates from driver uncertainty models
# These are summarized from cloglog_beta_driver_uncertainty models in 04_tidyEffectSizes.r
rho_core_in <- readRDS(here("data", "summary/rho_core_sd_effects.rds")) %>%
	filter(model_name != "all_covariates") %>%
	select(-pretty_name) %>%
	# Remove _beta_regression suffix to match abundance data format
	mutate(model_id = gsub("_beta_regression$", "", model_id))

# Check available model types
cat("Available model types:", paste(unique(rho_core_in$model_name), collapse=", "), "\n")
cat("Total rows in rho_core_in:", nrow(rho_core_in), "\n")

# Filter to driver uncertainty models (2013-2018 period with legacy covariate)
# Model IDs should match pattern: model_name_species_20130601_20180101_with_legacy_covariate
if("model_id" %in% colnames(rho_core_in)) {
	driver_uncertainty_pattern <- "20130601_20180101_with_legacy_covariate"
	n_driver_uncertainty <- sum(grepl(driver_uncertainty_pattern, rho_core_in$model_id))
	cat("Driver uncertainty models (2013-2018):", n_driver_uncertainty, "\n")
	
	# Optionally filter to only driver uncertainty models if other models exist
	if(n_driver_uncertainty > 0 && n_driver_uncertainty < nrow(rho_core_in)) {
		rho_core_in <- rho_core_in %>%
			filter(grepl(driver_uncertainty_pattern, model_id))
		cat("Filtered to driver uncertainty models only\n")
	}
}

# Load abundance data for normalization
in_list <- readRDS(here("data/summary/fcast_horizon_input.rds"))
fcast_horizon_null_site <-  in_list[[3]]

# Merge abundance data - ensure model_id format matches
# The abundance data may also have _beta_regression suffix that needs to be removed
if("model_id" %in% colnames(fcast_horizon_null_site)) {
	fcast_horizon_null_site$model_id <- gsub("_beta_regression$", "", fcast_horizon_null_site$model_id)
}

# fcast_horizon_null_site has one row per (model_id, siteID); merge would duplicate rho rows.
# Aggregate to one row per model_id so each model's rho appears once.
abundance_by_model <- fcast_horizon_null_site %>%
	group_by(model_id) %>%
	summarise(abundance = mean(abundance, na.rm = TRUE), .groups = "drop")
rho_core_in <- merge(rho_core_in, abundance_by_model, by = "model_id", all.x = TRUE)

cat("After merging abundance data:", nrow(rho_core_in), "rows\n")
cat("Models with abundance data:", sum(!is.na(rho_core_in$abundance)), "\n")

# Separate precision and rho parameters
precision_data = rho_core_in %>%
	filter(rowname == "precision") %>% 
	mutate(adj_sd = ifelse(is.na(abundance) | abundance == 0, Mean, Mean/abundance))

rho_data = rho_core_in %>%
	filter(rowname == "rho")

# Add confidence intervals and significance for rho
rho_data$hi <- rho_data$Mean + rho_data$SD*1.96
rho_data$lo <- rho_data$Mean - rho_data$SD*1.96
rho_data$significant <- microbialForecast:::is_significant(rho_data$lo, rho_data$hi)
rho_data$effSize <- abs(rho_data$Mean)

# Check available model types for rho data
cat("\n=== MODEL TYPE SUMMARY ===\n")
cat("Available model types in rho_data:\n")
if(nrow(rho_data) > 0) {
	model_counts <- rho_data %>% 
		group_by(model_name, fcast_type) %>%
		summarise(n = n(), .groups = "drop")
	print(model_counts)
} else {
	cat("No rho data available\n")
}

# Create comparison plots
# 1. Rho parameter comparison - combined kingdoms (env_cycl model from driver uncertainty)
if(nrow(rho_data %>% filter(model_name == "env_cycl")) > 0) {
	rho_plot_data <- rho_data %>% filter(model_name == "env_cycl")
	
	# Statistical tests (only if sufficient sample sizes) - no kingdom grouping
	rho_stats <- data.frame()
	group_counts <- rho_plot_data %>%
		group_by(fcast_type) %>%
		summarise(n = n(), .groups = "drop")
	
	# Check if all groups have at least 2 observations
	sufficient_data <- all(group_counts$n >= 2)
	
	if(sufficient_data) {
		rho_stats <- rho_plot_data %>%
			t_test(effSize ~ fcast_type) %>%
			add_significance() %>%
			add_xy_position(x = "fcast_type", dodge = 0.8)
	} else {
		cat("Insufficient data for statistical tests - some groups have < 2 observations\n")
	}
	
	# Create plot with violin plots - no kingdom grouping
	p1 <- ggplot(rho_plot_data,
				 aes(x = fcast_type, y = effSize)) +
		geom_violin(alpha=0.7, trim=FALSE, fill="lightblue", draw_quantiles = c(.5)) +
		geom_point(size=2, alpha=0.3, position=position_jitter(width = 0.1)) +
		xlab("Forecast Type") +
		theme_minimal(base_size = 14) +
		ylab("Temporal Stability (ρ)") +
		ggtitle("A) Temporal Stability Parameter Estimates") +
		scale_y_continuous(trans=ssqrt_trans, limits = c(0, 0.7)) +
		scale_x_discrete(labels = c("functional" = "Functional Groups", "taxon" = "Taxonomic Groups")) +
		theme(
			axis.text.x = element_text(angle = 45, vjust=1, hjust = 1, size = 12),
			axis.title = element_text(size = 14, face = "bold"),
			plot.title = element_text(size = 16, face = "bold", hjust = 0)
		)
	
	# Add statistical annotations
	if(nrow(rho_stats) > 0) {
		p1 <- p1 + stat_pvalue_manual(rho_stats, 
									  label = "{p.signif}", 
									  bracket.nudge.y = 0.1,
									  tip.length = 0.02,
									  size = 4)
	}
	
	# Display the plot
	print(p1)
	
	# Save the plot
	ggsave(here("figures", "figS5_rho_temporal_memory.png"), p1, width = 4, height = 6, dpi = 300)
	cat("Saved: figures/figS5_rho_temporal_memory.png\n")

	# Rho parameter comparison by kingdom (Bacteria vs Fungi), with separate statistical tests per kingdom
	rho_plot_data_kingdom <- rho_plot_data %>% filter(!is.na(pretty_group))
	if (nrow(rho_plot_data_kingdom) > 0) {
		rho_stats_by_kingdom_list <- lapply(split(rho_plot_data_kingdom, rho_plot_data_kingdom$pretty_group), function(d) {
			group_counts_kg <- d %>% group_by(fcast_type) %>% summarise(n = n(), .groups = "drop")
			if (any(group_counts_kg$n < 2)) return(NULL)
			out <- d %>%
				t_test(effSize ~ fcast_type) %>%
				add_significance() %>%
				add_xy_position(x = "fcast_type", dodge = 0.8)
			out$pretty_group <- d$pretty_group[1]
			out
		})
		rho_stats_by_kingdom <- bind_rows(rho_stats_by_kingdom_list[!sapply(rho_stats_by_kingdom_list, is.null)])
		# Set bracket y position per facet so it sits above the data in each panel
		if (nrow(rho_stats_by_kingdom) > 0) {
			ymax_by_kingdom <- rho_plot_data_kingdom %>%
				group_by(pretty_group) %>%
				summarise(y_max = max(effSize, na.rm = TRUE), .groups = "drop")
			rho_stats_by_kingdom <- rho_stats_by_kingdom %>%
				left_join(ymax_by_kingdom, by = "pretty_group") %>%
				mutate(y.position = ifelse(is.na(y_max), 0.5, y_max * 1.08))
		}
		p1_kingdom <- ggplot(rho_plot_data_kingdom, aes(x = fcast_type, y = effSize)) +
			geom_violin(alpha = 0.7, trim = FALSE, fill = "lightblue", draw_quantiles = c(.5)) +
			geom_point(size = 2, alpha = 0.3, position = position_jitter(width = 0.1)) +
			facet_wrap(~ pretty_group, ncol = 2) +
			xlab("Forecast Type") +
			theme_minimal(base_size = 14) +
			ylab("Temporal Stability (ρ)") +
			ggtitle("A) Temporal Stability Parameter Estimates by Kingdom") +
			scale_y_continuous(trans = ssqrt_trans, limits = c(0, 0.7)) +
			scale_x_discrete(labels = c(
				"functional" = "Functional Groups", "taxon" = "Taxonomic Groups",
				"Functional" = "Functional Groups", "Taxonomic" = "Taxonomic Groups")) +
			theme(
				axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
				axis.title = element_text(size = 14, face = "bold"),
				plot.title = element_text(size = 16, face = "bold", hjust = 0),
				strip.text = element_text(face = "bold", size = 12)
			)
		if (nrow(rho_stats_by_kingdom) > 0) {
			p1_kingdom <- p1_kingdom + stat_pvalue_manual(rho_stats_by_kingdom,
				label = "{p.signif}",
				bracket.nudge.y = 0.02,
				tip.length = 0.02,
				size = 4)
		}
		print(p1_kingdom)
		ggsave(here("figures", "rho_parameter_comparison_by_kingdom.png"), p1_kingdom, width = 7, height = 5, dpi = 300)
		cat("Rho parameter plot by kingdom saved to figures/rho_parameter_comparison_by_kingdom.png\n")
	}
} else {
	print("No rho data available for env_cycl model")
}

# 2. Precision parameter comparison - combined kingdoms
if(nrow(precision_data %>% filter(model_name == "env_cycl")) > 0) {
	precision_plot_data <- precision_data %>% filter(model_name == "env_cycl")
	
	# Statistical tests (only if sufficient sample sizes) - no kingdom grouping
	precision_stats <- data.frame()
	group_counts_precision <- precision_plot_data %>%
		group_by(fcast_type) %>%
		summarise(n = n(), .groups = "drop")
	
	# Check if all groups have at least 2 observations
	sufficient_data_precision <- all(group_counts_precision$n >= 2)
	
	if(sufficient_data_precision) {
		precision_stats <- precision_plot_data %>%
			t_test(adj_sd ~ fcast_type) %>%
			add_significance() %>%
			add_xy_position(x = "fcast_type", dodge = 0.8)
	} else {
		cat("Insufficient data for precision statistical tests - some groups have < 2 observations\n")
	}
	
	# Create plot with violin plots - no kingdom grouping
	p2 <- ggplot(precision_plot_data,
				 aes(x = fcast_type, y = adj_sd)) +
		geom_violin(alpha=0.7, trim=FALSE, fill="lightcoral", draw_quantiles = c(.5)) +
		geom_point(size=2, alpha=0.3, position=position_jitter(width = 0.1)) +
		theme_minimal(base_size = 14) +
		theme(
			axis.text.x = element_text(angle = 45, vjust=1, hjust = 1, size = 12),
			axis.title = element_text(size = 14, face = "bold"),
			plot.title = element_text(size = 16, face = "bold", hjust = 0)
		) +
		xlab("Forecast Type") +
		ylab("Precision / Mean Abundance") +
		ggtitle("B) Precision Parameter Estimates") +
		scale_y_continuous(trans = "log10", limits = c(10, 50000)) +
		scale_x_discrete(labels = c("functional" = "Functional Groups", "taxon" = "Taxonomic Groups"))
	
	# Add statistical annotations
	if(nrow(precision_stats) > 0) {
		p2 <- p2 + annotate("text", x = 1.5, y = 40000, 
							label = precision_stats$p.signif, 
							size = 4, hjust = 0.5, vjust = 0.5) +
					annotate("segment", x = 1, xend = 2, y = 35000, yend = 35000, 
							linewidth = 0.5) +
					annotate("segment", x = 1, xend = 1, y = 35000, yend = 30000, 
							linewidth = 0.5) +
					annotate("segment", x = 2, xend = 2, y = 35000, yend = 30000, 
							linewidth = 0.5)
	}
	
	# Display the plot
	print(p2)
	
	# Save the plot
	ggsave(here("figures", "precision_parameter_comparison.png"), p2, width = 6, height = 6, dpi = 300)
	cat("Precision parameter plot saved to figures/precision_parameter_comparison.png\n")
} else {
	print("No precision data available for env_cycl model")
}

# 3. Summary statistics and statistical tests
cat("\n=== SUMMARY STATISTICS ===")
cat("\nComparing parameters from driver uncertainty models (env_cycl)\n")
cat("Models: 2013-2018 period with legacy covariate\n")

# Rho parameter summary - combined kingdoms
if(nrow(rho_data %>% filter(model_name == "env_cycl")) > 0) {
	cat("\nA) Temporal Stability (Rho) Parameter Estimates (Combined Kingdoms):\n")
	rho_summary <- rho_data %>% 
		filter(model_name == "env_cycl") %>%
		group_by(fcast_type) %>%
		summarise(
			n = n(),
			mean_rho = round(mean(Mean, na.rm = TRUE), 4),
			sd_rho = round(sd(Mean, na.rm = TRUE), 4),
			mean_effSize = round(mean(effSize, na.rm = TRUE), 4),
			prop_significant = round(mean(significant, na.rm = TRUE), 3),
			.groups = "drop"
		)
	print(rho_summary)
	
	# Statistical test results for rho - combined kingdoms
	cat("\nStatistical Tests for Temporal Stability (Rho) - Combined Kingdoms:\n")
	group_counts_rho <- rho_data %>% 
		filter(model_name == "env_cycl") %>%
		group_by(fcast_type) %>%
		summarise(n = n(), .groups = "drop")
	
	if(all(group_counts_rho$n >= 2)) {
		rho_test_results <- rho_data %>%
			filter(model_name == "env_cycl") %>%
			t_test(effSize ~ fcast_type) %>%
			add_significance()
		print(rho_test_results)
	} else {
		cat("Insufficient data for statistical tests - some groups have < 2 observations\n")
		print(group_counts_rho)
	}

	# Rho by kingdom: summary and separate statistical tests per kingdom
	rho_data_kingdom <- rho_data %>% filter(model_name == "env_cycl", !is.na(pretty_group))
	if (nrow(rho_data_kingdom) > 0) {
		cat("\nA2) Temporal Stability (Rho) by Kingdom:\n")
		rho_summary_kingdom <- rho_data_kingdom %>%
			group_by(pretty_group, fcast_type) %>%
			summarise(
				n = n(),
				mean_rho = round(mean(Mean, na.rm = TRUE), 4),
				sd_rho = round(sd(Mean, na.rm = TRUE), 4),
				mean_effSize = round(mean(effSize, na.rm = TRUE), 4),
				prop_significant = round(mean(significant, na.rm = TRUE), 3),
				.groups = "drop"
			)
		print(rho_summary_kingdom)
		cat("\nStatistical Tests for Temporal Stability (Rho) - By Kingdom:\n")
		rho_test_by_kingdom <- rho_data_kingdom %>%
			group_by(pretty_group) %>%
			group_modify(function(d, key) {
				group_n <- d %>% group_by(fcast_type) %>% summarise(n = n(), .groups = "drop")
				if (any(group_n$n < 2)) return(tibble())
				d %>% t_test(effSize ~ fcast_type) %>% add_significance()
			}, .keep = TRUE)
		if (nrow(rho_test_by_kingdom) > 0) {
			print(rho_test_by_kingdom)
		} else {
			cat("Insufficient data in one or both kingdoms for t-test (need >= 2 per fcast_type)\n")
		}
	}
}

# Precision parameter summary - combined kingdoms
if(nrow(precision_data %>% filter(model_name == "env_cycl")) > 0) {
	cat("\nB) Precision Parameter Estimates (Combined Kingdoms):\n")
	precision_summary <- precision_data %>% 
		filter(model_name == "env_cycl") %>%
		group_by(fcast_type) %>%
		summarise(
			n = n(),
			mean_precision = round(mean(Mean, na.rm = TRUE), 2),
			sd_precision = round(sd(Mean, na.rm = TRUE), 2),
			mean_adj_sd = round(mean(adj_sd, na.rm = TRUE), 2),
			.groups = "drop"
		)
	print(precision_summary)
	
	# Statistical test results for precision - combined kingdoms
	cat("\nStatistical Tests for Precision - Combined Kingdoms:\n")
	group_counts_precision_summary <- precision_data %>% 
		filter(model_name == "env_cycl") %>%
		group_by(fcast_type) %>%
		summarise(n = n(), .groups = "drop")
	
	if(all(group_counts_precision_summary$n >= 2)) {
		precision_test_results <- precision_data %>% 
			filter(model_name == "env_cycl") %>%
			t_test(adj_sd ~ fcast_type) %>%
			add_significance()
		print(precision_test_results)
	} else {
		cat("Insufficient data for statistical tests - some groups have < 2 observations\n")
		print(group_counts_precision_summary)
	}
}

cat("\n=== COMPARISON COMPLETE ===\n")