
source("source.R")

# Check if scoring metrics data exists
if (!file.exists(here("data/summary/scoring_metrics_cv.rds"))) {
  stop("scoring_metrics_cv.rds not found. Please regenerate this file from the model analysis pipeline.")
}

scores_list <- readRDS(here("data/summary/scoring_metrics_cv.rds"))

if ("scoring_metrics" %in% names(scores_list)) {
  # Get scores by rank
  scores_by_rank <- scores_list$scoring_metrics %>%
    filter(model_name == "all_covariates") %>%
    group_by(pretty_name) %>%
    summarize(
      mean_RSQ = mean(RSQ, na.rm = TRUE),
      mean_CRPS = mean(CRPS, na.rm = TRUE),
      mean_RMSE = mean(RMSE, na.rm = TRUE),
      n_taxa = n()
    )
  
  if (nrow(scores_by_rank) > 0) {
    # Create scores by rank plot
    p1 <- ggplot(scores_by_rank, aes(x = pretty_name, y = mean_RSQ)) +
      geom_col(aes(fill = pretty_name)) +
      theme_bw() +
      labs(title = "Mean R-squared by Taxonomic Rank",
           x = "Taxonomic Rank",
           y = "Mean R-squared",
           fill = "Rank") +
      theme(text = element_text(size = 14),
            axis.text.x = element_text(angle = 45, hjust = 1))
    
    print(p1)
    cat("Scores by rank plot created successfully\n")
    
    # Create CRPS by rank plot
    p2 <- ggplot(scores_by_rank, aes(x = pretty_name, y = mean_CRPS)) +
      geom_col(aes(fill = pretty_name)) +
      theme_bw() +
      labs(title = "Mean CRPS by Taxonomic Rank",
           x = "Taxonomic Rank",
           y = "Mean CRPS",
           fill = "Rank") +
      theme(text = element_text(size = 14),
            axis.text.x = element_text(angle = 45, hjust = 1))
    
    print(p2)
    cat("CRPS by rank plot created successfully\n")
  } else {
    stop("No scores data available. Please check data generation.")
  }
} else {
  stop("scoring_metrics element not found in scores_list. Please check data structure.")
}

# Check if tax filter data exists
if (file.exists(here("data/summary/tax_filter_pass.rds"))) {
  tax_filter_pass <- readRDS(here("data/summary/tax_filter_pass.rds"))
  
  if (nrow(tax_filter_pass) > 0) {
    # Create tax filter summary
    tax_filter_summary <- tax_filter_pass %>%
      group_by(rank) %>%
      summarize(
        n_passed = n(),
        mean_abundance = mean(mean_abundance, na.rm = TRUE)
      )
    
    if (nrow(tax_filter_summary) > 0) {
      # Create tax filter summary plot
      p3 <- ggplot(tax_filter_summary, aes(x = rank, y = n_passed)) +
        geom_col(aes(fill = rank)) +
        theme_bw() +
        labs(title = "Number of Taxa Passing Filter by Rank",
             x = "Taxonomic Rank",
             y = "Number of Taxa",
             fill = "Rank") +
        theme(text = element_text(size = 14),
              axis.text.x = element_text(angle = 45, hjust = 1))
      
      print(p3)
      cat("Tax filter summary plot created successfully\n")
    } else {
      cat("No tax filter summary data available for plotting\n")
    }
  } else {
    cat("tax_filter_pass data has 0 rows\n")
  }
} else {
  cat("tax_filter_pass.rds not found. Data may need to be regenerated.\n")
}

