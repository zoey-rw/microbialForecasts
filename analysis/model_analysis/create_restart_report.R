#!/usr/bin/env Rscript

# Simple script to create restart report
source("batch_restart_models.R")

cat("Creating restart report...\n")
restart_report <- create_restart_report()

if (!is.null(restart_report)) {
  write.csv(restart_report, "restart_candidates_cycl_only_legacy.csv", row.names = FALSE)
  cat("Created restart report for", nrow(restart_report), "models\n")
  cat("\nFirst few models:\n")
  print(head(restart_report[, c("model_name", "rank.name", "species", "scenario", "n_chains")]))

  cat("\nModel type summary:\n")
  print(table(restart_report$model_name))

  cat("\nRank summary:\n")
  print(table(restart_report$rank.name))

  cat("\nScenario summary:\n")
  print(table(restart_report$scenario))

} else {
  cat("No models found for restart\n")
}
