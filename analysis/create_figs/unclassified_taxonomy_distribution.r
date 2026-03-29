# Distribution of unclassified reads across taxonomy (proxy for phylogeny).
# Tests whether unclassified ASVs are spread across many clades or concentrated in few.
# Single script: one table + one figure.
# Requires phyloseq objects in data/clean (or paths set via options).

source("source.R")

library(dplyr)
library(tidyr)
library(ggplot2)
library(here)

# Optional: set paths if phyloseq files live elsewhere
# options(phyloseq_16S_path = "/path/to/phyloseq_16S.rds", phyloseq_ITS_path = "/path/to/phyloseq_ITS.rds")
path_16S <- getOption("phyloseq_16S_path", here("data", "clean", "phyloseq_16S.rds"))
path_ITS <- getOption("phyloseq_ITS_path", here("data", "clean", "phyloseq_ITS.rds"))

if (!file.exists(path_16S) && !file.exists(path_ITS)) {
  stop(
    "No phyloseq RDS files found. Place phyloseq_16S.rds and/or phyloseq_ITS.rds in data/clean/ ",
    "or set options(phyloseq_16S_path = ..., phyloseq_ITS_path = ...)."
  )
}
if (!file.exists(path_16S)) message("phyloseq_16S.rds not found; skipping bacteria.")
if (!file.exists(path_ITS)) message("phyloseq_ITS.rds not found; skipping fungi.")

# Taxonomy column names (lowercase as in prep scripts)
tax_cols_16S <- c("kingdom", "phylum", "class", "order", "family", "genus")
tax_cols_ITS <- c("kingdom", "phylum", "class", "order", "family")

# Unclassified = all functional-group columns are "other" (no FG assignment)

label_unclassified <- function(tax_df, fg_cols) {
  if (length(fg_cols) == 0L) return(tax_df)
  other_mat <- tax_df[, fg_cols, drop = FALSE] == "other"
  tax_df$unclassified <- rowSums(other_mat, na.rm = TRUE) == ncol(other_mat)
  tax_df
}

summarize_by_rank <- function(tax_df, rank_col, kingdom_label) {
  rank_col <- as.name(rank_col)
  tax_df %>%
    mutate(rank_val = as.character(!!rank_col)) %>%
    filter(
      !is.na(.data$rank_val),
      .data$rank_val != "",
      !grepl("unassigned|uncultured|unknown", .data$rank_val, ignore.case = TRUE)
    ) %>%
    group_by(.data$rank_val) %>%
    summarise(
      n_unclassified = sum(.data$unclassified, na.rm = TRUE),
      n_classified = sum(!.data$unclassified, na.rm = TRUE),
      n_total = dplyr::n(),
      .groups = "drop"
    ) %>%
    mutate(
      pct_unclassified = 100 * .data$n_unclassified / .data$n_total,
      kingdom = kingdom_label
    )
}

process_phyloseq <- function(path, tax_cols, kingdom_label) {
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package phyloseq is required. Install with BiocManager::install('phyloseq').")
  }
  ps <- readRDS(path)
  tax_df <- as.data.frame(as.matrix(phyloseq::tax_table(ps)), stringsAsFactors = FALSE)
  colnames(tax_df) <- tolower(colnames(tax_df))
  fg_cols <- setdiff(colnames(tax_df), tax_cols)
  fg_cols <- intersect(fg_cols, colnames(tax_df)[sapply(tax_df, function(x) any(x == "other", na.rm = TRUE))])
  tax_df <- label_unclassified(tax_df, fg_cols)
  list(tax_df = tax_df, tax_cols = tax_cols, kingdom_label = kingdom_label, fg_cols = fg_cols)
}

# Load and process available kingdoms
proc_16S <- if (file.exists(path_16S)) process_phyloseq(path_16S, tax_cols_16S, "Bacteria") else NULL
proc_ITS <- if (file.exists(path_ITS)) process_phyloseq(path_ITS, tax_cols_ITS, "Fungi") else NULL

# By phylum (main test: broad vs concentrated)
by_phylum_16S <- if (!is.null(proc_16S)) summarize_by_rank(proc_16S$tax_df, "phylum", "Bacteria") else NULL
by_phylum_ITS <- if (!is.null(proc_ITS)) summarize_by_rank(proc_ITS$tax_df, "phylum", "Fungi") else NULL
by_phylum <- bind_rows(by_phylum_16S, by_phylum_ITS) %>%
  rename(phylum = rank_val)

# By class (finer resolution)
by_class_16S <- if (!is.null(proc_16S)) summarize_by_rank(proc_16S$tax_df, "class", "Bacteria") else NULL
by_class_ITS <- if (!is.null(proc_ITS)) summarize_by_rank(proc_ITS$tax_df, "class", "Fungi") else NULL
by_class <- bind_rows(by_class_16S, by_class_ITS) %>%
  rename(class = rank_val)

# Table: phylum-level summary (and save)
tab_phylum <- by_phylum %>%
  select(kingdom, phylum, n_unclassified, n_classified, n_total, pct_unclassified) %>%
  arrange(kingdom, desc(n_unclassified))

out_dir <- here("figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_table_path <- here("figures", "unclassified_by_taxonomy.csv")
readr::write_csv(tab_phylum, out_table_path)
cat("Table saved:", out_table_path, "\n")

# Summary stats
for (proc in list(proc_16S, proc_ITS)) {
  if (is.null(proc)) next
  n_unclass <- sum(proc$tax_df$unclassified, na.rm = TRUE)
  n_total <- nrow(proc$tax_df)
  by_phy <- by_phylum %>% filter(kingdom == proc$kingdom_label)
  n_phyla_with <- sum(by_phy$n_unclassified > 0)
  n_phyla_total <- nrow(by_phy)
  cat(proc$kingdom_label, ": unclassified ASVs in ", n_phyla_with, " of ", n_phyla_total, " phyla (",
      round(100 * n_unclass / n_total, 1), "%, ", n_unclass, "/", n_total, " ASVs)\n", sep = "")
}

# Figure: unclassified vs classified counts by phylum (shows spread across clades)
plot_df <- by_phylum %>%
  pivot_longer(
    cols = c(n_unclassified, n_classified),
    names_to = "classification",
    values_to = "n_ASVs"
  ) %>%
  mutate(classification = ifelse(classification == "n_unclassified", "Unclassified", "Classified"))

p <- ggplot(plot_df, aes(x = reorder(phylum, n_total), y = n_ASVs, fill = classification)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c(Classified = "gray70", Unclassified = "steelblue")) +
  facet_wrap(~kingdom, scales = "free_y", ncol = 1) +
  coord_flip() +
  labs(
    x = "Phylum",
    y = "Number of ASVs",
    fill = NULL,
    title = "Unclassified ASVs distributed across phyla (taxonomy as phylogeny proxy)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

out_fig_path <- here("figures", "unclassified_taxonomy_distribution.png")
ggsave(out_fig_path, p, width = 8, height = 8, dpi = 200)
cat("Figure saved:", out_fig_path, "\n")

# Optional: save class-level table for supplementary
readr::write_csv(
  by_class %>% select(kingdom, class, n_unclassified, n_classified, n_total, pct_unclassified) %>% arrange(kingdom, desc(n_unclassified)),
  here("figures", "unclassified_by_class.csv")
)
