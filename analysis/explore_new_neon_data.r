## Explore new NEON microbe soil data (DP1.10081.002) and compare with existing model data
## Run from project root: Rscript analysis/explore_new_neon_data.r
##
## Data structure:
##   mct_soilSampleMetadata_16S/ITS CSVs: sample metadata (siteID, plotID, collectDate, dnaSampleID)
##   *_16S_map.csv / *_ITS_map.csv: per-ASV taxonomy + read counts per dnaSampleID

source("source.R")
library(dplyr)
library(tidyr)
library(ggplot2)

new_data_dir <- here("data", "new_neon_microbe")
zip_path <- file.path(new_data_dir, "NEON_tax-microbe-soil.zip")
extracted_dir <- file.path(new_data_dir, "NEON_tax-microbe-soil")

# --- 1. Extract metadata files if needed ---
# Only extract the small metadata CSVs (not the large ASV map files)
meta_files_16S <- list.files(extracted_dir, pattern = "mct_soilSampleMetadata_16S",
                              recursive = TRUE, full.names = TRUE)
if (length(meta_files_16S) == 0) {
  cat("Extracting metadata files from zip...\n")
  system2("unzip", c("-o", zip_path,
                     "NEON_tax-microbe-soil/*mct_soilSampleMetadata*",
                     "-d", new_data_dir), stdout = TRUE)
  meta_files_16S <- list.files(extracted_dir, pattern = "mct_soilSampleMetadata_16S",
                                recursive = TRUE, full.names = TRUE)
}

meta_files_ITS <- list.files(extracted_dir, pattern = "mct_soilSampleMetadata_ITS",
                              recursive = TRUE, full.names = TRUE)

cat("Found", length(meta_files_16S), "16S metadata files and",
    length(meta_files_ITS), "ITS metadata files\n")

# --- 2. Load and combine all metadata ---
read_and_bind <- function(files) {
  dfs <- lapply(files, function(f) {
    tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
  })
  bind_rows(dfs[!sapply(dfs, is.null)])
}

cat("\nReading 16S metadata...\n")
meta_16S <- read_and_bind(meta_files_16S)
cat("Reading ITS metadata...\n")
meta_ITS <- read_and_bind(meta_files_ITS)

# Parse dates and create dateID
meta_16S$collectDate <- as.Date(substr(meta_16S$collectDate, 1, 10))
meta_16S$dateID <- format(meta_16S$collectDate, "%Y%m")
meta_ITS$collectDate <- as.Date(substr(meta_ITS$collectDate, 1, 10))
meta_ITS$dateID <- format(meta_ITS$collectDate, "%Y%m")

# Filter to samples that passed QC
if ("sequenceCountQF" %in% names(meta_16S)) {
  n_before <- nrow(meta_16S)
  meta_16S_qc <- meta_16S %>% filter(sequenceCountQF == "OK" | is.na(sequenceCountQF))
  cat(sprintf("\n16S QC filter: %d -> %d samples (dropped %d)\n",
              n_before, nrow(meta_16S_qc), n_before - nrow(meta_16S_qc)))
} else {
  meta_16S_qc <- meta_16S
}
if ("sequenceCountQF" %in% names(meta_ITS)) {
  n_before <- nrow(meta_ITS)
  meta_ITS_qc <- meta_ITS %>% filter(sequenceCountQF == "OK" | is.na(sequenceCountQF))
  cat(sprintf("ITS QC filter: %d -> %d samples (dropped %d)\n",
              n_before, nrow(meta_ITS_qc), n_before - nrow(meta_ITS_qc)))
} else {
  meta_ITS_qc <- meta_ITS
}

# --- 3. Overview of new data ---
cat("\n========================================\n")
cat("=== NEW DATA OVERVIEW ===\n")
cat("========================================\n")

cat("\n--- 16S (bacteria/archaea) ---\n")
cat("Total samples:", nrow(meta_16S_qc), "\n")
cat("Sites:", length(unique(meta_16S_qc$siteID)), "\n")
cat("Date range:", as.character(range(meta_16S_qc$collectDate, na.rm = TRUE)), "\n")
cat("dateID range:", range(meta_16S_qc$dateID, na.rm = TRUE), "\n")

cat("\n--- ITS (fungi) ---\n")
cat("Total samples:", nrow(meta_ITS_qc), "\n")
cat("Sites:", length(unique(meta_ITS_qc$siteID)), "\n")
cat("Date range:", as.character(range(meta_ITS_qc$collectDate, na.rm = TRUE)), "\n")
cat("dateID range:", range(meta_ITS_qc$dateID, na.rm = TRUE), "\n")

# --- 4. Compare with existing data ---
cat("\n========================================\n")
cat("=== COMPARISON WITH EXISTING DATA ===\n")
cat("========================================\n")

existing_16S <- readRDS(here("data", "clean", "groupAbundances_16S_2023.rds"))
existing_ITS <- readRDS(here("data", "clean", "groupAbundances_ITS_2023.rds"))

existing_16S_ph <- existing_16S$phylum_bac
existing_ITS_ph <- existing_ITS[[1]]

existing_sites_16S <- sort(unique(existing_16S_ph$siteID))
existing_sites_ITS <- sort(unique(existing_ITS_ph$siteID))
existing_dates_16S <- sort(unique(existing_16S_ph$dateID))
existing_dates_ITS <- sort(unique(existing_ITS_ph$dateID))

cat("\nExisting 16S: ", length(existing_sites_16S), " sites, dateID ",
    min(existing_dates_16S), "-", max(existing_dates_16S),
    " (", nrow(existing_16S_ph), " obs)\n", sep = "")
cat("Existing ITS: ", length(existing_sites_ITS), " sites, dateID ",
    min(existing_dates_ITS), "-", max(existing_dates_ITS),
    " (", nrow(existing_ITS_ph), " obs)\n", sep = "")

new_sites_16S <- sort(unique(meta_16S_qc$siteID))
new_sites_ITS <- sort(unique(meta_ITS_qc$siteID))

overlap_16S <- intersect(new_sites_16S, existing_sites_16S)
overlap_ITS <- intersect(new_sites_ITS, existing_sites_ITS)
new_only_16S <- setdiff(new_sites_16S, existing_sites_16S)
new_only_ITS <- setdiff(new_sites_ITS, existing_sites_ITS)

cat("\n--- Site overlap ---\n")
cat("16S: ", length(overlap_16S), "/", length(new_sites_16S),
    " new sites overlap with existing (", length(new_only_16S), " new-only)\n", sep = "")
cat("ITS: ", length(overlap_ITS), "/", length(new_sites_ITS),
    " new sites overlap with existing (", length(new_only_ITS), " new-only)\n", sep = "")

if (length(new_only_16S) > 0) {
  cat("\nNew-only 16S sites:", paste(new_only_16S, collapse = ", "), "\n")
}

# --- 5. Temporal overlap analysis ---
cat("\n--- Temporal overlap ---\n")
new_dates_16S <- sort(unique(meta_16S_qc$dateID))
new_dates_ITS <- sort(unique(meta_ITS_qc$dateID))

overlap_dates_16S <- intersect(new_dates_16S, existing_dates_16S)
post_hindcast_16S <- new_dates_16S[new_dates_16S > max(existing_dates_16S)]

cat("\n16S date-months:\n")
cat("  Overlapping with existing:", length(overlap_dates_16S), "\n")
cat("  Post-hindcast (>", max(existing_dates_16S), "):", length(post_hindcast_16S), "\n")
if (length(post_hindcast_16S) > 0) {
  cat("  Post-hindcast range:", min(post_hindcast_16S), "-", max(post_hindcast_16S), "\n")
}

# --- 6. Detailed per-site overlap ---
cat("\n========================================\n")
cat("=== PER-SITE OVERLAP (16S) ===\n")
cat("========================================\n")
cat(sprintf("%-6s %6s %6s %6s %6s %s\n",
            "Site", "N_new", "Shared", "NewDt", "Post", "Post-hindcast dates"))

for (site in sort(overlap_16S)) {
  new_site_dates <- sort(unique(meta_16S_qc$dateID[meta_16S_qc$siteID == site]))
  old_site_dates <- sort(unique(existing_16S_ph$dateID[existing_16S_ph$siteID == site]))

  shared <- intersect(new_site_dates, old_site_dates)
  new_only <- setdiff(new_site_dates, old_site_dates)
  post <- new_only[new_only > max(existing_dates_16S)]
  n_new <- sum(meta_16S_qc$siteID == site)

  cat(sprintf("%-6s %6d %6d %6d %6d %s\n",
              site, n_new, length(shared), length(new_only), length(post),
              paste(head(post, 6), collapse = ", ")))
}

# --- 7. Taxonomic overlap (sample a few 16S map files) ---
cat("\n========================================\n")
cat("=== TAXONOMIC OVERLAP (16S sample) ===\n")
cat("========================================\n")

# Extract a handful of 16S map files for one site to check taxonomy
# Use CLBJ as the example
clbj_map_files <- list.files(extracted_dir, pattern = "16S_map\\.csv$",
                              recursive = TRUE, full.names = TRUE)

if (length(clbj_map_files) == 0) {
  # Extract CLBJ 16S map files
  cat("Extracting CLBJ 16S map files for taxonomy check...\n")
  system2("unzip", c("-o", zip_path,
                     "NEON_tax-microbe-soil/NEON.D11.CLBJ*/*16S_map*",
                     "-d", new_data_dir), stdout = TRUE)
  clbj_map_files <- list.files(extracted_dir, pattern = "CLBJ.*16S_map\\.csv$",
                                recursive = TRUE, full.names = TRUE)
}

# Also try to find any already-extracted map files
if (length(clbj_map_files) == 0) {
  # Extract just a few from any site
  cat("Extracting sample 16S map files...\n")
  system2("unzip", c("-o", zip_path,
                     "NEON_tax-microbe-soil/NEON.D16.ABBY*/*16S_map*",
                     "-d", new_data_dir), stdout = TRUE)
  clbj_map_files <- list.files(extracted_dir, pattern = "16S_map\\.csv$",
                                recursive = TRUE, full.names = TRUE)
}

if (length(clbj_map_files) > 0) {
  # Read a few map files to get taxonomy
  sample_files <- head(clbj_map_files, 5)
  cat("Reading", length(sample_files), "ASV map files...\n")

  asv_sample <- read_and_bind(sample_files)
  cat("ASV rows:", nrow(asv_sample), "\n")
  cat("Columns:", paste(names(asv_sample), collapse = ", "), "\n")

  # Phylum-level comparison
  if ("phylum" %in% names(asv_sample)) {
    new_phyla <- sort(unique(asv_sample$phylum))
    new_phyla <- new_phyla[new_phyla != "" & !is.na(new_phyla)]

    existing_phyla <- names(existing_16S_ph)[!names(existing_16S_ph) %in%
                                                c("siteID", "plotID", "dateID", "sampleID",
                                                  "dates", "plot_date", "other")]

    cat("\nNew data phyla (", length(new_phyla), "):\n", sep = "")
    cat("  ", paste(head(new_phyla, 30), collapse = "\n  "), "\n")
    cat("\nExisting model phyla (", length(existing_phyla), "):\n", sep = "")
    cat("  ", paste(existing_phyla, collapse = "\n  "), "\n")

    # Fuzzy match: normalize names for comparison
    normalize_taxon <- function(x) {
      x <- tolower(x)
      x <- gsub("[^a-z0-9]", "", x)
      x
    }
    new_norm <- normalize_taxon(new_phyla)
    old_norm <- normalize_taxon(existing_phyla)

    matched <- intersect(new_norm, old_norm)
    cat("\nMatched phyla (normalized):", length(matched), "/", length(old_norm), "\n")
    cat("  ", paste(matched, collapse = ", "), "\n")

    unmatched_old <- old_norm[!old_norm %in% new_norm]
    if (length(unmatched_old) > 0) {
      cat("Existing phyla NOT found in new data:", paste(unmatched_old, collapse = ", "), "\n")
    }
  }

  # Check relative abundance computation feasibility
  if ("individualCount" %in% names(asv_sample)) {
    cat("\n--- Read count summary per sample ---\n")
    sample_totals <- asv_sample %>%
      group_by(dnaSampleID) %>%
      summarise(
        total_reads = sum(individualCount, na.rm = TRUE),
        n_asvs = n(),
        n_phyla = n_distinct(phylum[!is.na(phylum) & phylum != ""]),
        .groups = "drop"
      )
    cat("Samples:", nrow(sample_totals), "\n")
    cat("Reads per sample: median", median(sample_totals$total_reads),
        ", range", min(sample_totals$total_reads), "-", max(sample_totals$total_reads), "\n")
    cat("ASVs per sample: median", median(sample_totals$n_asvs),
        ", range", min(sample_totals$n_asvs), "-", max(sample_totals$n_asvs), "\n")
    cat("Phyla per sample: median", median(sample_totals$n_phyla),
        ", range", min(sample_totals$n_phyla), "-", max(sample_totals$n_phyla), "\n")

    # Compute phylum relative abundances for these samples
    phylum_abund <- asv_sample %>%
      filter(!is.na(phylum) & phylum != "") %>%
      group_by(dnaSampleID, phylum) %>%
      summarise(reads = sum(individualCount, na.rm = TRUE), .groups = "drop") %>%
      group_by(dnaSampleID) %>%
      mutate(rel_abund = reads / sum(reads)) %>%
      ungroup()

    cat("\nPhylum relative abundance ranges (across sampled files):\n")
    phylum_summary <- phylum_abund %>%
      group_by(phylum) %>%
      summarise(
        mean_rel = mean(rel_abund),
        min_rel = min(rel_abund),
        max_rel = max(rel_abund),
        .groups = "drop"
      ) %>%
      arrange(desc(mean_rel))
    print(head(phylum_summary, 15))
  }
} else {
  cat("No 16S map files found/extracted. Run unzip manually to check taxonomy.\n")
}

# --- 8. Samples per site-month (for hindcast evaluation feasibility) ---
cat("\n========================================\n")
cat("=== SAMPLING DENSITY ===\n")
cat("========================================\n")

site_month_counts <- meta_16S_qc %>%
  count(siteID, dateID, name = "n_samples") %>%
  arrange(siteID, dateID)

cat("\n16S samples per site-month:\n")
cat("  Median:", median(site_month_counts$n_samples), "\n")
cat("  Range:", min(site_month_counts$n_samples), "-", max(site_month_counts$n_samples), "\n")

# For existing sites, show post-hindcast density
post_hindcast <- site_month_counts %>%
  filter(siteID %in% existing_sites_16S, dateID > max(existing_dates_16S))

if (nrow(post_hindcast) > 0) {
  cat("\nPost-hindcast sampling density (existing sites only):\n")
  cat("  Site-months:", nrow(post_hindcast), "\n")
  cat("  Samples per site-month: median", median(post_hindcast$n_samples),
      ", range", min(post_hindcast$n_samples), "-", max(post_hindcast$n_samples), "\n")

  # Timeline plot
  p <- ggplot(post_hindcast, aes(x = dateID, y = siteID, fill = n_samples)) +
    geom_tile() +
    scale_fill_viridis_c() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
    labs(title = "Post-hindcast sampling density (16S)",
         subtitle = paste("Dates >", max(existing_dates_16S)),
         x = "Date (YYYYMM)", y = "Site", fill = "N samples")
  ggsave(here("figures", "new_data_sampling_density.png"), p,
         width = 14, height = 8, dpi = 150)
  cat("\nSaved: figures/new_data_sampling_density.png\n")
}

# --- 9. Summary ---
cat("\n========================================\n")
cat("=== ASSESSMENT SUMMARY ===\n")
cat("========================================\n")

n_eval_site_months <- nrow(post_hindcast)
n_eval_sites <- length(unique(post_hindcast$siteID))

cat(sprintf("
For hindcast evaluation:
- %d/%d existing 16S sites have new data
- %d site-months of post-hindcast data (after %s) across %d sites
- New data spans %s to %s
- Data product is DP1.10081.002 (marker gene taxonomy)
- ASV map files contain full taxonomy (kingdom -> species) + read counts
- Phylum-level relative abundances can be computed from individualCount

Next steps:
1. Extract all 16S/ITS map files (large: ~28k files)
2. Aggregate ASV reads to phylum-level relative abundances per sample
3. Match dnaSampleID -> siteID/plotID/dateID via metadata
4. Compare with hindcast predictions at overlapping site-dates
5. Focus on post-%s dates for true out-of-sample evaluation
",
length(overlap_16S), length(existing_sites_16S),
n_eval_site_months, max(existing_dates_16S), n_eval_sites,
as.character(min(meta_16S_qc$collectDate, na.rm = TRUE)),
as.character(max(meta_16S_qc$collectDate, na.rm = TRUE)),
max(existing_dates_16S)))

cat("\nDone.\n")
