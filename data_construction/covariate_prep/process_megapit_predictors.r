# Process NEON megapit data (DP1.00096) into site-level predictors
# Input: downloaded megapit files in data/clean/NEON_soil-megapit/
# Output: data/clean/site_effect_predictors.rds (updated with new variables)
#
# Run from project root: Rscript data_construction/covariate_prep/process_megapit_predictors.r

library(dplyr)
library(here)
here::i_am("source.R")

cat("=== Processing NEON megapit data ===\n")

# ---- 1. Read the downloaded megapit CSVs ----

megapit_dir <- here("data/clean/NEON_soil-megapit")
site_dirs <- list.dirs(megapit_dir, recursive = FALSE, full.names = TRUE)
site_dirs <- site_dirs[grepl("NEON\\.D\\d+\\.", basename(site_dirs))]

read_megapit_table <- function(dirs, pattern) {
  files <- unlist(lapply(dirs, function(d) list.files(d, pattern = pattern, full.names = TRUE)))
  lapply(files, read.csv, stringsAsFactors = FALSE) %>% bind_rows()
}

bgc_all <- read_megapit_table(site_dirs, "mgp_perbiogeosample")
bd_all  <- read_megapit_table(site_dirs, "mgp_perbulksample")

cat("Biogeosample rows:", nrow(bgc_all), " | Bulk density rows:", nrow(bd_all), "\n")
cat("Sites:", length(unique(bgc_all$siteID)), "\n")

# ---- 2. Extract biogeochemistry (shallow horizons only) ----

bgc <- bgc_all
cat("Biogeochem rows:", nrow(bgc), "\n")
cat("Sites:", length(unique(bgc$siteID)), "\n")

# Filter to shallow soil (top 25 cm) — matches original pipeline
bgc_shallow <- bgc %>% filter(biogeoTopDepth < 25)
cat("Shallow (<25cm) rows:", nrow(bgc_shallow), "\n")

# Aggregate to site means for all numeric columns
bgc_site <- bgc_shallow %>%
  group_by(siteID) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(totalP = pMjelm)

# ---- 3. Extract bulk density ----

if (nrow(bd_all) > 0) {
  depth_col <- intersect(c("bulkDensTopDepth", "biogeoTopDepth"), names(bd_all))
  if (length(depth_col) > 0) {
    bd_shallow <- bd_all %>% filter(.data[[depth_col[1]]] < 25)
  } else {
    bd_shallow <- bd_all
  }

  # Use bulkDensExclCoarseFrag (standard NEON measurement)
  if ("bulkDensExclCoarseFrag" %in% names(bd_shallow)) {
    bd_shallow$bulkDens <- bd_shallow$bulkDensExclCoarseFrag
  } else if ("bulkDensThirdBar" %in% names(bd_shallow)) {
    bd_shallow$bulkDens <- ifelse(is.na(bd_shallow$bulkDensThirdBar),
                                  bd_shallow$bulkDensFieldMoist,
                                  bd_shallow$bulkDensThirdBar)
  }

  bd_site <- bd_shallow %>%
    group_by(siteID) %>%
    summarise(bulkDens = mean(bulkDens, na.rm = TRUE), .groups = "drop")
} else {
  bd_site <- data.frame(siteID = character(0), bulkDens = numeric(0))
}

# ---- 4. Get site climate/location metadata from NEON ----

loc_cache <- here("data/clean/neon_fieldsites_loc.rds")
if (file.exists(loc_cache)) {
  fieldsites_loc <- readRDS(loc_cache)
} else {
  fieldsites_raw <- read.csv("https://www.neonscience.org/sites/default/files/NEON_Field_Site_Metadata_20220412.csv")
  fieldsites_loc <- fieldsites_raw %>%
    filter(!grepl("Aquatic", field_site_type)) %>%
    select(latitude = field_latitude, longitude = field_longitude, siteID = field_site_id)
  saveRDS(fieldsites_loc, loc_cache)
}

fieldsites_raw <- read.csv("https://www.neonscience.org/sites/default/files/NEON_Field_Site_Metadata_20220412.csv")
fieldsites_climate <- fieldsites_raw %>%
  filter(!grepl("Aquatic", field_site_type)) %>%
  select(siteID = field_site_id,
         MAT = field_mean_annual_temperature_C,
         MAP = field_mean_annual_precipitation_mm,
         latitude = field_latitude,
         longitude = field_longitude,
         nlcd = field_dominant_nlcd_classes,
         ecoregion = field_domain_id)

# ---- 5. Merge everything ----

df_predictors <- fieldsites_climate %>%
  left_join(bgc_site, by = "siteID") %>%
  left_join(bd_site, by = "siteID")

# Compute C:N ratio before scaling (needs raw values)
df_predictors$cnRatio <- df_predictors$carbonTot / df_predictors$nitrogenTot
df_predictors$cnRatio[!is.finite(df_predictors$cnRatio)] <- NA

cat("\nMerged predictor data:", nrow(df_predictors), "sites\n")

# ---- 6. Scale numeric predictors ----

# Keep latitude and longitude unscaled (raw values needed for variograms)
no_scale <- c("latitude", "longitude")
for (col in names(df_predictors)) {
  if (is.numeric(df_predictors[[col]]) && !(col %in% no_scale)) {
    df_predictors[[col]] <- as.vector(scale(df_predictors[[col]]))
  }
}

# ---- 7. Report completeness ----

cat("\n=== Predictor completeness ===\n")
for (col in sort(names(df_predictors))) {
  if (!is.numeric(df_predictors[[col]])) next
  n_valid <- sum(!is.na(df_predictors[[col]]))
  cat(sprintf("  %-18s %2d/%d sites\n", col, n_valid, nrow(df_predictors)))
}

# ---- 8. Save ----

saveRDS(df_predictors, here("data/clean/site_effect_predictors.rds"))
cat("\nSaved: data/clean/site_effect_predictors.rds\n")
cat("Columns:", paste(names(df_predictors), collapse = ", "), "\n")
