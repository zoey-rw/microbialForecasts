source("../../source.R")
library(kableExtra)

# Production calibration is 20130601-20180101; validation is the remaining
# data through 20200101. min.prev=3 used to match production model fits.
CAL_MIN <- "20130601"
CAL_MAX <- "20180101"
VAL_MAX <- "20200101"

fieldsites_raw <- read.csv("https://www.neonscience.org/sites/default/files/NEON_Field_Site_Metadata_20220412.csv") %>%
	filter(!grepl("Aquatic", field_site_type)) %>%
	select(siteID = field_site_id,
				 `site name` = field_site_name,
				 MAT = field_mean_annual_temperature_C,
				 MAP = field_mean_annual_precipitation_mm,
				 latitude = field_latitude,
				 longitude = field_longitude,
				 nlcd = field_dominant_nlcd_classes,
				 ecoregion = field_domain_id)

abun_its <- readRDS(here("data/clean/groupAbundances_ITS_2023.rds"))
abun_16s <- readRDS(here("data/clean/groupAbundances_16S_2023.rds"))

# Per-site counts: built from two prepBetaRegData calls.
# (1) "obs" call (default full_timeseries=FALSE, min.prev=3) gives the 30
#     observed sites with the same filtering as production model fits, so
#     cal+val totals here match the model-input counts in the Methods.
# (2) "new" call (full_timeseries=TRUE, min.prev=0) restores the 13
#     held-out new-sites that have no cal samples but contribute to
#     new-site validation.
per_site_counts <- function(rank.df, abundance_col) {
	obs <- prepBetaRegData(rank.df = rank.df, min.prev = 3,
												 min.date = CAL_MIN, max.date = VAL_MAX)
	new <- prepBetaRegData(rank.df = rank.df, min.prev = 0,
												 min.date = CAL_MIN, max.date = VAL_MAX,
												 full_timeseries = TRUE)

	jan2018 <- unique(obs$truth.plot.long$timepoint[
										obs$truth.plot.long$dateID == "201801"])[1]

	obs_counts <- obs$sample_values %>%
		filter(!is.na(.data[[abundance_col]])) %>%
		mutate(period = ifelse(timepoint <= jan2018, "cal", "val")) %>%
		group_by(siteID, period) %>% tally() %>%
		pivot_wider(names_from = period, values_from = n, values_fill = 0)

	new_site_ids <- setdiff(unique(new$sample_values$siteID),
													unique(obs$sample_values$siteID))
	new_counts <- new$sample_values %>%
		filter(siteID %in% new_site_ids, !is.na(.data[[abundance_col]])) %>%
		group_by(siteID) %>% tally(name = "val") %>%
		mutate(cal = 0L) %>% select(siteID, cal, val)

	bind_rows(obs_counts, new_counts) %>%
		mutate(cal = ifelse(cal == 0, NA_integer_, cal),
					 val = ifelse(val == 0, NA_integer_, val))
}

bac_counts <- per_site_counts(abun_16s$copiotroph, "copiotroph") %>%
	rename(`Calibration bacterial samples` = cal,
				 `Validation bacterial samples`  = val)
fun_counts <- per_site_counts(abun_its$ectomycorrhizal, "ectomycorrhizal") %>%
	rename(`Calibration fungal samples` = cal,
				 `Validation fungal samples`  = val)

bac_site_n     <- bac_counts %>% select(siteID, `Calibration bacterial samples`)
val_bac_site_n <- bac_counts %>% select(siteID, `Validation bacterial samples`)
fun_site_n     <- fun_counts %>% select(siteID, `Calibration fungal samples`)
val_fun_site_n <- fun_counts %>% select(siteID, `Validation fungal samples`)

fieldsites_samples <- merge(bac_site_n, fun_site_n, all = TRUE) %>%
	merge(merge(val_bac_site_n, val_fun_site_n, all = TRUE), all = TRUE) %>%
	merge(fieldsites_raw, all = TRUE) %>%
	filter(!(is.na(`Validation bacterial samples`) &
					 is.na(`Calibration bacterial samples`) &
					 is.na(`Validation fungal samples`) &
					 is.na(`Calibration fungal samples`))) %>%
	rename("Land cover (NLCD class)" = "nlcd")

# Reorder columns to match the existing Forecasting_TableS1.csv layout
fieldsites_samples <- fieldsites_samples[, c(
	"siteID",
	"Calibration bacterial samples", "Calibration fungal samples",
	"Validation bacterial samples", "Validation fungal samples",
	"site name", "MAT", "MAP", "latitude", "longitude",
	"Land cover (NLCD class)", "ecoregion")]

write.csv(fieldsites_samples, here("figures", "Forecasting_TableS1.csv"),
					row.names = FALSE)

kable(fieldsites_samples, "html") %>%
	kable_styling(bootstrap_options = c("striped", "hover")) %>%
	cat(., file = here("figures", "sample_size.html"))

cat("Wrote figures/Forecasting_TableS1.csv with", nrow(fieldsites_samples), "sites\n")
cat("Totals (calibration", CAL_MIN, "->", CAL_MAX, "):\n")
cat("  Bacterial:", sum(fieldsites_samples$`Calibration bacterial samples`, na.rm = TRUE), "\n")
cat("  Fungal:   ", sum(fieldsites_samples$`Calibration fungal samples`,    na.rm = TRUE), "\n")
cat("Totals (validation", CAL_MAX, "->", VAL_MAX, "):\n")
cat("  Bacterial:", sum(fieldsites_samples$`Validation bacterial samples`, na.rm = TRUE), "\n")
cat("  Fungal:   ", sum(fieldsites_samples$`Validation fungal samples`,    na.rm = TRUE), "\n")

# Raw ESV table sizes reported in the Methods ("Bioinformatics and Data
# Processing"). The DADA2/singleton-removal ESV tables are the phyloseq objects
# copied from the upstream HPC pipeline: ntaxa() is the number of ESVs and
# nsamples() the number of samples, per kingdom.
if (requireNamespace("phyloseq", quietly = TRUE)) {
	ps_16s <- readRDS(here("data/clean/phyloseq_16S.rds"))
	ps_its <- readRDS(here("data/clean/phyloseq_ITS.rds"))
	cat("Raw ESV tables (after DADA2 processing and singleton removal):\n")
	cat("  Bacterial:", phyloseq::ntaxa(ps_16s), "ESVs across",
			phyloseq::nsamples(ps_16s), "samples\n")
	cat("  Fungal:   ", phyloseq::ntaxa(ps_its), "ESVs across",
			phyloseq::nsamples(ps_its), "samples\n")
} else {
	cat("phyloseq not installed; skipping raw ESV table size report\n")
}
