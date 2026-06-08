#!/usr/bin/env Rscript
# Download the large input files that are NOT stored in git, needed to reproduce
# the downstream analyses and figures. Small inputs and all results that fit in
# git are committed in the repository; the larger inputs live on Zenodo and are
# fetched (and md5-verified) here.
#
# Most figures and the phylogenetic-signal analysis read these once present:
#   * Forecast scores / skill metrics (scoring_metrics_plsr2.rds) and the tidy
#     summaries (seasonal_amplitude, rho_core_sd_effects, site_effect_variograms,
#     fcast_horizon_*, site_effects*, *_converged_taxa_list, predictor/effect
#     summaries) -> most scoring, seasonality, and site-effect figures.
#   * phyloseq_16S.rds / phyloseq_ITS.rds -> tree-from-raw + unclassified-reads.
#   * monthly_soil_{temperature,moisture}.rds, soilCore_raw_measurements.rds
#     -> the soil-covariate supplementary figures.
#   * example_hindcasts_HARV_CPER.tar.gz -> per-site hindcasts for the example
#     forecast figure (extracted into data/hindcasts/driver_uncertainty/).
#
# NOTE: the aggregate hindcast parquet (all_hindcasts_plsr2.parquet, ~1 GB) used
# by load_hindcasts() for pooled-R2 figures is NOT yet in this manifest — it is
# being renamed/split by forecast type first; add it here once deposited.
#
# Setup: after the Zenodo deposit is published, set the record id below (or the
# MF_ZENODO_BASE environment variable) to the record's file base URL, e.g.
#   export MF_ZENODO_BASE="https://zenodo.org/records/1234567/files"

suppressMessages(library(here))

ZENODO_BASE <- Sys.getenv(
	"MF_ZENODO_BASE",
	"https://zenodo.org/records/REPLACE_WITH_RECORD_ID/files"   # <-- edit after upload
)

# Individual files: Zenodo file name -> local destination, with md5 to verify.
files <- list(
	# data/clean
	list(name = "phyloseq_16S.rds",            dest = here("data", "clean", "phyloseq_16S.rds"),
			 md5 = "0bdfa98e04eeec523fe6cfb253319899", note = "16S phyloseq; create_bacterial_phylogeny.r, unclassified reads"),
	list(name = "phyloseq_ITS.rds",            dest = here("data", "clean", "phyloseq_ITS.rds"),
			 md5 = "9a021f6a3ac3678f90e6c6e1695e8021", note = "ITS phyloseq; unclassified_taxonomy_distribution.r"),
	list(name = "monthly_soil_temperature.rds", dest = here("data", "clean", "monthly_soil_temperature.rds"),
			 md5 = "ddb4b9a607f6c6227e94d84fd15088cb", note = "soil temp gap-fill; figS14"),
	list(name = "monthly_soil_moisture.rds",   dest = here("data", "clean", "monthly_soil_moisture.rds"),
			 md5 = "6f3bc0c34269062887ba58297d6ae31d", note = "soil moisture gap-fill; figS15"),
	list(name = "soilCore_raw_measurements.rds", dest = here("data", "clean", "soilCore_raw_measurements.rds"),
			 md5 = "6e4e7581483806a6ff9f359d95449833", note = "soil pH/C raw; figS16"),
	# data/summary
	list(name = "scoring_metrics_plsr2.rds",   dest = here("data", "summary", "scoring_metrics_plsr2.rds"),
			 md5 = "2d712c9919f3b22fa15a6a8f6f72d044", note = "forecast scores; most scoring/skill figures + phylo_contribution_scores.r"),
	list(name = "seasonal_amplitude.rds",      dest = here("data", "summary", "seasonal_amplitude.rds"),
			 md5 = "0260eeb9f98411d760f566f0d7ee7682", note = "seasonal amplitude; seasonality/amplitude figures"),
	list(name = "rho_core_sd_effects.rds",     dest = here("data", "summary", "rho_core_sd_effects.rds"),
			 md5 = "07eb696684802e96c9fdeabfb6977b3f", note = "rho + core_sd; spatial/temporal trait figures"),
	list(name = "site_effect_variograms.rds",  dest = here("data", "summary", "site_effect_variograms.rds"),
			 md5 = "ed8c15af580de82676e23219b8e0268b", note = "variogram p-values; fig_variogram.r"),
	list(name = "fcast_horizon_df.rds",        dest = here("data", "summary", "fcast_horizon_df.rds"),
			 md5 = "eac958dc49b609459fdfcf540bf0b07b", note = "forecast horizon; horizon figures"),
	list(name = "fcast_horizon_input.rds",     dest = here("data", "summary", "fcast_horizon_input.rds"),
			 md5 = "3873f5ab1b48e57865e7564f6479f214", note = "forecast horizon input; skill-driver figures"),
	list(name = "site_effects_dredged.rds",    dest = here("data", "summary", "site_effects_dredged.rds"),
			 md5 = "4bf39f5c888b290fe0ad292cae351852", note = "PLSR site effects; fig4_site_eff_predictors / biplot"),
	list(name = "site_effects.rds",            dest = here("data", "summary", "site_effects.rds"),
			 md5 = "c8fc3a4d2c74e947e7fb8f4daa794609", note = "site effects; steps 05-06"),
	list(name = "weak_converged_taxa_list.rds", dest = here("data", "summary", "weak_converged_taxa_list.rds"),
			 md5 = "633d6b84b9495b1eeb73f08528d555b0", note = "weak-converged taxa; effect-size figures"),
	list(name = "converged_taxa_list.rds",     dest = here("data", "summary", "converged_taxa_list.rds"),
			 md5 = "531fa669b7c066ea585ab8ee4d6c3d95", note = "converged taxa; effect-size figures"),
	list(name = "logit_beta_fixed_priors_summaries.rds", dest = here("data", "summary", "logit_beta_fixed_priors_summaries.rds"),
			 md5 = "1de230a9df401ef07cd6c64524774557", note = "full parameter summary; reportModelParameters / legacy effect")
)

# Archives: tarballs whose members are repo-relative paths, extracted at the root.
archives <- list(
	list(name = "example_hindcasts_HARV_CPER.tar.gz",
			 md5  = "40628a46b89f045ca0ae16ab650d6553",
			 check = here("data", "hindcasts", "driver_uncertainty",
									 "hindcasts_env_cycl_cellulolytic_20130601_20180101_with_legacy_covariate_HARV_observed.rds"),
			 note = "per-site hindcasts (HARV/CPER, env_cycl) for fig_exampleHindcasts.r")
)

if (grepl("REPLACE_WITH_RECORD_ID", ZENODO_BASE))
	stop("Set MF_ZENODO_BASE (or edit ZENODO_BASE) to the published Zenodo record URL first.")

fetch <- function(name, dest, md5) {
	url <- file.path(ZENODO_BASE, paste0(name, "?download=1"))
	cat("Downloading", name, "->", dest, "\n")
	download.file(url, dest, mode = "wb")
	got <- unname(tools::md5sum(dest))
	if (got != md5)
		stop("md5 mismatch for ", name, ": got ", got, ", expected ", md5)
}

for (f in files) {
	dir.create(dirname(f$dest), showWarnings = FALSE, recursive = TRUE)
	if (file.exists(f$dest) && unname(tools::md5sum(f$dest)) == f$md5) {
		cat("OK (already present):", f$name, "\n"); next
	}
	fetch(f$name, f$dest, f$md5)
	cat("  verified --", f$note, "\n")
}

for (a in archives) {
	if (!is.null(a$check) && file.exists(a$check)) {
		cat("OK (already extracted):", a$name, "\n"); next
	}
	tmp <- tempfile(fileext = ".tar.gz")
	fetch(a$name, tmp, a$md5)
	utils::untar(tmp, exdir = here())   # members are repo-relative paths
	unlink(tmp)
	cat("  extracted --", a$note, "\n")
}

cat("All Zenodo-hosted input files present and verified.\n")
