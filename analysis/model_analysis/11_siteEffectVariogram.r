# Check for spatial autocorrelation in site effects vs PLSR cross-validated residuals
# Uses observed sites only (unobserved sites have predicted effects, no true residuals)

library(here)
source(here("source.R"))

site_eff_dredged_in <- readRDS(here("data/summary/site_effects_dredged.rds"))
scores_list <- readRDS(here("data/summary/scoring_metrics_plsr2.rds"))
converged_strict <- scores_list$converged_strict_list

library(gstat)
library(sp)
library(variosig)
library(foreach)
library(doParallel)

## Get site location data from NEON (cached locally if available)
cache_path <- here("data/clean/neon_fieldsites_loc.rds")
if (file.exists(cache_path)) {
  fieldsites_loc <- readRDS(cache_path)
} else {
  fieldsites_raw <- read.csv("https://www.neonscience.org/sites/default/files/NEON_Field_Site_Metadata_20220412.csv")
  fieldsites_loc <- fieldsites_raw %>%
    filter(!grepl("Aquatic", field_site_type)) %>%
    select(latitude = field_latitude, longitude = field_longitude, siteID = field_site_id)
  saveRDS(fieldsites_loc, cache_path)
}

# Observed sites: d[[3]] has TargetVar (raw site effect) and PLSR_CV_Residuals
obs_sites <- site_eff_dredged_in[[3]] %>%
  select(model_id, siteID, taxon, model_name,
         site_effect = TargetVar,
         resid = PLSR_CV_Residuals) %>%
  filter(!is.na(site_effect), !is.na(resid))

model_id_list <- obs_sites %>%
  filter(!grepl("other", taxon)) %>%
  distinct(model_id) %>%
  pull(model_id)

cat("Running variograms for", length(model_id_list), "models\n")

cl <- makeCluster(4, type = "FORK", outfile = "")
registerDoParallel(cl)

# 99 permutations is sufficient for detecting p < 0.05
N_PERM <- 99

sig_results_list <- foreach(mid = model_id_list, .errorhandling = "remove") %dopar% {
  dat <- obs_sites %>% filter(model_id == mid)
  points.df <- merge(fieldsites_loc, dat)
  if (nrow(points.df) < 5) return(NULL)

  TheData <- points.df
  coordinates(TheData) <- ~ longitude + latitude

  TheVariogram <- variogram(site_effect ~ 1, data = TheData)
  TheResidualVariogram <- variogram(resid ~ 1, data = TheData)

  eff_envelope <- envelope(TheVariogram, data = TheData, formula = site_effect ~ 1, nsim = N_PERM)
  eff_signif <- envsig(eff_envelope, method = "eb")

  resid_envelope <- envelope(TheResidualVariogram, data = TheData, formula = resid ~ 1, nsim = N_PERM)
  resid_signif <- envsig(resid_envelope, method = "eb")

  taxon_val <- unique(dat$taxon)
  model_name_val <- unique(dat$model_name)

  data.frame(
    taxon = taxon_val, rank = NA_character_, model_id = mid,
    model_name = model_name_val,
    site_effect_p = eff_signif$p.overall,
    residual_p = resid_signif$p.overall,
    stringsAsFactors = FALSE
  )
}

stopCluster(cl)

sig_results_list <- sig_results_list[!sapply(sig_results_list, is.null)]
sig_results <- do.call(rbind, sig_results_list)
colnames(sig_results) <- c("taxon", "rank", "model_id", "model_name",
                           "site effect", "site effect residuals")

# Save without plot objects (saves space, plots rarely used)
saveRDS(list(sig_results, list()), here("data/summary/site_effect_variograms.rds"))

cat("Saved variogram results for", nrow(sig_results), "models\n")
