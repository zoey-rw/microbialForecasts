# Check for spatial autocorrelation in random effects estimated from Bayesian hierarchical models

source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

df_predictors <- readRDS(here("data/summary/site_effect_predictors.rds"))

site_eff_dredged_in <- readRDS(here("data/summary/site_effects_dredged.rds"))
unobs_sites <- readRDS(here("data/summary/site_effects_unobserved.rds"))
converged_strict <- readRDS(here("data/summary/converged_taxa_list.rds"))

library(gstat)
library(sp)
library(variosig)


## Get site climate data from NEON
fieldsites_raw <- read.csv("https://www.neonscience.org/sites/default/files/NEON_Field_Site_Metadata_20220412.csv")
fieldsites_loc <- fieldsites_raw %>% filter(!grepl("Aquatic", field_site_type)) %>%
	select(latitude = field_latitude,
				 longitude = field_longitude,
				 siteID = field_site_id)

out_plots = list()
out_sig = list()

model_id_list = pred_sites %>% filter(!grepl("other", taxon)) %>%
	distinct(model_id) %>% unlist()

pacman::p_load(doParallel)
cl <- makeCluster(27, type="FORK", outfile="")
registerDoParallel(cl)


#Run for multiple groups, in parallel (via PSOCK)
sig_results_list = foreach(model_id=model_id_list, .errorhandling = 'remove') %dopar% {

site_effect_estimate = pred_sites %>%
	select(model_id, siteID, rank_only, taxon, model_name, site_effect = Median, pred_effect = pred) %>%
	mutate(resid = pred_effect - site_effect) %>%
	filter(model_id == !!model_id)
rank = unique(site_effect_estimate$rank_only)
taxon = unique(site_effect_estimate$taxon)
model_name = unique(site_effect_estimate$model_name)

points.df=merge(fieldsites_loc, site_effect_estimate)
TheData=points.df
coordinates(TheData)= ~ longitude+latitude

plot.new()
par(mfrow=c(2,1))
par(mar = c(4, 2, 3, 2))
TheVariogram=variogram(site_effect~1, data=TheData)
TheResidualVariogram=variogram(resid~1, data=TheData)

eff_envelope = envelope(TheVariogram, data=TheData, formula = site_effect~1)
envplot(eff_envelope, main = paste0("Site effects for: ", model_id), xlab = "Spatial distance bin")
eff_signif <- envsig(eff_envelope, method="eb")
mtext(text = paste0("P = ", round(eff_signif$p.overall, 4)),  outer = F, side=3)

resid_envelope = envelope(TheResidualVariogram, data=TheData, formula = resid~1)
envplot(resid_envelope, main = paste0("Residuals of site effect model for: ", model_id), xlab = "Spatial distance bin")
resid_signif <- envsig(resid_envelope, method="eb")
mtext(text = paste0("P = ", round(resid_signif$p.overall, 4)),  outer = F, side=3)

#out_plots[[taxon]] = recordPlot()
outplot = recordPlot()
out_sig = list(cbind(taxon, rank, model_id, model_name, eff_signif$p.overall, resid_signif$p.overall), outplot)
return(out_sig)
}

sig_results = lapply(sig_results_list, "[[", 1) %>% do.call(rbind, .) %>% as.data.frame()
plot_list =  lapply(sig_results_list, "[[", 2)

colnames(sig_results) = c("taxon", "rank","model_id","model_name","site effect", "site effect residuals")

saveRDS(list(sig_results, plot_list), here("data/summary/site_effect_variograms.rds"))

