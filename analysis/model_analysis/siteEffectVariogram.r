
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

df_predictors <- readRDS(here("data/summary/site_effect_predictors.rds"))

site_eff_dredged_in <- readRDS(here("data/summary/site_effects_dredged.rds"))
unobs_sites <- readRDS(here("data/summary/site_effects_unobserved.rds"))


library(gstat)
library(sp)
library(variosig)


## Get site climate data from NEON
fieldsites_raw <- read.csv("https://www.neonscience.org/sites/default/files/NEON_Field_Site_Metadata_20220412.csv")
fieldsites_loc <- fieldsites_raw %>% filter(!grepl("Aquatic", field_site_type)) %>% select(latitude = field_latitude,
																																											 longitude = field_longitude,
																																											 siteID = field_site_id)
pred_sites <- site_eff_dredged_in[[2]]


taxon = "acidimicrobiia"

out_plots = list()
out_sig = list()



taxon_list = pred_sites %>% filter(!grepl("other", taxon)) %>% distinct(taxon) %>% unlist()



pacman::p_load(doParallel)
cl <- makeCluster(27, type="FORK", outfile="")
registerDoParallel(cl)


#Run for multiple groups, in parallel (via PSOCK)
# Why does this work with "do" but not "dopar"!!!!!
sig_results_list = foreach(taxon=taxon_list, .errorhandling = 'remove') %dopar% {

site_effect_estimate = pred_sites %>%
	select(siteID, taxon_rank, taxon, site_effect = Mean, pred_effect = pred) %>%
	mutate(resid = pred_effect - site_effect) %>%
	filter(taxon == !!taxon)
rank = unique(site_effect_estimate$taxon_rank)

points.df=merge(fieldsites_loc, site_effect_estimate)
TheData=points.df
coordinates(TheData)= ~ longitude+latitude

plot.new()
par(mfrow=c(2,1))
par(mar = c(2, 2, 3, 2))
TheVariogram=variogram(site_effect~1, data=TheData)
TheResidualVariogram=variogram(resid~1, data=TheData)

eff_envelope = envelope(TheVariogram, data=TheData, formula = site_effect~1)
envplot(eff_envelope, main = paste0("Site effects for: ", taxon))
eff_signif <- envsig(eff_envelope, method="eb")
mtext(text = paste0("P = ", round(eff_signif$p.overall, 4)),  outer = F, side=3)

resid_envelope = envelope(TheResidualVariogram, data=TheData, formula = resid~1)
envplot(resid_envelope, main = paste0("Residuals of site effect model for: ", taxon))
resid_signif <- envsig(resid_envelope, method="eb")
mtext(text = paste0("P = ", round(resid_signif$p.overall, 4)),  outer = F, side=3)

#out_plots[[taxon]] = recordPlot()
outplot = recordPlot()
out_sig = list(cbind(taxon, rank, eff_signif$p.overall, resid_signif$p.overall), outplot)
return(out_sig)
}

sig_results = lapply(sig_results_list, "[[", 1) %>% do.call(rbind, .) %>% as.data.frame()
plot_list =  lapply(sig_results_list, "[[", 2)

colnames(sig_results) = c("taxon", "rank","site effect", "site effect residuals")
sig_results_long = sig_results %>% pivot_longer(cols=3:4) %>% mutate(value = as.numeric(value))
ggplot(sig_results_long, aes(x = name, y= value#, color = rank
														 )) +
	geom_violin(draw_quantiles = c(.5), alpha=.2) +
	#geom_point(position = position_jitterdodge(jitter.width = .05, jitter.height = 0)) +
	geom_point(position = position_jitter(width = .1, height = 0), alpha=.3, size=3) +
	geom_hline(yintercept = .05, color=2, linetype=2) +
	theme_minimal(base_size = 16)+
	ylab("Variogram p-value") +
	xlab(NULL) +
	ggtitle("Site effect autocorrelation is captured by linear model") +
	annotate("text", x = 1.5, y = .05, label = "p.value < .05 indicates \nspatial autocorrelation")

saveRDS(list(sig_results, plot_list), here("data/summary/site_effect_variograms.rds"))

variograms = readRDS(here("data/summary/site_effect_variograms.rds"))

variograms[[2]][4] + xlab("Spatial distance bin")
