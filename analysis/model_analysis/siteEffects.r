# Visualize site effects (random effects) from all model
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
pacman::p_load(stringr, forestplot, gridExtra, ggpubr) 


site_effects_all <- readRDS("./data/summary/site_effects.rds")
site_effects_refit <- site_effects_all %>% filter(time_period == "refit" & model_name == "all_covariates")

## Get micronutrient, climate, and land-cover data
fieldsites_raw <- read.csv("https://www.neonscience.org/sites/default/files/NEON_Field_Site_Metadata_20210226_0.csv")
fieldsites <- fieldsites_raw %>% filter(!grepl("Aquatic", field_site_type))

df_merged <- merge(site_effects_refit, fieldsites, by.x = "siteID", by.y = "field_site_id")


by_site <- df_merged %>% 
	pivot_wider(id_cols = c("rank","taxon"), names_from = "siteID", values_from = "Mean") %>% select(-c(1:2))
site_M <- cor(by_site)
corrplot(site_M, method="circle", sig.level = .05, type = "lower")




ggplot(df_merged) + geom_point(aes(x = field_mean_annual_precipitation_mm, y = Mean)) + facet_wrap(~taxon)
ggplot(df_merged) + geom_point(aes(x = field_mean_annual_temperature_C, y = Mean)) + facet_wrap(~taxon)



### Read in bulk density data
bd_raw <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/bulkDensity_allsites.rds")
bd <- bd_raw %>% filter(bulkDensTopDepth < 25) %>% group_by(siteID) %>%  
	summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

### Read in micronutrient data
bgc_raw <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/raw/dp1.10047.rds")
# remove deep horizons.
bgc <- bgc_raw %>% filter(biogeoTopDepth < 25)

#model_plotIDs <- unique(substr(model.dat$plotID, 6, 15))
model_plotIDs <- unique(model.dat$plotID)
model_siteIDs <- unique(model.dat$siteID)
length(intersect(model_plotIDs, unique(bgc$plotID)))
length(setdiff(model_plotIDs, unique(bgc$plotID)))

# bgc_site <- bgc %>% filter(siteID %in% model_siteIDs)
# ggplot(bgc_site) + geom_boxplot(aes(x= plotID, y = ))

relevant_cols <- c("siteID","plotID","biogeoTopDepth","biogeoBottomDepth","gypsumConc","caco3Conc","caNh4d","kNh4d", "mgNh4d", "naNh4d", "cecdNh4", "alSatCecd33","alKcl", "feKcl", "mnKcl","MehlichIIITotP","alOxalate", "feOxalate", "mnOxalate", "pOxalate", "siOxalate","estimatedOC","sulfurTot")

bgc_relevant <- bgc %>% select(!!!relevant_cols)

cor(bgc_relevant[,5:23])

library(corrplot)
most_relevant <- c("siteID","plotID",
									 "kSatx", "mgSatx","pSatx","naSatx",
									 "ececCecd33","mnKcl", 
									  "no2Satx", "no3Satx","alOxalate", "feOxalate", "so4Satx",
									 "mnOxalate", "pOxalate", "siOxalate","estimatedOC","sulfurTot","MehlichIIITotP","OlsenPExtractable","Bray1PExtractable"#,"caco3Conc"
									 )
bgc_relevant <- bgc %>% select(!!!most_relevant)
#corrplot(cor(bgc_relevant[,3:17], use = "pairwise.complete.obs"))
# not.na <- bgc_relevant[which(!is.na(bgc_relevant$mgSatx)),]
# ggplot(bgc_relevant) + geom_boxplot(aes(x = siteID, y = MehlichIIITotP, color = as.factor(plotID)), show.legend = F)

bgc_relevant$totalP <- ifelse(is.na(bgc_relevant$MehlichIIITotP),
															bgc_relevant$OlsenPExtractable,
															bgc_relevant$MehlichIIITotP)
bgc_relevant[,c("MehlichIIITotP","OlsenPExtractable","Bray1PExtractable")] <-NULL
bgc_relevant <- bgc_relevant %>% group_by(siteID) %>%  
	summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))


df_relevant <- df_merged %>% select(siteID,site_num,taxon_num,only_rank,Mean,taxon,field_mean_annual_temperature_C,field_mean_annual_precipitation_mm)

df_bgc <- merge(df_relevant, bgc_relevant)
df_bgc <- merge(df_bgc, bd[,c("siteID","bulkDensThirdBar")], all.y = F)

df_bgc[,c(7:26)] <- scale(df_bgc[,c(7:26)])[[1]]


library(corrplot)
corrplot(cor(df_bgc[,c(7:ncol(df_bgc))], use = "pairwise.complete.obs"))

library(MuMIn)

options(na.action = "na.fail")

dredge.out <- list()
weird <- list()
tax_list <- unique(df_bgc$taxon)
tax_list <- tax_list[!tax_list %in% c("other",NA)]#,"proteobacteria")]
s = "cyanobacteria"
for (s in tax_list){
	taxon_rank <- site_effects_refit %>% filter(taxon == s) %>% select(rank) %>% unique()
	
	species_dat <- df_bgc %>% filter(taxon==s) %>% select(-c(siteID,site_num,taxon_num,only_rank))
	species_dat <- species_dat[complete.cases(species_dat),]
	species_lm <-lm(Mean ~ field_mean_annual_temperature_C + 
										#bulkDensThirdBar +
										field_mean_annual_precipitation_mm + no3Satx + so4Satx +
										kSatx + mgSatx + #pSatx + 
										naSatx + #pOxalate + 
										ececCecd33 +
 + #mnKcl  + #estimatedOC  
 #	caco3Conc +
									#	sulfurTot + 
 	totalP, data = species_dat)
#results<-dredge(global.model = species_lm, m.lim = c(1,5))

models <- lapply(dredge(species_lm, evaluate = FALSE, m.lim = c(1,5)), eval)
ma <- model.avg(models,  subset = delta <= 2)
avg_results <- summary(ma)
avg <- as.data.frame(avg_results$coefmat.subset)
avg$taxon <- s
avg$taxon_rank <- taxon_rank[[1]]
#avg <- avg[avg$`Pr(>|z|)` < .1]
dredge.out[[s]] <- avg %>% rownames_to_column("predictor")
}
# 
# if(nrow(results[results$delta < 2,])==1){
# 	#results <- as.data.frame(results[1,])
# 	
# 	top_model <- get.models(results, subset = 1)[[1]]
# 	top_results <- summary(top_model)
# 	res <- as.data.frame(top_results$coefficients)
# 	res$taxon <- s
# 	res$taxon_rank <- taxon_rank[[1]]
# 	weird[[s]] <- res %>% rownames_to_column("predictor")
# } else {
# 	if(nrow(results[results$delta < 2,])> 200){
# 		weird[[s]] <- "idk"
# 		next()
# 	}
# 		
# 	subs <- results[results$delta < 2,]
# 	model.avg(subs)
# avg_results <- summary(model.avg(results, subset = delta <= 2))
# avg <- as.data.frame(avg_results$coefmat.subset)
# avg$taxon <- s
# avg$taxon_rank <- taxon_rank[[1]]
# #avg <- avg[avg$`Pr(>|z|)` < .1]
# dredge.out[[s]] <- avg %>% rownames_to_column("predictor")


dredge.full <- do.call(rbind, dredge.out)
weird.full <- do.call(rbind, weird)

all.out <- plyr::rbind.fill(dredge.full, weird.full)
all.out <- merge(all.out, unique(site_effects_refit[,c("rank","pretty_name","pretty_group")]), all.x=T, by.x="taxon_rank", by.y="rank")
all.out$`effect size` <- abs(all.out$Estimate)

all.out_sorted <- all.out %>% #group_by(pretty_group) %>% 
	arrange(`effect size`, predictor) 
all.out$predictor <- factor(all.out$predictor,levels=c("(Intercept)","naSatx","kSatx","mgSatx","so4Satx","ececCecd33","no3Satx","field_mean_annual_temperature_C","totalP","field_mean_annual_precipitation_mm"))


ggplot(all.out) + geom_point(aes(x = pretty_name, 
																 y = `effect size`, 
																 fill = pretty_group),
														 size=3, alpha = 1,
														 shape=21, position=position_jitterdodge(dodge.width = 1), color = "black") + 
	facet_wrap(~predictor, scale="free") +
	theme_bw() + theme(
		text = element_text(size = 16),
		axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
		axis.title=element_text(size=22,face="bold")
	) + scale_fill_manual(values = c("grey30","grey90"))


ggplot(df_bgc) + geom_point(aes(x = totalP, y = Mean)) + facet_wrap(~species)
ggplot(df_merged) + geom_point(aes(x = field_mean_annual_temperature_C, y = Mean)) + facet_wrap(~species)






# Visualize overall size of site effects
ggplot(data=site_effects_refit,
			 aes(x = pretty_name,y = abs(Mean)))+
	geom_jitter(aes(color = as.factor(siteID)), size = 4, width=.2) +
	labs(col = "Site", title = "Absolute site effect size") + 
	xlab("Rank")+ 	facet_grid(#rows = vars(only_rank), 
		rows = vars(pretty_group), drop = T,
		scales = "free", space = "free_x") +
	ylab(NULL)+
	theme_bw() + theme(
		text = element_text(size = 18),
		axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
		axis.title=element_text(size=22,face="bold")
	)




