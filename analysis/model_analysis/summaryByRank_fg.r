library(ggplot2)
library(scales)
library(viridis)
library(hrbrthemes)
library(coda)

source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepDirichletData.r")
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

summary_list <- list()
taxon_key_list <- list()

# Read in microbial abundances
d <- c(readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_16S_2021.rds"), 
			 readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_ITS_2021.rds"))

# Only keep functional groups
fg_names <- names(d)[!grepl("bac|fun$", names(d))]

# Create lists for outputs.
estimate_list <- list()
site_eff_list <- list()
beta_list <- list()
gelman_list <- list()

fg_names <- fg_names[!fg_names %in% c("d_galacturonicacid_simple")]
fg_names <- fg_names[!fg_names %in% c("lowcarbon_stress")]
fg_names <- fg_names[55:66]
fg_names <- fg_names[39:66]

rank.name <- "cellulolytic"
# Loop through all ranks
fg_names <- fg_names[!fg_names %in% c("copiotroph","benomyl_antibiotic","nonfsfeso4_anaerobic","gentamycin_antibiotic","nystatin_antibiotic","glycerol_simple","acetate_simple")]


fg_names <- c("endophyte", "plant_pathogen", "animal_pathogen", "ectomycorrhizal", 
							"lichenized", "wood_saprotroph", "soil_saprotroph", "litter_saprotroph", 
							"saprotroph")

scenario <- "no_uncertainty"
scenario <- "temporal_uncertainty"
scenario <- "spatial_uncertainty"
out_scenarios <- list()
missing_scenarios <- matrix(ncol = 2, dimnames = list(NULL, c("scenario", "rank.name")))
for (scenario in c("no_uncertainty","temporal_uncertainty","spatial_uncertainty","full_uncertainty")) {
#for (scenario in scenario.list){
for (rank.name in fg_names) {
	# for (rank.name in c("phylum_fun",
	# 	"class_fun","family_fun","order_fun",
	# 	"genus_fun")){
	
	
	# Subset to one rank.
	rank.df <- d[[rank.name]] 
	
	print(rank.name)
	fg_num <- which(fg_names==rank.name)


	# Prep model inputs/outputs.
	model.dat <- prepModelData(rank.df = rank.df, min.prev = 3)
	
	# taxon name-number key.
	taxon_key <- colnames(model.dat$y)
	names(taxon_key) <- seq(1, length(taxon_key))
	taxon_key2 <- 1:length(colnames(model.dat$y))
	names(taxon_key2) <- colnames(model.dat$y)
	
	# site name-number key.
	site_key <- unique(model.dat$siteID)
	names(site_key) <- 1:length(site_key)
	
	
	
	group.model.path <- paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/samples_", rank.name, "_", scenario,".rds")
	
	if (!file.exists(group.model.path)) {
		missing_scenarios <- rbind(missing_scenarios, c(rank.name, scenario))
		cat("missing: ",  c(rank.name, scenario))
		next()
	}
	
	read_in <- readRDS(group.model.path)
	
	# Read in samples for visualization
	full_summary <- read_in[[2]]
	means <- full_summary[[2]]
	
	# Calculate plot means and actual values per rank
	plot_rel_out <- read_in$plot_summary[[1]]
	pred.plot <- plot_rel_out %>% as.data.frame() %>%  rownames_to_column() %>%
		separate(rowname, sep=", ", into=c("plot_num","date_num")) %>%
		mutate(plot_num = as.numeric(gsub("plot_mu\\[", "", plot_num)),
					 date_num = as.numeric(gsub("\\]", "", date_num))) %>% 
		mutate(taxon = rank.name, scenario = scenario, rank = "functional_group")
	truth.plot.long <- model.dat$truth.plot.long %>% as.data.frame() %>%
		mutate(truth = as.numeric(as.character(truth)),
					 taxon = species)
	truth.plot.long$dates <- as.Date(paste0(truth.plot.long$dateID, "01"), "%Y%m%d")
	allplots <- merge(truth.plot.long, pred.plot, by = c("taxon","plot_num","date_num"), all=T)
	estimate_list[[rank.name]] <- allplots
	
	# Get site effect sizes per rank
	site_eff_out <- means %>% as.data.frame() %>% 
		rownames_to_column("rowname") %>% filter(grepl("site", rowname)) %>% 
		mutate(site_num = as.character(gsub("site_effect\\[|\\]", "", rowname))) %>% 
		mutate(siteID = recode(site_num, !!!site_key),
					 taxon = rank.name, rank = "functional_group")
	site_eff_list[[rank.name]] <- site_eff_out
	
	# Get beta sizes per rank
	beta_out <-  means %>% as.data.frame() %>% 
		rownames_to_column("rowname") %>% filter(grepl("beta", rowname)) %>% 
		mutate(beta_num = as.numeric(gsub("beta\\[|\\]", "", rowname))) %>% 
		mutate(beta = recode(beta_num,
												 "1" = "Temperature",
												 "2" = "Moisture",
												 "3" = "pH",
												 "4" = "pC",
												 "5" = "Plant species richness",
												 "6" = "% grasses",
												 "7" = "% invasive species"),
					 taxon = rank.name, rank = "functional_group")
	beta_list[[rank.name]] <- beta_out
	
	
	## Calculate gelman diagnostics to assess convergence
	samples <- read_in[[1]]
	gd <- gelman.diag(samples, multivariate = FALSE)
	gelman_list[[rank.name]] <- gd
	
}


beta_allranks <- do.call(rbind, beta_list)
site_allranks <- do.call(rbind, site_eff_list)
plot_est_allranks <- do.call(rbind, estimate_list)


beta_allranks$rank <- "functional_group"
beta_allranks$fg_cat <- NA
beta_allranks[grep("simple", beta_allranks$taxon),]$fg_cat <- "Simple substrates"
beta_allranks[grep("complex", beta_allranks$taxon),]$fg_cat <- "Complex substrates"
beta_allranks[grep("stress", beta_allranks$taxon),]$fg_cat <- "Stresses"
beta_allranks[grep("antibiotic", beta_allranks$taxon),]$fg_cat <- "Antibiotic resistance"
beta_allranks[grep("anaerobic", beta_allranks$taxon),]$fg_cat <- "Anaerobic"
beta_allranks[grep("nitr|fixa", beta_allranks$taxon),]$fg_cat <- "N-cycling"
beta_allranks[grep("sapr|path|arbusc|ecto|endo|lichen", beta_allranks$taxon),]$fg_cat <- "Trophic guild"
beta_allranks[grep("copio|oligo", beta_allranks$taxon),]$fg_cat <- "Life-history"
beta_allranks[is.na(beta_allranks$fg_cat),]$fg_cat <- "Other"

out_scenarios[[scenario]] <- list("plot_est" = plot_est_allranks, 
																 "site_eff" = site_allranks, 
																 "beta_eff" = beta_allranks,
																 "gelman" = gelman_list)
}

saveRDS(out_scenarios, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/fg_summaries.rds")



# saveRDS(list("plot_est" = plot_est_allranks, 
# 						 "site_eff" = site_allranks, 
# 						 "beta_eff" = beta_allranks,
# 						 "gelman" = gelman_list), "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/allrank_summaries_10tax.rds")

