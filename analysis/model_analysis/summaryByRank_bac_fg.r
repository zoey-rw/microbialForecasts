library(ggplot2)
library(scales)
library(viridis)
library(hrbrthemes)
library(coda)

source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/prepDirichletData.r")
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
rank.name <- "class_bac"

summary_list <- list()
taxon_key_list <- list()

# Have to read in *again* because I can't figure out how to pass data into function...
d1 <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_16S_2021.rds")
d2 <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/val_groupAbundances_16S_2021.rds")
ranks.keep <- names(d1)
ranks.keep <- ranks.keep[!grepl("bac", ranks.keep)]

in.path <- paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/samples_bac_fg_min3.rds")
read_in <- readRDS(in.path)

# Create lists for outputs.
estimate_list <- list()
site_eff_list <- list()
beta_list <- list()
gelman_list <- list()

rank.name <- "cellulolytic"
# Loop through all ranks
for (rank.name in ranks.keep) {
	# for (rank.name in c("phylum_fun",
	# 	"class_fun","family_fun","order_fun",
	# 	"genus_fun")){
	
	rank.df <- d[[rank.name]] 
	
	k <- which(ranks.keep==rank.name)
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
	
	# Read in samples for visualization
	group_in <- read_in[[k]]
	full_summary <- group_in[[2]]
	means <- full_summary[[2]]
	
	# Calculate plot means and actual values per rank
	plot_rel_out <- means[grep("plot_rel", rownames(means)),]
	pred.plot <- plot_rel_out %>% as.data.frame() %>%  rownames_to_column() %>%
		separate(rowname, sep=", ", into=c("plot_num","species_num","date_num")) %>%
		mutate(plot_num = as.numeric(gsub("plot_rel\\[", "", plot_num)),
					 date_num = as.numeric(gsub("\\]", "", date_num)),
					 species_num = as.numeric(species_num)) %>% 
		mutate(taxon = recode(species_num,
													!!!taxon_key),
					 rank = rank.name)
	truth.plot.long <- model.dat$plot.truth %>% as.data.frame() %>%
		mutate(truth = as.numeric(as.character(truth)),
					 taxon = species)
	truth.plot.long$dates <- as.Date(paste0(truth.plot.long$dateID, "01"), "%Y%m%d")
	allplots <- merge(truth.plot.long, pred.plot, by = c("taxon","plot_num","date_num"), all=T)
	estimate_list[[rank.name]] <- allplots
	
	# Get site effect sizes per rank
	site_eff_out <- means %>% as.data.frame() %>% 
		rownames_to_column("rowname") %>% filter(grepl("site", rowname)) %>% 
		separate(rowname, sep=", ", into=c("site_num","taxon_num")) %>% 
		mutate(site_num = as.character(gsub("site_effect\\[", "", site_num)),
					 taxon_num = as.character(gsub("\\]", "", taxon_num))) %>% 
		mutate(siteID = recode(site_num, !!!site_key),
					 taxon = recode(taxon_num, !!!taxon_key),
					 rank = rank.name)
	site_eff_list[[rank.name]] <- site_eff_out
	
	# Get beta sizes per rank
	beta_out <-  means %>% as.data.frame() %>% 
		rownames_to_column("rowname") %>% filter(grepl("beta", rowname)) %>% 
		separate(rowname, sep=", ", into=c("taxon_num","beta_num")) %>% 
		mutate(taxon_num = as.numeric(gsub("beta\\[", "", taxon_num)),
					 beta_num = as.numeric(gsub("\\]", "", beta_num))) %>% 
		mutate(beta = recode(beta_num,
												 "1" = "Temperature",
												 "2" = "Moisture",
												 "3" = "pH",
												 "4" = "pC",
												 "5" = "Plant species richness",
												 "6" = "% grasses",
												 "7" = "% invasive species"),
					 taxon_name = recode(taxon_num, !!!taxon_key),
					 rank = rank.name)
	beta_list[[rank.name]] <- beta_out
	
	
	## Calculate gelman diagnostics to assess convergence
	samples <- group_in[[1]]
	impt_cols <- !grepl("plot_mu|plot_rel|y|Ex", colnames(samples[[1]]))
	gd <- gelman.diag(samples[,impt_cols], multivariate = FALSE)
	gelman_list[[rank.name]] <- gd
	
}

beta_allranks <- do.call(rbind, beta_list)
site_allranks <- do.call(rbind, site_eff_list)
plot_est_allranks <- do.call(rbind, estimate_list)

saveRDS(list("plot_est" = plot_est_allranks, 
						 "site_eff" = site_allranks, 
						 "beta_eff" = beta_allranks,
						 "gelman" = gelman_list), "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/bac_fg_min3_summaries.rds")

