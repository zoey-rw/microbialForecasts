# Loop through output directories and summarize final biomasses by timepoint
library(tidyverse)
source("/projectnb/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/source.r")

simulation_types <- list.files(
	"/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_output", include.dirs = T,
	full.names = T)

sim_type = "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_output/n_only_refresh_nh4"


timestep <- 1

for (sim_type in simulation_types) {
	if (basename(sim_type) %in% c("spatial_o2","spatial_o2_sensormois")) next()

	print(sim_type)
	scenarios <- list.files(sim_type, pattern = "^total_biomass", recursive = T, full.names = T)

	all_biomass_files <- scenarios %>%
		setNames(., sub("\\.csv$", "", dirname(.))) %>%
		map(data.table::fread, nThread = 4)

	biomass <- map_df(all_biomass_files, rbind, .id = "id")
	biomass$id = basename(biomass$id)
	biomass$V1=NULL
	idkey = cbind.data.frame(new = c("Nitrosomonas_europaea_AOB",	"Nitrosomonas_eutropha_AOB",
																	 "Nitrosospira_multiformis_AOB",	"Nitrosococcus_oceani_AOB",
																	 "Nitrospira_defluvii_NOB",	"Nitrobacter_winogradskyi_NOB",
																	 "Nitrobacter_hamburgensis_NOB",	"Nitrospina_gracilis_NOB",
																	 "Geobacter_metallireducens",	"Methanosarcina_barkeri",
																	 "Salmonella_typhimurium",	"Escherichia_coli",	"Pseudomonas_stutzeri",
																	 "Pseudomonas_putida",	"Nitrospira_moscoviensis_NOB"),
													 original = c("not_metanet_N_europaea_AOB", "not_metanet_N_eutropha_AOB",
													 						 "not_metanet_N_multiformis_AOB", "not_metanet_N_oceani_AOB",
													 						 "not_metanet_N_defluvii_NOB", "not_metanet_N_winogradskyi_NOB",
													 						 "not_metanet_N_hamburgensis_NOB", "not_metanet_N_gracilis_NOB",
													 						 "not_metanet_iAF987",
													 						 "not_metanet_iAF692", "not_metanet_iRR1083", "not_metanet_iJO1366",
													 						 "not_metanet_iPB890", "not_metanet_iJN746",
													 						 "not_metanet_Nitrospira_NOB"))

	colnames(biomass) <- dplyr::recode(colnames(biomass), !!!setNames(as.character(idkey$new), idkey$original))



	biomass_long <- biomass %>% pivot_longer(cols = 3:ncol(biomass), names_to = "species")
	biomass_long$hour <- biomass_long$cycle*timestep

	biomass_long <- parseSimID(biomass_long)

	#biomass_long <- biomass_long %>%
	#	separate(col = "id", sep = "_n|_m|_o", into = c("diffusion","ammonium",'nitrate',"mois","o2"), remove=F)
	#biomass_long$diffusion <- as.numeric(gsub("sim_diff_", "", biomass_long$diffusion, fixed=T))
	#biomass_long$ammonium <- gsub("h4_", "", biomass_long$ammonium, fixed=T)
	#biomass_long$nitrate <- gsub("o3_", "", biomass_long$nitrate, fixed=T)
	#biomass_long$mois <- gsub("ois_", "", biomass_long$mois, fixed=T)
	#biomass_long$o2 <- gsub("2_", "", biomass_long$o2, fixed=T)
	#biomass_long$scenario_label = paste0("Moisture: ", biomass_long$mois, ", NH4: ",
	#																		 biomass_long$ammonium, ", NO3: ", biomass_long$nitrate, ", O2: ",
	#																		 biomass_long$o2)
	biomass_long$N_cycler <- ifelse(grepl("_NOB|_AOB", biomass_long$species), T, F)

	head(biomass_long)

	out_path = file.path("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary",
											 paste0(basename(sim_type), "_total_biomass.rds"))
	saveRDS(biomass_long, out_path)
	message("Saved to: ", out_path)
}


