# helper scripts for n cycle sims
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

parseSimID <- function(df) {

	df <- df %>%
		separate(col = "id", sep = "_n|_m|_o", into = c("diffusion","ammonium",'nitrate',"mois","o2"), remove=F)
	df$diffusion <- as.numeric(gsub("sim_diff_", "", df$diffusion, fixed=T))
	df$ammonium <- gsub("h4_", "", df$ammonium, fixed=T)
	df$nitrate <- gsub("o3_", "", df$nitrate, fixed=T)
	df$mois <- gsub("ois_", "", df$mois, fixed=T)
	df$o2 <- gsub("2_", "", df$o2, fixed=T)
	df$scenario_label = paste0("Moisture: ", df$mois, ", NH4: ", df$ammonium, ", NO3: ", df$nitrate, ", O2: ", df$o2)

	return(df)
	}

