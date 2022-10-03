extramod_biomass = read.table("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/biomass_0x2b598997eb20")


ggplot(extramod_biomass) + geom_line(aes(x = 0, y = V4))
library(tidyverse)


sim_directory <- "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_output/extramods/"
scenarios <- list.files(sim_directory, pattern = "^total_biomass", recursive = T, full.names = T)

scenarios <- scenarios[1:10]

scenarios <- "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_output/extramods/sim_diff_3.18_nh4_0.74_no3_0_mois_0.22_o2_0.144/total_biomass.csv"

# Spatial example file
scenarios <- "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_output/spatial_example/sim_diff_3.18_nh4_1.27_no3_0.05_mois_0.47_o2_0.001/total_biomass.csv"

all_biomass_files <- scenarios %>%
	setNames(., sub("\\.csv$", "", dirname(.))) %>%
	map(data.table::fread, nThread = 4)
biomass <- map_df(all_biomass_files, rbind, .id = "name")
biomass$name = basename(biomass$name)

colnames(biomass) <- c("id","NA","cycle","Nitrosomonas_europaea_AOB",	"Nitrosomonas_eutropha_AOB",	"Nitrosospira_multiformis_AOB",	"Nitrosococcus_oceani_AOB",	"Nitrospira_defluvii_NOB",	"Nitrobacter_winogradskyi_NOB",	"Nitrobacter_hamburgensis_NOB",	"Nitrospina_gracilis_NOB",	"Geobacter_metallireducens",	"Methanosarcina_barkeri",	"Salmonella_typhimurium",	"Escherichia_coli",	"Pseudomonas_stutzeri",	"Pseudomonas_putida",	"Nitrospira_moscoviensis_NOB")

timestep <- 1
biomass_long <- biomass %>% pivot_longer(cols = 4:ncol(biomass))
biomass_long$hour <- biomass_long$cycle*timestep

biomass_long <- biomass_long %>%
	separate(col = "id", sep = "_n|_m|_o", into = c("diffusion","ammonium",'nitrate',"mois","o2"), remove=F)
biomass_long$diffusion <- as.numeric(gsub("sim_diff_", "", biomass_long$diffusion, fixed=T))
biomass_long$ammonium <- gsub("h4_", "", biomass_long$ammonium, fixed=T)
biomass_long$nitrate <- gsub("o3_", "", biomass_long$nitrate, fixed=T)
biomass_long$mois <- gsub("ois_", "", biomass_long$mois, fixed=T)
biomass_long$o2 <- gsub("2_", "", biomass_long$o2, fixed=T)

biomass_long$scenario_label = paste0("Moisture: ", biomass_long$mois, ", NH4: ", biomass_long$ammonium, ", NO3: ", biomass_long$nitrate, ", O2: ", biomass_long$o2)
biomass_long$N_cycler <- ifelse(grepl("_NOB|_AOB", biomass_long$name), T, F)


saveRDS(biomass_long, "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/extramods_total_biomass.rds")
saveRDS(biomass_long, "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/spatial_example_biomass.rds")


library(ggrepel)

# VISUALIZE biomass over time
to_plot = biomass_long %>%
	group_by(id, name, N_cycler) %>%
	mutate(label = ifelse(value == max(value),
												 gsub("_|_NOB|_AOB", " ", as.character(name)),
												 NA_character_)) %>%
	ungroup()

to_plot %>%
	#filter(id %in% c("sim_diff_3.59_nh4_1.16_no3_0_mois_0.29_o2_0.131","sim_diff_3.59_nh4_0.83_no3_0_mois_0.22_o2_0.144")) %>%
	#filter(id %in% c("sim_diff_1.19_nh4_0.25_no3_0.02_mois_0.12_o2_0.122", "sim_diff_1.19_nh4_0.33_no3_0_mois_0.16_o2_0.149")) %>%
	filter(id %in% c(#"sim_diff_3.18_nh4_3.2_no3_0_mois_0.55_o2_0.013",
									 #"sim_diff_3.18_nh4_0.42_no3_0_mois_0.23_o2_0.173",
									 "sim_diff_3.18_nh4_6.08_no3_0_mois_0.47_o2_0.01",
									 "sim_diff_2.11_nh4_0.18_no3_0_mois_0.08_o2_0.01",
									 #"sim_diff_2.11_nh4_2.67_no3_1.93_mois_0.45_o2_0.007",
									 "sim_diff_2.11_nh4_2.65_no3_0.75_mois_0.2_o2_0.122")) %>%
	ggplot(aes(x = hour, y = value, color = as.factor(name))) + geom_line(size = 2, show.legend = F) +
	ylab("Biomass (g)") +
	theme_bw(base_size = 25) + #labs(color="Species") +
	geom_label_repel(aes(label = label), 									 nudge_x = 1,
									 nudge_y=.01,
									 force = .01,
									 na.rm = TRUE, show.legend = F, #x = 120,
									 size=6) +
	facet_grid(cols=vars(scenario_label), rows=vars(N_cycler), scales="free") +
	ggtitle("Competitive outcomes vary with soil conditions")



to_plot %>%
	filter(id %in% c("sim_diff_3.18_nh4_3.2_no3_0_mois_0.55_o2_0.013",
									 #"sim_diff_3.18_nh4_0.42_no3_0_mois_0.23_o2_0.173",
									 #"sim_diff_3.18_nh4_6.08_no3_0_mois_0.47_o2_0.01",
									 "sim_diff_2.11_nh4_0.18_no3_0_mois_0.08_o2_0.01",
									 #"sim_diff_2.11_nh4_2.67_no3_1.93_mois_0.45_o2_0.007",
									 "sim_diff_2.11_nh4_2.65_no3_0.75_mois_0.2_o2_0.122") #& N_cycler
				 ) %>%
	ggplot(aes(x = hour, y = value, color = as.factor(name))) + geom_line(size = 2, show.legend = F) +
	ylab("Biomass (g)") +
	theme_bw(base_size = 22) + #labs(color="Species") +
	geom_label_repel(aes(label = label),
									 nudge_x = 1,
									 #nudge_y=.1,
									 force = .1,
									 na.rm = TRUE,
									 show.legend = F, #x = 120,
									 size=6) +
	facet_wrap(~scenario_label, scales="free", nrow=3) +
	ggtitle("Competitive outcomes vary with soil conditions")



# Spatial example plot
to_plot %>%
	filter(id %in% c("sim_diff_3.18_nh4_1.27_no3_0.05_mois_0.47_o2_0.001") #& N_cycler
	) %>%
	ggplot(aes(x = hour, y = value, color = as.factor(name))) + geom_line(size = 2, show.legend = F) +
	ylab("Biomass (g)") +
	theme_bw(base_size = 22) + #labs(color="Species") +
	geom_label_repel(aes(label = label),
									 nudge_x = 1,
									 #nudge_y=.1,
									 force = .1,
									 na.rm = TRUE,
									 show.legend = F, #x = 120,
									 size=6) +
	facet_wrap(~scenario_label, scales="free", nrow=3)

