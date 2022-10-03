
# Compare horizon preferences

# First run "predict_bacterial_preferences"

library(tidyverse)

biomass_files = list.files("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/", pattern = "total", full.names = T)
biomass_files_in_list = biomass_files %>%
	setNames(., sub("\\.rds$", "", basename(.))) %>%
	map(readRDS)

all_biomass_scenarios <- map_df(biomass_files_in_list, rbind, .id = "simulation_type")
all_biomass_scenarios$N_cycler=NULL


all_sim_biomass = all_biomass_scenarios %>% group_by(simulation_type, id) %>%
	mutate(total_biomass = sum(value, na.rm=T)) %>% ungroup() %>%
	mutate(rel_abun = value/total_biomass)

sim_biomass_wide <- all_sim_biomass %>% pivot_wider(values_from = c(rel_abun, value), names_from = species)
sim_biomass_wide <- sim_biomass_wide %>% #filter(cycle==300) %>%
	group_by(simulation_type,id,hour) %>%
	mutate(rel_NOB =
				 	sum(rel_abun_Nitrospira_defluvii_NOB,
				 			rel_abun_Nitrospina_gracilis_NOB,
				 			rel_abun_Nitrobacter_hamburgensis_NOB,
				 			rel_abun_Nitrobacter_winogradskyi_NOB,
				 			rel_abun_Nitrospira_moscoviensis_NOB, na.rm = T),
				 rel_AOB = sum(rel_abun_Nitrosomonas_europaea_AOB,
				 							rel_abun_Nitrosomonas_eutropha_AOB,
				 							rel_abun_Nitrosospira_multiformis_AOB,
				 							rel_abun_Nitrosococcus_oceani_AOB, na.rm = T),
				 rel_AOB_NOB = sum(rel_AOB, rel_NOB, na.rm = T),
				 rel_Pseudomonas = sum(rel_abun_Pseudomonas_stutzeri,rel_abun_Pseudomonas_putida, na.rm = T),
				 rel_Nitrospira = sum(rel_abun_Nitrospira_moscoviensis_NOB,
				 										 rel_abun_Nitrospira_defluvii_NOB, na.rm = T),
				 rel_Nitrosomonas = sum(rel_abun_Nitrosomonas_eutropha_AOB,
				 											 rel_abun_Nitrosomonas_europaea_AOB, na.rm = T)
	) %>% as.data.frame()




eval_df <- read.csv("/projectnb/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/true_rates_df.csv")
sim_biomass_wide <- merge(sim_biomass_wide, eval_df[,c("sampleID","id")])
sim_biomass_wide <- sim_biomass_wide %>% mutate(siteID = substr(sampleID, 1, 4)) %>%
	mutate(dates = substr(sampleID, 20, 27)) %>%
	mutate(asDate = as.Date(as.character(dates), "%Y%m%d")) %>%
	mutate(dateID = substr(as.character(dates), 1, 6)) %>%
	mutate(plotID = substr(sampleID, 1, 8)) %>%
	mutate(site_date = paste0(siteID, "-", dateID)) %>%
	mutate(horizon = ifelse(grepl("-M-", sampleID), "M", "O")) %>%
	mutate(plot_date = paste0(plotID, "-", dateID)) %>%
	as.data.frame()



to_compare_horizons <- sim_biomass_wide %>%
	select(simulation_type, id, scenario_label, horizon, hour, plot_date,
				 rel_Nitrosomonas, rel_abun_Geobacter_metallireducens, rel_Pseudomonas, rel_Nitrospira) %>%
	pivot_longer(cols=
							 	c(rel_Nitrosomonas, rel_abun_Geobacter_metallireducens, rel_Pseudomonas, rel_Nitrospira),
							 names_to = "name") %>%
	group_by(simulation_type,hour,name) %>%
	mutate(abundance_scaled = scale(value)) %>%
	as.data.frame()



to_compare_horizons <- to_compare_horizons %>% mutate(genus = recode(name, "rel_Nitrospira" = "Nitrospira",
																																		 "rel_abun_Geobacter_metallireducens" = "Geobacter",
																																		 "rel_Pseudomonas" = "Pseudomonas",
																																		 "rel_Nitrosomonas" = "Nitrosomonas"),
																											data_type = "Predicted")

ggplot(to_compare_horizons %>% filter(hour == "200" & simulation_type %in%
															c("extramods_total_biomass","more_sugars_total_biomass","initial_o2_total_biomass","unlimited_total_biomass",NA))) +
	geom_boxplot(aes(x = abundance_scaled,
									 y = horizon, color = horizon)) +
	theme_minimal(base_size = 18) +
	facet_grid(genus~simulation_type)+
	#, scales="free") +
	#geom_hline(yintercept = 0) +
	geom_vline(xintercept = 0) #+ xlim(c(-2,))


meta_gen <- meta_gen %>% mutate(genus = recode(name, "g_Nitrospira" = "Nitrospira",
																							 "g_Geobacter" = "Geobacter",
																							 "g_Pseudomonas" = "Pseudomonas",
																							 "g_Nitrosomonas" = "Nitrosomonas"),
																data_type = "Observed",
																abundance_scaled = `scaled_percent`)


master_df <- data.table::rbindlist(list(to_compare_horizons %>%
																					filter(simulation_type %in%
																								 	c("initial_o2_total_biomass","more_sugars_total_biomass","unlimited_total_biomass",NA)),
																				meta_gen), fill = T)


ggplot(master_df %>% filter(hour == "200" & simulation_type %in%
															c("extramods_total_biomass","more_sugars_total_biomass","initial_o2_total_biomass","unlimited_total_biomass",NA))) +
	geom_boxplot(aes(x = abundance_scaled,
									 y = horizon,
									 color = horizon,
									 group=data_type),
							 position = "dodge") +
	theme_minimal(base_size = 18) +
	facet_grid(simulation_type~genus)+
	#, scales="free") +
	#geom_hline(yintercept = 0) +
	geom_vline(xintercept = 0) #+ xlim(c(-2,))


to_plot <- master_df %>% filter(hour == "200"| is.na(hour) & simulation_type %in% c("unlimited_total_biomass",NA))
ggplot(to_plot) +
	geom_boxplot(aes(x = abundance_scaled,
									 y = horizon,
									 fill = horizon),
							 position = "dodge",
							 show.legend = F) +
	geom_point(aes(x = abundance_scaled,
								 y = horizon,
								 fill=horizon),
						 position = position_jitterdodge(jitter.width = .3), alpha = .3,
						 show.legend = F) +
	theme_minimal(base_size = 18) +
	facet_grid(data_type~genus#, scales="free"
	) +
	xlab("Abundance (mean-centered)") +
	ylab("Soil horizon") +
	geom_vline(xintercept = 0) #+ xlim(c(-2,))
