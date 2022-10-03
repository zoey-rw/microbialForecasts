

sim_directory <- "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_output/no_static/"
biomass_in <- readRDS("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/no_static_total_biomass.rds")
flux_in = readRDS("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/no_static_fluxes.rds")

flux_in = readRDS("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/n_only_refresh_nh4_fluxes.rds")


flux_in = readRDS("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/static_o2_new_fluxes.rds")


flux_in = readRDS("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/refresh_o2_fluxes.rds")
biomass_in <- readRDS("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/refresh_o2_total_biomass.rds")

#flux_in = in_flux

to_plot = biomass_in %>%
	group_by(id, species, N_cycler) %>%
	mutate(label = ifelse(value == max(value),
												gsub("_|_NOB|_AOB", " ", as.character(species)),
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
	ggplot(aes(x = hour, y = value, color = as.factor(species))) + geom_line(size = 2, show.legend = F) +
	ylab("Biomass (g)") +
	theme_bw(base_size = 25) + #labs(color="Species") +
	geom_label_repel(aes(label = label), 									 nudge_x = 1,
									 nudge_y=.01,
									 force = .01,
									 na.rm = TRUE, show.legend = F, #x = 120,
									 size=6) +
	facet_grid(cols=vars(scenario_label), rows=vars(N_cycler), scales="free") +
	ggtitle("Competitive outcomes vary with soil conditions")





mets_to_plot = flux_in[[3]]

ggplot(mets_to_plot %>% filter(metabolite %in% c("o2_e","nh4_e","no2_e","no3_e")),
			 aes(y =total, x = cycle, color = name,
			 		group=name)) +
	geom_path(size = 1, alpha=.8, show.legend = F) +
	facet_wrap(~metabolite, scales="free") +
	theme_bw(base_size = 18)

ggplot(mets_to_plot %>% filter(metabolite %in% c("no_e")),
			 aes(y =total, x = cycle, color = name,
			 		group=name)) +
	geom_path(size = 1, alpha=.8, show.legend = F) +
	facet_wrap(~metabolite) +
	theme_bw(base_size = 18)



ggplot(mets_to_plot %>% filter(metabolite %in% c("o2_e")),
			 aes(y =total, x = cycle, color = name,
			 		group=name)) +
	geom_path(size = 1, alpha=.8, show.legend = F) +
	facet_wrap(~metabolite) +
	theme_bw(base_size = 18)


flux_in[[1]]$nitr_pred
