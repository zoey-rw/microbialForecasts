library(gganimate)
library(tidyverse)


spatial_biomass_file <- "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_output/spatial_example/sim_diff_3.18_nh4_1.27_no3_0.05_mois_0.47_o2_0.001/biomass.csv"

spatial_biomass_orig <- data.table::fread(spatial_biomass_file,  nThread = 16, data.table = F, col.names = c(NA,"timepoint","x","y","modelID","abundance"))
spatial_biomass_orig$V1 <- NULL




sim_directory <- "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_output/spatial_example/"
scenarios <- list.files(sim_directory, pattern = "^biomass", recursive = T, full.names = T)
all_spatial_biomass_files <- scenarios %>%
	setNames(., sub("\\.csv$", "", dirname(.))) %>%
	map(data.table::fread, nThread = 4, data.table = F, col.names = c(NA,"timepoint","x","y","species","abundance"))
spatial_biomass_orig <- map_df(all_spatial_biomass_files, rbind, .id = "name")
spatial_biomass_orig$V1 <- NULL
spatial_biomass_orig$name = basename(spatial_biomass_orig$name)



biomass_orig = spatial_biomass_orig


rock_locs <- vector("list", length = max(biomass_orig$x))
for (x in 1:max(biomass_orig$x)){
	for (y in 1:max(biomass_orig$y)) {
		loc <- biomass_orig %>% filter(x == !!x & y == !!y)
		if (nrow(loc) == 0){
			rock_locs[[x]][[y]] <- cbind(0, x, y, "rock", 0)
		}
	}
	rock_locs[[x]] <- do.call(rbind, rock_locs[[x]])
}
rock_locs_all <- as.data.frame(do.call(rbind, rock_locs))
colnames(rock_locs_all) <- c("timepoint", "x", "y", "modelID", "abundance")
rock_locs_all$microbe = "rocks (barriers to diffusion)"
rock_locs_all <- rock_locs_all %>% mutate(abundance = as.numeric(abundance),
																				x = as.numeric(x),
																				y = as.numeric(y),
																				timepoint = as.numeric(timepoint))


biomass_orig <- data.table::rbindlist(list(biomass_orig, rock_locs_all), fill=T)

timestep <- 1
biomass_orig$hour <- biomass_orig$cycle*timestep
biomass_orig$id <- biomass_orig$name
biomass_orig <- biomass_orig %>%
	separate(col = "id", sep = "_n|_m|_o", into = c("diffusion","ammonium",'nitrate',"mois","o2"), remove=F)
biomass_orig$diffusion <- as.numeric(gsub("sim_diff_", "", biomass_orig$diffusion, fixed=T))
biomass_orig$ammonium <- gsub("h4_", "", biomass_orig$ammonium, fixed=T)
biomass_orig$nitrate <- gsub("o3_", "", biomass_orig$nitrate, fixed=T)
biomass_orig$mois <- gsub("ois_", "", biomass_orig$mois, fixed=T)
biomass_orig$o2 <- gsub("2_", "", biomass_orig$o2, fixed=T)
biomass_orig$scenario_label = paste0("Moisture: ", biomass_orig$mois, ", NH4: ", biomass_orig$ammonium, ", NO3: ", biomass_orig$nitrate, ", O2: ", biomass_orig$o2)


biomass_orig <- biomass_orig %>% filter(scenario_label %in% c("Moisture: 0.16, NH4: 1.25, NO3: 0.06, O2: 0.169", "Moisture: 0.47, NH4: 1.27, NO3: 0.05, O2: 0.001"))

nitrospira <- biomass_orig %>% filter(species == "not_metanet_Nitrospira_NOB")  %>% mutate(microbe = "N. defluvii \n(nitrite-oxidizer)")
nitrosomonas <- biomass_orig %>% filter(species == "not_metanet_N_europaea_AOB")  %>% mutate(microbe = "N. europaea \n(ammonia-oxidizer)")
winogradskyi <- biomass_orig %>% filter(species == "not_metanet_N_winogradskyi_NOB")  %>% mutate(microbe = "N. winogradskyi \n(ammonia-oxidizer)")
geobacter <- biomass_orig %>% filter(species == "not_metanet_iAF987") %>% mutate(microbe = "G. metallireducens \n(metal-reducer)")
methano <- biomass_orig %>% filter(species == "not_metanet_iAF692") %>% mutate(microbe = "M. sarcina \n(methanogen)")



to_plot <- data.table::rbindlist(list(methano, nitrosomonas, nitrospira, winogradskyi, geobacter))
rock_locs_expanded <- data.table::rbindlist(list(
	cbind.data.frame(rock_locs_all, scenario_label = "Moisture: 0.16, NH4: 1.25, NO3: 0.06, O2: 0.169"),
	cbind.data.frame(rock_locs_all, scenario_label = "Moisture: 0.47, NH4: 1.27, NO3: 0.05, O2: 0.001")))
rock_locs_expanded <- rock_locs_expanded %>% group_by(scenario_label, x, y) %>%
	expand(timepoint = unique(to_plot$timepoint)) %>%
	expand(microbe = unique(to_plot$microbe))

to_plot <- data.table::rbindlist(list(to_plot, rock_locs_expanded), fill=T)

#to_plot <- data.table::rbindlist(list(methano, nitrosomonas, rock_locs_all), fill=T)

#to_plot$depth = -to_plot$x




ggplot(to_plot %>%
			 	filter(timepoint == 600)) +
	geom_point(data = rock_locs_expanded,
						 aes(x=as.numeric(x),y=as.numeric(y)),
						 fill = "grey",
						 color = "grey",
						 shape = 22, size = 1) +
	geom_point(#data = methano,
						 aes(x=x,y=y,
				 fill = abundance,
				 #fill = modelID,
				 #alpha = abundance,
				 color = abundance
				 ),
		 shape = 22,
		size=1) +
	scale_y_reverse(limits = c(100, 0)) +
	#, position = position_jitter(width=.0001)) +
	#theme_dark() +
	xlab("") + ylab("Soil depth (.1 cm)") +
	scale_fill_gradient(low = "black", high = "red") +
	scale_color_gradient(low = "black", high = "red") +
	guides(size="none") +
	theme(panel.background = element_rect(fill = "black", colour = "black"),
				panel.grid.major = element_line(colour = "black"),
				panel.grid.minor = element_line(colour = "black"),
				strip.text.x = element_text(size = 7)) +
	facet_grid(scenario_label~microbe)


methano_plot =
	ggplot(to_plot %>%
				 	filter(timepoint == 600 & species=="not_metanet_iAF692")) +
	geom_point(data = rock_locs_all,
						 aes(x=as.numeric(x),y=as.numeric(y)),
						 fill = "grey",
						 color = "grey",
						 shape = 22, size = 4, inherit.aes = F) +
	geom_point(#data = methano,
		aes(x=x,y=y,
				fill = abundance,
				color = abundance),
		shape = 22, alpha = .5,
		size=4) +
	scale_y_reverse(limits = c(100, 0)) +
	#, position = position_jitter(width=.0001)) +
	#theme_dark() +
	xlab("") + ylab("Soil depth (.1 cm)") +
	scale_fill_gradient(low = "black", high = "skyblue") +
	scale_color_gradient(low = "black", high = "skyblue") +
	guides(size="none") +
	theme(panel.background = element_rect(fill = "black", colour = "black"),
				panel.grid.major = element_line(colour = "black"),
				panel.grid.minor = element_line(colour = "black")) + facet_grid(~scenario_label)
methano_plot



neuro_plot =
	ggplot(to_plot %>%
				 	filter(timepoint == 600 & species=="not_metanet_N_europaea_AOB")) +
	geom_point(data = rock_locs_all,
						 aes(x=as.numeric(x),y=as.numeric(y)),
						 fill = "grey",
						 color = "grey",
						 shape = 22, size = 4, inherit.aes = F) +
	geom_point(#data = methano,
		aes(x=x,y=y,
				fill = abundance,
				color = abundance),
		shape = 22, #alpha = .5,
		size=4) +
	scale_y_reverse(limits = c(100, 0)) +
	#, position = position_jitter(width=.0001)) +
	#theme_dark() +
	xlab("") + ylab("Soil depth (.1 cm)") +
	scale_fill_gradient(low = "black", high = "red") +
	scale_color_gradient(low = "black", high = "red") +
	guides(size="none", fill="none") +
	theme(panel.background = element_rect(fill = "black", colour = "black"),
				panel.grid.major = element_line(colour = "black"),
				panel.grid.minor = element_line(colour = "black"))  + facet_grid(~name)
neuro_plot

library(ggpubr)
ggarrange(neuro_plot, methano_plot)

anim <- ggplot(nitrospira %>%
							 	#filter(timepoint == 300),
							 aes(x=x,y=y,
											 #size = abundance,
											 fill = abundance,
											 color = abundance
)) +
	geom_point(data = rock_locs_all,
						 aes(x=as.numeric(x),y=as.numeric(y)),
						 fill = "grey",
						 color = "grey",
						 shape = 22, size = 4, inherit.aes = F) +
	geom_point(#aes(x=x,y=y,
						# 		 #size = abundance,
						# 		 fill = abundance,
						# 		 color = abundance
						# 		 ),
						#  shape = 22,
						 size=4) +

	scale_y_reverse(limits = c(100, 0)) +
	#, position = position_jitter(width=.0001)) +
	#theme_dark() +
	xlab("") + ylab("Soil depth (.1 cm)") +
	scale_fill_gradient(low = "black", high = "red") +
	scale_color_gradient(low = "black", high = "red") +
	guides(size="none") +
	theme(panel.background = element_rect(fill = "black", colour = "black"),
				panel.grid.major = element_line(colour = "black"),
				panel.grid.minor = element_line(colour = "black")) +
	# Here comes the gganimate code
	transition_states(
		timepoint,
		transition_length = .1,
		state_length = .1
	) + labs(title = 'Hour: {closest_state}')




animate(anim, height = 8, width = 3, units = "in", res = 150)

anim_save("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/spatial_gif_v1.gif")






to_plot <- data.table::rbindlist(list(methano, nitrospira, rock_locs_all), fill=T)


nitro_plot =
	ggplot(to_plot %>%
				 	filter(timepoint == 600 & modelID=="not_metanet_Nitrospira_NOB")) +
	geom_point(data = rock_locs_all,
						 aes(x=as.numeric(x),y=as.numeric(y)),
						 fill = "grey",
						 color = "grey",
						 shape = 22, size = 4, inherit.aes = F) +
	geom_point(#data = methano,
		aes(x=x,y=y,
				fill = abundance,
				color = abundance),
		shape = 22, #alpha = .5,
		size=4) +
	scale_y_reverse(limits = c(100, 0)) +
	#, position = position_jitter(width=.0001)) +
	#theme_dark() +
	xlab("") + ylab("Soil depth (.1 cm)") +
	scale_fill_gradient(low = "black", high = "red") +
	scale_color_gradient(low = "black", high = "red") +
	guides(size="none", fill="none") +
	theme(panel.background = element_rect(fill = "black", colour = "black"),
				panel.grid.major = element_line(colour = "black"),
				panel.grid.minor = element_line(colour = "black")) + facet_grid(~scenario_label)
nitro_plot



nitro_plot =
	ggplot(winogradskyi %>%
				 	filter(timepoint == 600)) +
	geom_point(data = rock_locs_all,
						 aes(x=as.numeric(x),y=as.numeric(y)),
						 fill = "grey",
						 color = "grey",
						 shape = 22, size = 4, inherit.aes = F) +
	geom_point(#data = methano,
		aes(x=x,y=y,
				fill = abundance,
				color = abundance),
		shape = 22, #alpha = .5,
		size=4) +
	scale_y_reverse(limits = c(100, 0)) +
	#, position = position_jitter(width=.0001)) +
	#theme_dark() +
	xlab("") + ylab("Soil depth (.1 cm)") +
	scale_fill_gradient(low = "black", high = "red") +
	scale_color_gradient(low = "black", high = "red") +
	guides(size="none", fill="none") +
	theme(panel.background = element_rect(fill = "black", colour = "black"),
				panel.grid.major = element_line(colour = "black"),
				panel.grid.minor = element_line(colour = "black"))  + facet_grid(~scenario_label)
nitro_plot




media_spatial = data.table::fread("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_output/spatial_example/sim_diff_3.18_nh4_1.27_no3_0.05_mois_0.47_o2_0.001/media.csv", nThread = 16)
media_spatial$name = "sim_diff_3.18_nh4_1.27_no3_0.05_mois_0.47_o2_0.001"

ggplot(media_spatial %>%
			 	filter(cycle == 600 & metabolite=="o2_e")) +

	geom_point(
		aes(x=x,y=y,
				fill = conc_mmol,
				color = conc_mmol),
		shape = 22, #alpha = .5,
		size=4) +
	scale_y_reverse(limits = c(100, 0)) +
	#, position = position_jitter(width=.0001)) +
	#theme_dark() +
	xlab("") + ylab("Soil depth (.1 cm)") +
	scale_fill_gradient(low = "black", high = "blue") +
	scale_color_gradient(low = "black", high = "blue") +
	guides(size="none", fill="none") +
	theme(panel.background = element_rect(fill = "black", colour = "black"),
				panel.grid.major = element_line(colour = "black"),
				panel.grid.minor = element_line(colour = "black"))





anim_neuro = 	ggplot(to_plot %>%
										 	filter(species=="not_metanet_N_europaea_AOB")) +
	geom_point(#data = methano,
		aes(x=x,y=y,
				fill = abundance,
				color = abundance),
		shape = 22, #alpha = .5,
		size=4) +
	geom_point(data = rock_locs_all,
						 aes(x=as.numeric(x),y=as.numeric(y)),
						 fill = "grey",
						 color = "grey",
						 shape = 22, size = 4, inherit.aes = F) +

	scale_y_reverse(limits = c(100, 0)) +
	#, position = position_jitter(width=.0001)) +
	#theme_dark() +
	xlab("") + ylab("Soil depth (.1 cm)") +
	scale_fill_gradient(low = "black", high = "red") +
	scale_color_gradient(low = "black", high = "red") +
	guides(size="none", fill="none") +
	theme(panel.background = element_rect(fill = "black", colour = "black"),
				panel.grid.major = element_line(colour = "black"),
				panel.grid.minor = element_line(colour = "black"),
				strip.text.x = element_text(size = 10))  + facet_grid(~scenario_label) + 	# Here comes the gganimate code
	transition_states(
		timepoint,
		transition_length = .01,
		state_length = .01
	) + labs(title = "N. europaea abundance",
		subtitle = 'Hour: {closest_state}')
animate(anim_neuro, height = 7, width = 4, units = "in", res = 100)
anim_save("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/neuro_gif_v1.gif")





anim_methano = 	ggplot(to_plot %>%
										 	filter(species=="not_metanet_iAF692")) +
	geom_point(#data = methano,
		aes(x=x,y=y,
				fill = abundance,
				color = abundance),
		shape = 22, #alpha = .5,
		size=2) +
	geom_point(data = rock_locs_all,
						 aes(x=as.numeric(x),y=as.numeric(y)),
						 fill = "grey",
						 color = "grey",
						 shape = 22, size = 2, inherit.aes = F) +

	scale_y_reverse(limits = c(100, 0)) +
	#, position = position_jitter(width=.0001)) +
	#theme_dark() +
	xlab("") + ylab("Soil depth (.1 cm)") +
	scale_fill_gradient(low = "black", high = "skyblue") +
	scale_color_gradient(low = "black", high = "skyblue") +
	guides(size="none", fill="none") +
	theme(panel.background = element_rect(fill = "black", colour = "black"),
				panel.grid.major = element_line(colour = "black"),
				panel.grid.minor = element_line(colour = "black"),
				strip.text.x = element_text(size = 10))  + facet_grid(~scenario_label) + 	# Here comes the gganimate code
	transition_states(
		timepoint,
		transition_length = .01,
		state_length = .01
	) + labs(title = "M. sarcina abundance",
					 subtitle = 'Hour: {closest_state}')
animate(anim_methano, height = 7, width = 4, units = "in", res = 100)
anim_save("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/methano_gif_v1.gif")






full_anim = ggplot(to_plot) +
	geom_point(data = rock_locs_expanded,
						 aes(x=as.numeric(x),y=as.numeric(y)),
						 fill = "grey",
						 color = "grey",
						 shape = 22, size = 1) +
	geom_point(#data = methano,
		aes(x=x,y=y,
				fill = abundance,
				#fill = modelID,
				#alpha = abundance,
				color = abundance
		),
		shape = 22,
		size=1) +
	scale_y_reverse(limits = c(100, 0)) +
	#, position = position_jitter(width=.0001)) +
	#theme_dark() +
	xlab("") + ylab("Soil depth (.1 cm)") +
	scale_fill_gradient(low = "black", high = "red") +
	scale_color_gradient(low = "black", high = "red") +
	guides(size="none") +
	theme(panel.background = element_rect(fill = "black", colour = "black"),
				panel.grid.major = element_line(colour = "black"),
				panel.grid.minor = element_line(colour = "black"),
				strip.text.x = element_text(size = 7)) +
	facet_grid(scenario_label~microbe) + 	# Here comes the gganimate code
	transition_states(
		timepoint,
		transition_length = .01,
		state_length = .01
	) + labs(title = 'Hour: {closest_state}')
animate(full_anim, height = 8, width = 6, units = "in", res = 100)
anim_save("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/full_sim_gif_v1.gif")
