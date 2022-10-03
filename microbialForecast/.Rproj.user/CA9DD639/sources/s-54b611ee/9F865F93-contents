library(tidyverse)
library(ggrepel)
library(ggpubr)
options(scipen=999)

in_flux = readRDS("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/static_o2_fluxes.rds")
in_flux = readRDS("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/n_only_refresh_nh4_fluxes.rds")

in_flux = readRDS("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/n_only_fluxes.rds")

in_flux = readRDS("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/extramods_fluxes.rds")
in_flux = readRDS("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/more_sugars_fluxes.rds")
in_flux = readRDS("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/no_static_fluxes.rds")
in_flux = readRDS("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/static_o2_new_fluxes.rds")
in_flux = readRDS("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/refresh_o2_fluxes.rds")
in_flux = readRDS("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/initial_o2_fluxes.rds")
in_flux = readRDS("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/unlimited_fluxes.rds")

#for_eval = list(eval_out_df,lm_summary,met_by_timepoint,all_fluxes)

in_flux_lm_summary = in_flux[[2]]
eval_out_df = in_flux[[1]]

in_flux_lm_summary <- in_flux_lm_summary[,-c(7:8)]

best_times = in_flux_lm_summary %>% filter(nmin_slope > 0 & nitr_slope > 0) #%>% select(init_time, final_time)

if(nrow(best_times)==0) best_times = in_flux_lm_summary %>% filter(nitr_r.squared == max(nitr_r.squared))

best_times$time_set = paste0(best_times$init_time, "_", best_times$final_time)
best_times


eval_out_df$time_set = paste0(eval_out_df$init_time, "_", eval_out_df$final_time)
eval_out_df_best = eval_out_df %>% filter(eval_out_df$time_set %in% best_times$time_set)

ggplot(eval_out_df_best,
			 aes(y =nmin_obs, x = nmin_pred, color = time_set)) +
	geom_point(size = 1, alpha=.8) + geom_smooth(method = "lm") + theme_bw()

ggplot(eval_out_df_best,
			 aes(y =nitr_obs, x = nitr_pred, color = time_set)) +
	geom_point(size = 1, alpha=.8) + geom_smooth(method = "lm") + theme_bw()




eval_out_df_best %>% filter(time_set=="0_200") %>%
	ggplot(., aes(y =nitr_obs, x = nitr_pred)) +
	geom_point(size = 3, alpha=.5) + #geom_smooth(method = "lm") +
	theme_bw()  +
	stat_smooth(method = "lm") +
	stat_regline_equation(
		aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
		 label.x = 0, label.y = .4)


eval_out_df_best %>% filter(time_set=="0_200") %>%
	ggplot(., aes(y =nmin_obs, x = nmin_pred)) +
	geom_point(size = 3, alpha=.5) + #geom_smooth(method = "lm") +
	theme_bw()  +
	stat_smooth(method = "lm") +
	stat_regline_equation(
		aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
		label.x = 0)



mets_to_plot = in_flux[[3]]

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




flux_files = list.files("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/", pattern = "fluxes", full.names = T)
flux_files_in_list = flux_files %>%
	setNames(., sub("\\.rds$", "", basename(.))) %>%
	map(readRDS)

in_flux_preds <- lapply(flux_files_in_list, '[[', 1)

all_flux_scenarios <- map_df(in_flux_preds, rbind, .id = "name")

pos_predictions = all_flux_scenarios %>% filter(nmin_pred > 0.01)

