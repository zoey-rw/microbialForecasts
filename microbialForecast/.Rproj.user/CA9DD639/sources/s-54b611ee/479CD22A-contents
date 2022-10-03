library(tidyverse)
library(ggrepel)
options(scipen=999)

source("/projectnb/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/source.r")


eval_df <- read.csv("/projectnb/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/true_rates_df.csv")
#eval_df <- eval_df %>% filter(nmin_obs < 4)
#eval_df %>% select(nitr_obs, nmin_obs, id)

sim_directory <- "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_output/static_o2/"

sim_directory <- "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_output/no_static/"
sim_directory <- "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_output/n_only/"
sim_directory <- "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_output/n_only_refresh_nh4/"

sim_directory <- "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_output/more_sugars/"

sim_directory <- "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_output/static_o2_new/"

scenarios <- list.files(sim_directory, pattern = "^media", recursive = T, full.names = T)

scenarios = "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_output/extramods/sim_diff_3.18_nh4_0.74_no3_0_mois_0.22_o2_0.144/media.csv"

sim_directory <- "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_output/extramods/"
#scenarios <- scenarios[1:10]



simulation_types <- list.files(
	"/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_output", include.dirs = T,
	full.names = T)

sim_type = "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_output/n_only_refresh_nh4"
sim_type = "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_output/static_o2_new/"
sim_type = "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_output/refresh_o2/"
sim_type = "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_output/initial_o2/"
sim_type = "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_output/unlimited/"


for (sim_type in simulation_types) {
	if (basename(sim_type) %in% c("spatial_o2","spatial_o2_sensormois")) next()
	if (basename(sim_type) %in% c("static_o2_new","refresh_o2","initial_o2")) {
		timestep <- 2
	} else {
		timestep <- 1
	}

	scenarios <- list.files(sim_type, pattern = "^media", recursive = T, full.names = T)


all_media_files <- scenarios %>%
	setNames(., sub("\\.csv$", "", dirname(.))) %>%
	map(data.table::fread, nThread = 16)
all_media <- map_df(all_media_files, rbind, .id = "name")
all_media$name = basename(all_media$name)

# summarize by timepoint
all_mets <- all_media %>% group_by(name, metabolite, cycle) %>%
	summarize(total = sum(conc_mmol))
mets_keep =  c("nh4_e","no2_e", "no3_e","n2o_e","n2_e","no_e","h2o_e","h_e","h_c","co2_e", "nh3_e","biomass_e","o2_e","ch4_e")
select_mets <- all_mets %>% filter(metabolite %in% mets_keep)
met_by_timepoint <- select_mets %>%
	mutate(inorganic = ifelse(metabolite %in% c("nh4_e", #"nh3_e",
																							"no2_e", "no3_e"), T, F),
				 nitrate_nitrite = ifelse(metabolite %in% c("no2_e", "no3_e"), T, F),
				 gaseous = ifelse(metabolite %in% c("n2_e", "n2o_e", "no_e"), T, F),
				 n2o = ifelse(metabolite %in% c("n2o_e"), T, F),
				 methane = ifelse(metabolite %in% c("ch4_e"), T, F))

# summarize by flux type
fluxes = met_by_timepoint %>% group_by(name, cycle, inorganic) %>%
	mutate(inorganic_n = sum(total)) %>% ungroup() %>%
	group_by(name, cycle, nitrate_nitrite) %>%
	mutate(nitrate_nitrite_n = sum(total)) %>% ungroup() %>%
	group_by(name, cycle, gaseous) %>%
	mutate(gaseous_n = sum(total)) %>% ungroup() %>%
group_by(name, cycle, n2o) %>%
	mutate(n2o_n = sum(total)) %>% ungroup() %>%
group_by(name, cycle, methane) %>%
	mutate(methane_c = sum(total)) %>% ungroup()
nitrate_nitrite_n = fluxes %>%
	filter(nitrate_nitrite ==T) %>%
	select(name, cycle, nitrate_nitrite_n) %>%
	pivot_longer(cols = nitrate_nitrite_n, names_to = "flux")

gaseous_n = fluxes %>%
	filter(gaseous ==T) %>%
	select(name, cycle, gaseous_n) %>%
	pivot_longer(cols = gaseous_n, names_to = "flux")

inorganic_n = fluxes %>%
	filter(inorganic ==T) %>%
	select(name, cycle, inorganic_n) %>%
	pivot_longer(cols = inorganic_n, names_to = "flux")

n2o_n = fluxes %>%
	filter(n2o ==T) %>%
	select(name, cycle, n2o_n) %>%
	pivot_longer(cols = n2o_n, names_to = "flux")

methane_c = fluxes %>%
	filter(methane ==T) %>%
	select(name, cycle, methane_c) %>%
	pivot_longer(cols = methane_c, names_to = "flux")

all_fluxes <- do.call(rbind, list(nitrate_nitrite_n, gaseous_n, inorganic_n, n2o_n,methane_c)) %>%
	distinct()

all_fluxes$hour <- all_fluxes$cycle*timestep
all_fluxes$id <- all_fluxes$name
all_fluxes <- parseSimID(all_fluxes)

# Calculate flux rates at various time steps
df_calc <- all_fluxes
nitr_lm_out = list()
nmin_lm_out = list()
ammo_lm_out = list()

j =1
time_set = list(c(0, 50))
time_set = list(c(0, 25),
								c(0, 50),
								c(0, 75),
								c(0, 100),
								c(0, 200),
								c(25, 50),
								c(25, 75),
								c(25, 100),
								c(25, 200),
								c(25, 300),
								c(50, 75),
								c(50, 100),
								c(50, 125),
								c(50, 150),
								c(75, 100),
								c(100, 150),
								c(100, 200),
								c(100, 300),
								c(150, 200),
								c(150, 300),
								c(200, 250),
								c(200, 300),
								c(25, 300),
								c(0, 600),
								c(25, 600),
								c(200, 600),
								c(200, 400),
								c(50, 300),
								c(50, 275),
								c(10, 100),
								c(10, 200))

j = 4 # best nmin time
j = 17 # best nitr time


j = 23 # best overall time
#j = 24 # best overall time


i = c(25,300)

i = c(25,100)
i = c(0,25)

area = 400
eval_out = list()
for (j in 1:length(time_set)){
	i = time_set[[j]]
	duration = i[[2]] - i[[1]] / 24 # convert hourly to daily rate

	init_nitrate_col = paste0("nitrate_nitrite_n_",i[[1]])
	final_nitrate_col = paste0("nitrate_nitrite_n_",i[[2]])

	init_ammon_col = paste0("inorganic_n_",i[[1]])
	final_ammon_col = paste0("inorganic_n_",i[[2]])

	init_n2o_col = paste0("n2o_n_",i[[1]])
	final_n2o_col = paste0("n2o_n_",i[[2]])

	init_ch4_col = paste0("methane_c_",i[[1]])
	final_ch4_col = paste0("methane_c_",i[[2]])

	if (!i[[1]] %in% df_calc$hour | !i[[2]] %in% df_calc$hour) next()

	df_rates <- df_calc %>% filter(hour %in% i) %>% #group_by(id) %>%
		select(id, flux, hour, value) %>%
		pivot_wider(id_cols = c(id), names_from = c("flux", "hour"), values_from = c("value")) %>% as.data.frame()

	# df_rates$nmin = (df_rates[,..final_ammon_col] - df_rates[,..init_ammon_col]) / duration
	# df_rates$nitr = (df_rates[,..final_nitrate_col] - df_rates[,..init_nitrate_col]) / duration

	if (!init_ch4_col %in% colnames(df_rates)) df_rates[,init_ch4_col] <- 0
	if (!init_ammon_col %in% colnames(df_rates)) df_rates[,init_ammon_col] <- 0
	if (!init_nitrate_col %in% colnames(df_rates)) df_rates[,init_nitrate_col] <- 0

	# df_rates$net_inorganic = (df_rates[,final_ammon_col] + df_rates[,final_nitrate_col]) -
	# 	(df_rates[,init_ammon_col] + df_rates[,init_nitrate_col])

	df_rates$nmin_pred = (df_rates[,final_ammon_col] - df_rates[,init_ammon_col]) / duration / area
	df_rates$nitr_pred = (df_rates[,final_nitrate_col] - df_rates[,init_nitrate_col]) / duration / area
	df_rates$ammo_pred = df_rates$nmin_pred-df_rates$nitr_pred

	df_rates$n2o_pred = (df_rates[,final_n2o_col] - df_rates[,init_n2o_col]) / duration / area
	#df_rates$ch4_pred = (df_rates[,final_ch4_col] - df_rates[,init_ch4_col]) / duration / area

	eval <- merge(df_rates, eval_df, by="id")
	eval$ammo_obs = eval$nmin_obs - eval$nitr_obs
	nitr_lm_out[[j]] = summary(lm(nitr_obs ~ nitr_pred, eval))
	nmin_lm_out[[j]] = summary(lm(nmin_obs ~ nmin_pred, eval))
	ammo_lm_out[[j]] = summary(lm(ammo_obs ~ ammo_pred, eval))
	#	plot(nitr_obs ~ nitr_pred, eval)
#	plot(nmin_obs ~ nmin_pred, eval) #; abline(0,1)
eval_out[[j]] = eval %>% mutate(init_time = i[[1]], final_time = i[[2]])
}

#eval_out_df = do.call(rbind, eval_out)
eval_out_df = data.table::rbindlist(eval_out, fill = T)
times = eval_out_df[,c("init_time","final_time")] %>% unique()

# init_time = lapply(time_set, '[[', 1) %>% unlist()
# final_time = lapply(time_set, '[[', 2) %>% unlist()
nmin_lm_out <- nmin_lm_out %>% discard(is.null)
nitr_lm_out <- nitr_lm_out %>% discard(is.null)

nmin_lm_summary = lapply(nmin_lm_out, function(x) c(x$coefficients[2,1], x$r.squared)) %>%
	do.call(rbind, .) %>% cbind(., init_time = times[,1], final_time = times[,2])
colnames(nmin_lm_summary) = c("nmin_slope", "nmin_r.squared", "init_time", "final_time")
nitr_lm_summary = lapply(nitr_lm_out, function(x) c(x$coefficients[2,1], x$r.squared)) %>%
	do.call(rbind, .) %>% cbind(., init_time = times[,1], final_time = times[,2])
colnames(nitr_lm_summary) = c("nitr_slope", "nitr_r.squared", "init_time", "final_time")

lm_summary = cbind.data.frame(nmin_lm_summary, nitr_lm_summary)
lm_summary

for_eval = list(eval_out_df,lm_summary,met_by_timepoint,all_fluxes)

out_path = file.path("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary",
										 paste0(basename(sim_type), "_fluxes.rds"))
saveRDS(for_eval, out_path)
message("Saved to: ", out_path)
}


best_init = nitr_lm_summary[which.max(nitr_lm_summary[,"nitr_r.squared"]),]["init_time"]
best_final = nitr_lm_summary[which.max(nitr_lm_summary[,"nitr_r.squared"]),]["final_time"]

best_init = nmin_lm_summary[which.max(nmin_lm_summary[,"nmin_r.squared"]),]["init_time"]
best_final = nmin_lm_summary[which.max(nmin_lm_summary[,"nmin_r.squared"]),]["final_time"]
print(best_init)
print(best_final)
best_init = 25
best_final = 300

best_init = 25
best_final = 300
total_duration = best_final - best_init
init_nitrate_col = paste0("nitrateNitrite_",best_init)
final_nitrate_col = paste0("nitrateNitrite_",best_final)
init_ammon_col = paste0("inorganic_",best_init)
final_ammon_col = paste0("inorganic_",best_final)

summary(lm(nitr_obs ~ nitr_pred, eval))
plot(nitr_obs ~ nitr_pred, eval)
summary(lm(nmin_obs ~ nmin_pred, eval))
plot(nmin_obs ~ nmin_pred, eval); abline(0,1)

# We do predict concurrent methane and n2o production!
ggplot(eval) + geom_point(aes(x = n2o_pred, y = ch4_pred, color = moisture)) + theme_bw()

ggplot(eval) + geom_point(aes(x = moisture, y = ch4_pred))



for_eval = list(eval,lm_summary,met_by_timepoint,all_fluxes)

saveRDS(for_eval, "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/extramods_newdiffusion_fluxes_all.rds")


saveRDS(for_eval, "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/extramods_statico2_fluxes.rds")
saveRDS(for_eval, "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/extramods_newdiffusion_fluxes.rds")
saveRDS(for_eval, "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/statico2_fluxes.rds")
saveRDS(for_eval, "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/more_sugars_fluxes.rds")
saveRDS(for_eval, "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/no_static_fluxes.rds")

saveRDS(for_eval, "/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/n_only_fluxes.rds")


ggplot(df_calc %>% filter(!flux %in% c(#"o2_e",
	"h2o_e")), aes(y =value, x = cycle, color = diffusion, group=id)) +
	geom_line(size = 1, alpha=.3) +
	facet_grid(rows=vars(flux), scales="free") +
	xlim(c(0,125)) + theme_bw(base_size = 18)


ggplot(select_mets %>%
			 	filter(name %in% c("sim_diff_3.18_nh4_3.2_no3_0_mois_0.55_o2_0.013",
																				 #"sim_diff_3.18_nh4_0.42_no3_0_mois_0.23_o2_0.173",
																				 #"sim_diff_3.18_nh4_6.08_no3_0_mois_0.47_o2_0.01",
																				 "sim_diff_2.11_nh4_0.18_no3_0_mois_0.08_o2_0.01",
																				 #"sim_diff_2.11_nh4_2.67_no3_1.93_mois_0.45_o2_0.007",
																				 "sim_diff_2.11_nh4_2.65_no3_0.75_mois_0.2_o2_0.122") #& N_cycler
) %>%
	filter(!metabolite %in% c(#"o2_e",
	"h2o_e")), aes(y =total, x = cycle, color = metabolite, group=name)) +
	geom_path(size = 1, alpha=.3) +
	facet_grid(rows=vars(name), scales="free") +
	xlim(c(0,125)) +
	theme_bw(base_size = 18)





ggplot(select_mets %>%
			 	filter(name %in%  c("sim_diff_3.18_nh4_3.2_no3_0_mois_0.55_o2_0.013",
			 											#"sim_diff_3.18_nh4_0.42_no3_0_mois_0.23_o2_0.173",
			 											#"sim_diff_3.18_nh4_6.08_no3_0_mois_0.47_o2_0.01",
			 											"sim_diff_2.11_nh4_0.18_no3_0_mois_0.08_o2_0.01",
			 											#"sim_diff_2.11_nh4_2.67_no3_1.93_mois_0.45_o2_0.007",
			 											"sim_diff_2.11_nh4_2.65_no3_0.75_mois_0.2_o2_0.122") &
			 										 	!metabolite %in% c("co2_e", "h_e","biomass_e","o2_e","n2o_e",
			 		"h2o_e")),
			 		aes(y =total, x = cycle, color = metabolite, group=name)) +
	geom_path(size = 1, alpha=.8) +
	facet_grid(cols=vars(name), rows=vars(metabolite), scales="free_y") +
	xlim(c(0,125)) +
	theme_bw(base_size = 18)



eval_out_df$time_set = paste0(eval_out_df$init_time, "_", eval_out_df$final_time)
ggplot(eval_out_df,
			 aes(y =nmin_obs, x = nmin_pred, color = time_set)) +
	geom_point(size = 1, alpha=.8) + geom_smooth(method = "lm") + theme_bw()

ggplot(eval_out_df,
			 aes(y =nitr_obs, x = nitr_pred, color = time_set)) +
	geom_point(size = 1, alpha=.8) + geom_smooth(method = "lm") + theme_bw()




ggplot(select_mets %>% filter(metabolite %in% c("o2_e","nh4_e","no2_e","no3_e")),
			 aes(y =total, x = cycle, color = name, group=name)) +
	geom_path(size = 1, alpha=.8, show.legend = F) +
	facet_wrap(~metabolite) +
	theme_bw(base_size = 18)



ggplot(select_mets %>% filter(metabolite %in% c("no2_e","no3_e")),
			 aes(y =total, x = cycle, color = name, group=name)) +
	geom_path(size = 1, alpha=.8, show.legend = F) +
	facet_wrap(~metabolite) +
	theme_bw(base_size = 18)

ggplot(all_mets %>% filter(metabolite %in% c("glc__D_e","cit_e","cys__L_e")),
			 aes(y =total, x = cycle, color = name, group=name)) +
	geom_path(size = 1, alpha=.8, show.legend = F) +
	facet_wrap(~metabolite) +
	theme_bw(base_size = 18)


ggplot(all_mets %>% filter(metabolite %in% c("glc__D_e","cit_e","cys__L_e")),
			 aes(y =total, x = cycle, color = name, group=name)) +
	geom_path(size = 1, alpha=.8, show.legend = F) +
	facet_wrap(~metabolite) +
	theme_bw(base_size = 18)
