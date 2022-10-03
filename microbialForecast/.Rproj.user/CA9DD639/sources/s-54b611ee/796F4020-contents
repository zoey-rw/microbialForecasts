#Harvard forest parameters from Abramoff et al. 2017

parameters <- read.csv("/projectnb/talbot-lab-data/zrwerbin/interactions/damm/DAMM-MCNiPv0/parameters.csv")
p = parameters$fitted2009 #select either default parameters or fitted to 2009 flux data
bd <- p[25]                           #bulk density
pd <- p[26]                           #particle density
porosity = 1 - bd/pd            #calculate porosity: .697

# Liquid diffusion: Equation S2.45 from Sihi et al., who cites Davidson et al., 2012
diffusion_liquid = 1/(porosity^3)

# soilMoisture: volumetric
# availability of aqueous substrates to enzymatic active sites, micromol/L
# Equation S2.47 from Sihi et al., used for C, NH4, NO3
substrate_avail = substrate * diffusion_liquid * (soilMoisture/100)^3


calc_substrate_avail = function(substrate, soilMoisture, porosity = .697) {
	# Liquid diffusion: Equation S2.45 from Sihi et al., who cites Davidson et al., 2012
	diffusion_liquid = 1/(porosity^3)
	substrate_avail = substrate * diffusion_liquid * (soilMoisture/100)^3
	return(substrate_avail)
}

# tortuosity of diffusion pathway for gases
# Equation 8 from Sihi et al 2021
alpha_4_3 = (porosity - (soilMoisture/100))^(4/3) * (soilTemperature + 273.15/293.15)^1.75

# Air diffusivity constant from Eq S2.42 from Sihi et al.
Dgas = .139 # cm-2 s-1
# Equation 9 from Sihi et al 2021
gas_diff = Dgas * atmospheric_concentration * alpha_4_3



soilData_harv_all <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/interactions/data/HARV_DAMM_all_data2.rds")  %>% filter(soilInorganicNugPerGram_initial < 25) %>% distinct(.keep_all = T)

soilData_harv = soilData_harv_all %>% select(sampleID, horizon, plotID, soilInorganicNugPerGram_initial, soilAmmoniumNugPerGram_initial, soilNitrateNitriteNugPerGram_initial,
																		nitrogenPercent, organicCPercent, CNratio,
																		litterDepth,standingWaterDepth, soilMoisture, soilInCaClpH, soilTemp, vol_SoilMoisture, sensor_moisture, bd, netNitugPerGramPerDay, netNminugPerGramPerDay) %>% distinct(.keep_all = T)

eval_df <- soilData_harv_all %>% mutate(porosity = 1 - soilData_harv_all$bd/pd,
																					diffusion_liquid = 1/(porosity^3),
																					avail_ammonium =
														 	calc_substrate_avail(soilAmmoniumNugPerGram_initial,
														 											 soilMoisture = vol_SoilMoisture*100,
														 											 porosity = porosity),
														 avail_nitrate_nitrite =
														 	calc_substrate_avail(soilNitrateNitriteNugPerGram_initial,
														 											 soilMoisture = vol_SoilMoisture*100,
														 											 porosity = porosity),
														 moisture = vol_SoilMoisture)



simple_df <- soilData_harv_all %>% mutate(porosity = 1 - soilData_harv_all$bd/pd,
																					diffusion_liquid = 1/(porosity^3),
																					avail_ammonium =
																						calc_substrate_avail(soilAmmoniumNugPerGram_initial,
																																 soilMoisture = vol_SoilMoisture*100,
																																 porosity = porosity),
																					avail_nitrate_nitrite =
																						calc_substrate_avail(soilNitrateNitriteNugPerGram_initial,
																																 soilMoisture = vol_SoilMoisture*100,
																																 porosity = porosity),
																					moisture = vol_SoilMoisture) %>%
	select(sampleID,porosity,soilNitrateNitriteNugPerGram_initial,avail_nitrate_nitrite,soilAmmoniumNugPerGram_initial,avail_ammonium,diffusion_liquid,moisture)




soilTemperature = 20
soilMoisture = .3

CO2airfrac = 385.5 * 1e-6 * 1.01e+5 # ppm co2 concentration in air converted to mol / m3 dry air
CH4airfrac = 1.8 * 1e-6 * 1.01e+5 # ppm ch4 concentration in air converted to mol / m3 dry air
O2airfrac = 0.209*10^(6) * 1e-6 * 1.01e+5  # ppm o2 concentration in air converted to mol / m3 dry air
H2airfrac = 0.0000005*10^(6) * 1e-6 * 1.01e+5  # ppm h2 concentration in air converted to mol / m3 dry air

out_df_list = list()
another_list = list()

# Use Katie Atherton's BD values from Harvard Forest instead
soilData_harv$bd <- ifelse(soilData_harv$horizon=="M", .8, .55)

for (i in 1:nrow(soilData_harv)) {

	print(i)

	bd = soilData_harv$bd[[i]]
	porosity = 1 - (bd/pd)
	#soilMoist = soilData_harv$sensor_moisture[[i]]
	soilMoist = soilData_harv$vol_SoilMoisture[[i]]
	soilTemp = soilData_harv$soilTemp[[i]]
	init_ammonium = soilData_harv$soilAmmoniumNugPerGram_initial[[i]]
	init_nitrate_nitrite = soilData_harv$soilNitrateNitriteNugPerGram_initial[[i]]
	Dliq = 1/(porosity^3)

	nmin_obs = soilData_harv$netNminugPerGramPerDay[[i]]
	nitr_obs = soilData_harv$netNitugPerGramPerDay[[i]]

	a_4_3 <- ((porosity-soilMoist/100)^(4/3))*((soilTemp/273.15)^1.75)
a_4_3_max <- ((porosity-0/100)^(4/3))*(((25+273.15)/273.15)^1.75)
Dgas <- 1/a_4_3_max
CO2Dif <- Dgas * CO2airfrac * a_4_3
CH4Dif <- Dgas * CH4airfrac * a_4_3
O2Dif <- Dgas * O2airfrac * a_4_3
H2Dif <- Dgas * H2airfrac * a_4_3

# Convert back to mmol/cm3 for COMETS inputs
CO2Dif = CO2Dif/1000
CH4Dif = CH4Dif/1000
O2Dif = O2Dif/1000
H2Dif = H2Dif/1000

avail_ammonium = init_ammonium * Dliq * (soilMoist)
avail_nitrate_nitrite = init_nitrate_nitrite * Dliq * (soilMoist)


out_df = data.frame(sampleID = soilData_harv$sampleID[[i]],
										horizon = soilData_harv$horizon[[i]],
										mois = soilMoist,
										temp = soilTemp,
										bd = bd,
										Dgas = Dgas,
										Dliq = Dliq,
										CO2 = CO2Dif,
										CH4 = CH4Dif,
										O2 = O2Dif,
										H2 = H2Dif,
										avail_ammonium = avail_ammonium,
										avail_nitrate_nitrite = avail_nitrate_nitrite,
										nmin_obs = nmin_obs,
											nitr_obs = nitr_obs)
out_df_list[[i]] <- out_df

}


sim_gas = data.table::rbindlist(out_df_list)
plot(sim_gas$O2 ~ sim_gas$mois)
plot(sim_gas$O2 ~ sim_gas$CO2)

ghg_inputs <- read.csv("/projectnb/talbot-lab-data/zrwerbin/interactions/PR-model/flux.csv")
plot(ghg_inputs$O2 ~ ghg_inputs$SoilM)
plot(ghg_inputs$SoilT ~ ghg_inputs$SoilM)



eval_df <- sim_gas %>% mutate(moisture = round(mois, 2),
															liq_diffusion = round(Dliq, 2),
															gas_diffusion = round(Dgas, 2),
															ammonium = round(avail_ammonium, 2),
															nitrate = round(avail_nitrate_nitrite, 2),
															oxygen = round(O2, 3),
															co2 = round(CO2, 2),
															#diffusion = ifelse(horizon=="O",1e-6, 10e-6),
															id = paste0('sim_diff_', liq_diffusion, '_nh4_',
																					ammonium,'_no3_', nitrate,'_mois_',moisture,"_o2_",oxygen))

write.csv(eval_df, "/projectnb/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/true_rates_df.csv")



sim_categories <- eval_df %>% select(moisture, moisture,
																		 diffusion = liq_diffusion,
																		 ammonium, nitrate,oxygen,co2, id) %>%
	distinct(id, .keep_all=T)
rownames(sim_categories) <- NULL
sim_categories$rownames <- 1:nrow(sim_categories)

write.csv(sim_categories, "/projectnb/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_categories.csv")

