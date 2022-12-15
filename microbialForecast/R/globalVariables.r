# Global variables to make plotting easier
message("Loading global variables")
#' @importFrom nimble nimbleCode
#' @importFrom grDevices axisTicks
#' @importFrom graphics par
#' @importFrom utils head tail
#' @import coda
#' @import dplyr
#' @import tidyverse

utils::globalVariables(c("keep_fg_names", "tax_names", "all_covariates_key","rank_spec_names","date_recode"))


#### global variables ####
keep_fg_names <- c("cellulolytic", "assim_nitrite_reduction", "dissim_nitrite_reduction",
									 "assim_nitrate_reduction", "n_fixation", "dissim_nitrate_reduction",
									 "nitrification", "denitrification", "chitinolytic", "lignolytic",
									 "copiotroph", "oligotroph", "benomyl_antibiotic", "glucose_simple",  "pyruvate_simple",
									 "streptomycin_antibiotic", "sucrose_complex", "acetogen_anaerobic",
									 "chloramphenicol_antibiotic", "erythromycin_antibiotic",
									 "gentamycin_antibiotic", "glycerol_simple",
									 "acetate_simple",
									 "acidic_stress", "cellobiose_complex",
									 "cellulose_complex", "chitin_complex", "galactose_simple",
									 "xylose_simple", "salt_stress", "herbicide_stress", "osmotic_stress",
									 "heat_stress", "light_stress", "endophyte", "plant_pathogen",
									 "animal_pathogen", "ectomycorrhizal", "lichenized", "saprotroph")

fg_names <- c("cellulolytic", "assim_nitrite_reduction", "dissim_nitrite_reduction",
							"assim_nitrate_reduction", "n_fixation", "dissim_nitrate_reduction",
							"nitrification", "denitrification", "chitinolytic", "lignolytic",
							"methanotroph", "copiotroph", "oligotroph", "benomyl_antibiotic",
							"citrate_simple", "glucose_simple", "glycine_simple", "pyruvate_simple",
							"streptomycin_antibiotic", "sucrose_complex", "acetogen_anaerobic",
							"fsfeso4_anaerobic", "ironcitrate_anaerobic", "nonfsfeso4_anaerobic",
							"potassiumnitrate_anaerobic", "chloramphenicol_antibiotic", "erythromycin_antibiotic",
							"gentamycin_antibiotic", "nystatin_antibiotic", "glycerol_simple",
							"acetate_simple", "glutamate_simple", "late_stress", "propionate_simple",
							"acidic_stress", "alkaline_stress", "d_galacturonicacid_simple",
							"d_glucuronicacid_simple", "arabinose_simple", "cellobiose_complex",
							"cellulose_complex", "chitin_complex", "galactose_simple", "glucosamine_simple",
							"mannose_simple", "n_acetylglucosamine_simple", "pectin_complex",
							"rhamnose_simple", "trehalose_complex", "xylan_complex", "xylose_simple",
							"salt_stress", "herbicide_stress", "lowcarbon_stress", "osmotic_stress",
							"heat_stress", "light_stress", "endophyte", "plant_pathogen",
							"animal_pathogen", "ectomycorrhizal", "lichenized", "wood_saprotroph",
							"soil_saprotroph", "litter_saprotroph", "saprotroph")

tax_names <- c("phylum_bac", "class_bac", "order_bac", "family_bac", "genus_bac",
							 "phylum_fun", "class_fun", "order_fun", "family_fun", "genus_fun")

div_scenarios <- c("no_uncertainty_ITS", "spatial_uncertainty_ITS", "temporal_uncertainty_ITS",
									 "full_uncertainty_ITS", "no_uncertainty_16S", "spatial_uncertainty_16S",
									 "temporal_uncertainty_16S", "full_uncertainty_16S")


all_covariates_key <- c("1" = "Temperature",
												"2" = "Moisture",
												"3" = "pH",
												"4" = "pC",
												"5" = "Ectomycorrhizal trees",
												"6" = "LAI",
												"7" = "sin",
												"8" = "cos",
												"NA" = "NA")

cycl_only_key <- list("1" = "sin",
											"2" = "cos")

var_list <- c("tau_obs","tau_proc","plot_var","site_var", "plot_effect", "time_var", "beta", "site_effect","plot_mean_hat","alpha")
var_listSimple <- c("tau_obs","tau_proc","plot_var", "plot_effect",#"glob_mean",
										#"time_var",
										"plot_mean",
										"beta","plot_mean_hat","alpha")



pretty_names_old <- list("cellulolytic" = "Cellulolytic bacteria",
												 "assim_nitrite_reduction" = "Assimilatory nitrite reducing bacteria",
												 "dissim_nitrite_reduction" = "Dissimilatory nitrite reducing bacteria",
												 "assim_nitrate_reduction" = "Assimilatory nitrate reducing bacteria",
												 "n_fixation" = "Nitrogen fixing bacteria",
												 "dissim_nitrate_reduction" = "Dissimilatory nitrate reducing bacteria",
												 "nitrification" = "Nitrifying bacteria",
												 "denitrification" = "Denitrifying bacteria",
												 "chitinolytic" = "Chitinolytic bacteria",
												 "lignolytic" = "Ligninolytic bacteria",
												 "methanotroph" = "Methanotrophic bacteria",
												 "copiotroph" = "Copiotrophic bacteria",
												 "oligotroph" = "Oligotrophic bacteria",
												 "Arbuscular" = "Arbuscular mycorrhizal fungi",
												 "Animal_Pathogen" = "Animal pathogen fungi",
												 "Plant_Pathogen" = "Plant pathogen fungi",
												 "Saprotroph" = "Saprotrophs fungi",
												 "Wood_Saprotroph" = "Wood saprotrophs (fungi)",
												 "Ectomycorrhizal" = "Ectomycorrhizal fungi",
												 "Bdellovibrio" = "Genus Bdellovibrio",
												 "Pseudogymnoascus" = "Genus Pseudogymnoascus",
												 "Cyanobacteria" = "Phylum Cyanobacteria"
)

pretty_rank_names <- list("genus_bac" = "Genus",
													"family_bac" = "Family",
													"order_bac" = "Order",
													"class_bac" = "Class",
													"phylum_bac" = "Phylum",
													"genus_fun" = "Genus",
													"family_fun" = "Family",
													"order_fun" = "Order",
													"class_fun" = "Class",
													"phylum_fun" = "Phylum",
													"functional_group" = "Functional group",
													"diversity_16S" = "Diversity",
													"diversity_ITS" = "Diversity",
													"diversity" = "Diversity"
)

pretty_names <- list("cellulolytic" = "Cellulose degraders",
										 "assim_nitrite_reduction" = "Assimilatory nitrite reducers",
										 "dissim_nitrite_reduction" = "Dissimilatory nitrite reducers",
										 "assim_nitrate_reduction" = "Assimilatory nitrate reducers",
										 "n_fixation" = "Nitrogen fixers",
										 "dissim_nitrate_reduction" = "Dissimilatory nitrate reducers",
										 "nitrification" = "Nitrifiers",
										 "denitrification" = "Denitrifiers",
										 "chitinolytic" = "Chitin degraders",
										 "lignolytic" = "Lignin degraders",
										 "methanotroph" = "Methanotrophs",
										 "copiotroph" = "Copiotrophs",
										 "oligotroph" = "Oligotrophs",
										 "Arbuscular" = "Arbuscular mycorrhizae",
										 "Animal_Pathogen" = "Animal pathogens",
										 "Plant_Pathogen" = "Plant pathogens",
										 "Saprotroph" = "Saprotrophs",
										 "Wood_Saprotroph" = "Wood saprotrophs",
										 "Ectomycorrhizal" = "Ectomycorrhizae",
										 "Bdellovibrio" = "Genus Bdellovibrio",
										 "Pseudogymnoascus" = "Genus Pseudogymnoascus",
										 "Cyanobacteria" = "Phylum Cyanobacteria"
)


FG_kingdoms <- list("cellulolytic" = "Bacterial_functional_group",
										"assim_nitrite_reduction" = "Bacterial_functional_group",
										"dissim_nitrite_reduction" = "Bacterial_functional_group",
										"assim_nitrate_reduction" = "Bacterial_functional_group",
										"n_fixation" = "Bacterial_functional_group",
										"dissim_nitrate_reduction" = "Bacterial_functional_group",
										"nitrification" = "Bacterial_functional_group",
										"denitrification" = "Bacterial_functional_group",
										"chitinolytic" = "Bacterial_functional_group",
										"lignolytic" = "Bacterial_functional_group",
										"methanotroph" = "Bacterial_functional_group",
										"copiotroph" = "Bacterial_functional_group",
										"oligotroph" = "Bacterial_functional_group",
										"Arbuscular" = "Fungal_functional_group",
										"Animal_Pathogen" = "Fungal_functional_group",
										"Plant_Pathogen" = "Fungal_functional_group",
										"Saprotroph" = "Fungal_functional_group",
										"Wood_Saprotroph" = "Fungal_functional_group",
										"Ectomycorrhizal" = "Fungal_functional_group")


N_cyclers <- c("assim_nitrite_reduction", "dissim_nitrite_reduction","assim_nitrate_reduction","n_fixation","dissim_nitrate_reduction","nitrification","denitrification")

date_recode <- c("20151101_20180101" = "2015-11_2018-01",
								 "20160101_20200101" = "2016-01_2020-01",
								 "20130601_20170101" = "2013-06_2017-01",
								 "20130601_20150101" = "2013-06_2015-01",
								 "20130601_20200101" = "2013-06_2020-01",
								 "20151101_20200101" = "2015-11_2020-01")
calibration_label <- c("2015-11_2018-01"= "Excluding legacy data",
											 "2016-01_2020-01"= "Mistake",
											 "2013-06_2017-01"= "Including legacy data",
											 "2013-06_2015-01"= "Legacy only",
											 "2013-06_2020-01"= "Legacy + current (full dataset)")

rank_spec_names = list(phylum_bac = c("acidobacteriota", "actinobacteriota", "armatimonadota",
										"bacteroidota", "bdellovibrionota", "chloroflexi", "chloroflexota",
										"crenarchaeota", "cyanobacteria", "dependentiae", "desulfobacterota",
										"elusimicrobiota", "fibrobacterota", "firmicutes", "firmicutes_a",
										"firmicutes_f", "gal15", "gemmatimonadota", "mbnt15", "methylomirabilota",
										"myxococcota", "nb1.j", "nitrospirota", "patescibacteria", "planctomycetota",
										"proteobacteria", "rcp2.54", "spirochaetota", "verrucomicrobiota",
										"wps.2"), class_bac = c("acidimicrobiia", "acidobacteriae", "acidobacteriota_class",
																						"actinobacteria", "actinobacteriota_class", "ad3", "alphaproteobacteria",
																						"bacilli", "bacteroidia", "blastocatellia", "clostridia", "cyanobacteriia",
																						"gammaproteobacteria", "gemmatimonadetes", "kd4.96", "ktedonobacteria",
																						"mb.a2.108", "microgenomatia", "myxococcia", "nitrospiria", "phycisphaerae",
																						"planctomycetes", "polyangia", "rcp2.54_class", "rubrobacteria",
																						"subgroup.5", "thermoleophilia", "tk10", "verrucomicrobiae",
																						"vicinamibacteria"), order_bac = c("acetobacterales", "acidobacteriales",
																																							 "actinomycetales", "bacillales", "bryobacterales", "burkholderiales",
																																							 "caulobacterales", "chitinophagales", "chthoniobacterales", "corynebacteriales",
																																							 "cytophagales", "enterobacterales", "frankiales", "gaiellales",
																																							 "gammaproteobacteria.incertae.sedis", "gemmatimonadales", "imcc26256",
																																							 "lachnospirales", "lactobacillales", "pedosphaerales", "polyangiales",
																																							 "pyrinomonadales", "reyranellales", "rhizobiales", "solibacterales",
																																							 "solirubrobacterales", "sphingobacteriales", "sphingomonadales",
																																							 "tepidisphaerales", "xanthomonadales"), family_bac = c("acetobacteraceae",
																																							 																											 "bacillaceae", "beijerinckiaceae", "bryobacteraceae", "burkholderiaceae",
																																							 																											 "caulobacteraceae", "chitinophagaceae", "chthoniobacteraceae",
																																							 																											 "gaiellaceae", "gemmatimonadaceae", "hyphomicrobiaceae", "lachnospiraceae",
																																							 																											 "mycobacteriaceae", "nitrosomonadaceae", "nitrospiraceae", "nocardioidaceae",
																																							 																											 "pedosphaeraceae", "pirellulaceae", "pyrinomonadaceae", "reyranellaceae",
																																							 																											 "rhizobiaceae", "solibacteraceae", "solirubrobacteraceae", "sphingomonadaceae",
																																							 																											 "streptomycetaceae", "vibrionaceae", "wd2101.soil.group", "xanthobacteraceae",
																																							 																											 "xiphinematobacteraceae"), genus_bac = c("acidibacter", "adurb.bin063.1",
																																							 																											 																				 "bacillus", "bryobacter", "candidatus.koribacter", "candidatus.solibacter",
																																							 																											 																				 "candidatus.udaeobacter", "candidatus.xiphinematobacter", "caulobacter",
																																							 																											 																				 "chthoniobacter", "conexibacter", "dongia", "dorea", "ellin6067",
																																							 																											 																				 "gaiella", "haliangium", "mnd1", "mycobacterium", "nitrospira",
																																							 																											 																				 "nocardioides", "pedomicrobium", "puia", "rb41", "reyranella",
																																							 																											 																				 "rhodoplanes", "rickettsia", "sphingomonas", "steroidobacter",
																																							 																											 																				 "streptococcus", "streptomyces"), phylum_fun = c("ascomycota",
																																							 																											 																				 																								 "basidiomycota", "chytridiomycota", "glomeromycota", "mortierellomycota",
																																							 																											 																				 																								 "mucoromycota", "rozellomycota"), class_fun = c("agaricomycetes",
																																							 																											 																				 																								 																								"archaeorhizomycetes", "chytridiomycetes", "dothideomycetes",
																																							 																											 																				 																								 																								"eurotiomycetes", "geminibasidiomycetes", "glomeromycetes", "leotiomycetes",
																																							 																											 																				 																								 																								"microbotryomycetes", "mortierellomycetes", "mucoromycetes",
																																							 																											 																				 																								 																								"mucoromycotina_cls_incertae_sedis", "orbiliomycetes", "pezizomycetes",
																																							 																											 																				 																								 																								"rozellomycotina_cls_incertae_sedis", "saccharomycetes", "sordariomycetes",
																																							 																											 																				 																								 																								"spizellomycetes", "tremellomycetes", "umbelopsidomycetes"),
		 order_fun = c("agaricales", "archaeorhizomycetales", "cantharellales",
		 							"capnodiales", "chaetosphaeriales", "chaetothyriales", "eurotiales",
		 							"glomerales", "helotiales", "hypocreales", "mortierellales",
		 							"pezizales", "pleosporales", "russulales", "sordariales",
		 							"thelebolales", "thelephorales", "tremellales", "umbelopsidales",
		 							"xylariales"), family_fun = c("archaeorhizomycetaceae", "aspergillaceae",
		 																						"chaetomiaceae", "chaetosphaeriaceae", "cladosporiaceae",
		 																						"clavariaceae", "clavicipitaceae", "entolomataceae", "glomeraceae",
		 																						"helotiales_fam_incertae_sedis", "herpotrichiellaceae", "hypocreaceae",
		 																						"mortierellaceae", "myxotrichaceae", "nectriaceae", "pleosporaceae",
		 																						"russulaceae", "thelephoraceae", "trichocomaceae", "umbelopsidaceae"
		 							), genus_fun = c("archaeorhizomyces", "aspergillus", "cenococcum",
		 															 "chaetomium", "cladophialophora", "cladosporium", "coniochaeta",
		 															 "cortinarius", "entoloma", "exophiala", "fusarium", "glomus",
		 															 "metarhizium", "mortierella", "oidiodendron", "penicillium",
		 															 "russula", "talaromyces", "trichoderma", "umbelopsis"))

# usethis::use_data(cycl_only_key, all_covariates_key,
# 									keep_fg_names, fg_names, tax_names, pretty_rank_names, pretty_names,
# 									calibration_label,date_recode, N_cyclers, FG_kingdoms, rank_spec_names, internal = TRUE, overwrite=T)
