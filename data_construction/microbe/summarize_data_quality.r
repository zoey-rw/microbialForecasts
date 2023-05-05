# Summarize read depth
library(knitr)
library(kableExtra)
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
source("/projectnb/dietzelab/zrwerbin/microbialForecasts/source.R")

tax_16S = readRDS("/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/data/NEON_16S_tax.rds")

# Percent classified to genus
colSums(is.na(tax_16S))/nrow(tax_16S)

tax_ITS = readRDS("/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/data/NEON_ITS_tax.rds")
colSums(is.na(tax_ITS))/nrow(tax_ITS)

# Bacterial final functional groups
ps <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/phyloseq_16S.rds")

# Number of ESVs
dim(ps@otu_table); #58208

# Get average abundance per group
bacteria <- readRDS(here("data", "clean/groupAbundances_16S_2023.rds"))
fungi <- readRDS(here("data", "clean/groupAbundances_ITS_2023.rds"))
all_ranks = c(bacteria, fungi) 

mean_list=list()
sd_list=list()
for (k in 1:length(all_ranks)) {
	mean_list[[k]] = all_ranks[[k]] %>%
		summarise_if(is.numeric, mean, na.rm=T) %>% 
		stack  %>% rename("mean" = "values" , "fg" = "ind" )
	sd_list[[k]] = all_ranks[[k]] %>%
		summarise_if(is.numeric, sd, na.rm=T) %>% 
		stack %>% rename("sd" = "values" , "fg" = "ind" )
}
fg_mean = data.frame(do.call(rbind, mean_list)) %>% filter(fg != "other")
fg_sd = data.frame(do.call(rbind, sd_list)) %>% filter(fg != "other")
avg_abundance = merge(fg_mean, fg_sd)


fg_list = c("cellulolytic",
"assim_nitrite_reduction", "dissim_nitrite_reduction", "assim_nitrate_reduction",
"n_fixation", "dissim_nitrate_reduction", "nitrification", "denitrification",
"chitinolytic", "lignolytic", "methanotroph", "copiotroph", "oligotroph",
"benomyl_antibiotic", "citrate_simple", "glucose_simple", "glycine_simple",
"pyruvate_simple", "streptomycin_antibiotic", "sucrose_complex",
"acetogen_anaerobic", "fsfeso4_anaerobic", "ironcitrate_anaerobic",
"nonfsfeso4_anaerobic", "potassiumnitrate_anaerobic", "chloramphenicol_antibiotic",
"erythromycin_antibiotic", "gentamycin_antibiotic", "nystatin_antibiotic",
"glycerol_simple", "acetate_simple", "glutamate_simple", "late_stress",
"propionate_simple", "acidic_stress", "alkaline_stress", "d_galacturonicacid_simple",
"d_glucuronicacid_simple", "arabinose_simple", "cellobiose_complex",
"cellulose_complex", "chitin_complex", "galactose_simple", "glucosamine_simple",
"mannose_simple", "n_acetylglucosamine_simple", "pectin_complex",
"rhamnose_simple", "trehalose_complex", "xylan_complex", "xylose_simple",
"salt_stress", "herbicide_stress", "lowcarbon_stress", "osmotic_stress",
"heat_stress", "light_stress", "endophyte")

tax_df = as.data.frame(as.matrix(ps@tax_table))

out <- list()
fg = fg_list[[1]]
for (fg in fg_list){

out_list <- list()
print(fg)
fg_tax_df = tax_df %>% filter(get({{fg}}) != "other")
n_esv = nrow(fg_tax_df)
n_genus = length(unique(fg_tax_df$genus))

clean_fg_names = recode(fg, !!!microbialForecast:::pretty_names)

fg_kingdom = recode(fg, !!!microbialForecast:::FG_kingdoms)

out[[fg]] = cbind(fg, n_esv, n_genus, #mean_abun, sd_abun, 
									clean_fg_names, fg_kingdom)
}

fg_summary_table = data.frame(do.call(rbind, out))




fg_sources = c("Naylor et al. 2010",
							 "Albright et al. 2018",
							 "Berlemont & Martiny 2013",
							 "Ho et al. 2017")
#fg_source = recode(fg, !!!fg_sources)


# What percent were not included in any group?

fg_df <- tax_df[,c(7:ncol(tax_df))]
unclassified = fg_df[rowSums(fg_df == "other") == ncol(fg_df),]
dim(unclassified)
nrow(unclassified) / nrow(tax_df)




# Bacterial final functional groups
ps_ITS <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/phyloseq_ITS.rds")

# Number of ESVs
dim(ps_ITS@otu_table); #58208

tax_df_its = as.data.frame(as.matrix(ps_ITS@tax_table))

out <- list()
fg_list_its = c("endophyte", 
"plant_pathogen", "animal_pathogen", "ectomycorrhizal", "lichenized", 
"wood_saprotroph", "soil_saprotroph", "litter_saprotroph", "saprotroph"
)
for (fg in fg_list_its){
	
	out_list <- list()
	print(fg)
	fg_tax_df_ITS = tax_df_its %>% filter(get({{fg}}) != "other")
	n_esv = nrow(fg_tax_df)
	n_genus = length(unique(fg_tax_df_ITS$genus))
	
	clean_fg_names = recode(fg, !!!microbialForecast:::pretty_names)
	
	fg_kingdom = recode(fg, !!!microbialForecast:::FG_kingdoms)
	
	out[[fg]] = cbind(fg, n_esv, n_genus, #mean_abun, sd_abun, 
										clean_fg_names, fg_kingdom)
}
fg_summary_table_its = data.frame(do.call(rbind, out))


# What percent were not included in any group?

fg_df_its <- tax_df_its[,c(6:ncol(tax_df_its))]
unclassified_its = fg_df_its[rowSums(fg_df_its == "other") == ncol(fg_df_its),]
unclassified_its_pct = nrow(unclassified_its) / nrow(tax_df_its)
1 - unclassified_its_pct # 0.2326402



fg_summary_out <- rbind(fg_summary_table, fg_summary_table_its)

fg_summary_out <- merge(avg_abundance, fg_summary_out, all.x=F, by="fg") 
colnames(fg_summary_out) <- c("Group ID","Mean abundance", "SD abundance", "Unique ESVs","Unique genera","Functional group name", "Functional group kingdom")

fg_summary_out <- fg_summary_out %>% filter(`Group ID` %in% microbialForecast:::keep_fg_names)
saveRDS(fg_summary_out, here("figures", "fg_summary_table.rds"))




kable(fg_summary_out, "html") %>%
	kable_styling(bootstrap_options = c("striped", "hover")) %>%
	cat(., file = here("figures", "fg_summary_table.html"))



harv_ps <- ps_ITS %>%
	speedyseq::filter_sample_data(siteID %in% c("HARV"))
