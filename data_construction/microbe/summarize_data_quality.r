# Summarize read depth
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

tax_16S = readRDS("/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/data/NEON_16S_tax.rds")

# Percent classified to genus
colSums(is.na(tax_16S))/nrow(tax_16S)


tax_ITS = readRDS("/projectnb/dietzelab/zrwerbin/NEON_soil_microbe_processing/data/NEON_ITS_tax.rds")
colSums(is.na(tax_ITS))/nrow(tax_ITS)


# Bacterial final functional groups
ps <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/phyloseq_16S.rds")

# Number of ESVs
dim(ps@otu_table); #58208


dim(ps@tax_table); #58208

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
"heat_stress", "light_stress")

tax_df = as.data.frame(as.matrix(ps@tax_table))

out <- list()
fg = fg_list[[1]]
for (fg in fg_list[1:5]){

out_list <- list()
print(fg)
fg_tax_df = tax_df %>% filter(get({{fg}}) != "other")
n_esv = nrow(fg_tax_df)
n_genus = length(unique(fg_tax_df$genus))
glom <- speedyseq::tax_glom(ps, taxrank=({{fg}}))
glom_melt <- speedyseq::psmelt(glom)
form <- as.formula(paste0("sampleID ~ ", fg))
glom_wide <- reshape2::dcast(glom_melt, form, value.var = "Abundance", fun.aggregate = sum)
out_abun <- transform(glom_wide, row.names=sampleID, sampleID=NULL)
mean_abun = mean(out_abun[,1], na.rm = T)
sd_abun = sd(out_abun[,1], na.rm = T)

clean_fg_names = recode(fg, !!!microbialForecast:::pretty_names)

fg_kingdom = recode(fg, !!!microbialForecast:::FG_kingdoms)
fg_kingdom = recode(fg, !!!fg_sources)

out[[fg]] = cbind(n_esv, n_genus, mean_abun, sd_abun, clean_fg_names, fg_kingdom)
}

fg_summary_table = data.frame(do.call(rbind, out))
fg_summary_table

fg_sources = c("Naylor et al. 2010",
							 "Albright et al. 2018",
							 "Berlemont & Martiny 2013",
							 "Ho et al. 2017")

unclassified = tax_df[all_equal(tax_df[7:ncol(tax_df)])]


# What percent fell into any group at all? What percent were not included in any group?
