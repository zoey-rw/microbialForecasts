# Looking at coefficients of variation (volatility) for microbial time series

#install.packages("cvequality")
library(cvequality)

# Read in microbial abundances
d <- c(readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_16S_2021.rds"), 
			 readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_ITS_2021.rds"))



scores <- readRDS("./data/summary/CRPS_hindcasts.rds")

# Wide format taxonomic
scored_casts_tax <- scores$scored_hindcasts %>% filter(fcast_type=="Taxonomic" & model_name == "all_covariates" & taxon != "other")
scored_casts_tax_wide <- scored_casts_tax %>% 
	pivot_wider(id_cols = c(timepoint, plotID), names_from = taxon, values_from = truth) %>% 
	select(-c(timepoint, plotID))
tax_cv <- sapply(scored_casts_tax_wide, 
								 function(x) sd(x, na.rm = T) / mean(x, na.rm = T) * 100)
tax_cv_long <- tax_cv %>% stack()
colnames(tax_cv_long) <- c("CV", "taxon")


# Wide format functional group
scored_casts_fg <- scores$scored_hindcasts %>% 
	filter(fcast_type=="Functional group" & 
				 	model_name == "all_covariates" & taxon != "other") %>%  
	distinct()
scored_casts_fg_wide <- scored_casts_fg %>% 
	pivot_wider(id_cols = c(timepoint, plotID), names_from = species, values_from = truth) %>%
	select(-c(timepoint, plotID))
fg_cv <- sapply(scored_casts_fg_wide, 
								 function(x) sd(x, na.rm = T) / mean(x, na.rm = T) * 100)
fg_cv_long <- fg_cv %>% stack()
colnames(fg_cv_long) <- c("CV", "taxon")

cv <- rbind(tax_cv_long, fg_cv_long)

tax_scores <- scores$scored_hindcasts_taxon 
with_scores <- left_join(tax_scores, cv)


scores_in <- with_scores %>% filter(model_name == "all_covariates") %>%
	group_by(pretty_group, pretty_name) %>% 
	mutate(CV_scale = scale(CV)[,1],
				 crps_scale = scale(crps_mean)[,1])
#	mutate(outlier = ifelse(crps_mean == min(crps_mean), taxon, as.numeric(NA))) 
	
ggplot(scores_in %>% filter(newsite=="Observed site"), 
			 aes(x = CV_scale, y = crps_scale)) + 
	#coord_trans(y = "log10") +
	geom_point(aes(shape = pretty_group, color = pretty_name), #width=.2, #height = 0, 
							alpha = .4, size=4, 
							position=position_jitterdodge(dodge.width = 1)) + 
	#facet_wrap(~newsite) +  
	ylab("Mean continuous ranked probability score (CRPS)") + 
	xlab("Coefficient of Variation (scaled by rank)") + 
	ggtitle("Predictability ~ variability") +
	geom_smooth(method="lm") +
	theme_minimal(base_size=18) #+ 
#	geom_text(aes(label = outlier), na.rm = TRUE, hjust = -0.3) 


ggplot(scores_in, aes(x = CV, y = crps_mean)) + 
	coord_trans(y = "log10") +
	geom_point(aes(shape = pretty_group, color = pretty_group), #width=.2, #height = 0, 
						 alpha = .4, size=4, 
						 position=position_jitterdodge(dodge.width = 1)) + 
	ylab("Mean continuous ranked probability score (CRPS)") + xlab("Coefficient of Variation") + 
	ggtitle("Predictability ~ variability") +
	theme_minimal(base_size=18) + geom_text(aes(label = outlier), na.rm = TRUE, hjust = -0.3) + 
	facet_grid(pretty_name ~pretty_group, scales="free")+ geom_text(aes(label = outlier), na.rm = TRUE, hjust = -0.3) 


#### prev attempts #####


# why aren't these working??
# if you are reading this. yes you. tell me what I am doing wrong here?
asymptotic_test(seed = 1, cellulolytic, acidobacteriota)
mslr_test(nr = 10000, x = cellulolytic, y = acidobacteriota)

# tibble::glimpse(cellulolytic)
# num [1:2810] 0.0521 0.0587 0.0534 0.0392 0.0853 ... # all numeric, no NAs
# 
# tibble::glimpse(acidobacteriota)
# num [1:2810] 0.189 0.179 0.162 0.163 0.175 ... # all numeric, no NAs

# cvequality::mslr_test() function code

x <- cellulolytic
y <- acidobacteriota

if (!is.numeric(x) && !is.numeric(y) && !is.character(y)) {
	warning("x is not numeric or y is not numeric or character: returning NA")
	return(NA_real_)
}
if (anyNA(x)) {
	warning("x cannot contain any NA values: returning NA")
	return(NA_real_)
}
if (anyNA(y)) {
	warning("y cannot contain any NA values: returning NA")
	return(NA_real_)
}
n <- data.frame(table(y))$Freq
s <- aggregate(x, by = list(y), FUN = sd)$x
x <- aggregate(x, by = list(y), FUN = mean)$x
k <- length(x)
gv <- as.vector(nr)
df <- n - 1
xst0 <- LRT_STAT(n, x, s)
uh0 <- xst0[1:k]
tauh0 <- xst0[k + 1]
stat0 <- xst0[k + 2]
