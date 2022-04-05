library(runjags)
library(coda)
library(forestplot)
library(ggplot2)
library(gridExtra)
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/helperFunctions.r")


sum.all.tax <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/allrank_summaries.rds")
sum.all.fg <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/fg_summaries.rds")
sum.all.tax$beta_eff$fg_cat <- NA
sum.all.tax$beta_eff$fg_cat <- ifelse(grepl("fun$", sum.all.tax$beta_eff$rank), "Fungal taxa", "Bacterial taxa")

df <- plyr::rbind.fill(sum.all.tax$beta_eff, sum.all.fg$full_uncertainty$beta_eff)

df <- df[which(!df$taxon_name %in% c("other")),]



df$pretty_name <- recode(df$rank,
															"genus_bac" = "Genus",
															"family_bac" = "Family",
															"order_bac" = "Order", 
															"class_bac" = "Class", 
															"phylum_bac" = "Phylum",
												 "genus_fun" = "Genus",
												 "family_fun" = "Family",
												 "order_fun" = "Order", 
												 "class_fun" = "Class", 
												 "phylum_fun" = "Phylum",
															"functional_group" = "Functional group")

df$only_rank <- sapply(str_split(df$rank, "_",  n = 2), `[`, 1)
df$only_rank <- ordered(df$only_rank, levels = c("genus",
																			 "family",
																			 "order", 
																			 "class",
																			 "phylum", "functional"))
df$rank <- ordered(df$rank, levels = c("genus_bac","genus_fun",
																			 "family_bac","family_fun",
																			 "order_bac", "order_fun",
																			 "class_bac", "class_fun",
																			 "phylum_bac","phylum_fun",
																			 "functional_group"))

#sum.all <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/effectSizeSummary_bac_fg.rds")


df$sig_strict <- ifelse(df$`2.5%` < 0 & df$`97.5%` < 0 |
                   df$`2.5%` > -0 & df$`97.5%` > -0,
                      1, 0)
# df$Lower95 <- df$Mean - 1.96*df$SD
# df$Upper95 <- df$Mean + 1.96*df$SD
df$effSize <- abs(df$`50%`)
# ## GGPLOT
sum.all <- df

fg_summary <- sum.all[sum.all$rank=="functional_group",]
fg_summary[grep("lytic", fg_summary$taxon),]$fg_cat <- "Decomposers"
fg_summary[grep("methan", fg_summary$taxon),]$fg_cat <- "Methanotroph"

saveRDS(sum.all, "/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/for_plotting_ranks.rds")

## NO FUNCTIONAL GROUPS


### DOT PLOTS ###
# By rank with every taxon - cluttered
ggplot(data=sum.all[c(sum.all$fg_cat %in% c("Bacterial taxa","Fungal taxa")),],
       aes(x = taxon_name,y = effSize)) +
  geom_point(aes(shape = as.factor(sig_strict), color = beta), size = 3) +
  labs(col = "Parameter", title = "Absolute effect size") + 
  xlab("Taxon")+ 
  ylab(NULL) +
  facet_grid(rows = vars(fg_cat), cols = vars(only_rank), scales = "free", space = "free_x") + #,strip.position="bottom",nrow=2) +
 theme_bw() + theme(#axis.ticks.x=element_blank(),
 	text = element_text(size = 16),
 	 axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
 	 	angle = 320, vjust=1, hjust = -0.05),
 	axis.title=element_text(size=22,face="bold")
 	#strip.text.y = element_text(size=24,hjust=0,vjust = 1,angle=180,face="bold")
 ) + scale_shape_manual(values = c(21, 16), name = NULL, 
 											 labels = c("Not significant","Significant")) 


# By rank with every taxon - cluttered - only significant effects
ggplot(data=sum.all[which(sum.all$fg_cat %in% c("Bacterial taxa","Fungal taxa") & sum.all$sig_strict == 1),],
			 aes(x = taxon_name,y = effSize)) +
	geom_point(aes(color = beta), size = 3) +
	labs(col = "Parameter", title = "Absolute effect size (significant effects only)") + 
	xlab("Taxon")+ 
	ylab(NULL) +
	facet_grid(rows = vars(fg_cat), cols = vars(only_rank), scales = "free", space = "free_x") + #,strip.position="bottom",nrow=2) +
	theme_bw() + theme(#axis.ticks.x=element_blank(),
		text = element_text(size = 16),
		axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
			angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=22,face="bold")
		#strip.text.y = element_text(size=24,hjust=0,vjust = 1,angle=180,face="bold")
	) 

# Fungal taxa - all taxa - not absolute effects
ggplot(data=sum.all[which(sum.all$fg_cat %in% c("Fungal taxa")),],
			 aes(x = taxon_name,y = `50%`)) +
	geom_point(aes(color = beta), size = 3) +
	labs(col = "Parameter", title = "Effect size") + 
	xlab("Taxon")+ 
	ylab(NULL) +
	facet_grid(rows = vars(fg_cat), cols = vars(only_rank), scales = "free", space = "free_x") + #,strip.position="bottom",nrow=2) +
	theme_bw() + theme(#axis.ticks.x=element_blank(),
		text = element_text(size = 16),
		axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
			angle = 320, vjust=1, hjust = -0.05),
		axis.title=element_text(size=26,face="bold"),
		strip.text.x = element_text(size=18),
strip.text.y = element_text(size=18)) 


pretty_names <- c("Genus","Family","Order","Class","Phylum","Functional group")

### BOX PLOTS ###
# Summarized by rank - only significant effects
ggplot(data=sum.all[which(sum.all$fg_cat %in% c("Bacterial taxa","Fungal taxa") & sum.all$sig_strict == 1),],
		aes(x = only_rank,y = effSize))+
	geom_boxplot(aes(fill = as.factor(beta)), show.legend = F) +
	geom_jitter(width = .2, size = 3,) +
	labs(col = "Parameter", title = "Absolute effect size (significant effects only)") + 
	xlab("Taxonomic rank")+ 
	ylab(NULL)+
	#facet_grid(rows=vars(beta), scales = "free") + #,strip.position="bottom",nrow=2) +
	#facet_wrap(~beta, scales = "fixed",nrow=3) + #,strip.position="bottom",nrow=2) +
	facet_grid(rows = vars(fg_cat), cols = vars(beta), scales = "free", space = "free_x") + #,strip.position="bottom",nrow=2) +
			 	
	theme_bw() + theme(#axis.ticks.x=element_blank(),
		text = element_text(size = 16),
		axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
		#	angle = 320, vjust=-2, hjust = 1),
		axis.title=element_text(size=22,face="bold")
		#strip.text.y = element_text(size=24,hjust=0,vjust = 1,angle=180,face="bold")
		) + scale_x_discrete(labels = pretty_names)

# Summarized by rank - all effects
ggplot(data=sum.all[which(sum.all$fg_cat %in% c("Bacterial taxa","Fungal taxa")),],
			  aes(x = only_rank,y = effSize))+
			 	geom_boxplot(aes(fill = as.factor(beta)), show.legend = F) +
	geom_jitter(aes(shape = as.factor(sig_strict)), size = 3, width=.2) +
	labs(col = "Parameter", title = "Absolute effect size") + 
			 	xlab("Taxonomic rank")+ 
			 	ylab(NULL)+
			 	#facet_grid(rows=vars(beta), scales = "free") + #,strip.position="bottom",nrow=2) +
			 	#facet_wrap(~beta, scales = "fixed",nrow=3) + #,strip.position="bottom",nrow=2) +
			 	facet_grid(rows = vars(fg_cat), cols = vars(beta), scales = "free", space = "free_x") + #,strip.position="bottom",nrow=2) +
			 	
			 	theme_bw() + theme(#axis.ticks.x=element_blank(),
			 		text = element_text(size = 16),
			 		axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
			 		#	angle = 320, vjust=-2, hjust = 1),
			 		axis.title=element_text(size=22,face="bold")
			 		#strip.text.y = element_text(size=24,hjust=0,vjust = 1,angle=180,face="bold")
			 	) + scale_x_discrete(labels = pretty_names) + scale_shape_manual(values = c(21, 16), name = NULL, 
			 																																	 labels = c("Not significant","Significant")) 
			 
			 
			 
			 
			 
# View by parameter instead
ggplot(data=sum.all,
			 aes(x = beta,y = effSize)) +
	geom_jitter(aes(color = beta,shape = as.factor(sig_strict)), size = 3, width=.2) +
	labs(col = "Parameter", title = "Absolute effect size") + 
	xlab("Taxon")+ 
	ylab(NULL) +
	facet_wrap(~ pretty_name, nrow = 2) + 
	theme_bw() + theme(
		text = element_text(size = 16),
		axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
		axis.title=element_text(size=22,face="bold")
	) + scale_shape_manual(values = c(21, 16), name = NULL, 
												 labels = c("Not significant","Significant")) 


# Functional groups only
ggplot(data=fg_summary,
			 aes(x = fg_cat,y = effSize))+
	geom_boxplot(aes(fill = as.factor(beta)), show.legend = F) +
	geom_jitter(width = .2) +
	labs(col = "Parameter", title = "Absolute effect size (significant effects only)") + 
	xlab("Taxonomic rank")+ 
	ylab(NULL)+
	#facet_grid(rows=vars(beta), scales = "free") + #,strip.position="bottom",nrow=2) +
	facet_wrap(~beta, scales = "fixed", nrow=3) + #,strip.position="bottom",nrow=2) +
	theme_bw() + theme(#axis.ticks.x=element_blank(),
		text = element_text(size = 16),
		axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
		#	angle = 320, vjust=-2, hjust = 1),
		axis.title=element_text(size=22,face="bold")
		#strip.text.y = element_text(size=24,hjust=0,vjust = 1,angle=180,face="bold")
	)


# Only decomposers
ggplot(data=fg_summary[fg_summary$fg_cat=="Decomposers",],
			 aes(x = reorder(beta,effSize),y = effSize))+
	geom_point(aes(color = as.factor(beta), shape=taxon), show.legend = T, size = 5) +
	labs(col = "Parameter", title = "Decomposers: absolute effect size") + 
	xlab("Parameter")+ 
	ylab(NULL)+
	theme_bw() + theme(#axis.ticks.x=element_blank(),
		text = element_text(size = 22),
		axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
		axis.title=element_text(size=22,face="bold")
	)

ggplot(data=sum.all[sum.all$sig_lax == 1,],
       aes(x = cat,y = effSize))+
  #geom_point(aes(color = as.factor(sig_strict))) +
  geom_point(aes(color = as.factor(kingdom))) +
  labs(col = "Parameter", title = "Absolute effect size") + 
  xlab("Taxonomic rank")+ 
  ylab(NULL)+
  facet_wrap(~param,strip.position="bottom",nrow=2) +
  theme(axis.ticks.x=element_blank(),text = element_text(size = 24),
        axis.text.x=element_text(angle = 320,vjust=-1),
        axis.title=element_text(size=22,face="bold"),
        strip.text.y = element_text(size=24,hjust=0,vjust = 0,angle=180,face="bold")) 




sig.count <- sum.all %>% group_by(fg_cat, param) %>% count(sig_lax)
ggplot(sig.count[sig.count$sig_lax==1,], aes(cat, n)) + geom_point()+ facet_grid(~param)


temp <- sum.all[sum.all$param=="Temperature",]
temp <- sum.all[sum.all$param=="Precipitation",]
temp <- sum.all[sum.all$beta=="pH",]
temp <- sum.all[sum.all$param=="C:N",]

temp_sig <- temp[temp$sig_lax == 1,]
temp_sig <- temp[temp$sig_strict == 1,]
p <- ggplot(data=temp[temp$sig_strict == 1,],
       aes(x = rank,y = effSize))+
  geom_point(aes(color = as.factor(sig_strict))) +
  labs(col = "Parameter", title = "Absolute effect size") + 
  xlab("Taxonomic rank")+ 
  ylab(NULL)+
  theme(axis.ticks.x=element_blank(),text = element_text(size = 24),
        axis.text.x=element_text(angle = 320),
        axis.title=element_text(size=22,face="bold"),
        strip.text.y = element_text(size=24,hjust=0,vjust = 1,angle=180,face="bold")) 

library(agricolae)
# TUKEY PLOT
abs_max <- max(temp$effSize)
abs_max <- max(temp_sig$effSize)
# get the highest point for each class
maxs <- temp %>%
  group_by(rank) %>%
  summarise(tot=max(effSize) + 0.3 * abs_max)
Tukey_test_waic <- aov(effSize~rank, data=temp) %>%
  agricolae::HSD.test("rank", group=TRUE) %>%
  .$groups %>%
  as_tibble(rownames="rank") %>%
  rename("Letters_Tukey"="groups") %>% 
  dplyr::select(-"effSize") %>%
  left_join(maxs, by="rank") #%>% mutate(cat = "bacteria")
p + geom_text(data=Tukey_test_waic, aes(y = tot, label=Letters_Tukey))


tukey_out <- list()
for (b in unique(sum.all$beta)){
	temp <- sum.all[sum.all$beta==b,]
	maxs <- temp %>%
		group_by(rank) %>%
		summarise(tot=max(effSize) + 0.3 * abs_max)
	Tukey_test_waic <- aov(effSize~rank, data=temp) %>%
		agricolae::HSD.test("rank", group=TRUE) %>%
		.$groups %>%
		as_tibble(rownames="rank") %>%
		rename("Letters_Tukey"="groups") %>% 
		dplyr::select(-"effSize") %>%
		left_join(maxs, by="rank") 
	Tukey_test_waic$beta <- b
	tukey_out[[b]] <- Tukey_test_waic
}
tukey <- do.call(rbind, tukey_out)


p1 + geom_text(data=tukey, aes(y = tot, label=Letters_Tukey), size = 6)




tax.all <- sum.all[is.na(sum.all$fg_cat),]
p2 <- ggplot(data=tax.all,
						 aes(x = rank,y = effSize))+
	geom_boxplot(aes(fill = as.factor(beta)), show.legend = F) +
	geom_jitter(width = .2) +
	labs(col = "Parameter", title = "Absolute effect size") + 
	xlab("Taxonomic rank")+ 
	ylab(NULL)+
	#facet_grid(rows=vars(beta), scales = "free") + #,strip.position="bottom",nrow=2) +
	facet_wrap(~beta, scales = "fixed",nrow=3) + #,strip.position="bottom",nrow=2) +
	theme_bw() + theme(#axis.ticks.x=element_blank(),
		text = element_text(size = 16),
		axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
		#	angle = 320, vjust=-2, hjust = 1),
		axis.title=element_text(size=22,face="bold")
		#strip.text.y = element_text(size=24,hjust=0,vjust = 1,angle=180,face="bold")
	) + scale_x_discrete(labels = pretty_names) + geom_text(data=tukey, aes(y = tot, label=Letters_Tukey), size = 6)





gen <- sum.all[sum.all$rank=="genus_bac",]
ggplot(data=sum.all,
			 aes(x = taxon_name,y = effSize))+
	facet_grid(cols = vars(rank), scales="free_x") +
	geom_point(aes(color = as.factor(beta)), size = 4) +
	labs(col = "Parameter", title = "Absolute effect size") + 
	xlab("Genus")+ 
	ylab(NULL)+
	theme_bw() + theme(#axis.ticks.x=element_blank(),
		text = element_text(size = 16),
		axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
		#	angle = 320, vjust=-2, hjust = 1),
		axis.title=element_text(size=22,face="bold")
		#strip.text.y = element_text(size=24,hjust=0,vjust = 1,angle=180,face="bold")
	)

# BiocManager::install("metagenomeFeatures")
library(metagenomeFeatures)






dada_res
dada_res$seqtab.nochim 
