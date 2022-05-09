# Calculate percentage of points within 95% credible intervals

crps_in <- readRDS("./data/summary/CRPS_hindcasts.rds")
observed <- crps_in$scored_hindcasts %>% filter(!is.na(crps) & newsite == "Observed site")

observed$in_95 <- NA
observed$in_95 <- ifelse(observed$truth < observed$hi & 
												 	observed$truth > observed$lo,
																1, 0)



overdispersion <-  observed %>%
	group_by(fcast_type, pretty_group, model_name, newsite, pretty_name, taxon) %>% 
	mutate(countT= n()) %>%
	mutate(in_95_count = 
				 	sum(in_95, na.rm=T)) %>% 
	mutate(in_95_pct = in_95_count/countT) 

overdispersion <- overdispersion %>% 
	select(fcast_type, pretty_group, model_name, newsite, pretty_name, taxon, in_95_pct) %>% distinct()



dispersion_ribbons=data.frame(x1=c(0,0), x2=c(8,8), y1=c(0,.95), y2=c(.95,1), t=c('overdispersion','underdispersion'), r=c(1,2))

ggplot(overdispersion) + 
	geom_point(aes(x = pretty_name, y = in_95_pct, color = pretty_group), 
						 alpha = 1, size=4,
						 position=position_jitterdodge(dodge.width = 1)) + 
	facet_grid(rows=vars(model_name), drop=T, scales = "free") +  
	ylab("% observations in 95% CI") + xlab(NULL) + 
	theme_minimal(base_size=20) + ggtitle("Overdispersion at observed sites") + 
	theme(text = element_text(size = 22),
				axis.text.x=element_text(angle = 320, vjust=1, hjust = -0.05),
				axis.title=element_text(size=24), legend.position = c(.9,1)) + 
	guides(color=guide_legend(title=NULL)) +
	geom_rect(data=dispersion_ribbons, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=t), alpha=0.5) +
	geom_text(data=dispersion_ribbons, aes(x=8, y=1.1, label=r), size=4)
	# 
	# 
	# annotate("rect", ymin = 0, ymax =.95, xmin=0, xmax=8,
	# 				 alpha = .1,fill = "blue") +
	# annotate("rect", ymin = .95, ymax =1, xmin=0, xmax=8,
	# 				 alpha = .1,fill = "red") + annotate()
	# 
