# Create distance matrix of all NEON soil samples
source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
options(scipen=999)
library(ecodist)
library(metagMisc)

# devtools::install_github("vmikk/metagMisc")

# Path to soil sample info (from cohesion analysis)
phys.path = "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/old/data/soil_phys_allsites.rds"
require(stringr)
require(vegan)
require(spatstat)
require(fields)

soil_phys <- readRDS(phys.path)

	# dec.degrees.mat<-as.matrix(cbind(soil_phys$adjDecimalLongitude, soil_phys$adjDecimalLatitude))
	# rownames(dec.degrees.mat)<-soil_phys$sampleID
	# #create great circle matrix
	# distance.matrix<-rdist.earth(dec.degrees.mat[,1:2],miles=FALSE)
	# distance.mat <- reshape2::melt(distance.matrix)
	# colnames(distance.mat) <- c("soilSample1","soilSample2","spatial_distance")
	# distance.mat$soilSample1 <- as.character(distance.mat$soilSample1)
	# distance.mat$soilSample2 <- as.character(distance.mat$soilSample2)
	# distance.mat <- distance.mat[which(!is.na(distance.mat$spatial_distance)),]
	# distance.mat <- distance.mat[!duplicated(distance.mat),]
	# distance.mat <- distance.mat[distance.mat$spatial_distance != 0,]
	# distance.mat$siteID_samp1 <- substr(distance.mat$soilSample1, 1, 4)
	# distance.mat$siteID_samp2 <- substr(distance.mat$soilSample2, 1, 4)
	#
	#
	# saveRDS(distance.mat, here("data/clean", "soil_phys_allsites.rds"))





	# Read in microbial abundances
	cal <- c(readRDS(here("data", "clean/cal_groupAbundances_16S_2021.rds")),
					 readRDS(here("data", "clean/cal_groupAbundances_ITS_2021.rds")))

	k <- 1
	i <- 8
	moran.stat_rank <- list()
	for (k in 1:length(cal)){
		rank.name <- names(cal)[[k]]
		print(rank.name)
		species_df <- cal[[k]]

		n_spec <- ncol(cal[[k]])
			# Add lat/lon to taxon df
		species_df <- species_df[which(!is.na(species_df$sampleID)),]
		species_df$without_horizon <- parseNEONsampleIDs(species_df$sampleID)$without_horizon
			soil_phys$without_horizon <- parseNEONsampleIDs(as.character(soil_phys$sampleID))[,"without_horizon"]
			species_df$lat <- soil_phys[match(species_df$without_horizon, soil_phys$without_horizon),]$adjDecimalLatitude
			species_df$lon <- soil_phys[match(species_df$without_horizon, soil_phys$without_horizon),]$adjDecimalLongitude
			species_df$lat_lon <- paste(species_df$lat, species_df$lon)
			species_df <- species_df[!duplicated(species_df$lat_lon),]
			species_df <- species_df %>% drop_na()

			timepoints <- unique(species_df$dateID)
			moran.date.stat <- list()

			for (dateID in timepoints) {
				print(dateID)

d <- species_df[species_df$dateID==dateID,]
			# Create lat/lon df, check for matchup
			dec.degrees.mat<-as.matrix(cbind(d$lon, d$lat))
			rownames(dec.degrees.mat)<-d$sampleID
			identical(d$sampleID, rownames(dec.degrees.mat))

			# Approach from colin's Moran's I calculation
			tax.dist <- geosphere::distm(dec.degrees.mat)
			inv.tax.dist <- 1/tax.dist
			diag(inv.tax.dist) <- 0

			moran.stat <- list()


			keep_list <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/unconverged_taxa_list.rds")
			spec_names <- keep_list[[rank.name]]$taxon.name


			for(i in 1:length(spec_names)){
				if (sd(d[,i])==0 | nrow(d)==1) next()
			moran.stat[[i]]  <- ape::Moran.I(d[,i], inv.tax.dist) %>% rbind %>% as.data.frame() %>%
				mutate(taxon = names(d)[i])
			}
			# bind taxa for each date
			moran.date.stat[[dateID]] <-  do.call(rbind, moran.stat) %>% as.data.frame() %>%
				mutate(dateID = !!dateID)
			}
			# bind dates for each rank
			moran.stat_rank[[k]] <- do.call(rbind, moran.date.stat) %>% as.data.frame() %>%
				mutate(rank = names(cal)[[k]])
	}
	# bind ranks
	moran.stat_all_rank <- do.call(rbind, moran.stat_rank) %>% as.data.frame()

	moran.stat_all_rank <- apply(moran.stat_all_rank, 2, unlist) %>% cbind.data.frame() %>%  mutate_at(c('observed', 'expected', 'sd', 'p.value'), as.numeric) %>% rename(spatial_autocorrelation = observed)

	moran.stat_all_rank <- moran.stat_all_rank %>% group_by(rank, taxon) %>%
		summarize(mean_morans= mean(spatial_autocorrelation, na.rm=T),
							sd_morans= sd(spatial_autocorrelation, na.rm=T),)

	saveRDS(moran.stat_all_rank, here("data/clean", "moran_stat.rds"))



	moran.stat_all_rank = readRDS(here("data/clean", "moran_stat.rds"))



	cl <- makeCluster(36, outfile="")
	registerDoParallel(cl)

	parms_rank = foreach(k = 1:length(cal), .errorhandling = "pass") %dopar% {
		source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
		source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/microbialForecast/R/helperFunctions.r")
		library(ecodist)

	#for (k in 1:length(cal)){

		rank.name <- names(cal)[[k]]
		print(rank.name)
		species_df <- cal[[k]]

		n_spec <- ncol(cal[[k]])

		d <- species_df
		d <- d[!duplicated(d$sampleID),]
		d <- d[!is.na(d$sampleID),]


		# SPATIAL DISTANCE MATRIX
		# Create lat/lon df, check for matchup
		# Add lat/lon to taxon df
		d$without_horizon <- parseNEONsampleIDs(d$sampleID)$without_horizon
		soil_phys$without_horizon <- parseNEONsampleIDs(as.character(soil_phys$sampleID))[,"without_horizon"]
		d$lat <- soil_phys[match(d$without_horizon, soil_phys$without_horizon),]$adjDecimalLatitude
		d$lon <- soil_phys[match(d$without_horizon, soil_phys$without_horizon),]$adjDecimalLongitude
		d$lat_lon <- paste(d$lat, d$lon)
		d <- d[!duplicated(d$lat_lon),]
		d <- d %>% drop_na()


		dec.degrees.mat<-as.matrix(cbind(d$lon, d$lat))
		rownames(dec.degrees.mat)<-d$sampleID
		identical(d$sampleID, rownames(dec.degrees.mat))

		# Approach from colin's Moran's I calculation
		tax.dist <- geosphere::distm(dec.degrees.mat)
		inv.tax.dist <- 1/tax.dist
		diag(inv.tax.dist) <- 0
		colnames(inv.tax.dist) <- d$sampleID
		rownames(inv.tax.dist) <- d$sampleID
		inv.space.melted <- reshape2::melt(inv.tax.dist) %>% rename(space = "value")


		# isolate sampling dates, put in order
		sampleID <- rownames(d)
		siteID <- substr(sampleID, 1, 4)
		YYYYMMDD <- str_extract(sampleID, "\\d\\d\\d\\d\\d\\d\\d\\d")
		date <- as.Date(YYYYMMDD, "%Y%m%d")
		metadata <- cbind.data.frame(siteID, sampleID, YYYYMMDD, date)
		metadata <- metadata[order(metadata$date),]


		# distance.mat <- distance.mat[!duplicated(distance.mat),]

		# create matrix of date differences
		# INTER-ANNUAL DISTANCE MATRIX
		days <- as.numeric(metadata$date)
		i <- outer(days, days, '-') # get days since previous sample
		i[upper.tri(i)] <- NA
		date.df <- as.data.frame(i)
		date.df$time2 <- metadata$date
		colnames(i) <- d$sampleID
		rownames(i) <- d$sampleID
		days.melted <- reshape2::melt(i)  %>% rename(interAnnual = "value")

		# INTRA-ANNUAL DISTANCE MATRIX
		doy <- lubridate::yday(date)
		i <- outer(doy, doy, '-') # get days since previous sample
		i[upper.tri(i)] <- NA
		doy.df <- as.data.frame(i)
		doy.df$time2 <- metadata$date
		colnames(i) <- d$sampleID
		rownames(i) <- d$sampleID
		doy.melted <- reshape2::melt(i) %>% rename(intraAnnual = "value")


		parms_taxon <- list()


		#keep_list <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/summary/unconverged_taxa_list.rds")

		if (grepl("_fun", rank.name) | grepl("_bac", rank.name) ) {
			#spec_names <- keep_list[[rank.name]]$taxon.name

			spec_names <- microbialForecast:::rank_spec_names[[rank.name]]
			} else spec_names <- rank.name

		for(j in 1:length(spec_names)){
		#for(i in 7:n_spec){
taxon = spec_names[[j]]
print(taxon)

			spec_dist <- dist(d %>% select(!!taxon))
			spec_dist_mat <- as.matrix(spec_dist)
			spec_dist_mat[upper.tri(spec_dist_mat)] <- NA
			colnames(spec_dist_mat) <- d$sampleID
			rownames(spec_dist_mat) <- d$sampleID
			spec_dist_melt <- reshape2::melt(spec_dist_mat) %>% rename(taxon_dist = "value")

			# Something in this line is slowing things down
			mat_melted <- list(spec_dist_melt, doy.melted, days.melted, inv.space.melted) %>%
				reduce(inner_join) %>% drop_na %>% filter(!(intraAnnual == 0 & interAnnual == 0 & space == 0 & taxon_dist == 0))

			m <- MRM((taxon_dist) ~ (intraAnnual) + (interAnnual) + (space), data = mat_melted); m
			parms <- as.data.frame(m$coef); parms
			#use lm to get standard errors of parameter estimates.
			lm.m <- lm((taxon_dist) ~ (intraAnnual) + (interAnnual) + (space), data = mat_melted)
			#grab standard error, multiply by sqrt N
			parms$sd <- summary(lm.m)$coefficients[,2] * sqrt(nrow(d))
			parms$param <- rownames(parms)
			parms <- rbind(parms, c(m$r.squared[1], m$r.squared[2], NA, "overall_Rsq"))
			parms$taxon =  taxon

			parms_taxon[[j]] = parms
		}

		return(do.call(rbind, parms_taxon) %>% as.data.frame() %>%
					 	mutate(rank = names(cal)[[k]]))
		# parms_rank[[k]] = do.call(rbind, parms_taxon) %>% as.data.frame() %>%
		# 	mutate(rank = names(cal)[[k]])
		# return(do.call(rbind, parms_rank) %>% as.data.frame())

	}

	# parms_rank_all = c(parms_rank_fun, parms_rank)
	# parms_MRM = map(parms_rank_all, possibly(as.data.frame, otherwise = NULL)) %>% rbindlist(fill = T)


	parms_MRM = map(parms_rank, possibly(as.data.frame, otherwise = NULL)) %>% rbindlist(fill = T)
	parms_MRM_wide = parms_MRM %>%  pivot_wider(id_cols = c("rank","taxon"), names_from = param, values_from = `(taxon_dist)`)

	moran.stat_all_rank <- moran.stat_all_rank %>% filter(taxon != "other")
	moran.MRM <- merge(moran.stat_all_rank, parms_MRM_wide, all=T)
	saveRDS(moran.MRM, here("data/clean", "moran_stat_MRM.rds"))

