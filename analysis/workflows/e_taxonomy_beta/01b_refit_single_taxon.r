# Fit taxa one-by-one to evaluate convergence.

source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

#Get arguments from the command line
argv <- commandArgs(TRUE)
# Check if the command line is not empty and convert values to numerical values
if (length(argv) > 0){
	k <- as.numeric( argv[1] )
}

# Running dirichlet model
# Read in microbial data

d <- c(readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_16S_2021.rds"),
			 readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_ITS_2021.rds"))
# spec_names_out <- list()
# for(k in 1:10){
# # Subset to one rank.
# rank.name <- microbialForecast:::tax_names[k]
# rank.df <- d[[rank.name]]
# spec_names <- colnames(rank.df)[!colnames(rank.df) %in% c("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date","other")]
# spec_names_out[[k]] <- spec_names
# }
# names()


# Create function that calls run_MCMC for each uncertainty scenario
run_scenarios <- function(k, ...) {
	pacman::p_load(doParallel, reshape2)

	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
	#	source(here("functions", "MCMC_single_taxon.r"))

	# Read in microbial data
	d <- c(readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_16S_2021.rds"),
				 readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_ITS_2021.rds"))
	# Subset to one rank.
	rank.name <- microbialForecast:::tax_names[k]
	rank.df <- d[[rank.name]]
	spec_names <- colnames(rank.df)[!colnames(rank.df) %in% c("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date","other")]

	#start_no <- ifelse(k == 5, 23, 29)
#
#
# 	cal.file.list = intersect(list.files(here("data/model_outputs/single_taxon/all_covariates/"),recursive = T,
# 																			 pattern = "20151101_20180101"),
# 														list.files(here("data/model_outputs/single_taxon/all_covariates/"), recursive = T,
# 																			 pattern = "samples"))
# 	modeled_already <- list()
# 	for (i in seq_along(cal.file.list)) {
# 		info <- basename(cal.file.list[[i]]) %>% str_split("_") %>% unlist()
# 		rank.name <- info %>% head(-2) %>% tail(-1) #%>% paste0(collapse = "_")
# 		modeled_already[[i]] <- info %>% head(-2) %>% tail(-3)
# 	}
# 	modeled_already <- unlist(modeled_already)
# 	spec_names <- spec_names[!spec_names %in% modeled_already]

	cl <- makeCluster(28, type="PSOCK", outfile=paste0("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/analysis/workflows/qsub/", rank.name, "_201501_202001.txt"))
	registerDoParallel(cl)

	output.list = foreach(s=spec_names,
												.errorhandling = 'pass') %dopar% {
													source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
													print(s)
													out <- run_MCMC_single_taxon_bychain(k = k,
																											 #iter = 500000, burnin = 300000, thin = 10,
																											 #iter = 1001, burnin = 502, thin = 1,
																											 test=F,

																											 iter_per_chunk = 50000,
																											 init_iter = 100000, burnin = 50000,
																											 temporalDriverUncertainty = TRUE,
																											 spatialDriverUncertainty = TRUE,
																											 scenario ="201511_202001",
																											 model_name = "all_covariates",
																											 max.date = "20200101",
																											 chain_no = chain,
																											 s = s)
													return(out)
												}
	return(output.list)
	stopCluster(cl)
}


test_out <- run_scenarios(k)

# cl <- makeCluster(10, type="PSOCK", outfile="")
# registerDoParallel(cl)
# #Run for all 10 taxonomic ranks, in parallel (via PSOCK)
# output.list = foreach(k=c(1:10),
# 											.errorhandling = 'pass') %dopar% {
# 												run_scenarios(k)
# 											}
# stopCluster(cl)


