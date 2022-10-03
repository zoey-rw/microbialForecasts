# Fit taxa one-by-one to evaluate convergence.


# # For testing
iter <- 500
burnin <- 200
thin <- 1
test =F
k = 1
temporalDriverUncertainty <- TRUE
spatialDriverUncertainty <- TRUE
scenario <- "full_uncertainty"
s = "acidobacteriota"



#Get arguments from the command line
argv <- commandArgs(TRUE)
# Check if the command line is not empty and convert values to numerical values
if (length(argv) > 0){
	k <- as.numeric( argv[1] )
}

# Running dirichlet model 

	
# Create function that calls run_MCMC for each uncertainty scenario
run_scenarios <- function(k, ...) {
	pacman::p_load(doParallel, reshape2) 
	
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")
	source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/functions/MCMC_single_taxon.r")
#	source(here("functions", "MCMC_single_taxon.r"))
	
	# Read in microbial data
	d <- c(readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_16S_2021.rds"), 
				 readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_ITS_2021.rds"))
	# Subset to one rank.
	rank.name <- tax_names[k]
	rank.df <- d[[rank.name]] 
	spec_names <- colnames(rank.df)[!colnames(rank.df) %in% c("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date","other")] 
	
	#start_no <- ifelse(k == 5, 23, 29)
	
	
	cl <- makeCluster(16, type="PSOCK", outfile="")
	registerDoParallel(cl)
	
	output.list = foreach(s=spec_names,.export = c("run_MCMC", "k"),
	.errorhandling = 'pass') %dopar% {

		print(s)
	out <- run_MCMC(k = k, 
									iter = 1000000, burnin = 500000, thin = 10, 
									#iter = 50001, burnin = 20002, thin = 3, 
									test=F, 
									temporalDriverUncertainty = TRUE, 
									spatialDriverUncertainty = TRUE,
									scenario ="full_uncertainty",
									s = s)
	return(out)
	}
	return(output.list)
	stopCluster(cl)
}


test_out <- run_scenarios(k, run_MCMC)

# cl <- makeCluster(10, type="PSOCK", outfile="")
# registerDoParallel(cl)
# #Run for all 10 taxonomic ranks, in parallel (via PSOCK)
# output.list = foreach(k=c(1:10),
# 											.errorhandling = 'pass') %dopar% {
# 												run_scenarios(k)
# 											}
# stopCluster(cl)


