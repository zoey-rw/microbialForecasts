#' @title run_MCMC_single_taxon_bychain
#' @description run_MCMC_single_taxon_bychain
#' @export
#
# # For testing
# iter <- 500
# burnin <- 200
# thin <- 1
# test =F
# k = 1
# temporalDriverUncertainty <- TRUE
# spatialDriverUncertainty <- TRUE
# scenario <- "full_uncertainty"
# s = "acidobacteriota"
# model_name = "all_covariates"
# chain_no = 1
# min.date = "20151101"
# max.date = "20200101"
# max.date = "20180101"
# iter_per_chunk = 500
# init_iter = 1000
# s = "acidobacteriota"
# model_name = "cycl_only"

run_MCMC_single_taxon_bychain <- function(k = 1,
																					iter = 1000,  burnin = 500, thin = 1,
                                  test = F, chain_no = 1,
                                  temporalDriverUncertainty = TRUE, spatialDriverUncertainty = TRUE,
                                  scenario=NULL,
                                  s = "acidobacteriota",
                                  min.date = "20151101",
                                  max.date = "20180101",
                                  model_name = "cycl_only",

																					iter_per_chunk = 5000,
																					init_iter = 10000,
                                  ...) {
  pacman::p_load(reshape2, parallel, nimble, coda, tidyverse)

  # Subset to one rank.
  rank.name <- microbialForecast:::tax_names[k]

  # Read in microbial abundances
  cal <- c(readRDS(here("data", "clean/cal_groupAbundances_16S_2021.rds")),
           readRDS(here("data", "clean/cal_groupAbundances_ITS_2021.rds")))
  val <- c(readRDS(here("data", "clean/val_groupAbundances_16S_2021.rds")),
           readRDS(here("data", "clean/val_groupAbundances_ITS_2021.rds")))

  cal.rank.df <- cal[[rank.name]]
  val.rank.df <- val[[rank.name]]
  rank.df <- rbind(cal.rank.df, val.rank.df)

  # Reduce size for testing
  if (test == T) {
    rank.df <- rank.df %>% arrange(siteID, plotID, dateID)
    rank.df = rank.df[1:500,]
  }

  # Prep model inputs/outputs.
  print(paste0("Preparing model data for ", rank.name))

  spec_names <- colnames(rank.df)[!colnames(rank.df) %in% c("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date","other")]
  rank.df_spec <- rank.df %>%
    select("siteID", "plotID", "dateID", "sampleID", "dates", "plot_date", !!s)
  rank.df_spec$other <- 1-rank.df_spec[[s]]

  model.dat <- prepTaxonomicData(rank.df = rank.df_spec,
                                 min.prev = 3,
                                 min.date = min.date,
                                 max.date = max.date)

  out.path <- here("data",paste0("model_outputs/single_taxon/",model_name,"/samples_", rank.name, "_", s, "_",min.date,"_",max.date, ".rds"))


  print(paste("Completed model data for", rank.name, ", group:", s))
  constants <- model.dat[c("plotID",  "timepoint","plot_site", "site_start", "plot_start", "plot_index",
                           "plot_num", "plot_site_num",
                           "N.plot", "N.spp", "N.core", "N.site", "N.date",
                           "mois", "mois_sd", "temp", "temp_sd", "pH", "pH_sd",
                           "pC", "pC_sd", "LAI", "relEM", "sin_mo", "cos_mo")]


  if (model_name == "cycl_only"){
    constants$N.beta = 2
    Nimble_model = nimbleModTaxa_cycl_only
    constants <- constants[c("plotID",  "timepoint","plot_site", "site_start", "plot_start", "plot_index",
                             "plot_num", "plot_site_num",
                             "N.plot", "N.spp", "N.core", "N.site", "N.date", "N.beta","sin_mo", "cos_mo")]
  } else if(model_name == "all_covariates"){
    constants$N.beta = 8
    Nimble_model = nimbleModTaxa
  } else message("Missing specification of Nimble model.")

  constants$omega <- 0.0001 * diag(constants$N.spp)
  constants$zeros = rep(0, constants$N.spp)
  if (constants$N.spp < 8) {
    constants$omega <- 0.0001 * diag(8)
    constants$zeros = rep(0, 8)
  }
  if (constants$N.spp < 8) {
    constants$omega <- 0.0001 * diag(8)
    constants$zeros = rep(0, 8)
  }



  # for output.
  metadata <- list("rank.name" = rank.name,
                   "niter" = iter,
                   "nburnin" = burnin,
                   "thin" = thin,
                   "model_data" = model.dat$truth.plot.long)

  inits <- initsFun(constants, type = "tax")
  #
  # inits <- inits[c("y", "plot_mu", "intercept", "sig", "beta", "rho", "sigma",
  # 								 "plot_rel", "site_effect")]
  ## Configure & compile model
  Rmodel <- nimbleModel(code = Nimble_model,
                        constants = constants, data = list(y=model.dat$y),
                        inits = inits)
  #Rmodel$Ex <- inits$Ex
  #Rmodel$site_effect <- rep(0, constants$N.site)
  #Rmodel$y <- inits$y
  # Compile model
  cModel <- compileNimble(Rmodel)
  nimbleOptions(multivariateNodesAsScalars = TRUE)
  # Configure & compile MCMC
  mcmcConf <- configureMCMC(cModel, monitors = c("beta","sigma","site_effect",
                                                 "sig","intercept",
                                                 "rho"),
                            monitors2 = c("plot_rel"), thin2 = 25,
                            useConjugacy = T)

  myMCMC <- buildMCMC(mcmcConf)
  compiled <- compileNimble(myMCMC, project = Rmodel, resetFunctions = TRUE)


  compiled$run(niter=init_iter, thin=thin, thin2 = 25, nburnin = burnin)
  cat(paste0("\nInitial run finished for chain", chain_no))

  # Sample from MCMC

  out.run<-as.mcmc(as.matrix(compiled$mvSamples))
  out.run2<-as.mcmc(as.matrix(compiled$mvSamples2))

  out.path2 <- gsub(".rds", paste0("_chain", chain_no, ".rds"), out.path)

  saveRDS(list(samples = out.run, samples2 = out.run2, metadata = metadata),
          out.path2)

  continue <- check_continue(out.run, min_eff_size = 5)
  loop_counter = 0
  while (continue){
    message("\nEffective sample size too low; running for another ", iter_per_chunk, " iterations\n")
  	if (loop_counter < 100) {

    compiled$run(niter=iter_per_chunk, thin=thin, reset=F)

    # Shorten if more than 10k samples have accumulated
    out.run = window_chain(compiled$mvSamples, max_size = 20000)
    out.run2 <- window_chain(compiled$mvSamples2, max_size = 20000)
    out.run2 <- na.omit(out.run2)

    # Remove timepoints that had no sampling at all (only zeros)
    out.run2 <- mcmc(out.run2[, which(colSums(out.run2) != 0)])

    saveRDS(list(samples = out.run, samples2 = out.run2, metadata = metadata),
            out.path2)
    continue <- check_continue(out.run, min_eff_size = 5)
    loop_counter = loop_counter + 1
  	} else {
  		message("Exceeded 100 loops. No more sampling!")
  		continue <- FALSE
  	}
  }

  return("ok")
}
