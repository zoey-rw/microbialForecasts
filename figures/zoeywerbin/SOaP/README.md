# SOaP
GE585 Spring 2019: "Soil Organism chAnging Predictably" microbial functional group forecast

Name: Yetianjian Wang
Email: wangytj@bu.edu

Name: Steve Gougherty
Email: gougher@bu.edu

Name: Ryan Quinn 
Email: rkq@bu.edu 

Name: Zoey Werbin
Email: zrwerbin@gmail.com

### Updated 2019/05/07

### This repo contains a graduate student forecasting project created for EE585 (Ecological Forecasting, Dietze) at Boston University.

### Brief project description:
We created hindcasts for the ratios of bacteria/archaea to fungi in the soil ("B:F ratios"). Our Bayesian dynamic linear models were fit to 5 sites from the National Ecological Observatory Network (NEON), using data from 2013-2018 and the JAGS software. We aimed to include soil covariates such as pH and litter depth, but NEON's very early (2013-2014) data was a bit too sparse for that (code to download it is still in the repo, though). We used precipitation and air temperature ensembles from the ERA5 climate reanalysis dataset to model B:F ratios, but model selection via DIC indicated that air temperature was not a covariate worth keeping, so precipitation is the one covariate in the model. Using 2013-2014 as calibration data, we then forecasted to later dates; as we incorporated new data, we re-fit the model (no advanced data-assimilation happening here, though we aimed to implement a particle filter). Our scripts create gifs for the forecasts for the 5 sites.

### How to run our scripts

1. Run all scripts within the `data_construction/` directory, starting with the `01` scripts and then the `02` scripts. `01_download_ERA5_ensembles.R` cannot be run from any computer other than BU's SCC, but the outputs are already in the `data/` directory.
2. Model scripts can be run in any order (a lot of the forecasting code repeats between them). `runForecast.R` is the main one - this will fit the calibration model to the earlier data, and forecast forward, saving GIFs of the output (in which every frame is a different timestep). It will probably take 10-20 minutes to run, but you can change the loop to only look at one site. `variancePartitioning.R` creates forecasts for STER (the NEON site with the most data), and partitions forecast uncertainty into parameter, process, driver, and initial condition uncertainty. `modelSelection.R` and `predictedVsObserved.R` are also diagnostic plots. The functionality of `calibrationModel.R` is subsumed by the forecast script and by the function `fitModel.R` (in the `functions/` directory).

Note: the `old/` subdirectory has a lot of files that we didn't end up using for our final analysis, but might someday want to revisit. This includes code to download and visualize Canopy Height Model data (a remote sensing product from NEON sites), daymet weather data, worldClim historical climate data, and NEON ITS taxonomic group abundance data.
