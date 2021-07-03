# Function written by Clara Qin as part of the neonMicrobe R package


#' Download NEON Soil Data Associated with Marker Gene Sequencing Data
#'
#' Downloads the following NEON data product:
#' - DP1.10086.001: "Soil physical and chemical properties, periodic", tables sls_soilCoreCollection, sls_soilMoisture, sls_soilpH, and sls_soilChemistry.
#' This function uses \code{\link[neonUtilities]{loadByProduct}} to conduct the downloads.
#'
#' @param sites Either the string 'all' (default), meaning all available sites, or a character vector of 4-letter NEON site codes, e.g. c('ONAQ','RMNP'). Defaults to PRESET_SITES parameter in params.R.
#' @param startYrMo,endYrMo Either NA (default), meaning all available dates, or a character vector in the form YYYY-MM, e.g. 2017-01. Defaults to PRESET_START_YR_MO in params.R.
#' @param dpID NEON data product(s) of interest. Default is DP1.10086.001 ("Soil physical and chemical properties, periodic").
#' @param outDir Default NEONMICROBE_DIR_SOIL() If a local copy of the filtered metadata is desired, provide path to output directory.
#' @param rmSamplingImpractical Default TRUE. Whether to remove soil data records when sampling did not actually occur.
#' @param rmNTransBouts Default TRUE. Whether to remove soil data records from bouts to collect N-transformation incubation tubes. These are only useful for calculating N-transformation rates and aren't associated with microbial data.
#' @param rmFailedCNDataQF Default TRUE. Whether to remove soil data records where cnPercentQF indicates failure. While other QF fields exist and are simply passed to output, this particular check may be desirable because this function later aggregates nitrogenPercent and organicCPercent values for cnSampleIDs with analytical replicates.
#'
#' @return If return_data==TRUE, returns a dataframe consisting of joined soil data records from DP1.10086 ("Soil physical and chemical properties, periodic"). Otherwise, no value is returned.
#' @export
#'
#' @importFrom magrittr "%>%"
downloadSoilData <- function(sites='all', startYrMo = NA, endYrMo = NA,
														 dpID = c("DP1.10086.001"), outDir=NEONMICROBE_DIR_SOIL(),
														 rmSamplingImpractical=TRUE, rmNTransBouts=TRUE,
														 rmFailedCNDataQF=TRUE, token_user = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJ6cndlcmJpbkBidS5lZHUiLCJzY29wZSI6InJhdGU6cHVibGljIiwiaXNzIjoiaHR0cHM6Ly9kYXRhLm5lb25zY2llbmNlLm9yZy8iLCJleHAiOjE3NTc4ODg2NzAsImlhdCI6MTYwMDIwODY3MCwiZW1haWwiOiJ6cndlcmJpbkBidS5lZHUifQ.8eW8vxUOiton-kQ_Xyvva0QSHD_BDd2E5IGeNKW3WHib-m7UpTnEhGFAUlAHGdsUyz-dKE1jMOAGS5A_NRYXGg") {
	if(!dir.exists(outDir)) {
		message("Output directory does not exist. Returning NULL.")
		return(NULL)
	}
	
	library(dplyr)
	library(neonUtilities)
	
	# check valid data values entered
	## validate dpID ##
	if(!all(grepl("DP1", dpID) & grepl('\\.001', dpID) & grepl('10086', dpID))) {
		message("Invalid Data Product ID: must follow convention 'DP1.[5-digit value].001' and must be a distributed periodic soil data product ID. (DP1.10078.001 is no longer supported.)")
		return(NULL)
	} else {
		dpID <- dpID
	}
	
	# validate site(s)
	terrSiteList <- c("all","HARV","SCBI","OSBS","GUAN","UNDE","KONZ","ORNL","TALL","WOOD","CPER","CLBJ","YELL","NIWO",
										"SRER","ONAQ","WREF","SJER","TOOL","BONA","PUUM","BART","BLAN","SERC","SCBI","DSNY","JERC","LAJA",
										"TREE","STEI","KONA","UKFS","MLBS","GRSM","LENO","DELA","NOGP","DCFS","STER","RMNP","OAES","MOAB",
										"JORN","ABBY","TEAK","SOAP","BARR","DEJU","HEAL")
	if(!any(sites %in% terrSiteList)){
		message("Invalid site(s): must be a valid NEON site or 'all'")
		return(NULL)
	} else {
		sites <- sites
	}
	
	slsL1 <- list()
	message("loading soil data...")
	for(i in 1:length(dpID)) {
		slsL1[[dpID[i]]] <- tryCatch({
			neonUtilities::loadByProduct(dpID[i], sites, package = 'expanded', check.size = F, startdate = startYrMo, enddate = endYrMo, token = token_user) # output is a list of lists of each soil data file
		}, error = function(e) {
			warning("No data was found for data product ", dpID[i], " at the specified sites and dates.")
			NA
		})
	}
	
	# Columns by which to join sls tables
	joining_cols <- c("domainID", "siteID", "plotID", "sampleID")
	
	tables_available <- data.frame(available=rep(FALSE, 4), row.names=c("sls_soilCoreCollection", "sls_soilMoisture", "sls_soilpH", "sls_soilChemistry"))
	
	warnIfColsMissing <- function(table, cols, table_name) {
		if(!all(keep_cols %in% names(table))) {
			not_available <- which(!(keep_cols %in% names(table)))
			warning("NEON soil data products have changed, and the following columns are no longer available in '", table_name, "': ",
							paste(keep_cols[not_available], collapse=", "),
							". Please file an Issue on the GitHub page: https://github.com/claraqin/neonMicrobe/issues")
		}
	}
	
	# If DP1.10086.001 is found
	if(any(grepl('10086', dpID)) & !all(is.na(slsL1[["DP1.10086.001"]]))) {
		tables_available$available[1] <- TRUE
		
		# start with soilCoreCollection data...
		keep_cols <- c("domainID", "siteID", "plotID", "namedLocation", "plotType", "nlcdClass", "coreCoordinateX", "coreCoordinateY", "geodeticDatum",
									 "decimalLatitude", "decimalLongitude", "elevation", "samplingProtocolVersion", "collectDate", "sampleTiming", "standingWaterDepth",
									 "nTransBoutType", "boutType", "samplingImpractical", "sampleID", "horizon", "soilTemp", "litterDepth", "sampleTopDepth",
									 "sampleBottomDepth", "soilSamplingDevice", "geneticSampleID", "dataQF")
		warnIfColsMissing(slsL1[["DP1.10086.001"]]$"sls_soilCoreCollection", keep_cols, "sls_soilCoreCollection")
		
		dat_soil <- dplyr::select(slsL1[["DP1.10086.001"]]$"sls_soilCoreCollection", all_of(keep_cols)) %>%
			dplyr::rename(sccSamplingProtocolVersion = samplingProtocolVersion,
										sccDataQF = dataQF)
		
		# Remove records that were not actually sampled, if desired
		if(rmSamplingImpractical) {
			dat_soil <- dplyr::filter(dat_soil, samplingImpractical == "OK" | is.na(samplingImpractical) )
		}
		# Remove records that were only used to measure N transformations, if desired
		if(rmNTransBouts) {
			dat_soil <- dplyr::filter(dat_soil, boutType != "fieldOnly")
		}
		
		# merge with soilMoisture data...
		if(!is.null(slsL1[["DP1.10086.001"]]$"sls_soilMoisture")) {
			tables_available$available[2] <- TRUE
			
			keep_cols <- c("moistureSampleID", "samplingProtocolVersion", "soilMoisture", "smDataQF")
			warnIfColsMissing(slsL1[["DP1.10086.001"]]$"sls_soilMoisture", keep_cols, "sls_soilMoisture")
			
			dat_soil <- merge(
				dat_soil,
				dplyr::select(slsL1[["DP1.10086.001"]]$"sls_soilMoisture", all_of(c(joining_cols, keep_cols))) %>%
					dplyr::rename(smSamplingProtocolVersion = samplingProtocolVersion),
				by = joining_cols,
				all.x = TRUE
			)
		}
		
		# merge with soilpH data...
		if(!is.null(slsL1[["DP1.10086.001"]]$"sls_soilpH")) {
			tables_available$available[3] <- TRUE
			
			keep_cols <- c("pHSampleID", "samplingProtocolVersion", "soilInWaterpH", "soilInCaClpH", "pHDataQF")
			warnIfColsMissing(slsL1[["DP1.10086.001"]]$"sls_soilpH", keep_cols, "sls_soilpH")
			
			dat_soil <- merge(
				dat_soil,
				dplyr::select(slsL1[["DP1.10086.001"]]$"sls_soilpH", all_of(c(joining_cols, keep_cols))) %>%
					dplyr::rename(pHSamplingProtocolVersion=samplingProtocolVersion),
				by = joining_cols,
				all.x = TRUE
			)
		}
		
		# if available, merge with soil chemistry data
		if(!is.null(slsL1[["DP1.10086.001"]]$"sls_soilChemistry")) {
			tables_available$available[4] <- TRUE
			
			keep_cols <- c("cnSampleID", "nitrogenPercent", "organicCPercent", "CNratio",
										 "testMethod", "instrument", "cnPercentQF")
			warnIfColsMissing(slsL1[["DP1.10086.001"]]$"sls_soilChemistry", keep_cols, "sls_soilChemistry")
			
			# Need to collapse nitrogenPercent and organicCPercent values into the same rows. (They come in separate rows.)
			soilchem <- dplyr::select(slsL1[["DP1.10086.001"]]$"sls_soilChemistry", all_of(c(joining_cols, keep_cols))) %>%
				dplyr::rename(cnTestMethod=testMethod, cnInstrument=instrument) %>%
				tidyr::pivot_longer(c(nitrogenPercent, organicCPercent)) %>%
				dplyr::filter(!is.na(value))
			if(rmFailedCNDataQF==TRUE) {
				soilchem <- dplyr::filter(soilchem, cnPercentQF == "OK" | is.na(cnPercentQF))
			}
			# Aggregate measurements of analytical replicates
			soilchem_measurements <- tidyr::pivot_wider(soilchem, id_cols = domainID:cnSampleID,
																									names_from = name, names_sort = TRUE,
																									values_from = value, values_fn = mean)
			# Aggregation function to use for next pivot_wider calls
			indicateIfMixedAggregation <- function(x, value_if_mixed) {
				if(length(unique(x)) == 1) {
					return(x[1])
				} else {
					return(value_if_mixed)
				}
			}
			soilchem_methods <- tidyr::pivot_wider(soilchem, id_cols = domainID:cnSampleID, names_from = name,
																						 names_sort = TRUE, names_glue="{name}TestMethod",
																						 values_from = cnTestMethod, values_fn = function(x) {
																						 	indicateIfMixedAggregation(x, "Aggregated from mixed methods")
																						 })
			soilchem_QF <- tidyr::pivot_wider(soilchem, id_cols = domainID:cnSampleID, names_from = name,
																				names_sort = TRUE, names_glue="{name}QF",
																				values_from = cnPercentQF, values_fn = function(x) {
																					indicateIfMixedAggregation(x, "Aggregated from mixed quality flags")
																				})
			soilchem <- merge(merge(soilchem_measurements, soilchem_methods, all.x=TRUE), soilchem_QF, all.x=TRUE)
			# Collapsing complete
			# Now merge with the rest of soil data
			dat_soil <- merge(dat_soil, soilchem, by = joining_cols, all.x = TRUE
			)
			
		}
		# If DP1.10086.001 is not found
	} else {
		dat_soil <- NULL
	}
	
	if(!is.null(dat_soil)) {
		message("Soil data availability for specified sites and dates:")
		print(tables_available)
		message("Returning soil data. Note that this function does not return the entire data product(s);\n",
						"for specialized analyses, downloading directly from the NEON Data Portal, the NEON Data API,\n",
						"or from neonUtilities may be necessary. Lastly, check the 'QF' (quality flag) columns to ensure\n",
						"you are only keeping records of sufficient quality for your purposes.")
	} else {
		warning("No soil data available at the specified sites and dates. Returning NULL.")
	}
	
	# download local copy if user provided output dir path
	if(outDir != "") {
		if(!dir.exists(outDir)) {
			dir.create(outDir)
		}
		write.csv(dat_soil, paste0(outDir, "/sls_soilData_", Sys.Date(), ".csv"),
							row.names=F)
		message(paste0("Soil data downloaded to: ", outDir, "/sls_soilData_", Sys.Date(), ".csv") )
	}
	
	return(dat_soil)
}


#' Alias for downloadSoilData()
#'
#' This function name is deprecated.
#'
#' @seealso \code{\link{downloadSoilData}}
#'
#' @export
downloadRawSoilData <- function(sites='all', startYrMo = NA, endYrMo = NA,
																dpID = c("DP1.10086.001"), outDir=NEONMICROBE_DIR_SOIL(),
																rmSamplingImpractical=TRUE, rmNTransBouts=TRUE,
																rmFailedCNDataQF=TRUE) {
	warning("downloadRawSoilData() is now downloadSoilData().")
}



##############################################################################################
#' @title Calculate more precise geolocations for specific NEON data products

#' @author 
#' Claire Lunch \email{clunch@battelleecology.org}

#' @description 
#' Calculation Function. Refine the geolocation data associated with NEON data products, based on product-specific rules and spatial designs.
#' 
#' @param data A data frame containing NEON named locations and other sampling information. For reliable results, use data tables as downloaded from the NEON data portal or API.
#' @param dataProd The table name of the NEON data product table to find locations for. Must be one of: ltr_pertrap, hbp_perbout, sls_soilCoreCollection, brd_perpoint or brd_countdata, mam_pertrapnight, div_1m2Data or div_10m2Data100m2Data, vst_mappingandtagging.
#' @param token User specific API token (generated within neon.datascience user accounts). Optional.

#' @return A data frame of geolocations for the input product and data

#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

#' @keywords Currently none

#' @examples 
#' d <- data.frame(namedLocation="GUAN_044.basePlot.ltr", subplotID=23, trapID="GUAN_044_385")
#' getLocTOS(d, "ltr_pertrap")

#' @seealso Currently none

#' @export

# changelog and author contributions / copyrights
#   Claire Lunch (2019-09-03)
#     adapted from and replaces def.calc.geo.os()
##############################################################################################
getLocTOS <- function(
	data,
	dataProd,
	token=NA_character_
){
	
	# convert format for safety
	data <- data.frame(data)
		
		data$rowid <- 1:nrow(data)
		
		# Use the getLocByName function to pull the plot geolocations from the API
		locCol <- "namedLocation"
		plot.all <- geoNEON::getLocByName(data, locCol=locCol, locOnly=T, token=token)
		
		# Use relevant columns
		plot.merg <- plot.all[,c("namedLocation","utmZone",
														 "northing","easting","namedLocationCoordUncertainty",
														 "decimalLatitude","decimalLongitude",
														 "elevation","namedLocationElevUncertainty")]
		colnames(plot.merg) <- c(locCol, 'utmZone',"adjNorthing","adjEasting",
														 "adjCoordinateUncertainty","adjDecimalLatitude",
														 "adjDecimalLongitude","adjElevation",
														 "adjElevationUncertainty")
		plot.loc <- base::merge(data, plot.merg, by=locCol, all.x=T)
		plot.loc <- plot.loc[order(plot.loc$rowid),]
		
		# Subtract 20 meters from the easting and northing values to get the 
		# location of the southwest corner
		plot.loc$adjEasting <- as.numeric(plot.loc$adjEasting) - 20
		plot.loc$adjNorthing <- as.numeric(plot.loc$adjNorthing) - 20
		
		# Add coreCoordinateX to the easting value and coreCoordinateY to the northing value
		plot.loc$adjEasting <- plot.loc$adjEasting + data$coreCoordinateX
		plot.loc$adjNorthing <- plot.loc$adjNorthing + data$coreCoordinateY
		
		# Set the coordinate uncertainty to 0.5 meter
		plot.loc$adjCoordinateUncertainty <- 0.5
		
		# SOMEthing in this command is causing crashes!!!
		# calculate latitude and longitude from the corrected northing and easting
		adjLatLong <- geoNEON::calcLatLong(easting=plot.loc$adjEasting, 
																			 northing=plot.loc$adjNorthing,
																			 utmZone=plot.loc$utmZone)
		plot.loc$adjDecimalLatitude <- adjLatLong$decimalLatitude
		plot.loc$adjDecimalLongitude <- adjLatLong$decimalLongitude
		
		# reorder to original order
		all.return <- plot.loc[order(plot.loc$rowid),]
		all.return <- all.return[,!names(all.return) %in% c('rowid')]
		
		return(all.return)
	}
