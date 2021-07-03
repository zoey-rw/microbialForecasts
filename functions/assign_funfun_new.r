# access funfun database, as well as funguild, and use both to assign fungal guilds
# modified from a function writted by Colin Averill, 2018

# funguild uses higher-up classifications (i.e. at the phylum level) 
# whereas funfun cleans up funguild and uses genus/species classifications
# this approach is quite messy, would love it if someone were to code a cleaner version...

#devtools::install_github("ropenscilabs/datastorr")
#devtools::install_github("traitecoevo/fungaltraits")
#devtools::install_github("brendanf/FUNGuildR")

assign_fungal_guilds <- function(tax_table, url = "http://www.stbates.org/funguild_db.php", fungaltraits_path = "/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/data/reference_data/FungalTraitsDB.csv", n.cores = NA){
  #check if dependencies are installed. If not, stop.----
	if (!require('fungaltraits',character.only = TRUE)){
		stop("please install the fungaltraits package from Github (traitecoevo/fungaltraits)")
	}
	if (!require('rvest',character.only = TRUE)){
    stop("please install the rvest package.")
  }
  if (!require('jsonlite',character.only = TRUE)){
    stop("please install the jsonlite package.")
  }
  if(!require('doParallel', character.only = TRUE)){
    stop('please install the doParallel package.')
  }
  #check that the input is formatted right. If not, stop, throw an error.
  if (!is.data.frame(tax_table)){
    stop('Your taxonomy table needs to be a data.frame. Try again.')
  }
	
	require(FUNGuildR)
	#setup parallel.----
	library(doParallel)
	if(is.na(n.cores)){
		n.cores <- detectCores()
	}
	#if number of cores is still NA for some reason set to 1.
	if(is.na(n.cores)){
		n.cores <- 1
		cat('detectCores() returned NA. Setting n.cores to 1.\n')
	}
	registerDoParallel(n.cores)

	tax_table_orig = tax_table
	
	# for testing only!
	tax_table = tax_table[1:500,]
	

	## Format taxonomy table
  #make sure tax column names are lower case.
  colnames(tax_table) <- tolower(colnames(tax_table))
  rows <- rownames(tax_table) # keep rownames for later
  tax_table <- as.data.frame(apply(tax_table,2,tolower))
  tax_table$species <- gsub("s__","",tax_table$species)
  # match tax_table species format to the format of funfun  
  tax_table$funfun_species <- paste0(tax_table$genus,"_",tax_table$species)
 
  # FUNGuild database
  tax_table_assign <- tax_table
  tax_table_assign <- cbind.data.frame(lapply(tax_table_assign, tools::toTitleCase))
  tax_table_assign <- tax_table_assign %>% unite("Taxonomy", 1:7, sep = ";")
  funguild_assigned <- funguild_assign(otu_table = tax_table_assign)
  
  # test_assign$notes <- NULL
  # test_assign$citationSource <- NULL
  
  #FUNFUN database
  # read in and reformat funfun database
  db <- fungaltraits::fungal_traits()
  funfun <- db
  colnames(funfun) <- tolower(colnames(funfun))
  funfun$genus <- tolower(funfun$genus)
  funfun$speciesmatched <- stringr::word(tolower(funfun$speciesmatched), 1, 2) 
  funfun$speciesmatched <- gsub(" ", "_", funfun$speciesmatched)
  funfun <- funfun[which(!is.na(funfun$guild_fg)),]
  
  #start with highest level of taxonomy and go down.
  #setup output list.
  out <- list() 
  out <-
  	foreach(i = 1:nrow(tax_table)) %dopar% {
  		to_return <- NA

			# FUNFUN   		
      #genus level match.
      if(tax_table$genus[i] %in% funfun$genus){
        to_return <- funfun[match(tax_table$genus[i], funfun$genus),c(2,3,52,58,59,61,72,99,101)]
      }
      #species level match.
      if(tax_table$funfun_species[i] %in% funfun$species){
      	to_return <- funfun[match(tax_table$funfun_species[i], funfun$species),c(2,3,52,58,59,61,72,99,101)]
      }
      if(tax_table$funfun_species[i] %in% funfun$speciesmatched){
      	to_return <- funfun[match(tax_table$funfun_species[i], funfun$speciesmatched),c(2,3,52,58,59,61,72,99,101)]
      }
      #return output.
      return(data.frame(to_return))
    } #end parallel loop.
  
  #bind up output
  out.all <- plyr::rbind.fill(out)
  

  
  
  
  # Now try adding the fungaltraits database
  fungaltraits <- read.csv(fungaltraits_path)
  fungaltraits$genus <- tolower(fungaltraits$GENUS)
  fungaltraits$guild <- paste(fungaltraits$primary_lifestyle, fungaltraits$Secondary_lifestyle)
  fungaltraits$guild_traits <- gsub("_", " ", fungaltraits$guild)
  
  out_traits <- foreach(i = 1:nrow(tax_table)) %dopar% {
  	to_return <- NA
  	#genus level match.
  	if(tax_table$genus[i] %in% fungaltraits$genus){
  		to_return <- fungaltraits[match(tax_table$genus[i], fungaltraits$genus),c(25, 27)]
  	}
  	#return output.
  	return(data.frame(to_return))
  } #end parallel loop.
  
  out_traits.all <- plyr::rbind.fill(out_traits)
  out_traits.all$to_return <- NULL
  
  
  # append output to taxonomy table.
  out.all.sources <- cbind(tax_table,
  												 funguild = funguild_assigned$guild, 
  												 funfun = out.all$guild_fg,
  												 funtraits = out_traits.all$guild_traits)

  # # Prioritize funfun over funguild
  # out.all$guild_assign <- ifelse(is.na(out.all$guild_fg), out.all$guild, out.all$guild_fg)
  
  #report and return output.
  cat(sum(!is.na(tax_table_out$guild))/(nrow(tax_table_out))*100,'% of taxa assigned a functional guild.\n', sep = '')
  return(tax_table_out)
}
