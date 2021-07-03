# access funfun database, as well as funguild, and use both to assign fungal guilds
# modified from a function writted by Colin Averill, 2018

# funguild uses higher-up classifications (i.e. at the phylum level) 
# whereas funfun cleans up funguild and uses genus/species classifications
# this approach is quite messy, would love it if someone were to code a cleaner version...

#devtools::install_github("ropenscilabs/datastorr")
#devtools::install_github("traitecoevo/fungaltraits")
#devtools::install_github("brendanf/FUNGuildR")

assign_fungal_guilds <- function(tax_table, url = "http://www.stbates.org/funguild_db.php", n.cores = NA){
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
	
	# read in funfun database
	db <- fungaltraits::fungal_traits()
	funfun <- db
  
	
	## Format taxonomy table
  #make sure tax column names are lower case.
  colnames(tax_table) <- tolower(colnames(tax_table))
  rows <- rownames(tax_table) # keep rownames for later
  tax_table <- as.data.frame(apply(tax_table,2,tolower))
  tax_table$species <- gsub("s__","",tax_table$species)
  # match tax_table species format to the format of funfun  
  tax_table$funfun_species <- paste0(tax_table$genus,"_",tax_table$species)
  
  #download FUNGuild database, convert it to something R interpretable.----
  # fg <- url %>% 
  #   xml2::read_html() %>%
  #   rvest::html_text() 
  # fg <- jsonlite::fromJSON(gsub("funguild_db", "", fg))
  # fg$taxon <- tolower(fg$taxon)
  
  # Alternative to JSON method using metagMisc package
  fg <- metagMisc::parse_funguild(url = url, tax_name = TRUE)
  fg$taxon <- tolower(fg$taxon)
  
  
  tax_table_assign <- tax_df
  tax_table_assign <- cbind.data.frame(lapply(tax_table_assign, tools::toTitleCase))
  tax_table_assign <- tax_table_assign %>% unite("Taxonomy", 1:7, sep = ";")
  test_assign <- funguild_assign(otu_table = tax_table_assign)
  
  
  #start with highest level of taxonomy and go down.
  colnames(funfun) <- tolower(colnames(funfun))
  funfun$genus <- tolower(funfun$genus)
  funfun$speciesmatched <- stringr::word(tolower(funfun$speciesmatched), 1, 2) 
  funfun$speciesmatched <- gsub(" ", "_", funfun$speciesmatched)
  funfun <- funfun[which(!is.na(funfun$guild_fg)),]
  
  #setup output list.
  out <- list() 
  
  
  
  out <-
  	foreach(i = 1:nrow(tax_table)) %dopar% {
  		to_return <- NA
  		#phylum level match.
  		if(tax_table$phylum[i] %in% fg$taxon){
  			to_return <- fg[match(tax_table$phylum[i], fg$taxon),4:10]
  		}
  		#class level match.
  		if(tax_table$class[i] %in% fg$taxon){
  			to_return <- fg[match(tax_table$class[i], fg$taxon),4:10]
  		}
  		#order level match.
  		if(tax_table$order[i] %in% fg$taxon){
  			to_return <- fg[match(tax_table$order[i], fg$taxon),4:10]
  		}
  		#family level match.
  		if(tax_table$family[i] %in% fg$taxon){
  			to_return <- fg[match(tax_table$family[i], fg$taxon),4:10]
  		}
  		#genus level match.
  		if(tax_table$genus[i] %in% fg$taxon){
  			to_return <- fg[match(tax_table$genus[i], fg$taxon),4:10]
  		}
  		#species level match.
  		if(tax_table$species[i] %in% fg$taxon){
  			to_return <- fg[match(tax_table$species[i], fg$taxon),4:10]
  		}
  #run parallel fg assign loop.
  # out <-
  #   foreach(i = 1:nrow(tax_table)) %dopar% {
  #     to_return <- NA
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
  
  #bind up output, append to taxonomy table.
  out.all <- plyr::rbind.fill(out)
  out.all$guild <- ifelse(is.na(out.all$guild_fg), out.all$guild, out.all$guild_fg)
  # out <- do.call(plyr::rbind.fill,out)
  # out <- do.call(rbind,out)
  tax_table_out <- cbind.data.frame(tax_table, guild = out.all$guild)
  #report and return output.
  cat(sum(!is.na(tax_table_out$guild))/(nrow(tax_table_out))*100,'% of taxa assigned a functional guild.\n', sep = '')
  #cat(sum(!is.na(tax_table$guild_fg))/(nrow(tax_table))*100,'% of taxa assigned a functional guild.\n', sep = '')
  return(tax_table_out)
}