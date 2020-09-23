####### Prep data for Boston Trees Shiny app #######
library(googlesheets)
  
  #### prepare tree-characteristic data ####
  gs <- gs_title("Ecological Considerations")
  app_df <- gs %>% gs_read("Application.input.2")
  app_df <- data.frame(app_df)
  
  # format canopy spread column
  spread <- app_df$Spread..ft.
  spread[spread %in% c("50-60","50-70","50-75")] <- "Large"
  spread[spread %in% c("35-50","40-50","40-60","45-60","40-70")] <- "Med-Large"
  spread[spread %in% c("30-35","30-40","35-40")] <- "Med-Small"
  spread[spread %in% c("15-20","15-25","20-25","25-30","25-35")] <- "Small"
  app_df$spread <- spread
  
  # format heat-reduction column
  good_heat_reduc <- c("Good (top 20%)", "Good (top 10%)", "Good (top 30%)")
  app_df$heat_reduc <- 0
  app_df$heat_reduc[which(app_df$UHI %in% good_heat_reduc)] <- 1
  app_df$heat_reduc_icon[app_df$heat_reduc == 1] <- '<i class="fas fa-temperature-low" style="font-size:20px;color:blue;"></i> <br> Heat-reducing tree'

  # create "more information" column
  app_df$more.info <- paste0("<a href='", app_df$Link, "'>Fact sheet</a>")
  app_df$old.Description <- app_df$Description
  app_df$Description <- gsub("\\n\\n", "<br>", app_df$old.Description)
  
  # create bitmap images since direct html isn't working for local images
  img_uri <- function(x) { sprintf('<img src="%s" style="width:100px;height:100px;"/>', knitr::image_uri(x)) }
  filename <- paste0('www/', app_df$Sci.name,'.PNG')
  for (i in 1:length(app_df$Sci.name)){ 
    if(file.exists(filename[i])) {
      app_df$Scientific.name[i] <- paste(img_uri(filename[i]), app_df$Sci.name[i])
    } else {
      app_df$Scientific.name[i] <- paste(filename[i], app_df$Sci.name[i])
    }
  }
  # save output
  saveRDS(app_df, "data/treeData.rds")

