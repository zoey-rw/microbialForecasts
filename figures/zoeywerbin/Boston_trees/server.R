server <- function(input, output, session) {
  
  #### TAB 1: MAPS OF TEMPERATURE OR HVI, WITH PRIORITY REGIONS OVERLAID, AND CURRENT/POTENTIAL TREE MAP
  
  observeEvent(input$goToTab2, {
    updateTabsetPanel(session, "inTabset",
                      selected = "other_considerations")
  })
  observeEvent(input$goToTab3, {
    updateTabsetPanel(session, "inTabset",
                      selected = "choose_tree")
  })
  observeEvent(input$goToTab4, {
    updateTabsetPanel(session, "inTabset",
                      selected = "keep_healthy")
  })
  
  output$temp_hvi <- renderLeaflet({
    tracts %>%  
      leaflet() %>%
      addTiles() %>%  
      setView(-71.08, 42.3, zoom = 11) %>%
      addProviderTiles("CartoDB.Positron")  %>% 
      clearControls()
  })
  
  observeEvent(input$choose_temp_hvi, {
    if(input$choose_temp_hvi == "temp"){
      var <- var1
      pal <- pal1
      highlight <- "priority_heat_trees"
      markers <- gCentroid(tract_data[tract_data[[highlight]]==T,],byid=TRUE)
      title <- "Mean morning <br>summer <br>temperature (C)"
    } else if(input$choose_temp_hvi == "HVI"){
      var <- var2
      pal <- pal2
      highlight <- "priority_hvi_trees"
      markers <- gCentroid(tract_data[tract_data[[highlight]]==T,],byid=TRUE)
      title <- "Heat Vulnerability Index"
    }
    
    leafletProxy("temp_hvi") %>% 
      clearControls() %>%
      clearShapes() %>% 
      clearMarkers() %>% 
      addPolygons(data = tract_data , 
                  fillColor = ~pal(var), 
                  fillOpacity = 0.7, 
                  weight = 0.2, 
                  smoothFactor = 0.2,
                  popup = ~tract_data$neighborhood) %>%
      addLegend(pal = pal, 
                values = var, 
                position = "bottomright", 
                title = title)  %>%  
      addMarkers(lng = markers$x, 
                 lat = markers$y,
                 icon = exclamation.point,
                 label = "Priority region")
  })
  
  
  output$current_trees <- renderLeaflet({
    tracts %>%
      leaflet() %>% addTiles() %>%
      setView(-71.08, 42.3, zoom = 12) %>%
      addProviderTiles("CartoDB.Positron") %>%
      addPolygons(data = tract_data ,
                  fillColor = ~pal3(var3),
                  fillOpacity = 0.7,
                  weight = 0.2,
                  smoothFactor = 0.2,
                  popup = ~tract_data$popups,
                  group = "Current tree cover (%)") %>%
      addPolygons(data = tract_data ,
                  fillColor = ~pal4(var4),
                  fillOpacity = 0.7,
                  weight = 0.2,
                  smoothFactor = 0.2,
                  popup = ~tract_data$popups,
                  group = "Potential tree cover (%)") %>%
      addLayersControl(
        baseGroups = c("Current tree cover (%)", "Potential tree cover (%)"),
        options = layersControlOptions(collapsed = FALSE)) %>%
      hideGroup("Potential tree cover (%)")
  })
  
  observeEvent(input$current_trees_groups, {
    current <- leafletProxy("current_trees") %>% clearControls()
    if (input$current_trees_groups == 'Current tree cover (%)'){
      current <- current %>% addLegend(pal = pal3,
                                       values = var3,
                                       opacity = .7,
                                       position = "bottomright",
                                       title = "Current <br>tree cover (%)")
    } else if (input$current_trees_groups == 'Potential tree cover (%)'){
      current <- current %>% addLegend(pal = pal4,
                                       values = var4,
                                       opacity = .7,
                                       labels = labelFormat(6),
                                       position = "bottomright",
                                       title = "Potential <br>tree cover (%)")
    }
  })
  
  ### TAB 2: OTHER CONSIDERATIONS
  
  output$other_considerations <- renderLeaflet({
    tracts %>%  
      leaflet() %>% addTiles() %>%  
      setView(-71.08, 42.3, zoom = 12) %>%
      addProviderTiles("CartoDB.Positron")
  })
  
  output$other_considerations_overlay <- renderPlot({
    leafletProxy("other_considerations") %>%
      clearControls() %>%
      clearShapes() %>% 
      addPolygons(data = tract_data , 
                  fillColor = ~pal1(var1), 
                  fillOpacity = 0.7, 
                  weight = 0.2, 
                  smoothFactor = 0.2, 
                  popup = ~tract_data$popups) %>%
      addPolygons(data = parcel[parcel$PTYPE %in% c(902),],
                  color = "#444444", weight = 1, smoothFactor = 0.5,
                  opacity = 1.0, fillOpacity = 1,
                  fillColor = "green",
                  group = "City-owned property") %>% 
      addPolygons(data = parcel[parcel$PTYPE %in% c(901, 910:929),],
                  color = "#444444", weight = 1, smoothFactor = 0.5,
                  opacity = 1.0, fillOpacity = 1,
                  fillColor = "green",
                  group = "State-owned property") %>% 
      addPolygons(data = parcel[parcel$PTYPE %in% c(900),],
                  color = "#444444", weight = 1, smoothFactor = 0.5,
                  opacity = 1.0, fillOpacity = 1,
                  fillColor = "green",
                  group = "Federal-owned property") %>% 
      addPolygons(data = parcel[parcel$PTYPE %in% c(903, 908,973, 986, 965, 978, 984),],
                  color = "#444444", weight = 1, smoothFactor = 0.5,
                  opacity = 1.0, fillOpacity = 1,
                  fillColor = "green",
                  #group = "BHA") %>% 
                  group = "Boston Housing Authority, <br>Boston Redevelopment Authority, <br>or other government-associated or public land") %>% 
      addLayersControl(overlayGroups = 
                         c("City-owned property", "State-owned property", 
                           "Federal-owned property", #"BHA"),
                           "Boston Housing Authority, <br>Boston Redevelopment Authority, <br>or other government-associated or public land"),
                       options = layersControlOptions(collapsed = FALSE)) %>% 
      addLegend(pal = pal1, 
                values = var1, 
                position = "bottomright", 
                title = "Mean summer <br>temperature (C)") 
  })
  
  
  #### TAB 3: TREE-SELECTION ####
  
  # create one reactivevalue that stores the output from each question
  approved_trees <- reactiveValues(size_trees=0, 
                                   allergen_trees=0, 
                                   site_trees=0, 
                                   light_trees=0,
                                   resistant_trees=0)
  
  #### CANOPY SIZE ####
  
  observeEvent(input$canopySize, {
    approved_size <- list()
    if ("Small" %in% input$canopySize) {
      approved_small <- app_df[which(app_df$spread=="Small"|
                                       app_df$spread=="Med-Small"),]
      approved_size <- append(approved_size, approved_small$Common.Name)
    }
    if ("Medium" %in% input$canopySize) {
      approved_medium <- app_df[which(app_df$spread=="Med-Large"|
                                        app_df$spread=="Med-Small"),]
      approved_size <- append(approved_size, approved_medium$Common.Name)
    }
    if ("Large" %in% input$canopySize) {
      approved_large <- app_df[which(app_df$spread=="Large"|
                                       app_df$spread=="Med-Large"),]
      approved_size <- append(approved_size, approved_large$Common.Name)
    }
    if ("NoInput" %in% input$canopySize) {
      approved_size <- append(approved_size, app_df$Common.Name)
    }
    approved_trees$size_trees <- unique(approved_size)
  })
  
  #### ALLERGENS ####  
  
  observeEvent(input$allergens, {
    approved_allergens <- list()
    if ("Low" %in% input$allergens) {
      low <- app_df[which(app_df$Allergens=="low"),]
      approved_allergens <- append(approved_allergens, low$Common.Name)
    }
    if ("Moderate" %in% input$allergens) {
      moderate <- app_df[which(app_df$Allergens=="low"|
                                 app_df$Allergens=="moderate"),]
      approved_allergens <- append(approved_allergens, moderate$Common.Name)
    }
    if ("NoInput" %in% input$allergens) {
      approved_allergens <- append(approved_allergens, app_df$Common.Name)
    }
    approved_trees$allergen_trees <- unique(approved_allergens)
  })
  
  #### SITE TYPE ####
  
  observeEvent(input$siteType, {
    approved_site <- list()
    if ("Park_yard" %in% input$siteType) {
      park <- app_df[which(app_df$Site.type=="park/yard" | app_df$Site.type=="either"),]
      approved_site <- append(approved_site, park$Common.Name)
    }
    if ("Street" %in% input$siteType) {
      street <- app_df[which(app_df$Site.type=="street" | app_df$Site.type=="either"),]
      approved_site <- append(approved_site, street$Common.Name)
    }
    if ("NoInput" %in% input$siteType) {
      approved_site <- append(approved_site, app_df$Common.Name)
    }
    approved_trees$site_trees <- unique(approved_site)
  })
  
  #### LIGHT AVAILABILITY ####
  
  observeEvent(input$lightReqs, {
    approved_light <- list()
    if ("FullShade" %in% input$lightReqs) {
      shade <- app_df[which(app_df$Light.requirement %in% 
                              c("Full shade", "Full sun/full shade", "Any")),]
      approved_light <- append(approved_light, shade$Common.Name)
    }
    if ("PartShade" %in% input$lightReqs) {
      part <- app_df[which(app_df$Light.requirement %in% 
                             c("Full sun/full shade", "Part shade/part sun/full sun", "Any")),]
      approved_light <- append(approved_light, part$Common.Name)
    }
    if ("FullSun" %in% input$lightReqs) {
      sun <- app_df[which(app_df$Light.requirement %in% 
                            c("Full sun", "Full sun/full shade", "Part shade/part sun/full sun", "Any")),]
      approved_light <- append(approved_light, sun$Common.Name)
    }
    if ("NoInput" %in% input$lightReqs) {
      approved_light <- append(approved_light, app_df$Common.Name)
    }
    approved_trees$light_trees <- unique(approved_light)
  })
  
  #### BREAKAGE ####
  
  observeEvent(input$resistant, {
    if (input$resistant) {
      approved_breakage <-
        unique(app_df[which(app_df$Resistant.to.Breakage == "yes"), ]$Common.Name)
    } else {
      approved_breakage <- unique(app_df$Common.Name)
    }
  })
  
  output$treeRecommendations <- renderUI({
    recommended_trees <-Reduce(intersect, list(approved_trees$size_trees,
                                               approved_trees$allergen_trees,
                                               approved_trees$site_trees,
                                               approved_trees$light_trees
    ))
    df <- app_df[which(app_df$Common.Name %in% recommended_trees),]
    to_display <- df[,c("heat_reduc_icon", "Scientific.name", "Common.Name", "Description", "more.info")]
    dt <- DT::datatable(to_display, escape = FALSE, filter="none", rownames = FALSE, 
                        selection = 'none', colnames=NULL,
                        options = list( columnDefs = list(
                            list(width = '8%', targets = 0), # heat-reduction icon
                            list(className = 'dt-center', targets = 0), # heat-reduction icon
                            list(width = '8%', targets = 1), # scientific.name column 
                            list(width = '10%', targets = 2), # common.name column 
                            list(width = '72%', targets = 3), # description column 
                            list(width = '10%', targets = 4) # source column 
                          ), #list(visible=FALSE, targets=4)), # hide heat-reduction column
                          scrollX=T, sDom  = '<"top">lrt<"bottom">ip', # remove search options
                          pageLength = 3, lengthChange = FALSE, bSort=FALSE)) 
    output$recommendations <- DT::renderDT({ dt
    })
    DTOutput("recommendations")
  }) # close tree recommendations renderUI
}
