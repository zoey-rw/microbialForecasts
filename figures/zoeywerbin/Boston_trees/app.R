# Shiny web application for choosing among trees that the City of Boston has approved.
# Written by Zoey Werbin, Mar 2019
# May or may not be maintained, email zoeywerbin@gmail.com with any issues (or use Github issues)

library(shiny)
library(DT)
suppressMessages(library(dplyr))
library(tidyverse)
library(bsplus)
library(leaflet) # for map base
library(htmlwidgets) # css fix for map
library(shinyWidgets)
library(shinythemes)
library(shinydashboard)

#redownloadData = TRUE
source("data_prep.R")


# Application title
header <-  dashboardHeader(title = "Right Place, Right Tree | Boston, MA",
                           titleWidth = 530
) # close header panel

sidebar <- dashboardSidebar(
  
  sidebarMenu(
    menuItem("Step 1: Choose region", tabName = "priority", icon = icon("map-marked-alt")),
    menuItem(tags$div("Step 2: Learn about ",tags$br(), "regional considerations", style = "display: inline-block; vertical-align: middle;"), icon = icon("question-circle"), tabName = "other_considerations"),
  
    menuItem(tags$div("Step 3: Choose ",tags$br(), "the right tree", style = "display: inline-block; vertical-align: middle;"), icon = icon("tree"), tabName = "choose_tree"),
    menuItem(tags$div("Step 4: Keep your tree ", tags$br(), "healthy!", style = "display: inline-block; vertical-align: middle;"), icon = icon("check-circle"), tabName = "keep_healthy")
    ) # close sidebarmenu
) # close sidebar

body <-   dashboardBody(
  # tags$head(
  #   tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
  # ),
  h3("Did you know that Boston's urban canopy has the potential to shape our local climate? Street trees can keep us cooler in the summer, reducing the 'urban heat island' effect, and they can keep us warmer in the winter by moderating wind speeds. But planting the wrong tree in the wrong place can lead to roots lifting up sidewalks, tree branches falling into streets, and other unintended consequences. Use this tool to help decide which species of tree to plant and where to plant them.", style = "margin-top: 0;
                    margin-block-start: 0px; "),
  use_bs_tooltip(),
  use_bs_popover(),
  tags$head(
    HTML('<link rel="stylesheet" href="https://use.fontawesome.com/releases/v5.8.1/css/all.css" integrity="sha384-50oBUHEmvpQ+1lW4y57PTFmhCaXp0ML5d60M1M7uH2+nqUivzIebhndOJK28anvf" crossorigin="anonymous">'),
    tags$style(HTML("@import url('//fonts.googleapis.com/css?family=Josefin+Sans');
                    /* logo */
                    .skin-blue .main-header .logo {
                    background-color: #1366bf;
                    font-family: Josefin Sans; 
                    /* font-weight: bold; */
                    font-size: 30px;
                    }
                    
                    /* logo when hovered */
                    .skin-blue .main-header .logo:hover {
                    background-color: #1366bf;
                    }
                    
                    /* navbar (rest of the header) */
                    .skin-blue .main-header .navbar {
                    background-color: #1366bf;
                    }        ")),
    tags$style(type='text/css', "* { font-family: Josefin Sans; font-size: 16px;}" ),
    tags$style("td {padding:4px 4px !important;}"),
    tags$style("div.info.legend.leaflet-control br {clear: both;}"),
    tags$style("h2 {margin-top: 5px;}")
    ),
  
  tabItems(
    tabItem(tabName = "priority",
             fluidRow(
               box(width = 12,
                   h2("Step 1. Choose a region for planting trees."),
                   "Where can new trees provide the largest benefits? Skip this step if you already know where you're planting!")),
      fluidRow(
        box(width=12, height = 800,
                        "Use the maps below to explore regions that may benefit from trees.",
                        tags$b("'Priority' regions (bordered in black) also have high potential for canopy expansion. Scroll down to see regions where other heat interventions may be necessary."),
                        fluidRow(
                          column(6,
                                 h3("High summer temperatures"),
                                 "Regions with darker red colors have higher summer land surface temperatures. The 'priority' regions by temperature are bordered in black. Expanding tree canopy in hot regions can help lower temperatures, reducing heat-related morbitidity and mortality."),
                          column(6,
                                 h3("High vulnerability to heat"),
                                 "Regions with darker red colors have higher scores on the Heat Vulnerability Index (HVI), which incorporates factors such as old housing, proportion of seniors and children, and income. The 'priority' regions by heat vulnerability are bordered in black. Expanding tree canopy in these regions can provide health benefits and reduce heat-related morbitidity and mortality.")),
                        fluidRow(
                          column(6,
                          leafletOutput("temp_trees", height = 500)), # temperature with priority areas highlighted
                          column(6,
                          leafletOutput("priority_trees", height = 500)))
            ) # close box of tree-priority areas 
      ), # close row of tree-priority areas
        fluidRow(
          box(width=12, height = 800,
              h2("Alternatives to tree-planting recommended"),
                        
              "These maps highlight regions that have low potential for canopy expansion, but still have high summer temperatures and high heat-vulnerability. Areas with low potential for canopy expansion may benefit from other interventions, such as green roofing, retrofitting houses with cooling units, and creating emergency cooling centers.",
              fluidRow(
                        column(6,
                               h3("High summer temperatures"),
                               "Regions with darker red colors have higher summer land surface temperatures. The hottest regions are bordered in black. Expanding tree canopy in hot regions can help lower temperatures, reducing heat-related morbitidity and mortality."),
                        column(6,
                               h3("High vulnerability to heat"),
                               "Regions with darker red colors have higher scores on the Heat Vulnerability Index (HVI), which takes into account factors such as old housing, proportion of seniors and children, and income. The most vulnerable regions are bordered in black. Expanding tree canopy in these regions can provide significant health benefits and reduce heat-related morbitidity and mortality.")),
          fluidRow(
                        column(6,
                               leafletOutput("temp_other", height = 500)), # temperature with priority areas highlighted
                        column(6,
                               leafletOutput("priority_other", height = 500))))) # HVIwith priority areas highlighted
      ), # close priority tab
#pests and disease map: https://usfs.maps.arcgis.com/apps/webappviewer/index.html?id=52cb2bcc3c2b4868ac87b66f622062ab
    tabItem(tabName = "other_considerations",
            fluidRow(
              box(width = 12,
                  h2("Step 2. Learn about regional considerations"),
                  "When planting a tree, you should carefully consider its context, such as community characteristics, land ownership, and local pests or diseases. Is asthma prevalent in the region? Maybe consider planting a low-allergen tree. Are there community organizations you could partner with to ensure maintenance? Contact them! " )),
            fluidRow(
              tabBox(width=12, height = 650,
                     tabPanel("Existing canopy cover in Boston",
                              "This map displays the existing canopy cover in each Boston census tract, as of 2016. Contrast this with the map of potential tree canopy. Click on the map to learn more about a region.",
                              leafletOutput("current_trees", height = 500)),
                     tabPanel("Potential canopy cover in Boston",
                              "This map displays the potential space in each Boston census tract available for immediately planting additional trees, as of 2016. This includes land that is currently lawn or open park, but does not include impervious surfaces such as concrete, which require additional work before planting. Certain regions on this map, such as Logan Airport, may not be reasonable targets for additional trees due to other considerations.",
                              leafletOutput("potential_trees", height = 500))
              ) # close tabBox
            ) # close canopy map row
            ), # close tabItem

    tabItem(tabName = "choose_tree",
            fluidRow(
              box(width = 12,
                  h2("Step 3. Choose the right tree species"),
                  "Enter your site characteristics to view a list of trees that may be perfect for your site. If planting multiple trees, plant a mixture of species - this makes our urban forest more resilient to future pests and diseases." )),
  fluidRow(
    bs_modal(id = "about_sizing", title = "Determining size requirements for trees", body = "When planting a tree, you should take into account the expected canopy spread at the tree's maturity, as well as its height and the extent of its roots. Trees can upheave sidewalks or building foundations if planted too closely. Mature trees should be at least X feet from the sides of buildings, X feet from the corners of buildings, and X feet from sidewalks. Choose a canopy size that allows for at least this amount of space."),
    bs_modal(id = "about_allergies", title = "Planting with allergies in mind", body = HTML("<p>Allergic reactions (also known as 'hay fever') are often more frequent and more severe among sufferers of chronic asthma. If planting in a census tract with high chronic asthma prevalence, consider planting trees that are less likely to trigger allergic reactions. Click on a census tract in this application's interactive map to see the chronic asthma prevalence in your neighborhood. </p> <p>Individuals can have reactions to various types of plants, but allergies to airborne pollen from pine trees (rather than insect-pollinated trees) are common. Firs, cedars, and many fruit trees produce very low amounts of pollen, while elm, oak, zelkova and liquidambar produce higher amounts of pollen.</p> <p>Other resources for tree allergies:</p> <p> <a href='https://weather.com/forecast/allergy/l/USMA0046:1:US'> Allergy forecast from Weather.com </a></p>")),
    bs_modal(id = "about_site_type", title = "What are the requirements of different site types?", body = "Street trees have different requirements than those planted in parks or spacious yards. Street trees can create significant hazards when not maintained regularly, and they must be more tolerant to pollutants from cars (such as ozone). If power lines are present, smaller ornamental trees should be planted rather than larger shade trees."),
    bs_modal(id = "about_light", title = "How much light does a tree need?", body = "Here are some guidelines:"),
    box(width=12,
        #box(
        column(width = 3,          
               h4("Filter by tree characteristics:", style = "margin-top: 0px;"),
               radioButtons(inputId = "canopySize", label = "Canopy spread",
                            choices = c("Small (15-35 ft)" = "Small",
                                        "Medium (30-70 ft)" = "Medium",
                                        "Large (50-75 ft)" = "Large",
                                        "Doesn't matter to me" = "NoInput"),
                            selected = "NoInput") %>%
                 shinyInput_label_embed(
                   shiny_iconlink() %>%
                     bs_attach_modal(paste0("about_sizing"))),
               #content = "<a href='https://www.arborday.org/trees/righttreeandplace/size.cfm'> Basic sizing guide </a>")),
               radioButtons(inputId = "allergens", label = "Allergen preferences", 
                            choices = c("Low allergen trees" = "Low",
                                        "Moderate or low allergen trees" = "Moderate",
                                        "Doesn't matter to me" = "NoInput"), 
                            selected = "NoInput") %>%
                 shinyInput_label_embed(
                   shiny_iconlink() %>%
                     bs_attach_modal(paste0("about_allergies"))),
               # bs_embed_popover(placement = "right", title = "Seasonal allergies?",
               #                  content = tags$a(href="https://www.nysinuscenter.com/2017/06/best-and-worst-trees-for-allergies/", "Best and worst trees for allergies"))),
               radioButtons(inputId = "siteType", label = "Site type", 
                            choices = c("Park or yard" = "Park_yard",
                                        "Street" = "Street",
                                        "Doesn't matter to me" = "NoInput"),
                            selected = "NoInput")%>%
                 shinyInput_label_embed(
                   shiny_iconlink() %>%
                     bs_attach_modal(paste0("about_site_type"))),
               #bs_embed_tooltip(placement = "right", title =  "here's info about site type",
               #                trigger = "hover")),
               radioButtons(inputId = "lightReqs", label = "Light availability at site", 
                            choices = c("Full shade" = "FullShade",
                                        "Part shade/part sun" = "PartShade",
                                        "Full sun" = "FullSun",
                                        "Doesn't matter to me" = "NoInput"),
                            selected = "NoInput") %>%
                 shinyInput_label_embed(
                   shiny_iconlink() %>%
                     bs_attach_modal(paste0("about_light"))),
               # bs_embed_tooltip(placement = "right", title =  "here's info about light availability",
               #                  trigger = "hover")),
               checkboxInput(inputId = "resistant", label = "Resistant to breakage", value = FALSE)
        ), 
        column(width = 9,
               uiOutput("treeRecommendations"),
               br(),
               "These tree species are recommended by the City of Boston.", 
               a("Read about Boston's urban forest here.", href="https://www.boston.gov/departments/parks-and-recreation/caring-bostons-urban-forest"))
    )
  ) # close tree species row 
  ), # close choose_tree tab
  # fluidRow(
  #   box(width = 12,
  #       h4("Click on a census tract to see neighborhood-specific considerations."),
  #       column(9, h5("Change map layers to visualize relevant spatial patterns across Boston. Priority areas for tree planting take the following into account: Heat Vulnerability Index (which describes vulnerable social groups), the current tree cover, and the mean daytime temperature. Clicking on a census tract will report other factors to consider, such as asthma prevalence, pest sightings, and neighborhood groups that work to support Boston's urban forest.")),
  #       column(3,
  # bs_modal(id = "about_data", title = "Calculations and data sources", body = HTML("<p>How values for each map layer were derived:</p> <p> The Heat Vulnerability Index represents the intersection of multiple census-derived characteristics that can increase the vulnerability of a region to high temperatures. The following variables were included:</p>")),
  # bs_button(paste0("About this data")) %>%
  #   bs_attach_modal(paste0("about_data"))),
  #       br(),
  #       leafletOutput("bostonMap")
  #       )
  tabItem(tabName = "keep_healthy",
          fluidRow(
            box(width = 12,
                h2("Step 4. Keep your tree healthy!"),
                "Here are resources for ensuring that your tree thrives in its new home." )),
  fluidRow(
    box(width = 12,
        column(2,
               tags$a(img(src = "https://pbs.twimg.com/profile_images/631076340280635396/nH93RnFb_400x400.png",
                          width = "80px"),
                      href="https://www.cityofboston.gov/311/",  target="_blank")),
      column(width = 4,
        h3("Request tree maintenance from the City of Boston.", style = "margin-top: 0;
                    margin-block-start: 0px; ")),
        column(6,
               "Boston's 311 is a telephone number, web service, and smartphone application that allows you to request maintenance (as well as new trees!) on public land."))),
  fluidRow(
    box(width = 12,
        column(2,
               tags$a(img(src = "https://upload.wikimedia.org/wikipedia/commons/3/37/USDA_logo.png",
                          width = "80px"),
                      href="https://www.aphis.usda.gov/aphis/resources/pests-diseases/hungry-pests/pest-tracker/states/massachusetts",  target="_blank")),
        column(4,
        h3("Stay aware of pests in your area.", style = "margin-top: 0;
                    margin-block-start: 0px; ")),
        column(6,
               "Check out USDA resources on the current ranges of common MA pests."))),
    fluidRow(
      box(width = 12,
          column(2,
                 tags$a(img(src = "https://ww1.prweb.com/prfiles/2017/10/08/14781858/COB_B_Blue_wName-01.png",
                            height = "90px"),
                        href="https://www.boston.gov/departments/parks-and-recreation/caring-bostons-urban-forest#tree-care-tips",  target="_blank")),
        column(4,
      h3("Read the City of Boston's tips on caring for street trees.", style = "margin-top: 0;
                    margin-block-start: 0px; ")),
        column(6,
               "Check out resources describing how you can help care for the city's urban forest."))),
  fluidRow(
    box(width = 12,
        column(5,
               h5("This application was created by the Boston University Terrestrial Biogeoscience Practicum, Spring 2019.")), 
        column(5,             
               bs_modal(id = "about_data", title = "Calculations and data sources", body = HTML("<p>How values for each map layer were derived:</p> <p> The Heat Vulnerability Index represents the intersection of multiple census-derived characteristics that can increase the vulnerability of a region to high temperatures. The following variables were included:</p>")),
               bs_button(paste0("Data sources and methods")) %>%
                 bs_attach_modal(paste0("about_data")),
               bs_modal(id = "exec_summary", title = "Project Executive Summary", body = HTML("<p>Executive summary text will be pasted here, or linked as PDF.</p>")),
               bs_button(paste0("Project Executive Summary")) %>%
                 bs_attach_modal(paste0("about_data"))),
        column(2,
               img(src='https://upload.wikimedia.org/wikipedia/commons/thumb/3/31/Boston_University_wordmark.svg/1280px-Boston_University_wordmark.svg.png', align = "right", width = "100px", height="50px"))
    ))
  ) # close get_planted tab
  ) # close tabItems
    ) # close dashboardBody



ui <- dashboardPage(#theme = "sandstone",
  header, sidebar, body)


server <- function(input, output, session) {
  
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
    # if ("High" %in% input$allergens) {
    #   high <- app_df[which(app_df$Allergens=="high"),]$Common.Name
    #   approved_allergens <- append(approved_allergens, high)
    # }
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
    dt <- DT::datatable(to_display, 
                        escape = FALSE, 
                        filter="none", 
                        rownames = FALSE, 
                        selection = 'none',
                        colnames=NULL,
                        options = list(
                          #autoWidth = TRUE,
                          columnDefs = list(
                            list(width = '8%', targets = 0), # heat-reduction icon
                            list(className = 'dt-center', targets = 0), # heat-reduction icon
                            list(width = '8%', targets = 1), # scientific.name column 
                            list(width = '10%', targets = 2), # common.name column 
                            list(width = '72%', targets = 3), # description column 
                            list(width = '10%', targets = 4) # source column 
                          ),#list(visible=FALSE, targets=4)), # hide heat-reduction column
                          scrollX=T,
                          sDom  = '<"top">lrt<"bottom">ip', # remove search options
                          pageLength = 3,
                          lengthChange = FALSE,
                          bSort=FALSE)
    ) 
    output$recommendations <- DT::renderDT({ dt
    })
    
    DTOutput("recommendations")
    
  })
  
  output$layerInfo <- renderUI({
    input <- input$map_layer
    description <- switch(input,
                          "HVI" = "This layer shows a Boston-specific heat vulnerability index (HVI), which incorporates information on poverty, race, English proficiency, disability, unemployment, housing density, and building age. The HVI creates a ranking of 12 (least vulnerable to heat) to 22 (most vulnerable to heat) for each census tract in the City of Boston.",
                          "existing_canopy" = "existing canopy description",
                          "potential_canopy" = "potential canopy description")
    description
    
    bs_button(paste0("More about these values")) %>%
      bs_attach_modal(paste0("about_", input))
    #)
    #)
  })
  
  output$temp_trees <- renderLeaflet({
    leafletMap(var = var1, pal = pal1, title = "Mean summer <br>temperature", highlight = "priority_heat_trees")
  })
  
  output$priority_trees <- renderLeaflet({
    leafletMap(var = var2, pal = pal2, title = "Heat Vulnerability <br>Index", highlight = "priority_hvi_trees")
  })
  
  output$temp_other <- renderLeaflet({
    leafletMap(var = var1, pal = pal1, title = "Mean summer <br>temperature", highlight = "priority_heat_other")
  })
  
  output$priority_other <- renderLeaflet({
    leafletMap(var = var2, pal = pal2, title = "Heat Vulnerability <br>Index", highlight = "priority_hvi_other")
  })
  
  output$current_trees <- renderLeaflet({
    m <-  tracts %>%  
      leaflet() %>% addTiles() %>%  
      setView(-71.08, 42.3, zoom = 12) %>%
      addProviderTiles("CartoDB.Positron") %>%
      addPolygons(data = tract_data , 
                  fillColor = ~pal3(var3), 
                  fillOpacity = 0.7, 
                  weight = 0.2, 
                  smoothFactor = 0.2, 
                  popup = ~tract_data$popups,
                  group = "Current tree cover (percentage)") %>%
      addPolygons(data = tract_data , 
                  fillColor = ~pal4(var4), 
                  fillOpacity = 0.7, 
                  weight = 0.2, 
                  smoothFactor = 0.2, 
                  popup = ~tract_data$popups,
                  group = "Current tree cover (sq ft)") %>% 
      addLayersControl(
        baseGroups = c("Current tree cover (percentage)", "Current tree cover (sq ft)"),
        options = layersControlOptions(collapsed = FALSE)) %>% 
      hideGroup("Current tree cover (sq ft)")
    m
  })
  
  observeEvent(input$current_trees_groups, {
    current <- leafletProxy("current_trees") %>% clearControls()
    if (input$current_trees_groups == 'Current tree cover (percentage)'){
      current <- current %>% addLegend(pal = pal3,
                                       values = var3,
                                       position = "bottomright",
                                       title = "Current tree cover <br>(percentage)")
    } else if (input$current_trees_groups == 'Current tree cover (sq ft)'){
      current <- current %>% addLegend(pal = pal4,
                                       values = var4,
                                       labels = labelFormat(6),
                                       position = "bottomright",
                                       title = "Current tree cover <br>(sq ft)")
    }
  })
  
  output$potential_trees <- renderLeaflet({
    m <-  tracts %>%  
      leaflet() %>% addTiles() %>%  
      setView(-71.08, 42.3, zoom = 12) %>%
      addProviderTiles("CartoDB.Positron") %>%
      addPolygons(data = tract_data , 
                  fillColor = ~pal5(var5), 
                  fillOpacity = 0.7, 
                  weight = 0.2, 
                  smoothFactor = 0.2, 
                  popup = ~tract_data$popups,
                  group = "Potential tree cover (percentage)") %>%
      addPolygons(data = tract_data , 
                  fillColor = ~pal6(var6), 
                  fillOpacity = 0.7, 
                  weight = 0.2, 
                  smoothFactor = 0.2, 
                  popup = ~tract_data$popups,
                  group = "Potential tree cover (sq ft)") %>% 
      addLayersControl(
        baseGroups = c("Potential tree cover (percentage)", "Potential tree cover (sq ft)"),
        options = layersControlOptions(collapsed = FALSE)) %>% 
      hideGroup("Potential tree cover (sq ft)")
    m
  })
  
  observeEvent(input$potential_trees_groups, {
    potential <- leafletProxy("potential_trees") %>% clearControls()
    if (input$potential_trees_groups == 'Potential tree cover (percentage)'){
      potential <- potential %>% addLegend(pal = pal5,
                                           values = var5,
                                           position = "bottomright",
                                           title = "Potential tree cover <br>(percentage)")
    } else if (input$potential_trees_groups == 'Potential tree cover (sq ft)'){
      potential <- potential %>% addLegend(pal = pal6,
                                           values = var6,
                                           labels = labelFormat(6),
                                           position = "bottomright",
                                           title = "Potential tree cover <br>(sq ft)")
    }
  })    
  
}

# Run the application 
shinyApp(ui = ui, server = server)

