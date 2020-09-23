# Shiny web application for facilitating successful tree-planting in Boston, MA.
# Written by Zoey Werbin, 2019
# May or may not be maintained, email zoeywerbin@gmail.com with any issues (or use Github issues)

library(shiny)
library(shinyBS)
library(DT)
suppressMessages(library(dplyr))
library(tidyverse)
library(bsplus)
library(leaflet) # for map base
library(htmlwidgets) # css fix for map
library(shinyWidgets)
library(shinythemes)
library(shinydashboard)

# call source file, which loads in data and sets some variable names.
# text for modals (pop-ups) is all in the source file as well.
source("source.R")

# title
header <-  dashboardHeader(title = "Right Place, Right Tree | Boston, MA", titleWidth = 530)

# sidebar
sidebar <- dashboardSidebar(#id = "inTabset",
                            
  sidebarMenu(id = "inTabset",
    menuItem("Did you know that Boston's street trees can change the local climate? This tool can help you decide where to plant trees, what species will thrive at your site, and how to take care of new trees."),
    menuItem("Step 1: Choose region", icon = icon("map-marked-alt"), tabName = "priority"),
    menuItem("Step 2: Learn about regional considerations", icon = icon("question-circle"), tabName = "other_considerations"),
    menuItem("Step 3: Choose the right tree", icon = icon("tree"), tabName = "choose_tree"),
    menuItem("Step 4: Keep your tree healthy!", icon = icon("check-circle"), tabName = "keep_healthy"),  
    menuItem(bs_button("Data sources and methods")) %>%
               bs_attach_modal("about_data"),
    helpText("This application was created by the Boston University Biogeoscience and URBAN graduate programs.")
  )) 

# body (split into tabs)
body <- dashboardBody(
  
  bs_modal(id = "about_data", 
           title = "Data sources and methods", 
           body = HTML(modal_data_text)),
  
  tags$head(
    # include stylesheets for the "fontAwesome" icons, as well as the Josefin Sans font.
    HTML('<link rel="stylesheet" href="https://use.fontawesome.com/releases/v5.8.1/css/all.css" integrity="sha384-50oBUHEmvpQ+1lW4y57PTFmhCaXp0ML5d60M1M7uH2+nqUivzIebhndOJK28anvf" crossorigin="anonymous">
         <link rel="stylesheet" href="//fonts.googleapis.com/css?family=Josefin+Sans" crossorigin="anonymous">')),
  
  tabItems(
    tabItem(tabName = "priority",
            fluidRow(
              box(width = 12,
                  h2("Step 1. Choose a region for planting trees."), "Skip this step if you already know where you're planting!")),
            fluidRow(
              box(width=12, 
                  column(6, "Street trees can cool the city down in the summer, reducing the 'urban heat islands.' Heat can be especially dangerous for certain groups, like those who can't leave their homes or don't have air conditioning.",
                         tags$div(HTML('<i class="fas fa-exclamation"></i> marks "priority" regions that have lots of room for new trees.'))),
                  column(6,
                         radioButtons("choose_temp_hvi", "View map of Boston by:", 
                                       c("Summer morning temperatures - how hot does this region get, on average?" = "temp",
                                          "Heat Vulnerability Index - how dangerous are heat waves for this community?" = "HVI"))),
                    leafletOutput("temp_hvi", height = 500)) # output of map with temperature or HVI highlighted
              ),
            fluidRow(
              bs_modal(id = "about_alternatives", 
                       title = "Alternatives to tree-planting", 
                       body = HTML(modal_alternatives_text)),
              box(width=10, 
            bsCollapse(id = "collapseExample", 
                       bsCollapsePanel("Not sure if there is enough space for trees?", 
                                           "Explore the map below to see the current and potential tree canopy cover for each Boston census tract. In regions with low potential for tree canopy expansion, it may make sense to explore other greening options.",
                                       bs_button("Learn about alternatives to tree-planting") %>%
                                         bs_attach_modal("about_alternatives"),
                                           leafletOutput("current_trees", height = 500))
                                       )),
            actionButton('goToTab2', 'Step 2', icon("chevron-circle-right"), 
                         style="color: #fff; background-color: #009933;")
            
            )
            ), # close tabItem "priority"
    
    
    tabItem(tabName = "other_considerations",
            fluidRow(
              box(width = 12,
                  h2("Step 2. Learn about regional considerations"),
                  "When planting a tree, you should carefully consider its context, such as community characteristics, land ownership, and local pests or diseases. Is asthma prevalent in the region? Consider planting a low-allergen tree. Are there community organizations you could partner with to ensure maintenance? Contact them!", 
                  br(), br(),
                  "You can request that the city plant trees on your own land, or on public land. Use the map layers below to find potential tree sites." ,
              leafletOutput("other_considerations", height = 500),
              plotOutput("other_considerations_overlay", height = 10) # hacky way of getting the map to update
              )
            ), # close map row
            column(10),
            column(2,
            actionButton('goToTab3', 'Step 3', icon("chevron-circle-right"), 
                         style="color: #fff; background-color: #009933;"))
    ), # close tabItem "other_considerations"
    
    
    tabItem(tabName = "choose_tree",
            fluidRow(
              box(width = 12,
                  h2("Step 3. Choose the right tree species"),
                  "Enter your site characteristics to view a list of trees that may be perfect for your site. If planting multiple trees, plant a mixture of species - this makes our urban forest more resilient to future pests and diseases." )),
            fluidRow( 
              # modals below are for informational pop-ups; their text is in source.R
              bs_modal(id = "about_sizing", 
                       title = "Determining size requirements for trees", 
                       body = HTML(modal_sizing_text)),
              bs_modal(id = "about_allergies", 
                       title = "Planting with allergies in mind", 
                       body = HTML(modal_allergen_text)),
              bs_modal(id = "about_site_type", 
                       title = "What are the requirements of different site types?", 
                       body = modal_site_text),
              bs_modal(id = "about_light", 
                       title = "What counts as 'full sun,' or 'part sun', or 'shade'?", 
                       body = HTML(modal_light_text)),
              
              box(width=12,
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
                         radioButtons(inputId = "allergens", label = "Allergen preferences", 
                                      choices = c("Low allergen trees" = "Low",
                                                  "Moderate or low allergen trees" = "Moderate",
                                                  "Doesn't matter to me" = "NoInput"), 
                                      selected = "NoInput") %>%
                           shinyInput_label_embed(
                             shiny_iconlink() %>%
                               bs_attach_modal(paste0("about_allergies"))),
                         radioButtons(inputId = "siteType", label = "Site type", 
                                      choices = c("Park or yard" = "Park_yard",
                                                  "Street" = "Street",
                                                  "Doesn't matter to me" = "NoInput"),
                                      selected = "NoInput")%>%
                           shinyInput_label_embed(
                             shiny_iconlink() %>%
                               bs_attach_modal(paste0("about_site_type"))),
                         radioButtons(inputId = "lightReqs", label = "Light availability at site", 
                                      choices = c("Full shade" = "FullShade",
                                                  "Part shade/part sun" = "PartShade",
                                                  "Full sun" = "FullSun",
                                                  "Doesn't matter to me" = "NoInput"),
                                      selected = "NoInput") %>%
                           shinyInput_label_embed(
                             shiny_iconlink() %>%
                               bs_attach_modal(paste0("about_light"))),
                         checkboxInput(inputId = "resistant", label = "Resistant to breakage", value = FALSE)
                  ), 
                  column(9, uiOutput("treeRecommendations"),
                         br(),
                         "These tree species are recommended by the City of Boston.", 
                         a("Read about Boston's urban forest here.", 
                           href="https://www.boston.gov/departments/parks-and-recreation/caring-bostons-urban-forest"))
              )
            ), # close tree species row 
            fluidRow(
            column(2, offset = 10,
                   actionButton('goToTab4', 'Step 4', icon("chevron-circle-right"), 
                                style="color: #fff; background-color: #009933;")))
    ), # close choose_tree tab
   
    
    tabItem(tabName = "keep_healthy",
            fluidRow(
              box(width = 12,
                  h2("Step 4. Keep your tree healthy!"),
                  "Here are resources for ensuring that your tree thrives in its new home." )),
            fluidRow(
              box(width = 12,
                  column(2, tags$a(img(src = "https://pbs.twimg.com/profile_images/631076340280635396/nH93RnFb_400x400.png", width = "80px"),
                                href="https://www.cityofboston.gov/311/",  target="_blank"), style="height: 100%;"),
                  column(width = 4,
                         h3("Request tree maintenance from the City of Boston.")),
                  column(6, "Boston's 311 is a telephone number, smartphone application, and ",
                         tags$a("web service", 
                                href="https://www.cityofboston.gov/311/",  
                                target="_blank"), 
                         " that allows you to request maintenance (as well as new trees!) on public land.", style="margin-top: 10px;"
                         ))),
            fluidRow(
              box(width = 12,
                  column(2, tags$a(img(src = "https://upload.wikimedia.org/wikipedia/commons/3/37/USDA_logo.png",
                                    width = "80px"),
                                href="https://www.aphis.usda.gov/aphis/resources/pests-diseases/hungry-pests/pest-tracker/states/massachusetts",  
                                target="_blank")),
                  column(4, h3("Stay aware of pests in your area.")),
                  column(6, tags$a("Check out USDA resources on the current ranges of common MA pests.", 
                                   href="https://www.aphis.usda.gov/aphis/resources/pests-diseases/hungry-pests/pest-tracker/states/massachusetts",
                                   target="_blank"),
                         br(), br(),
                         tags$a("Explore the National Insect & Disease Risk Map", 
                                href="https://usfs.maps.arcgis.com/apps/webappviewer/index.html?id=52cb2bcc3c2b4868ac87b66f622062ab",  
                                target="_blank")))
              ), # close row 
            fluidRow(
              box(width = 12,
                  column(2, tags$a(img(src = "https://ww1.prweb.com/prfiles/2017/10/08/14781858/COB_B_Blue_wName-01.png",
                                    height = "90px"),
                                href="https://www.boston.gov/departments/parks-and-recreation/caring-bostons-urban-forest#tree-care-tips",  
                                target="_blank")),
                  column(4, h3("Read the City of Boston's tips on caring for street trees.")),
                  column(6, tags$a("Check out resources describing how you can help care for the city's urban forest.", 
                                   href="https://www.boston.gov/departments/parks-and-recreation/caring-bostons-urban-forest#tree-care-tips",  
                                   target="_blank"),
                         br(), br(),
                         tags$a("Caring for trees on private property? Contact a certified arborist for assistance!", href = "https://massarbor.org/directory.php",  
                                target="_blank")))),
            fluidRow(
                  box(width = 12,
                  column(2, tags$a(img(src = "http://www.digsafe.com/img/logo.png",width = "80px"),
                                   href="http://www.digsafe.com/",  target="_blank")),
                  column(4, h3("Planting on private property? Notify DigSafe to ensure utility lines are not disrupted.")),
                  column(6, "MA state law requires that DigSafe is notified of at least 72 hours in advance of any projects that involve digging, including planting shrubs or trees. Call 811 or go to ",
                         tags$a("DigSafe.com", 
                                href="http://www.digsafe.com/",  
                                target="_blank"), 
                         " to submit a notification."
                  )))
              ) # close get_planted tab
              ) # close tabItems
            ) # close dashboardBody

ui <- dashboardPage(header, sidebar, body, skin = "green",
    tags$head(includeCSS("www/custom.css"))
  )
