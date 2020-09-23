#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)

# Define UI for application that draws a histogram
ui <- dashboardPage(
  skin = "green",
  dashboardHeader(title = "Buy some green from Dan and Zoey",
                  titleWidth = 420,
                  dropdownMenuOutput("messageMenu")
                  ), # close dashboard Header 
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("house")),
      menuItem("Design a terrarium", icon = icon("stone"), tabName = "terrarium"),
      menuItem("Plant shop", icon = icon("plant"), tabName = "plants"),
      menuItem("Art shop", icon = icon("paintbrush"), tabName = "art")
    )
  ), # close dashboard Sidebar
  dashboardBody(
    tabItems(
    tabItem(tabName = "home",
            box("welcome to our plant-stuff shop")
            ), # close "home" tab
    
    tabItem(tabName = "terrarium",
            
                fluidRow(
                box(title = "Pick a container size and type:",
                   radioButtons("container", "Container type/size",
                                  c("small glass orb ($4)" = "small_orb",
                                    "large glass orb ($6)" = "large_orb",
                                    "small glass teardrop ($6)" = "small+teardrop",
                                    "aquarium ($25)" = "aquarium",
                                    "centerpiece ($15)" = "centerpiece"))
                ),
                box(title = "Pick large plant(s):",
                    checkboxGroupInput("plants", "Choose between 0-3 large plants",
                                  c("arrowhead plant (syngonium sp.) ($4)" = "syngonium",
                                    "fittonia ($4)"= "fittonia",
                                    "lemon button fern ($9)" = "lemon_button",
                                    "prayer plant ($8)" = "prayer_plant",
                                    "coleus ($8)" = "coleus",
                                    "heart-leaf philodendron (p. scandens) ($12)" = "micans"))),
                box(title = "Pick small plants(s):",
                      checkboxGroupInput("small_plants", label = "Choose between 0-3 small plants", 
                                          choices = c(
                                            "tall air plant ($4)" = "tall_airplant",
                                           "round air plant ($4)" = "round_airplant",
                                           "mini oakleaf creeping fig ($3)" = "creeping_fig")))
                ), #close fluidRow #1
            fluidRow(
                box(title = "Choose substrate:",
                    checkboxGroupInput("substrate", label = "Choose at least 1 substrate (but up to 4 layers in aquarium/centerpiece):", 
                                       choices = c(
                                         "pebbles ($2)" = "pebbles",
                                         "decorative ('raindeer') moss ($2)" = "decorative_moss",
                                         "live green moss ($3)" = "live_moss",
                                         "white sand ($2)" = "white_sand",
                                         "soil (recommended, but not necessary for airplants) ($1)" = "soil")
                    )
                    ),
                box(title = "Add extras:",
                      checkboxGroupInput("extras", "Choose between 0-3 additional things:",
                                         c("unicorn ($1)" = "unicorn",
                                           "dinosaur ($1)" = "dinosaur",
                                           "crystal ($3)" = "crystal"
                                           )))
                ), #close fluidRow #2
            fluidRow(
              box(title = "How much light do you have to work with?",
                  sliderInput("light", "Available light (0 = far from a window; 10 = south-facing window or outside:", 1, 10, 5)
              ),
              box(
                textOutput(outputId = "lightOpinion")
              ),
              box("Total cost:",
                  br(),
                  "1 billion dollars",
                  br(),
                  "also",
                  checkboxInput("pre_assembled", "we assemble it for you (+$6)",
                                value = FALSE))
            )
    ), # close "terrarium" tab
    tabItem("plants",
            box("we have cute plants for $1 trillion")
    ),
    tabItem("art",
            box("look at these plant-inspired prints and originals from various artists. each painting is $175k"))    
    ) # close "tabItems" 
) # close dashboardBody

) # close dashboardPage

# 
# tags$head(tags$style(HTML('
#       .main-header .logo {
#                           font-family: "Georgia", Times, "Times New Roman", serif;
#                           font-weight: bold;
#                           font-size: 19px;
#                           }
#                           ')))



# Define server logic required to draw a histogram
server <- function(input, output) {
  
  enoughLight <- reactive({
    lightOpinion <- "ok"
    if (input$light < 3) {
      lightOpinion <- "uh oh, low-light may only be tolerable for the ferns and philodendron"
    } else if (input$light > 8) {
      lightOpinion <-  "uh oh, bright light may not be tolerable for many of these plants"
    }
      lightOpinion
  })
  
  output$lightOpinion <- renderPrint({
    enoughLight()
  })
}



# Run the application 
shinyApp(ui = ui, server = server)

