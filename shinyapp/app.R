source("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/source.R")

library(shiny)
library(shinydashboard)
library(hrbrthemes)


effect_sizes <- readRDS("./data/summary/all_fcast_effects.rds")

#hindcasts <-  readRDS("./data/summary/all_hindcasts.rds")
hindcasts_orig <-  readRDS("./data/summary/hindcast_div.rds")
hindcasts <- hindcasts_orig %>% 
	mutate(calibration_label = recode(as.character(time_period), !!!calibration_label),
				 truth = as.numeric(truth))
div_hindcasts <- hindcasts %>% filter(fcast_type=="Diversity")


# samples <- data_in$samples
# uncert_samples <- samples[[1]]

ui <- dashboardPage(
	dashboardHeader(title = "Soil microbiome forecasts"),
	## Sidebar content
	dashboardSidebar(
		sidebarMenu(
			menuItem("Diversity", tabName = "diversity", icon = icon("tachometer-alt")),
			menuItem("Taxa", tabName = "taxa", icon = icon("th")),
			menuItem("Functional groups", tabName = "function", icon = icon("th"))
		)
	),
	
	## Body content
	dashboardBody(
		tabItems(
			# First tab content
			tabItem(tabName = "diversity",
							h2("Diversity forecasts"),
							
							fluidRow(
								box(plotOutput("plot1", height = 450)),
								
								box(
									title = "Controls",
									radioButtons("group",
															 "Select microbial group",
															 choices = c("Fungi" = "ITS",
															 						"Bacteria" = "16S"), 
															 selected = "Fungi"),
									radioButtons("model_name",
															 "Select model type",
															 choices = c("environmental" = "all_covariates",
															 						"phenology" = "cycl_only"), 
															 selected = "environmental"),
									
									selectInput("siteID", h3("Select site"), unique(hindcasts$siteID)),
															
									selectInput("plotID", h3("Select plot"), unique(hindcasts$plotID))
															# c("HARV_001" = "HARV_001", "OSBS_001" = "OSBS_001",
															# 	"CPER_001" = "CPER_001")
									) # Close box
								), # Close fluidRow
							fluidRow(
								verbatimTextOutput("summary")
								
							)
								
								
	) # Close tab
) # Close tabItems
) # Close dashboardBody
) # Close dashboardPage




server <- function(input, output, session) {
	
	# run every time data is updated
	observe({
		site_data <- div_hindcasts %>% 
			filter(siteID %in% input$siteID) %>% 
			filter(#model_name==input$model_name & 
						 	model_name=="all_covariates" &
						 	#group == input$group & 
						 	group == "ITS" & 
						 	calibration_label == "Excluding legacy data")
		site_plots <- unique(site_data$plotID)
	
		
		updateSelectInput(session, "plotID",
								 choices = c(site_plots, "all"), # update choices
								 selected = "all") # remove selection
	})


	# observeEvent(input$group, {
	# 	
	# 	uncert_samples$data <- samples[[1]]
	# 	
	# 
	# # View fcast w/ confidence intervals for no and full uncertainties
	# plot_est$observed <- ifelse(is.na(plot_est$truth), "Estimated", "Observed")
	# # Fcast
	# }
	
	plot_data <- reactive({

		plotting_data <- div_hindcasts %>% 
			filter(siteID %in% input$siteID) %>% 
			filter(model_name=="all_covariates" & 
						 	group == "ITS" & 
						 	calibration_label == "Excluding legacy data")
		
		if (input$plotID != "all") {
			plotting_data <- plotting_data %>% filter(plotID %in% input$plotID)
		}
		plotting_data
	})
	
	output$summary <- renderPrint({
		print(head(plot_data()))
		print(tail(plot_data()))
					})
	
	output$plot1 <- renderPlot({
		

		# input <- list()
		# input$plotID <- "HARV_001"
			
		# plotting_data <- div_hindcasts %>% 
		# 		dplyr::filter(plotID %in% input$plotID) %>% 
		# 		filter(model_name=="all_covariates" & group == "ITS" & 
		# 					 	calibration_label == "Excluding legacy data")
			
			output.plot <-  
				# View actual models & hindcasts for a single site
				ggplot(plot_data()) + 
				facet_grid(rows=vars(plotID), 
									 #cols = vars(calibration_label), 
									 drop=T, scales="free") +
				geom_line(aes(x = dates, y = med), show.legend = F, linetype=2) +
				geom_line(aes(x = dates, y = `50%`), show.legend = F) +
				geom_ribbon(aes(x = dates, ymin = lo, ymax = hi), alpha=0.6, fill="blue") +
				geom_ribbon(aes(x = dates, ymin = `2.5%`, ymax = `97.5%`),fill="red", alpha=0.6) +
				theme_bw()+
				scale_fill_brewer(palette = "Paired") +
				theme(text = element_text(size = 14), panel.spacing = unit(.2, "cm"),
							legend.position = "bottom",legend.title = element_text(NULL),
							plot.margin = unit(c(.2, .2, 2, .2), "cm")) + ylab(NULL) +
				geom_point(aes(x = dates, y = as.numeric(truth))) + xlab(NULL) + labs(fill='') + 
				xlim(c(as.Date("2013-06-01"), as.Date("2020-01-01"))) + 
				ggtitle(paste0("Diversity hindcasts at NEON plot: ", input$plotID))
		
		output.plot
	
	})
	
	# output$parameter_estimates <- renderPlot({
	# 	
	# 	# By rank with every taxon - cluttered
	# 	ggplot(data=beta_out,
	# 				 aes(x = reorder(beta, effSize),y = effSize)) +
	# 		facet_grid(~scenario) +
	# 		geom_point(aes(shape = as.factor(significant), color = beta), size = 3) +
	# 		labs(col = "Parameter", title = "Absolute effect size") + 
	# 		xlab("Parameter")+ 
	# 		ylab(NULL) +
	# 		theme_bw() + theme(#axis.ticks.x=element_blank(),
	# 			text = element_text(size = 16),
	# 			axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
	# 				angle = 320, vjust=1, hjust = -0.05),
	# 			axis.title=element_text(size=22,face="bold")
	# 			#strip.text.y = element_text(size=24,hjust=0,vjust = 1,angle=180,face="bold")
	# 		) + scale_shape_manual(values = c(21, 16), name = NULL, 
	# 													 labels = c("Not significant","Significant")) 
	# })
}

shinyApp(ui, server)
