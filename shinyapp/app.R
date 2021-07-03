
library(shiny)
library(shinydashboard)
library(dplyr)
library(ggplot2)
library(hrbrthemes)
library(coda)


data_in <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/model_outputs/div_summaries.rds")

samples <- data_in$samples
uncert_samples <- samples[[1]]
 
ui <- dashboardPage(
	dashboardHeader(title = "Soil microbiome forecasts"),
	## Sidebar content
	dashboardSidebar(
		sidebarMenu(
			menuItem("Diversity", tabName = "diversity", icon = icon("dashboard")),
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
									radioButtons("uncert","Driver uncertainty sources",
																		 choices = c("all" = "all",
																		 						"spatial" = "spatial", 
																		 						"none" = "none"), 
																		 selected = "all"),
									selectInput("plotID", h3("Select plot"), unique(data_in$plot_est$plotID)
															# c("HARV_001" = "HARV_001", "OSBS_001" = "OSBS_001",
															# 	"CPER_001" = "CPER_001")
															)),
							),
							fluidRow(
								box(plotOutput("parameter_estimates"))
							),
							fluidRow(
							#	box(
									plotOutput("traceplots_all", height = 600),
									#),
								box(actionButton("traceplots", "Generate parameter traceplots"))
							)
			)
		,
		
		# Taxonomy tab content
		tabItem(tabName = "taxa",
						h2("Taxa tab forecasts")
		),
		# Fg tab content
		tabItem(tabName = "function",
						h2("Functional group forecasts")
		)
	)
)
)




server <- function(input, output, session) {
	
	summary_df <- data_in$summary_df
	plot_est <- data_in$plot_est
	
	uncert_samples <- reactiveValues(data = NULL)
	
	observeEvent(input$traceplots, {

		uncert_samples$data <- samples[[1]]
		
		if ("spatial" %in% input$uncert) {
			uncert_samples$data <- samples$full_uncertainty_ITS
		} else if ("all" %in% input$uncert) {
			uncert_samples$data <- samples$spatial_uncertainty_ITS
		} else if ("none" %in% input$uncert){
			uncert_samples$data <- samples$no_uncertainty_ITS
		}
	})
	
	output$traceplots_all <- renderPlot({
		if (is.null(uncert_samples$data)) return()
		plot(uncert_samples$data)
	})
	
	# View fcast w/ confidence intervals for no and full uncertainties
	plot_est$observed <- ifelse(is.na(plot_est$truth), "Estimated", "Observed")
	# Fcast
	
	
	plot_data <- reactive({
		plot_est %>% dplyr::filter(plotID %in% input$plotID)
		
	})
	
	output$plot1 <- renderPlot({
		
		full_uncert <- plot_data() %>% dplyr::filter(scenario == "full_uncertainty_ITS")
		no_uncert <- plot_data() %>% dplyr::filter(scenario == "no_uncertainty_ITS")
		spat_uncert <- plot_data() %>% dplyr::filter(scenario == "spatial_uncertainty_ITS")
		temp_uncert <- plot_data() %>% dplyr::filter(scenario == "temporal_uncertainty_ITS")
		
		if ("spatial" %in% input$uncert){
			output.plot <- ggplot() +
				#facet_wrap(~scenario) + 
				geom_line(data = no_uncert,
									aes(x = dates, y = `50%`), show.legend = F) + 
				geom_ribbon(data = no_uncert, 
										aes(x = dates, ymin = `25%`, ymax = `75%`), fill = "lightblue", alpha=0.6, show.legend = F) +
				geom_point(data = no_uncert,
									 aes(x = dates, y = as.numeric(no_uncert$truth))) + 
				theme_ipsum(base_size = 14, strip_text_size = 22) + 
				ggtitle(paste0("Shannon diversity at ", input$plotID)) + scale_x_date() +	
				scale_fill_brewer(palette = "Paired") + 
				theme(panel.spacing = unit(0, "lines"),
							plot.margin = unit(c(.2, .2, .2, .2), "cm")) + ylim(c(4,6)) +
				geom_ribbon(data = spat_uncert,
										aes(x = dates, ymin = `25%`, ymax = `75%`), fill = "darkgreen", alpha=0.4, show.legend = F) +
				scale_fill_brewer(palette = "Paired")
			
		}

		if ("none" %in% input$uncert) {
			
		
		output.plot <- ggplot() +
			#facet_wrap(~scenario) + 
			geom_line(data = no_uncert,
								aes(x = dates, y = `50%`), show.legend = F) + 
			geom_ribbon(data = no_uncert, 
									aes(x = dates, ymin = `25%`, ymax = `75%`), fill = "lightblue", alpha=0.4, show.legend = F) +
			geom_point(data = no_uncert,
								 aes(x = dates, y = as.numeric(no_uncert$truth))) + 
			theme_ipsum(base_size = 14, strip_text_size = 22) + 
			ggtitle(paste0("Shannon diversity at ", input$plotID)) + scale_x_date() +	
			scale_fill_brewer(palette = "Paired") + 
			theme(panel.spacing = unit(0, "lines"),
						plot.margin = unit(c(.2, .2, .2, .2), "cm")) + ylim(c(4,6))
		}
		
		if ("all" %in% input$uncert) {

		output.plot <-  ggplot() +
			geom_line(data = no_uncert,
								aes(x = dates, y = `50%`), show.legend = F) + 
			geom_ribbon(data = no_uncert, 
									aes(x = dates, ymin = `25%`, ymax = `75%`), fill = "lightblue", alpha=0.4, show.legend = F) +
			geom_point(data = no_uncert,
								 aes(x = dates, y = as.numeric(no_uncert$truth))) + 
			theme_ipsum(base_size = 14, strip_text_size = 22) + 
			ggtitle(paste0("Shannon diversity at ", input$plotID)) + scale_x_date() +	
			scale_fill_brewer(palette = "Paired") + 
			theme(panel.spacing = unit(0, "lines"),
						plot.margin = unit(c(.2, .2, .2, .2), "cm")) + ylim(c(4,6)) +
			geom_ribbon(data = full_uncert,
				aes(x = dates, ymin = `25%`, ymax = `75%`), fill = "darkblue", alpha=0.4, show.legend = F) +
			scale_fill_brewer(palette = "Paired")
			
		}
		
		output.plot
	})
	
	output$parameter_estimates <- renderPlot({
		
		full_uncert <- summary_df %>% dplyr::filter(scenario == "full_uncertainty_ITS")
		no_uncert <- summary_df %>% dplyr::filter(scenario == "no_uncertainty_ITS")
		spat_uncert <- summary_df %>% dplyr::filter(scenario == "spatial_uncertainty_ITS")
		temp_uncert <- summary_df %>% dplyr::filter(scenario == "temporal_uncertainty_ITS")
		
		if ("spatial" %in% input$uncert) {
			beta_out <- spat_uncert[which(!is.na(spat_uncert$beta)),]
		} else if ("all" %in% input$uncert) {
			beta_out <- full_uncert[which(!is.na(full_uncert$beta)),]
		} else if ("none" %in% input$uncert){
			beta_out <- no_uncert[which(!is.na(no_uncert$beta)),]
		}
			# By rank with every taxon - cluttered
			ggplot(data=beta_out,
						 aes(x = reorder(beta, effSize),y = effSize)) +
				facet_grid(~scenario) +
				geom_point(aes(shape = as.factor(significant), color = beta), size = 3) +
				labs(col = "Parameter", title = "Absolute effect size") + 
				xlab("Parameter")+ 
				ylab(NULL) +
				theme_bw() + theme(#axis.ticks.x=element_blank(),
					text = element_text(size = 16),
					axis.text.x=element_text(#angle = 45, hjust = 1, vjust = 1),
						angle = 320, vjust=1, hjust = -0.05),
					axis.title=element_text(size=22,face="bold")
					#strip.text.y = element_text(size=24,hjust=0,vjust = 1,angle=180,face="bold")
				) + scale_shape_manual(values = c(21, 16), name = NULL, 
															 labels = c("Not significant","Significant")) 
		})
}

shinyApp(ui, server)
