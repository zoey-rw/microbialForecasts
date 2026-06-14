# Per-taxon microbial forecast explorer (MVP)
#
# Exploratory Shiny app for browsing single-taxon results from the cloglog beta
# driver-uncertainty pipeline. Two pages (left sidebar):
#   - Available models: browsable catalog of every converged taxon/group
#   - Explore taxon:     forecast time series, predictor effects, scoring metrics
#
# Only converged models are shown (weak convergence, Rhat < 1.2, is sufficient);
# unconverged taxon/model combinations are omitted entirely.
#
# The per-model hindcast parquet files are queried per-taxon on demand via duckdb
# (unioned) so the app does not have to hold all 19M rows in memory. Small summary
# tables are loaded once at startup.
#
# Run from the project root or this directory:
#   shiny::runApp("shinyapp", launch.browser = TRUE)

# --- Environment ---------------------------------------------------------------
# source.R defines the canonical color palettes (kingdom_colors, model_colors)
# and prepends R_library/ to .libPaths(). Locate it whether the app is launched
# from the project root or from shinyapp/.
if (file.exists("source.R")) {
	source("source.R")
} else if (file.exists("../source.R")) {
	source("../source.R")
} else {
	stop("Could not find source.R (run from project root or shinyapp/).")
}

suppressPackageStartupMessages({
	library(shiny)
	library(shinydashboard)
	library(DT)
	library(duckdb)
	library(DBI)
})

# --- Data paths ----------------------------------------------------------------
parquet_dir         <- here::here("data/summary/parquet")
parquet_files       <- file.path(parquet_dir,
                                  sprintf("hindcasts_%s.parquet", c("env_cycl", "cycl_only", "env_cov")))
parquet_files       <- parquet_files[file.exists(parquet_files)]
predictor_eff_path  <- here::here("data/summary/predictor_effects.rds")
scoring_path        <- here::here("data/summary/scoring_metrics_plsr2.rds")

if (length(parquet_files) == 0) {
	stop("Hindcast parquet files not found in ", parquet_dir,
			 "\nRun analysis/model_analysis/07_tidyHindcasts.r or download_data.R.")
}

# duckdb read expression unioning the available per-model hindcast files
hindcast_read <- sprintf("read_parquet([%s], union_by_name=true)",
												 paste(sprintf("'%s'", parquet_files), collapse = ", "))

# Model types (label = value). env_cycl is the preferred default when available.
model_choices <- c("Seasonality + environment" = "env_cycl",
									 "Seasonality only"          = "cycl_only",
									 "Environment only"          = "env_cov")
model_label <- function(m) names(model_choices)[match(m, model_choices)]

# Site-prediction modes shown to users. The new-site random-effect mode exists in
# the data but is intentionally omitted (modeled effect is the meaningful one).
site_pred_choices <- c("Observed site"      = "New time (observed site)",
											 "New site (modeled)" = "New time x site (modeled effect)")

forecast_start <- as.Date("2018-01-01")

# --- Load small summary tables (once) -----------------------------------------
predictor_effects <- tryCatch(
	as.data.frame(readRDS(predictor_eff_path)),
	error = function(e) NULL
)

# Convergence set: weak convergence (Rhat < 1.2) is treated as sufficient.
converged_set <- tryCatch({
	sm <- readRDS(scoring_path)
	list(metrics   = as.data.frame(sm$scoring_metrics),
			 converged = unname(as.character(sm$converged_list)))
}, error = function(e) list(metrics = NULL, converged = character(0)))

# Is (taxon, model) a converged model?
is_converged <- function(taxon, model) {
	set <- converged_set$converged
	length(set) > 0 &&
		any(grepl(model, set, fixed = TRUE) &
					grepl(paste0("_", taxon, "_"), set, fixed = TRUE))
}

# --- NEON site lookup ----------------------------------------------------------
# NEON terrestrial field-site names (codes are otherwise opaque). Stable public
# identifiers; match() tolerates any code not in this list (falls back to code).
neon_site_names <- c(
	ABBY = "Abby Road, WA",            BARR = "Utqiagvik (Barrow), AK",
	BART = "Bartlett Experimental Forest, NH", BLAN = "Blandy Experimental Farm, VA",
	BONA = "Caribou-Poker Creeks, AK", CLBJ = "LBJ National Grassland, TX",
	CPER = "Central Plains Experimental Range, CO", DEJU = "Delta Junction, AK",
	DELA = "Dead Lake, AL",            DSNY = "Disney Wilderness Preserve, FL",
	GRSM = "Great Smoky Mountains, TN",GUAN = "Guanica Forest, PR",
	HARV = "Harvard Forest, MA",       HEAL = "Healy, AK",
	JERC = "Jones Ecological Research Center, GA", JORN = "Jornada, NM",
	KONA = "Konza Prairie (Agricultural), KS", KONZ = "Konza Prairie Biological Station, KS",
	LAJA = "Lajas Experimental Station, PR", LENO = "Lenoir Landing, AL",
	MLBS = "Mountain Lake Biological Station, VA", MOAB = "Moab, UT",
	NIWO = "Niwot Ridge, CO",          NOGP = "Northern Great Plains, ND",
	OAES = "Klemme Range Research Station, OK", ONAQ = "Onaqui, UT",
	ORNL = "Oak Ridge, TN",            OSBS = "Ordway-Swisher Biological Station, FL",
	RMNP = "Rocky Mountain National Park, CO", SCBI = "Smithsonian Conservation Biology Institute, VA",
	SERC = "Smithsonian Environmental Research Center, MD", SJER = "San Joaquin Experimental Range, CA",
	SOAP = "Soaproot Saddle, CA",      SRER = "Santa Rita Experimental Range, AZ",
	STEI = "Steigerwaldt Land Services, WI", STER = "Sterling, CO",
	TALL = "Talladega National Forest, AL", TOOL = "Toolik Lake, AK",
	TREE = "Treehaven, WI",            UKFS = "University of Kansas Field Station, KS",
	UNDE = "UNDERC, MI",               WOOD = "Woodworth, ND",
	WREF = "Wind River Experimental Forest, WA", YELL = "Yellowstone Northern Range, WY")

# Site metadata table (latitude/longitude, NEON domain, land cover) used to give
# context to the opaque site codes. MAT/MAP in this file are standardized, so
# they are intentionally not shown as physical units.
site_meta <- tryCatch({
	sp <- as.data.frame(readRDS(here::here("data/clean/site_effect_predictors.rds")))
	sp <- sp[, c("siteID", "latitude", "longitude", "nlcd", "ecoregion")]
	sp$siteName <- unname(neon_site_names[sp$siteID])
	sp
}, error = function(e) NULL)

# Label a site code with its full name, e.g. "HARV - Harvard Forest, MA".
site_label <- function(ids) {
	if (is.null(site_meta)) return(ids)
	nm <- site_meta$siteName[match(ids, site_meta$siteID)]
	ifelse(is.na(nm), ids, paste0(ids, " - ", nm))
}

# --- Taxon catalog from the parquet (distinct species, once) -------------------
taxon_catalog <- local({
	con <- dbConnect(duckdb::duckdb())
	on.exit(dbDisconnect(con, shutdown = TRUE))
	sql <- paste0(
		"SELECT DISTINCT species, rank_name, pretty_group ",
		"FROM ", hindcast_read, " ",
		"WHERE species IS NOT NULL ORDER BY pretty_group, rank_name, species")
	tc <- dbGetQuery(con, sql)
	tc$pretty_group[is.na(tc$pretty_group)] <- "Unknown"
	tc
})

# --- Convergence availability + catalog (once) ---------------------------------
# Per-taxon: which model types converged. Taxa with no converged model are
# dropped entirely (we only show converged models).
avail <- taxon_catalog
for (m in model_choices) {
	avail[[m]] <- vapply(avail$species, is_converged, logical(1), model = m)
}
avail$any <- avail$env_cycl | avail$cycl_only | avail$env_cov
avail <- avail[avail$any, , drop = FALSE]

# Converged model types for a taxon, in preference order (env_cycl first).
converged_models <- function(taxon) {
	row <- avail[avail$species == taxon, , drop = FALSE]
	if (!nrow(row)) return(character(0))
	model_choices[as.logical(row[1, unname(model_choices)])]
}

# Grouped named list for the taxon selector: Kingdom > "taxon (rank)".
taxon_choices <- local({
	groups <- split(avail, avail$pretty_group)
	lapply(groups, function(g) {
		labels <- paste0(g$species, "  (", sub("_(bac|fun)$", "", g$rank_name), ")")
		stats::setNames(as.list(g$species), labels)
	})
})

default_taxon <- if ("acidobacteriota" %in% avail$species) "acidobacteriota" else avail$species[1]

# Display catalog for the "Available models" page: one row per converged taxon,
# with a check mark for each model type that converged.
model_catalog <- local({
	out <- data.frame(
		Taxon   = avail$species,
		Kingdom = avail$pretty_group,
		Type    = ifelse(grepl("_(bac|fun)$", avail$rank_name), "Taxonomic", "Functional"),
		Rank    = sub("_(bac|fun)$", "", avail$rank_name),
		check.names = FALSE, stringsAsFactors = FALSE)
	for (m in model_choices) {
		out[[names(model_choices)[model_choices == m]]] <-
			ifelse(avail[[m]], "✓", "")
	}
	out
})

# --- Per-taxon hindcast query --------------------------------------------------
# Returns the time series for one taxon + model, all site_prediction types.
# Built with paste0 + dbQuoteString so the literal "%" in column names like
# "2.5%" is never fed to sprintf.
#
# Rows before a plot's first sampling date are filler (median pinned at exactly
# 0) and would otherwise draw a flat line from 2013 up to the plot start. This
# leading block is trimmed by finding each plot's first non-filler date (the
# earliest date with med > 0) and keeping rows on/after it. Trimming by date
# (not by timepoint) keeps the forecast period, whose rows have a NULL timepoint.
query_taxon <- function(con, taxon, model_name) {
	# The complete quantile ladder is lo(2.5%) < lo_25(25%) < med(50%) <
	# hi_75(75%) < hi(97.5%); these are populated for every period and site mode.
	# The percent-named columns ("2.5%", "25%", ...) are only populated for
	# observed-site calibration, so they are not used here.
	inner <- paste0(
		"SELECT siteID, plotID, dates, truth, med, mean, ",
		"lo, lo_25, hi_75, hi, ",
		"fcast_period, new_site, site_prediction, ",
		"MIN(CASE WHEN med > 0 THEN dates END) ",
		"OVER (PARTITION BY plotID, site_prediction) AS start_date ",
		"FROM ", hindcast_read, " ",
		"WHERE species = ", dbQuoteString(con, taxon),
		" AND model_name = ", dbQuoteString(con, model_name))
	sql <- paste0("SELECT * FROM (", inner, ") t ",
								"WHERE start_date IS NULL OR dates >= start_date")
	out <- dbGetQuery(con, sql)
	if (nrow(out)) {
		out$dates <- as.Date(out$dates)
		out$truth <- suppressWarnings(as.numeric(out$truth))
	}
	out
}

# --- UI ------------------------------------------------------------------------
ui <- dashboardPage(
	dashboardHeader(title = "Microbial forecasts", titleWidth = 320),

	dashboardSidebar(
		width = 320,
		sidebarMenu(
			id = "menu",
			menuItem("Available models", tabName = "catalog", icon = icon("table")),
			menuItem("Explore taxon",   tabName = "explore", icon = icon("chart-line"))
		),
		# Explorer controls appear only on the Explore page.
		conditionalPanel(
			condition = "input.menu == 'explore'",
			selectInput("taxon", "Taxon", choices = taxon_choices, selected = default_taxon),
			radioButtons("model_name", "Model type",
									 choices = converged_models(default_taxon),
									 selected = converged_models(default_taxon)[1]),
			radioButtons("site_pred", "Site prediction",
									 choices = site_pred_choices, selected = site_pred_choices[1]),
			selectInput("siteID", "Site", choices = NULL),
			selectInput("plotID", "Plot", choices = NULL),
			tags$div(style = "padding:8px 15px;",
							 downloadButton("dl_forecast", "Download forecasts (CSV)"),
							 helpText(style = "margin-top:6px;",
							 				 "All sites and plots for the selected taxon and model."))
		)
	),

	dashboardBody(
		tabItems(
			tabItem(
				tabName = "catalog",
				h3("Available models"),
				helpText("All taxa and functional groups with at least one converged model.",
								 "A check mark indicates that model type converged for the taxon.",
								 "Use the search and column filters to browse; click a row to explore it."),
				DT::dataTableOutput("catalog_table")
			),
			tabItem(
				tabName = "explore",
				fluidRow(
					tabBox(
						width = 12, id = "tabs",
						tabPanel("Forecast",
										 uiOutput("site_info"),
										 plotOutput("forecast_plot", height = 480),
										 helpText("Dashed line = posterior mean; solid = posterior median.",
										 				 "Orange band = 95% interval, blue band = 50% interval.",
										 				 "Points = observed values. Dotted vertical line = start of forecast period.")),
						tabPanel("Predictor effects",
										 plotOutput("effects_plot", height = 480),
										 helpText("Signed posterior mean effect on relative abundance; bars are an",
										 				 "approximate 95% interval (mean +/- 1.96 SD).",
										 				 "Filled points are significant (95% credible interval excludes zero);",
										 				 "open points are not.")),
						tabPanel("Scores",
										 helpText("Forecast scoring metrics for this taxon across model types",
										 				 "and site-prediction modes (hindcast period)."),
										 DT::dataTableOutput("scores_table"))
					)
				)
			)
		)
	)
)

# --- Server --------------------------------------------------------------------
server <- function(input, output, session) {

	# One duckdb connection per session.
	con <- dbConnect(duckdb::duckdb())
	session$onSessionEnded(function() dbDisconnect(con, shutdown = TRUE))

	# Restrict the model-type choices to those that converged for the taxon,
	# keeping the current pick if still valid (else the preferred default).
	observeEvent(input$taxon, {
		mods <- converged_models(input$taxon)
		sel <- if (!is.null(input$model_name) && input$model_name %in% mods)
			input$model_name else mods[1]
		updateRadioButtons(session, "model_name", choices = mods, selected = sel)
	})

	# Full per-taxon + model time series (all shown site_prediction types).
	taxon_data <- reactive({
		req(input$taxon, input$model_name)
		query_taxon(con, input$taxon, input$model_name)
	})

	# Subset to the chosen site_prediction mode.
	mode_data <- reactive({
		d <- taxon_data()
		if (!nrow(d)) return(d)
		d[d$site_prediction == input$site_pred, , drop = FALSE]
	})

	# Update site choices when the taxon/model/mode changes.
	observeEvent(list(input$taxon, input$model_name, input$site_pred), {
		d <- mode_data()
		sites <- sort(unique(d$siteID))
		sel <- if (length(sites)) sites[1] else character(0)
		updateSelectInput(session, "siteID",
											choices = stats::setNames(sites, site_label(sites)),
											selected = sel)
	})

	# Update plot choices when the site changes.
	observeEvent(list(input$siteID, input$site_pred, input$taxon, input$model_name), {
		d <- mode_data()
		req(input$siteID)
		plots <- sort(unique(d$plotID[d$siteID == input$siteID]))
		updateSelectInput(session, "plotID",
											choices = c("All plots at site" = "all", plots),
											selected = "all")
	}, ignoreInit = TRUE)

	# Data for the forecast plot (site, optionally one plot).
	plot_data <- reactive({
		d <- mode_data()
		req(input$siteID, nrow(d) > 0)
		d <- d[d$siteID == input$siteID, , drop = FALSE]
		if (!is.null(input$plotID) && input$plotID != "all") {
			d <- d[d$plotID == input$plotID, , drop = FALSE]
		}
		d
	})

	# Tidy, clearly-named export of the selected taxon + model: all sites/plots,
	# both shown site-prediction modes (observed and new-site modeled).
	export_data <- reactive({
		d <- taxon_data()
		req(nrow(d) > 0)
		d <- d[d$site_prediction %in% site_pred_choices, , drop = FALSE]
		mode_lab <- stats::setNames(names(site_pred_choices), unname(site_pred_choices))
		out <- data.frame(
			taxon           = input$taxon,
			model           = model_label(input$model_name),
			siteID          = d$siteID,
			plotID          = d$plotID,
			date            = d$dates,
			period          = d$fcast_period,
			site_prediction = unname(mode_lab[d$site_prediction]),
			observed        = d$truth,
			mean            = d$mean,
			median          = d$med,
			lower_95        = d$lo,
			lower_50        = d$lo_25,
			upper_50        = d$hi_75,
			upper_95        = d$hi,
			stringsAsFactors = FALSE)
		out <- out[order(out$site_prediction, out$plotID, out$date), , drop = FALSE]
		num <- c("observed", "mean", "median", "lower_95", "lower_50", "upper_50", "upper_95")
		out[num] <- lapply(out[num], round, 5)
		out
	})

	output$dl_forecast <- downloadHandler(
		filename = function() {
			paste0(input$taxon, "_", input$model_name, "_forecast.csv")
		},
		content = function(file) {
			utils::write.csv(export_data(), file, row.names = FALSE, na = "")
		}
	)

	output$forecast_plot <- renderPlot({
		d <- plot_data()
		validate(need(nrow(d) > 0,
									"No forecasts available for this taxon / model / site combination."))

		ttl <- paste0(input$taxon, " - ", model_label(input$model_name),
									" (", input$siteID, ")")

		ggplot(d, aes(x = dates)) +
			facet_wrap(~ plotID, scales = "free_y") +
			geom_vline(xintercept = forecast_start, linetype = 3, color = "grey50") +
			geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#D55E00", alpha = 0.25) +
			geom_ribbon(aes(ymin = lo_25, ymax = hi_75), fill = "#0072B2", alpha = 0.35) +
			geom_line(aes(y = mean), linetype = 2, na.rm = TRUE) +
			geom_line(aes(y = med), na.rm = TRUE) +
			geom_point(aes(y = truth), size = 1.6, na.rm = TRUE) +
			theme_bw(base_size = 14) +
			labs(x = NULL, y = "Relative abundance", title = ttl) +
			theme(panel.spacing = unit(0.3, "cm"))
	})

	output$effects_plot <- renderPlot({
		validate(need(!is.null(predictor_effects),
									"predictor_effects.rds not available."))
		pe <- predictor_effects
		d <- pe[pe$taxon == input$taxon & pe$model_name == input$model_name, , drop = FALSE]
		validate(need(nrow(d) > 0,
									"No predictor effects for this taxon / model combination."))

		d$beta <- factor(as.character(d$beta))
		d$is_sig <- factor(ifelse(d$significant == 1, "Significant", "Not significant"),
											 levels = c("Significant", "Not significant"))

		# Plot the signed posterior mean with an approximate 95% interval
		# (mean +/- 1.96 SD). Significance is the stored credible-interval test
		# (interval excludes 0); drawing the bar at the same ~95% width makes the
		# open/filled coding consistent with whether the bar crosses zero.
		zCI <- 1.96
		ggplot(d, aes(x = stats::reorder(beta, effSize), y = Mean)) +
			geom_hline(yintercept = 0, color = "grey60") +
			geom_errorbar(aes(ymin = Mean - zCI * SD, ymax = Mean + zCI * SD),
										width = 0.2, color = "grey50") +
			geom_point(aes(shape = is_sig), size = 4, fill = "#0072B2", color = "#0072B2") +
			scale_shape_manual(values = c("Significant" = 21, "Not significant" = 1),
												 drop = FALSE, name = NULL) +
			coord_flip() +
			theme_bw(base_size = 14) +
			labs(x = NULL, y = "Effect size (posterior mean, bars ~ 95% interval)",
					 title = paste0("Predictor effects: ", input$taxon))
	})

	output$scores_table <- DT::renderDataTable({
		validate(need(!is.null(converged_set$metrics),
									"scoring_metrics_plsr2.rds not available."))
		sm <- converged_set$metrics
		d <- sm[sm$species == input$taxon &
							sm$site_prediction %in% site_pred_choices, , drop = FALSE]
		validate(need(nrow(d) > 0, "No scoring metrics for this taxon."))

		keep <- c("model_name", "site_prediction", "RMSE", "RMSE.norm",
							"CRPS", "RSQ", "BIAS", "MAE")
		keep <- intersect(keep, names(d))
		d <- d[, keep, drop = FALSE]
		num <- vapply(d, is.numeric, logical(1))
		d[num] <- lapply(d[num], round, 3)
		DT::datatable(d, rownames = FALSE, options = list(dom = "t", pageLength = 25))
	})

	# Human-readable context for the selected site code.
	output$site_info <- renderUI({
		req(input$siteID)
		if (is.null(site_meta)) return(NULL)
		r <- site_meta[site_meta$siteID == input$siteID, , drop = FALSE]
		if (!nrow(r)) return(NULL)
		nm <- if (is.na(r$siteName[1])) input$siteID else r$siteName[1]
		div(style = "padding:8px 12px; margin-bottom:8px; background:#f3f6f8; border-radius:4px; font-size:13px;",
				HTML(paste0(
					"<b>", input$siteID, " - ", nm, "</b><br>",
					"NEON domain: ", r$ecoregion[1],
					" &nbsp;|&nbsp; Lat/Lon: ", round(r$latitude[1], 3), ", ", round(r$longitude[1], 3),
					"<br>Land cover: ", r$nlcd[1])))
	})

	# Full browsable catalog of converged taxa/groups.
	output$catalog_table <- DT::renderDataTable({
		DT::datatable(model_catalog, rownames = FALSE, selection = "single",
									filter = "top",
									options = list(pageLength = 15,
																 order = list(list(1, "asc"), list(0, "asc"))))
	})

	# Clicking a catalog row loads that taxon and switches to the Explore page.
	observeEvent(input$catalog_table_rows_selected, {
		i <- input$catalog_table_rows_selected
		req(length(i) == 1)
		updateSelectInput(session, "taxon", selected = model_catalog$Taxon[i])
		updateTabItems(session, "menu", "explore")
	})
}

shinyApp(ui, server)
