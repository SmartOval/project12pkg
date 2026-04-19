#' Launch the marker gene Shiny app
#'
#' Opens an interactive Shiny application for filtering genes,
#' identifying marker genes, and visualizing results from the
#' example single-cell dataset.
#'
#' @return A shiny application object.
#' @export
run_app <- function() {
  ui <- shiny::fluidPage(
    shiny::titlePanel("Marker Gene Explorer"),

    shiny::sidebarPanel(
      shiny::sliderInput(
        "min_detect_rate",
        "Minimum Detection Rate",
        min = 0.01,
        max = 0.20,
        value = 0.05,
        step = 0.01
      ),

      shiny::tags$div(                                     # I used AI (ChatGPT) here to help me create a textbox to make a Warning Label
        style = "
      margin-top: 10px;
      padding: 10px;
      background-color: #fff3cd;
      border: 1px solid #ffeeba;
      color: #856404;
      border-radius: 4px;
      font-size: 13px;
    ",
        shiny::strong("Warning: "),
        "Marker analysis can take up to 8 minutes to complete. ",
        "Please ensure that your parameters are correct before waiting for the plot to display."
      ),

      shiny::numericInput(
        "n_markers",
        "Number of Markers per Cell Type",
        value = 5,
        min = 1,
        max = 20,
        step = 1
      )
    ),
      shiny::mainPanel(
        shiny::tabsetPanel(
          shiny::tabPanel(
            "Dot Plot",
            shiny::plotOutput("marker_plot", height = "700px")
          ),
          shiny::tabPanel(
            "Summary",
            shiny::tableOutput("marker_summary")
          ),
          shiny::tabPanel(
            "Top Markers",
            shiny::tableOutput("top_markers")
          )
        )
      )
    )
  server <- function(input, output, session) {
    example_sce_obj <- get("example_sce", envir = asNamespace("project12pkg"))

    analysis_results <- shiny::reactive({
      shiny::validate(
        shiny::need(input$n_markers >= 1, "Number of Markers must be at least 1")
      )

      withProgress(message = "Running marker analysis... (This may take a bit)", value = 0, {
        incProgress(0.2)
        sce_filt <- filter_genes(example_sce_obj, min_detect_rate = input$min_detect_rate)

        incProgress(0.2)
        mat_norm <- normalize_counts(sce_filt)

        incProgress(0.3)
        all_markers <- find_markers(sce_filt, mat_norm, group_col = "label")

        incProgress(0.1)
        top_markers <- select_top_markers(all_markers, n_markers = input$n_markers)
        marker_summary <- summarize_markers(all_markers)

        incProgress(0.2)
        cell_type <- SummarizedExperiment::colData(sce_filt)[["label"]]
        plot_df <- build_dotplot(all_markers, top_markers, mat_norm, cell_type)

        list(
          sce_filt = sce_filt,
          mat_norm = mat_norm,
          all_markers = all_markers,
          top_markers = top_markers,
          marker_summary = marker_summary,
          plot_df = plot_df
        )
      })
    })

    output$marker_plot <- shiny::renderPlot({
      res <- analysis_results()

      shiny::validate(
        shiny::need(!is.null(res$plot_df), "Plot data was not created."),
        shiny::need(nrow(res$plot_df) > 0, "No rows available to plot.")
      )

      print(plot_markers(res$plot_df))
    })

    output$marker_summary <- shiny::renderTable({
      analysis_results()$marker_summary
    })

    output$top_markers <- shiny::renderTable({
      analysis_results()$top_markers
    })
  }
  shiny::shinyApp(ui = ui, server = server)
}
