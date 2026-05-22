# Section: Multi-Plot Module ----

## Subsection: UI Component ----

#' @title Combined Pathway Plotting UI
#' @keywords internal

mod_multi_plot_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::div(
      style = "margin-top: 20px; padding: 15px; background-color: #f8f9fa; border-radius: 8px;",
      shiny::h4("Combined Pathway Plotting"),
      shiny::uiOutput(ns("selection_info")),
      shiny::selectizeInput(
        ns("pathway_selector"),
        label = "Selected Pathways (click x to remove):",
        choices = character(0),
        multiple = TRUE,
        options = list(plugins = list("remove_button"), placeholder = "Select pathways..."),
        width = "100%"
      ),
      shiny::actionButton(
        ns("clear_btn"),
        "Clear Selection",
        class = "btn-warning",
        style = "margin-bottom: 15px;"
      ),
      shiny::plotOutput(ns("multi_plot"), height = "800px", width = "100%")
    )
  )
}

## Subsection: Server Component ----

#' @title Combined Pathway Plotting Server
#' @keywords internal

mod_multi_plot_server <- function(id, data_prep, table_controller) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    selected_ids <- table_controller$selected_pathways
    updating <- shiny::reactiveVal(FALSE)

    # Sync selection to selectize
    shiny::observe({
      sel <- selected_ids()
      updating(TRUE)
      shiny::updateSelectizeInput(session, "pathway_selector", choices = sel, selected = sel)
      updating(FALSE)
    })

    # Listen for selectize removal
    shiny::observeEvent(input$pathway_selector,
      {
        if (updating()) {
          return()
        }
        current <- selected_ids()
        new_sel <- input$pathway_selector
        if (is.null(new_sel)) new_sel <- character(0)

        removed <- setdiff(current, new_sel)
        if (length(removed) > 0) {
          table_controller$remove_pathways(removed)
        }
      },
      ignoreInit = TRUE
    )

    # Clear button
    shiny::observeEvent(input$clear_btn, {
      table_controller$clear_selection()
    })

    # Render selection info
    output$selection_info <- shiny::renderUI({
      n <- length(selected_ids())
      if (n == 0) {
        shiny::tags$div(
          style = "color: #856404; background-color: #fff3cd; padding: 10px;",
          "Check the 'Combined Display' column in the pathway table to select pathways"
        )
      } else {
        shiny::tags$div(
          style = "color: #155724; background-color: #d4edda; padding: 10px;",
          sprintf("%d pathway(s) selected", n)
        )
      }
    })

    # Combined plot rendering
    output$multi_plot <- shiny::renderPlot({
      sel <- selected_ids()
      shiny::req(sel, length(sel) > 0)

      data_list <- data_prep()
      shiny::req(data_list)

      colors <- trimws(strsplit(data_list$custom_colors, ",")[[1]])
      if (length(colors) < length(sel)) {
        colors <- rep(colors, length.out = length(sel))
      }

      tryCatch(
        {
          # Single pathway: show pathway name; multiple: show combined count
          if (length(sel) == 1) {
            res_df <- as.data.frame(data_list$gsea_res@result)
            desc <- res_df$Description[res_df$ID == sel[1]]
            if (length(desc) == 0 || is.na(desc[1])) desc <- sel[1]
            plot_title <- desc[1]
          } else {
            plot_title <- sprintf("Combined Display: %d Pathway(s)", length(sel))
          }
          print(plot_directional_gsea(
            directional_gsea_obj = list(
              gsea_res = data_list$gsea_res,
              meta = list(
                left_group = data_list$left_group,
                right_group = data_list$right_group
              )
            ),
            target_pathways = sel,
            subPlot = data_list$plot_subtype,
            curveCol = colors,
            main_title = plot_title
          ))
        },
        error = function(e) {
          graphics::plot(1, type = "n", axes = FALSE)
          graphics::text(1, 1, sprintf("Plotting Error: %s", e$message), col = "red")
        }
      )
    })
  })
}
