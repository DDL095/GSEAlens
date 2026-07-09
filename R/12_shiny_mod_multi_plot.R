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

            shiny::plotOutput(ns("multi_plot"), height = "800px", width = "100%"),

            shiny::hr(),

            shiny::actionButton(

        ns("open_combined_export"),

        label = "Export Publication Plot",

        class = "btn-success",

        icon  = shiny::icon("file-image"),

        style = "width: 100%;"

            )

    )

    )

}



## Subsection: Server Component ----



#' @title Combined Pathway Plotting Server

#' @keywords internal



mod_multi_plot_server <- function(id, data_prep, table_controller, gsea_res = NULL) {

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



    # Combined plot reactive (extracted so the export modal can reuse it).
    # Returns a patchwork object from plot_directional_gsea, or NULL on error.
    .combined_plot_error <- shiny::reactiveVal(NULL)

    .combined_plot <- shiny::reactive({
            sel <- selected_ids()
            shiny::req(sel, length(sel) > 0)

            data_list <- data_prep()
            shiny::req(data_list)

            colors <- trimws(strsplit(data_list$custom_colors, ",")[[1]])
            if (length(colors) < length(sel)) {
        colors <- rep(colors, length.out = length(sel))
            }

            .combined_plot_error(NULL)
            tryCatch(
        {
                    if (!is.null(gsea_res)) {
            contrast_id <- data_list$contrast_id
            full_task <- extract_gsea_task(gsea_res, contrast_id, "ALL")
            plot_directional_gsea(
                            directional_gsea_obj = full_task,
                            target_pathways = sel,
                            subPlot = data_list$plot_subtype,
                            curveCol = colors
            )
                    } else {
            plot_directional_gsea(
                            directional_gsea_obj = list(
                gsea_res = data_list$gsea_res,
                meta = list(
                                    left_group = data_list$left_group,
                                    right_group = data_list$right_group
                )
                            ),
                            target_pathways = sel,
                            subPlot = data_list$plot_subtype,
                            curveCol = colors
            )
                    }
        },
        error = function(e) {
                    .combined_plot_error(e$message)
                    NULL
        }
            )
    })

    # Combined plot rendering
    output$multi_plot <- shiny::renderPlot({
            p <- .combined_plot()
            if (is.null(p)) {
        err <- .combined_plot_error()
        graphics::plot(1, type = "n", axes = FALSE)
        graphics::text(1, 1,
                    sprintf("Plotting Error: %s", if (is.null(err)) "unknown" else err),
                    col = "red")
        return(invisible(NULL))
            }
            p
    })

    # ============================================================
    # Export Image Modal (mirrors Joint Canvas / module 14)
    # ============================================================

    .build_combined_export_modal <- function() {
            shiny::modalDialog(
        title = "Export Publication Plot - Combined Pathway Plotting",
        size = "l", footer = NULL, easyClose = TRUE,
        shiny::fluidRow(
                    shiny::column(5,
            shiny::h5("Dimensions"),
            shiny::fluidRow(
                            shiny::column(6,
                shiny::numericInput(ns("cp_width"),  "Width (inch)",  value = 12, min = 4, max = 40)
                            ),
                            shiny::column(6,
                shiny::numericInput(ns("cp_height"), "Height (inch)", value = 10, min = 4, max = 40)
                            )
            ),
            shiny::numericInput(ns("cp_dpi"), "DPI", value = 300, min = 72, max = 600),
            shiny::hr(),
            shiny::h5("Download"),
            shiny::fluidRow(
                            shiny::column(6,
                shiny::downloadButton(ns("cp_download_pdf"),
                                    "Download PDF",
                                    class = "btn-danger btn-block",
                                    icon  = shiny::icon("file-pdf"))
                            ),
                            shiny::column(6,
                shiny::downloadButton(ns("cp_download_png"),
                                    "Download PNG",
                                    class = "btn-primary btn-block",
                                    icon  = shiny::icon("file-image"))
                            )
            ),
            shiny::hr(),
            shiny::fluidRow(
                            shiny::column(6,
                shiny::selectInput(ns("cp_format"), "Other format",
                                    choices = c("SVG" = "svg", "TIFF" = "tiff"),
                                    selected = "svg")
                            ),
                            shiny::column(6,
                shiny::div(style = "margin-top: 22px;",
                                    shiny::downloadButton(ns("cp_download_other"),
                    "Download", class = "btn-default btn-block")
                )
                            )
            ),
            shiny::hr(),
            shiny::fluidRow(
                            shiny::column(6,
                shiny::actionButton(ns("cp_copy_code"),
                                    "Copy R Code",
                                    class = "btn-info btn-block",
                                    icon = shiny::icon("clipboard"))
                            ),
                            shiny::column(6,
                shiny::actionButton(ns("cp_dismiss"), "Close",
                                    class = "btn-default btn-block")
                            )
            )
                    ),
                    shiny::column(7,
            shiny::h5("Live Preview"),
            shiny::div(style = "background:#f5f5f5; border:1px solid #ddd; border-radius:4px; padding:8px; width:100%; height:520px; display:flex; align-items:center; justify-content:center; overflow:auto;",
                            shiny::plotOutput(ns("cp_preview")) |>
                shinycssloaders::withSpinner(type = 6, color = "#28a745")
            ),
            shiny::tags$small(style = "color: #666; display:block; margin-top:6px;",
                            "Combined Pathway Plotting uses plot_directional_gsea (patchwork). Preview auto-scaled to fit while keeping the export aspect ratio.")
                    )
        )
            )
    }

    shiny::observeEvent(input$open_combined_export, {
            if (is.null(.combined_plot())) {
        shiny::showNotification("Select pathways first.", type = "warning"); return()
            }
            shiny::showModal(.build_combined_export_modal())
    })
    shiny::observeEvent(input$cp_dismiss, shiny::removeModal())

    # Preview device matching export physical size (res fixed at 72).
    .cp_preview_dims <- function(w_in, h_in, res = 72, max_px = 1600) {
            w_in <- suppressWarnings(as.numeric(w_in[1]))
            h_in <- suppressWarnings(as.numeric(h_in[1]))
            if (length(w_in) == 0 || is.na(w_in) || w_in <= 0) w_in <- 12
            if (length(h_in) == 0 || is.na(h_in) || h_in <= 0) h_in <- 10
            scale <- min(1, max_px / (max(w_in, h_in) * res))
            list(width  = round(w_in * res * scale),
                        height = round(h_in * res * scale))
    }

    output$cp_preview <- shiny::renderPlot(
            {
        .combined_plot()
            },
            res  = 72,
            width  = function() {
        d <- .cp_preview_dims(input$cp_width, input$cp_height); d$width
            },
            height = function() {
        d <- .cp_preview_dims(input$cp_width, input$cp_height); d$height
            }
    )

    .cp_render_to_file <- function(file, fmt) {
            p <- .combined_plot()
            if (is.null(p)) {
        shiny::showNotification("No plot to export.", type = "warning"); return()
            }
            ggplot2::ggsave(
        file, p,
        width  = if (is.null(input$cp_width))  12 else input$cp_width,
        height = if (is.null(input$cp_height)) 10 else input$cp_height,
        dpi    = if (is.null(input$cp_dpi))   300 else input$cp_dpi,
        device = fmt,
        limitsize = FALSE
            )
    }

    output$cp_download_pdf <- shiny::downloadHandler(
            filename = function() sprintf("GSEAlens_combined_plot_%s.pdf", Sys.Date()),
            content  = function(file) .cp_render_to_file(file, "pdf")
    )
    output$cp_download_png <- shiny::downloadHandler(
            filename = function() sprintf("GSEAlens_combined_plot_%s.png", Sys.Date()),
            content  = function(file) .cp_render_to_file(file, "png")
    )
    output$cp_download_other <- shiny::downloadHandler(
            filename = function() sprintf("GSEAlens_combined_plot_%s.%s", Sys.Date(),
                                    if (is.null(input$cp_format)) "svg" else input$cp_format),
            content  = function(file) .cp_render_to_file(file, input$cp_format)
    )

    # ----- Code export (mirrors module 13/14 pattern) -----
    # Generates self-contained R code reproducing the current
    # plot_directional_gsea() call so the user can reproduce the figure.
    .cp_export_code <- function() {
            sel <- selected_ids()
            if (is.null(sel) || length(sel) == 0) return("")
            data_list <- data_prep()
            if (is.null(data_list)) return("")

            colors <- trimws(strsplit(data_list$custom_colors, ",")[[1]])
            if (length(colors) == 0) colors <- NULL

            lg <- if (is.null(data_list$left_group))  "Left"  else data_list$left_group
            rg <- if (is.null(data_list$right_group)) "Right" else data_list$right_group
            sub <- if (is.null(data_list$plot_subtype)) 3 else data_list$plot_subtype

            generate_combined_plot_code(
        gsea_res_var   = "gsea_res",
        contrast_id    = data_list$contrast_id,
        target_pathways = sel,
        subPlot        = sub,
        curve_colors   = colors,
        left_group     = lg,
        right_group    = rg
            )
    }

    shiny::observeEvent(input$cp_copy_code, {
            code <- .cp_export_code()
            if (!nzchar(code)) {
        shiny::showNotification("Nothing to copy.", type = "warning"); return()
            }
            ok <- tryCatch({ clipr::write_clip(code); TRUE }, error = function(e) FALSE)
            if (isTRUE(ok)) {
        shiny::showNotification("R code copied to clipboard.", type = "message")
            } else {
        # Fallback: show code in a new modal so the user can manually copy
        shiny::showModal(shiny::modalDialog(
                    title = "Reproducible R Code (copy manually)",
                    size = "l", easyClose = TRUE,
                    shiny::div(
            style = "background:#f5f5f5; padding:15px; border-radius:5px; max-height:500px; overflow:auto;",
            shiny::tags$pre(shiny::code(code),
                            style = "font-size:11px; white-space:pre-wrap; word-break:break-all;")
                    ),
                    footer = shiny::modalButton("Close")
        ))
            }
    })
    })
}

