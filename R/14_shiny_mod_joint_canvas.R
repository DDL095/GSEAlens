# Section: Joint Canvas Module ----

#' @title Joint GSEA Fill Canvas Module UI (Display-Only Version)
#' @description
#'   All controls have been moved to the sidebar; only canvas results are displayed here.
#'   This module provides a display interface for joint GSEA visualization results
#'   with export functionality for reproducible R code generation.
#' @keywords internal

mod_joint_canvas_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::div(
      class = "white-box", style = "min-height: 900px;",
      shiny::h4("GSEA Joint Canvas"),
      shiny::uiOutput(ns("canvas_info")),
      shiny::plotOutput(ns("canvas_plot"), height = "auto", width = "100%")
    ),
    shiny::hr(),
    shiny::fluidRow(
      shiny::column(12, shiny::div(
        class = "white-box",
        shiny::h4("Export Image / Code"),
        shiny::fluidRow(
          shiny::column(3,
            shiny::actionButton(
              ns("open_export_modal"),
              label = "Export Publication Plot",
              class = "btn-success",
              icon  = shiny::icon("file-image"),
              style = "width: 100%;"
            )
          ),
          shiny::column(3,
            shiny::actionButton(
              ns("export_code_btn"),
              label = "Export Joint Canvas Code",
              class = "btn-secondary",
              icon  = shiny::icon("code"),
              style = "width: 100%;"
            )
          ),
          shiny::column(6,
            shiny::helpText("Image export uses ggplot2::ggsave (PDF/PNG/SVG/TIFF). ",
                            "Code export generates a self-contained R script for reproduction.")
          )
        )
      ))
    )
  )
}

# Section: Joint Canvas Server Module ----

#' @title Joint GSEA Fill Canvas Module Server (Sidebar Control Version)
#' @description
#'   Server-side logic for the joint GSEA canvas display module.
#'   This module:
#'   \enumerate{
#'     \item Receives sidebar control parameters (arrangement mode)
#'     \item Applies no maximum row limit
#'     \item Maintains aesthetic consistency with multi_plot
#'     \item Generates combined visualization with multiple contrasts and pathways
#'     \item Provides R code export functionality for reproducibility
#'   }
#' @param id Character string identifying the module namespace.
#' @param gsea_res Reactive expression containing GSEA results.
#' @param data_prep_list Reactive list containing data preparation parameters
#'   including contrasts, ncol, and custom colors.
#' @param table_result Reactive list containing pathway selection results.
#' @keywords internal

mod_joint_canvas_server <- function(id, gsea_res, data_prep_list, table_result) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # 存储当前画布结果
    canvas_result <- shiny::reactiveVal(NULL)

    shiny::observeEvent(data_prep_list$joint_generate(), {
      contrasts <- data_prep_list$joint_contrasts()
      ncol_val <- data_prep_list$joint_ncol()

      shiny::req(contrasts, ncol_val)
      if (length(contrasts) == 0) {
        shiny::showNotification("Please select contrast groups in the left control panel first", type = "error")
        return()
      }

      data_list <- data_prep_list$data()
      shiny::req(data_list)

      plot_subtype <- data_list$plot_subtype
      custom_colors <- data_list$custom_colors

      pathways <- table_result$selected_pathways()
      shiny::req(pathways)
      if (length(pathways) == 0) {
        shiny::showNotification("Please check the pathways for 'Joint Display' in the main workspace first", type = "error")
        return()
      }

      # 验证ncol
      if (is.null(ncol_val) || ncol_val < 1) ncol_val <- 2

      # 解析颜色
      colors <- trimws(strsplit(custom_colors, ",")[[1]])
      if (length(colors) < length(pathways)) {
        colors <- rep(colors, length.out = length(pathways))
      }

      shiny::incProgress(0.1, detail = "Extracting GSEA tasks...")

      plot_list <- list()

      for (i in seq_along(contrasts)) {
        contrast_id <- contrasts[i]

        task_obj <- tryCatch(
          {
            extract_gsea_task(gsea_res, contrast_id, "ALL")
          },
          error = function(e) {
            message(sprintf("Extraction failed for %s: %s", contrast_id, e$message))
            NULL
          }
        )

        if (is.null(task_obj)) {
          p <- ggplot2::ggplot() +
            ggplot2::annotate("text",
              x = 0.5, y = 0.5,
              label = sprintf("Not found:\n%s", contrast_id),
              size = 3, color = "red"
            ) +
            ggplot2::theme_void()
          plot_list[[i]] <- p
          next
        }

        main_title <- sprintf(
          "%s [%d pathways]",
          gsub("_vs_", " vs ", contrast_id),
          length(pathways)
        )

        p <- tryCatch(
          {
            plot_directional_gsea(
              directional_gsea_obj = task_obj,
              target_pathways = pathways,
              subPlot = as.numeric(plot_subtype), # 使用主控制栏的subPlot
              curveCol = colors,
              main_title = main_title,
              add_pval = FALSE,
              show_contrast_in_axis = TRUE
            )
          },
          error = function(e) {
            ggplot2::ggplot() +
              ggplot2::annotate("text",
                x = 0.5, y = 0.5,
                label = sprintf("Error:\n%s", substr(e$message, 1, 80)),
                size = 3, color = "red"
              ) +
              ggplot2::theme_void()
          }
        )

        plot_list[[i]] <- p

        if (i %% 2 == 0 || i == length(contrasts)) {
          shiny::incProgress(0.6 * i / length(contrasts),
            detail = sprintf("Plotting %d/%d ...", i, length(contrasts))
          )
        }
      }

      if (length(plot_list) == 0) {
        shiny::showNotification("No valid plots can be generated", type = "error")
        return()
      }

      n_plots <- length(plot_list)
      actual_nrow <- ceiling(n_plots / ncol_val)

      combined_plot <- patchwork::wrap_plots(
        plot_list,
        ncol = ncol_val,
        nrow = actual_nrow,
        byrow = TRUE,
        guides = "collect"
      ) + patchwork::plot_annotation(
        title = sprintf("Joint GSEA Canvas: %d contrast groups x %d pathways", n_plots, length(pathways)),
        subtitle = sprintf("Arrangement: %s", paste(contrasts, collapse = " -> ")),
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
          plot.subtitle = ggplot2::element_text(size = 8, color = "gray50", hjust = 0.5),
          plot.margin = ggplot2::margin(10, 10, 10, 10)
        )
      )

      canvas_result(list(
        plot = combined_plot,
        n_plots = n_plots,
        ncol = ncol_val,
        nrow = actual_nrow,
        contrasts = contrasts,
        pathways = pathways
      ))

      shiny::incProgress(1.0, detail = "Done!")
    })


    # Subsection: Export Image Modal ----
    # ------------------------------------------------------------
    # Joint Canvas is a patchwork object (canvas_result()$plot), so we can
    # ggsave it directly without any code generation. This mirrors the
    # export centers in modules 10/13/16 but is simpler because the plot
    # object already exists (no eval() needed).
    # ------------------------------------------------------------
    .build_canvas_export_modal <- function() {
      shiny::modalDialog(
        title = "Export Publication Plot - Joint Canvas",
        size = "l", footer = NULL, easyClose = TRUE,
        shiny::fluidRow(
          shiny::column(5,
            shiny::h5("Dimensions"),
            shiny::fluidRow(
              shiny::column(6,
                shiny::numericInput(ns("jc_width"),  "Width (inch)",  value = 16, min = 4, max = 40)
              ),
              shiny::column(6,
                shiny::numericInput(ns("jc_height"), "Height (inch)", value = 12, min = 4, max = 40)
              )
            ),
            shiny::numericInput(ns("jc_dpi"), "DPI", value = 300, min = 72, max = 600),
            shiny::hr(),
            shiny::h5("Download"),
            shiny::fluidRow(
              shiny::column(6,
                shiny::downloadButton(ns("jc_download_pdf"),
                  "Download PDF",
                  class = "btn-danger btn-block",
                  icon  = shiny::icon("file-pdf"))
              ),
              shiny::column(6,
                shiny::downloadButton(ns("jc_download_png"),
                  "Download PNG",
                  class = "btn-primary btn-block",
                  icon  = shiny::icon("file-image"))
              )
            ),
            shiny::hr(),
            shiny::fluidRow(
              shiny::column(6,
                shiny::selectInput(ns("jc_format"), "Other format",
                  choices = c("SVG" = "svg", "TIFF" = "tiff"),
                  selected = "svg")
              ),
              shiny::column(6,
                shiny::div(style = "margin-top: 22px;",
                  shiny::downloadButton(ns("jc_download_other"),
                    "Download", class = "btn-default btn-block")
                )
              )
            ),
            shiny::hr(),
            shiny::actionButton(ns("jc_dismiss"), "Close",
                                class = "btn-default btn-block")
          ),
          shiny::column(7,
            shiny::h5("Live Preview"),
            shiny::div(style = "background:#f5f5f5; border:1px solid #ddd; border-radius:4px; padding:8px;",
              shiny::plotOutput(ns("jc_preview"), height = "460px") |>
                shinycssloaders::withSpinner(type = 6, color = "#28a745")
            ),
            shiny::tags$small(style = "color: #666; display:block; margin-top:6px;",
              "The Joint Canvas is a patchwork of GSEA sub-plots. ",
              "Larger widths are recommended for multi-contrast layouts.")
          )
        )
      )
    }

    shiny::observeEvent(input$open_export_modal, {
      if (is.null(canvas_result()) || is.null(canvas_result()$plot)) {
        shiny::showNotification("Generate a canvas first.", type = "warning"); return()
      }
      shiny::showModal(.build_canvas_export_modal())
    })
    shiny::observeEvent(input$jc_dismiss, shiny::removeModal())

    output$jc_preview <- shiny::renderPlot({
      cr <- canvas_result()
      shiny::req(cr, cr$plot)
      cr$plot
    })

    .jc_render_to_file <- function(file, fmt) {
      cr <- canvas_result()
      if (is.null(cr) || is.null(cr$plot)) {
        shiny::showNotification("No canvas to export.", type = "warning"); return()
      }
      ggplot2::ggsave(
        file, cr$plot,
        width  = if (is.null(input$jc_width))  16 else input$jc_width,
        height = if (is.null(input$jc_height)) 12 else input$jc_height,
        dpi    = if (is.null(input$jc_dpi))   300 else input$jc_dpi,
        device = fmt,
        limitsize = FALSE
      )
    }

    output$jc_download_pdf <- shiny::downloadHandler(
      filename = function() sprintf("GSEAlens_joint_canvas_%s.pdf", Sys.Date()),
      content  = function(file) .jc_render_to_file(file, "pdf")
    )
    output$jc_download_png <- shiny::downloadHandler(
      filename = function() sprintf("GSEAlens_joint_canvas_%s.png", Sys.Date()),
      content  = function(file) .jc_render_to_file(file, "png")
    )
    output$jc_download_other <- shiny::downloadHandler(
      filename = function() sprintf("GSEAlens_joint_canvas_%s.%s", Sys.Date(),
                                    if (is.null(input$jc_format)) "svg" else input$jc_format),
      content  = function(file) .jc_render_to_file(file, input$jc_format)
    )

    # Subsection: Export Code Modal ----

    shiny::observeEvent(input$export_code_btn, {
      shiny::showModal(shiny::modalDialog(
        title = "Generated R Code",
        size = "l",
        easyClose = TRUE,
        shiny::fluidRow(
          shiny::column(
            12,
            shiny::div(
              style = "background: #f5f5f5; padding: 15px; border-radius: 5px; max-height: 500px; overflow: auto;",
              shiny::tags$pre(
                shiny::code(
                  generate_joint_canvas_code(
                    GSEAlens_res = gsea_res,
                    contrast_ids = data_prep_list$joint_contrasts(),
                    target_pathways = table_result$selected_pathways(),
                    ncol = data_prep_list$joint_ncol()
                  )
                ),
                style = "font-size: 11px; white-space: pre-wrap; word-break: break-all;"
              )
            )
          )
        ),
        footer = shiny::modalButton("Close")
      ))
    })


    # Subsection: Render Canvas ----

    # 渲染画布
    output$canvas_plot <- shiny::renderPlot(
      {
        shiny::req(canvas_result())
        canvas_result()$plot
      },
      height = function() {
        if (is.null(canvas_result())) {
          return(900)
        }
        nrow <- canvas_result()$nrow
        return(max(900, nrow * 400))
      },
      width = 1600,
      res = 72
    )

    # 画布信息
    output$canvas_info <- shiny::renderUI({
      shiny::req(canvas_result())
      info <- canvas_result()

      shiny::tags$div(
        style = "margin-bottom: 15px; padding: 10px; background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; border-radius: 8px;",
        shiny::HTML(sprintf(
          "<strong>Canvas Info:</strong> %d contrast groups | Layout: %d cols x %d rows | %d pathways<br>
           <small>Arrangement order: %s</small>",
          info$n_plots, info$ncol, info$nrow, length(info$pathways),
          paste(info$contrasts, collapse = " -> ")
        ))
      )
    })
  })
}
