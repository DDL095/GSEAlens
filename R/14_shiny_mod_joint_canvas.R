#' @title Joint GSEA Fill Canvas Module UI (Display-Only Version)
#' @description All controls have been moved to the sidebar; only canvas results are displayed here.
#' @keywords internal

mod_joint_canvas_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    # 🔧 Modification: Remove all control elements, keep only canvas display (full 12 columns)
    shiny::div(class = "white-box", style = "min-height: 900px;",
               shiny::h4("GSEA Joint Canvas"),
               shiny::uiOutput(ns("canvas_info")),
               shiny::plotOutput(ns("canvas_plot"), height = "auto", width = "100%")
    )
  )
}

#' @title Joint GSEA Fill Canvas Module Server (Sidebar Control Version)
#' @description
#'   1. Receive sidebar control parameters (arrangement mode)
#'   2. No maximum row limit
#'   3. Aesthetic consistency with multi_plot
#' @keywords internal

mod_joint_canvas_server <- function(id, gsea_res, data_prep_list, table_result) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # 存储当前画布结果
    canvas_result <- shiny::reactiveVal(NULL)

    # 🔧 Core: Monitor the generate button in sidebar (instead of local button)
    shiny::observeEvent(data_prep_list$joint_generate(), {

      # 🔧 Get parameters from data_prep_list (instead of input)
      contrasts <- data_prep_list$joint_contrasts()
      ncol_val <- data_prep_list$joint_ncol()

      shiny::req(contrasts, ncol_val)
      if (length(contrasts) == 0) {
        shiny::showNotification("Please select contrast groups (arrangement mode) in the left control panel first", type = "error")
        return()
      }

      # 🔧 Get plotting parameters from data_prep_list (consistent with multi_plot)
      data_list <- data_prep_list$data()
      shiny::req(data_list)

      plot_subtype <- data_list$plot_subtype
      custom_colors <- data_list$custom_colors

      # 🔧 Pathways: share selection with multi_plot
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

      # 🔧 Arrangement mode: strictly follow user-selected order (support both A_vs_B and B_vs_A coexisting)
      for (i in seq_along(contrasts)) {
        contrast_id <- contrasts[i]

        task_obj <- tryCatch({
          extract_gsea_task(gsea_res, contrast_id, "ALL")
        }, error = function(e) {
          message(sprintf("Extraction failed for %s: %s", contrast_id, e$message))
          NULL
        })

        if (is.null(task_obj)) {
          p <- ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0.5, y = 0.5,
                              label = sprintf("Not found:\n%s", contrast_id),
                              size = 3, color = "red") +
            ggplot2::theme_void()
          plot_list[[i]] <- p
          next
        }

        # 🔧 Aesthetic unification: consistent with multi_plot
        main_title <- sprintf(
          "%s [%d pathways]",
          gsub("_vs_", " vs ", contrast_id),
          length(pathways)
        )

        p <- tryCatch({
          plot_directional_gsea(
            directional_gsea_obj = task_obj,
            target_pathways = pathways,
            subPlot = as.numeric(plot_subtype),  # 使用主控制栏的subPlot
            curveCol = colors,
            main_title = main_title,
            add_pval = FALSE,
            show_contrast_in_axis = TRUE
          )
        }, error = function(e) {
          ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0.5, y = 0.5,
                              label = sprintf("Error:\n%s", substr(e$message, 1, 80)),
                              size = 3, color = "red") +
            ggplot2::theme_void()
        })

        plot_list[[i]] <- p

        if (i %% 2 == 0 || i == length(contrasts)) {
          shiny::incProgress(0.6 * i / length(contrasts),
                             detail = sprintf("Plotting %d/%d ...", i, length(contrasts)))
        }
      }

      if (length(plot_list) == 0) {
        shiny::showNotification("No valid plots can be generated", type = "error")
        return()
      }

      n_plots <- length(plot_list)
      # 🔧 Remove max row limit: auto-calculate
      actual_nrow <- ceiling(n_plots / ncol_val)

      # 🔧 Aesthetic optimization: consistent with multi_plot
      combined_plot <- patchwork::wrap_plots(
        plot_list,
        ncol = ncol_val,
        nrow = actual_nrow,
        byrow = TRUE,
        guides = "collect"
      ) + patchwork::plot_annotation(
        title = sprintf("Joint GSEA Canvas: %d contrast groups × %d pathways", n_plots, length(pathways)),
        subtitle = sprintf("Arrangement: %s", paste(contrasts, collapse = " → ")),
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

    # 渲染画布
    output$canvas_plot <- shiny::renderPlot({
      shiny::req(canvas_result())
      canvas_result()$plot
    }, height = function() {
      if (is.null(canvas_result())) return(900)
      nrow <- canvas_result()$nrow
      # 🔧 400px per row, no upper limit
      return(max(900, nrow * 400))
    }, width = 1600, res = 72)

    # 画布信息
    output$canvas_info <- shiny::renderUI({
      shiny::req(canvas_result())
      info <- canvas_result()

      shiny::tags$div(
        style = "margin-bottom: 15px; padding: 10px; background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; border-radius: 8px;",
        shiny::HTML(sprintf(
          "<strong>📊 Canvas Info:</strong> %d contrast groups | Layout: %d cols × %d rows | %d pathways<br>
           <small>Arrangement order: %s</small>",
          info$n_plots, info$ncol, info$nrow, length(info$pathways),
          paste(info$contrasts, collapse = " → ")
        ))
      )
    })

    # 🔧 Deletion: Save canvas function has been removed (per user request)

  })
}
