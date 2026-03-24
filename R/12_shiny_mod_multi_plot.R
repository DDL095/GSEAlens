#' @title 联合绘图 UI
#' @keywords internal
#' @noRd
mod_multi_plot_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::div(
      style = "margin-top: 20px; padding: 15px; background-color: #f8f9fa; border-radius: 8px;",
      shiny::h4("🎨 联合通路绘图"),
      shiny::uiOutput(ns("selection_info")),
      shiny::selectizeInput(
        ns("pathway_selector"),
        label = "当前已选择通路（可点击 × 删除）:",
        choices = character(0), multiple = TRUE,
        options = list(plugins = list("remove_button"), placeholder = "当前勾选的通路..."),width = "100%",
      ),
      shiny::actionButton(ns("clear_btn"), "🗑️ 清空选择", class = "btn-warning",
                          style = "margin-bottom: 15px;"),
      shiny::plotOutput(ns("multi_plot"), height = "800px",width = "100%")
    )
  )
}

#' @title 联合绘图 Server
#' @keywords internal
#' @noRd
mod_multi_plot_server <- function(id, data_prep, table_controller) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    selected_ids <- table_controller$selected_pathways
    updating <- shiny::reactiveVal(FALSE)

    # 同步选择到 selectize
    shiny::observe({
      sel <- selected_ids()
      updating(TRUE)
      shiny::updateSelectizeInput(session, "pathway_selector", choices = sel, selected = sel)
      updating(FALSE)
    })

    # 监听 selectize 删除
    shiny::observeEvent(input$pathway_selector, {
      if (updating()) return()
      current <- selected_ids()
      new_sel <- input$pathway_selector
      if (is.null(new_sel)) new_sel <- character(0)

      removed <- setdiff(current, new_sel)
      if (length(removed) > 0) {
        table_controller$remove_pathways(removed)
      }
    }, ignoreInit = TRUE)

    # 清空按钮
    shiny::observeEvent(input$clear_btn, {
      table_controller$clear_selection()
    })

    # 显示选择信息
    output$selection_info <- shiny::renderUI({
      n <- length(selected_ids())
      if (n == 0) {
        shiny::tags$div(style = "color: #856404; background-color: #fff3cd; padding: 10px;",
                        "⚠️ 请在工作台表格中勾选\"联合展示\"列")
      } else {
        shiny::tags$div(style = "color: #155724; background-color: #d4edda; padding: 10px;",
                        sprintf("✅ 已选择 %d 条通路", n))
      }
    })

    # 联合绘图
    output$multi_plot <- shiny::renderPlot({
      sel <- selected_ids()
      shiny::req(sel, length(sel) > 0)

      data_list <- data_prep()
      shiny::req(data_list)

      colors <- trimws(strsplit(data_list$custom_colors, ",")[[1]])
      if (length(colors) < length(sel)) {
        colors <- rep(colors, length.out = length(sel))
      }

      tryCatch({
        print(plot_directional_gsea(
          directional_gsea_obj = list(
            gsea_res = data_list$gsea_res,
            meta = list(left_group = data_list$left_group, right_group = data_list$right_group)
          ),
          target_pathways = sel,
          subPlot = data_list$plot_subtype,
          curveCol = colors,
          main_title = sprintf("联合展示: %d 条通路", length(sel))
        ))
      }, error = function(e) {
        graphics::plot(1, type = "n", axes = FALSE)
        graphics::text(1, 1, sprintf("绘图错误: %s", e$message), col = "red")
      })
    })
  })
}
