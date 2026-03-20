#' @title 主工作台表格 UI
#' @keywords internal
mod_master_table_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::div(
      class = "master-table-container",
      style = "width: 100%; overflow-x: auto;",
      shiny::tags$div(
        style = "margin-bottom: 10px; color: #666; font-size: 12px;",
        shiny::HTML(.tr("table.tooltip_drag"))
      ),
      DT::dataTableOutput(ns("table"))
    ),
    shiny::tags$div(
      style = "display: none;",
      shiny::textInput(ns("joint_selection_store"), label = NULL, value = "")
    )
  )
}

#' @title 主工作台表格 Server
#' @param id 模块 ID
#' @param data_prep 来自数据预处理模块的响应式数据
#' @return 列表，包含 selected_pathways (reactive) 和 show_modal (reactive)
#' @keywords internal
mod_master_table_server <- function(id, data_prep) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    joint_selected <- shiny::reactiveVal(character(0))
    has_interaction <- shiny::reactiveVal(FALSE)

    to_safe <- function(x) gsub("'", "\\\\'", x, fixed = TRUE)
    to_original <- function(x) gsub("\\\\'", "'", x, fixed = TRUE)

    shiny::observeEvent(input$joint_plot_toggle, {
      has_interaction(TRUE)
      toggle <- input$joint_plot_toggle
      if (is.null(toggle) || !is.list(toggle)) return()

      current <- joint_selected()
      if (isTRUE(toggle$checked)) {
        if (!(toggle$id %in% current)) joint_selected(c(current, toggle$id))
      } else {
        joint_selected(setdiff(current, toggle$id))
      }
    })

    output$table <- DT::renderDataTable({
      data_list <- data_prep()
      shiny::validate(shiny::need(data_list, "等待数据加载..."))
      shiny::validate(shiny::need(nrow(data_list$df) > 0, "当前筛选条件下无数据"))

      df <- data_list$df

      # 计算 Enriched_In 列
      left_grp <- data_list$left_group
      right_grp <- data_list$right_group
      df$Enriched_In <- ifelse(df$NES > 0, left_grp, right_grp)
      df$Enriched_In <- factor(df$Enriched_In, levels = c(left_grp, right_grp))

      # 合并 Description
      if (!is.null(data_list$task_obj$meta$meta_dict)) {
        meta_dict <- data_list$task_obj$meta$meta_dict
        if ("Brief_Description" %in% colnames(meta_dict)) {
          desc_map <- setNames(meta_dict$Brief_Description, meta_dict$ID)
          df$Description <- desc_map[df$ID]
        } else if ("Description" %in% colnames(meta_dict)) {
          desc_map <- setNames(meta_dict$Description, meta_dict$ID)
          df$Description <- desc_map[df$ID]
        } else {
          df$Description <- df$ID
        }
      } else {
        df$Description <- df$ID
      }
      df$Description[is.na(df$Description)] <- df$ID[is.na(df$Description)]

      # 格式化数值
      df$NES_display <- round(df$NES, 2)
      df$pvalue_display <- format(df$pvalue, digits = 2, scientific = TRUE)
      df$padj_display <- format(df$p.adjust, digits = 2, scientific = TRUE)

      # 生成交互列（HTML）
      df$Select_for_Plot <- sprintf(
        '<input type="checkbox" class="joint-plot-checkbox" data-id="%s" %s onclick="Shiny.setInputValue(&#39;%s&#39;, {id: &#39;%s&#39;, checked: this.checked}, {priority: &#39;event&#39;});"/>',
        df$Safe_ID,
        ifelse(df$Safe_ID %in% joint_selected(), 'checked="checked"', ''),
        ns("joint_plot_toggle"),
        df$Safe_ID
      )

      df$Detail_Btn <- sprintf(
        '<button class="btn btn-sm btn-success" onclick="Shiny.setInputValue(&#39;%s&#39;, &#39;%s&#39;, {priority: &#39;event&#39;})">🔍 %s</button>',
        ns("show_modal"),
        df$Safe_ID,
        "Dashboard"
      )

      # 🔧 关键修复：直接重命名数据框列（而非使用colnames参数）
      # 这样escape=TRUE也能正常工作，因为列名就是实际名称
      display_cols <- c("Rank", "Select_for_Plot", "Detail_Btn", "ID", "Enriched_In",
                        "NES_display", "pvalue_display", "padj_display", "setSize", "Description")
      display_cols <- intersect(display_cols, colnames(df))

      # 提取数据
      dt_data <- df[, display_cols]

      # 🔧 关键修复：直接修改列名（而非使用colnames参数）
      names(dt_data) <- c(
        "Rank", "联合展示", "Dashboard", "Pathway", "Enriched In",
        "NES", "P-value", "FDR", "Size", "Description"
      )[1:length(display_cols)]

      # 构建DT（不使用colnames参数，避免escape检查冲突）
      dt <- DT::datatable(
        dt_data,
        escape = FALSE,  # 因为包含HTML按钮
        selection = "multiple",
        rownames = FALSE,
        extensions = c('Buttons', 'Scroller', 'ColReorder'),
        options = list(
          scrollX = TRUE,
          scrollY = "60vh",
          scroller = TRUE,
          deferRender = TRUE,
          pageLength = -1,
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel'),
          colReorder = TRUE,
          columnDefs = list(
            list(width = '60px', targets = 1, className = 'dt-center', orderable = FALSE),
            list(width = '90px', targets = 2, orderable = FALSE),
            list(width = '100px', targets = 4),  # Enriched In
            list(width = '300px', targets = 9, className = 'description-cell')  # Description
          ),
          order = list(list(0, 'asc'))
        )
      )

      # formatStyle使用重命名后的列名（中文）
      dt <- dt %>% DT::formatStyle(
        columns = "Enriched In",
        backgroundColor = DT::styleEqual(c(left_grp, right_grp), c('#fee0d2', '#deebf7')),
        color = DT::styleEqual(c(left_grp, right_grp), c('#cc212f', '#1052bd')),
        fontWeight = 'bold'
      )

      dt <- dt %>% DT::formatStyle(
        columns = "NES",
        color = DT::styleInterval(0, c('#1052bd', '#cc212f')),
        fontWeight = 'bold'
      )

      dt <- dt %>% DT::formatStyle(
        columns = "FDR",
        backgroundColor = DT::styleInterval(c(0.01, 0.05), c('#fc9272', '#fee0d2', 'transparent')),
        fontWeight = DT::styleInterval(0.05, c('bold', 'normal'))
      )

      dt
    }, server = TRUE)

    selected_pathways <- shiny::reactive({
      if (isTRUE(has_interaction())) return(to_original(joint_selected()))
      data_list <- data_prep()
      if (is.null(data_list) || is.null(input$table_rows_selected)) return(character(0))
      return(data_list$df$ID[input$table_rows_selected])
    })

    show_modal_trigger <- shiny::reactive({ input$show_modal })

    remove_pathways <- function(ids) {
      ids <- ids[!is.na(ids)]
      if (length(ids) == 0) return()
      has_interaction(TRUE)
      safe_ids <- to_safe(ids)
      new_sel <- setdiff(joint_selected(), safe_ids)
      joint_selected(new_sel)
    }

    clear_selection <- function() {
      has_interaction(TRUE)
      joint_selected(character(0))
    }

    return(list(
      selected_pathways = selected_pathways,
      show_modal = show_modal_trigger,
      remove_pathways = remove_pathways,
      clear_selection = clear_selection
    ))
  })
}
