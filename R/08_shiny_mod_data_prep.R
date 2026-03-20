#' @title 数据预处理模块 UI
#' @keywords internal
mod_data_prep_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::selectInput(
      ns("selected_contrast"),
      label = .tr("ui.select_contrast"),
      choices = NULL
    ),
    shiny::hr(),
    shiny::h4(.tr("ui.plot_control")),
    shiny::selectInput(
      ns("plot_subtype"),
      label = .tr("ui.gseavis_style"),
      choices = c(
        "1: 仅经典富集" = "1",
        "2: 富集+热图带" = "2",
        "3: 完整带Rank" = "3"
      ),
      selected = "3"
    ),
    shiny::textInput(
      ns("custom_colors"),
      label = .tr("ui.custom_colors"),
      value = "#E41A1C, #377EB8, #4DAF4A, #984EA3",
      placeholder = "例如: #FF0000, #00FF00, #0000FF"
    ),
    shiny::hr(),
    shiny::h4(.tr("ui.expression_metric")),
    shiny::selectInput(
      ns("expression_type"),
      label = .tr("ui.expr_type"),
      choices = c(
        "log2(CPM)" = "logcpm",
        "CPM (原始)" = "cpm",
        "VST (DESeq2)" = "vst",
        "log2(FPKM)" = "logfpkm",
        "FPKM (原始)" = "fpkm"
      ),
      selected = "logcpm"
    ),
    shiny::hr(),
    shiny::h4(.tr("ui.data_slice")),
    shiny::selectizeInput(
      ns("selected_collections"),
      label = .tr("ui.select_collections"),
      choices = NULL,
      multiple = TRUE,
      options = list(
        plugins = list("remove_button"),
        placeholder = "选择亚组或保留ALL"
      )
    ),
    shiny::selectInput(
      ns("sort_by"),
      label = .tr("ui.sort_strategy"),
      choices = c(
        "按 NES (降序)" = "nes_desc",
        "按 NES (升序)" = "nes_asc",
        "按 NES 绝对值 (降序)" = "abs_nes_desc",
        "按 P-value (升序)" = "pval_asc",
        "按 FDR (升序)" = "fdr_asc"
      ),
      selected = "abs_nes_desc"
    ),
    shiny::actionButton(
      ns("run_btn"),
      label = .tr("ui.btn_confirm"),
      class = "btn-success",
      style = "width: 100%; font-weight: bold; margin-top: 15px;"
    )
  )
}

#' @title 数据预处理模块 Server
#' @description 处理对比组选择、数据切片、标准化数据流
#' @param id 模块 ID
#' @param gsea_res GseaRes 对象
#' @return 响应式列表，包含处理后的数据
#' @keywords internal
mod_data_prep_server <- function(id, gsea_res) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    result_data <- shiny::reactiveVal(NULL)
    has_initialized <- shiny::reactiveVal(FALSE)

    # 🔧 修改：动态填充对比组 - 包含正向和反向（共12个）
    shiny::observe({
      registry <- gsea_res$contrast_registry

      # 🆕 新增：生成所有对比（正向6个 + 反向6个 = 12个）
      all_contrasts <- list()

      for (i in 1:nrow(registry)) {
        row <- registry[i, ]

        # 正向对比（如 IBAA vs PREA）
        all_contrasts[[row$contrast_id]] <- paste(row$left_group, "vs", row$right_group)

        # 🆕 新增：反向对比（如 PREA vs IBAA）
        # 生成反向ID：交换left和right
        rev_id <- paste(row$right_group, row$left_group, sep = "_vs_")
        rev_display <- paste(row$right_group, "vs", row$left_group)
        all_contrasts[[rev_id]] <- rev_display
      }

      # 转换为命名向量（显示名 = 值）
      choices <- setNames(names(all_contrasts), unlist(all_contrasts))

      shiny::updateSelectInput(
        session,
        "selected_contrast",
        choices = choices,
        selected = registry$contrast_id[1]  # 默认选第一个正向
      )
    })





    # 核心数据处理函数 (适配新版架构)
    process_data_core <- function(contrast_id, collections, sort_by,
                                  expr_type, plot_subtype, colors, is_auto = FALSE) {

      if (is.null(contrast_id)) return(NULL)

      # 使用新版提取函数
      task_obj <- tryCatch({
        extract_gsea_task(gsea_res, contrast_id, target_collection = collections)
      }, error = function(e) {
        message(sprintf("提取任务失败: %s", e$message))
        return(NULL)
      })

      if (is.null(task_obj)) return(NULL)

      # 获取对比组信息
      meta <- task_obj$meta
      gsea_res_obj <- task_obj$gsea_res
      df <- as.data.frame(gsea_res_obj@result)

      if (nrow(df) == 0) return(NULL)

      # 处理排序
      df$abs_NES <- abs(df$NES)
      sort_config <- switch(sort_by,
                            "nes_desc" = list(col = "NES", desc = TRUE),
                            "nes_asc" = list(col = "NES", desc = FALSE),
                            "abs_nes_desc" = list(col = "abs_NES", desc = TRUE),
                            "pval_asc" = list(col = "pvalue", desc = FALSE),
                            "fdr_asc" = list(col = "p.adjust", desc = FALSE),
                            list(col = "abs_NES", desc = TRUE)
      )

      sort_order <- order(df[[sort_config$col]], decreasing = sort_config$desc)
      df <- df[sort_order, ]
      df$Rank <- seq_len(nrow(df))

      # 生成 Safe_ID (用于前端交互)
      df$Safe_ID <- gsub("'", "\\\\'", df$ID)

      # 获取当前对比组在 registry 中的信息
      reg_row <- gsea_res$contrast_registry[gsea_res$contrast_registry$contrast_id == contrast_id, ]

      if (nrow(reg_row) == 0) {
        # 可能是反向对比，尝试解析
        parts <- strsplit(contrast_id, "_vs_")[[1]]
        if (length(parts) == 2) {
          left <- parts[1]
          right <- parts[2]
        } else {
          left <- contrast_id
          right <- "Background"
        }
      } else {
        left <- reg_row$left_group[1]
        right <- reg_row$right_group[1]
      }

      # 构建返回数据结构
      list(
        df = df,
        gsea_res = gsea_res_obj,
        task_obj = task_obj,
        left_group = left,
        right_group = right,
        contrast_id = contrast_id,
        contrast_name = sprintf("%s vs %s", left, right),  # 简化显示名，用于图表标题
        expression_type = expr_type,
        plot_subtype = as.numeric(plot_subtype),
        custom_colors = colors,
        is_preview = is_auto,
        backend = gsea_res$backend_info$backend
      )
    }

    # 自动初始化 (选择最小基因集)
    shiny::observeEvent(input$selected_contrast, {
      if (has_initialized()) return()

      shiny::invalidateLater(100, session)
      shiny::isolate({
        if (is.null(input$selected_contrast)) return()

        # 获取可用集合
        task_temp <- tryCatch(
          extract_gsea_task(gsea_res, input$selected_contrast, "ALL"),
          error = function(e) NULL
        )
        if (is.null(task_temp)) return()

        df_temp <- as.data.frame(task_temp$gsea_res@result)
        available <- unique(c(df_temp$Collection, df_temp$Combo_Name))
        available <- setdiff(available, c(NA, "Unknown"))

        if (length(available) == 0) return()

        # 选择最小的集合
        counts <- sapply(available, function(x) {
          sum(startsWith(as.character(df_temp$Combo_Name), x))
        })
        min_col <- available[which.min(counts)[1]]

        shiny::updateSelectizeInput(
          session,
          "selected_collections",
          choices = c("ALL", sort(available)),
          selected = min_col
        )

        # 执行初始化
        auto_data <- process_data_core(
          input$selected_contrast,
          min_col,
          "abs_nes_desc",
          "logcpm",
          "3",
          "#E41A1C, #377EB8, #4DAF4A, #984EA3",
          TRUE
        )

        if (!is.null(auto_data)) {
          result_data(auto_data)
          has_initialized(TRUE)
          shiny::showNotification(
            sprintf("已自动加载: %s", min_col),
            type = "default",
            duration = 3
          )
        }
      })
    }, ignoreNULL = TRUE, ignoreInit = FALSE)

    # 手动确认按钮
    shiny::observeEvent(input$run_btn, {
      shiny::req(input$selected_contrast)

      cols <- input$selected_collections
      if (is.null(cols) || length(cols) == 0) cols <- "ALL"

      manual_data <- process_data_core(
        input$selected_contrast,
        cols,
        input$sort_by,
        input$expression_type,
        input$plot_subtype,
        input$custom_colors,
        FALSE
      )

      if (is.null(manual_data)) {
        shiny::showNotification("所选配置无数据", type = "error", duration = 5)
        return()
      }

      result_data(manual_data)
      shiny::showNotification("工作台已更新", type = "message", duration = 3)
    })

    return(shiny::reactive({ result_data() }))
  })
}
