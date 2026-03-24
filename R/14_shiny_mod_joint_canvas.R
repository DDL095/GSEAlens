#' @title 联合GSEA填充画布模块 UI（纯展示版）
#' @description 所有控制已移至侧边栏，此处仅显示画布结果
#' @keywords internal

mod_joint_canvas_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    # 🔧 修改：移除所有控制元素，仅保留画布展示（占满12列）
    shiny::div(class = "white-box", style = "min-height: 900px;",
               shiny::h4("GSEA 联合画布"),
               shiny::uiOutput(ns("canvas_info")),
               shiny::plotOutput(ns("canvas_plot"), height = "auto", width = "100%")
    )
  )
}

#' @title 联合GSEA填充画布模块 Server（侧边栏控制版）
#' @description
#'   1. 接收侧边栏控制参数（排列模式）
#'   2. 无最大行数限制
#'   3. 美学与multi_plot完全一致
#' @keywords internal

mod_joint_canvas_server <- function(id, gsea_res, data_prep_list, table_result) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # 存储当前画布结果
    canvas_result <- shiny::reactiveVal(NULL)

    # 🔧 核心：监听侧边栏的生成按钮（而非本地按钮）
    shiny::observeEvent(data_prep_list$joint_generate(), {

      # 🔧 从data_prep_list获取参数（而非input）
      contrasts <- data_prep_list$joint_contrasts()
      ncol_val <- data_prep_list$joint_ncol()

      shiny::req(contrasts, ncol_val)
      if (length(contrasts) == 0) {
        shiny::showNotification("请先在左侧控制栏选择对比组（排列模式）", type = "error")
        return()
      }

      # 🔧 从data_prep_list获取绘图参数（与multi_plot一致）
      data_list <- data_prep_list$data()
      shiny::req(data_list)

      plot_subtype <- data_list$plot_subtype
      custom_colors <- data_list$custom_colors

      # 🔧 通路：与multi_plot共用选择
      pathways <- table_result$selected_pathways()
      shiny::req(pathways)
      if (length(pathways) == 0) {
        shiny::showNotification("请先在主工作台勾选\"联合展示\"的通路", type = "error")
        return()
      }

      # 验证ncol
      if (is.null(ncol_val) || ncol_val < 1) ncol_val <- 2

      # 解析颜色
      colors <- trimws(strsplit(custom_colors, ",")[[1]])
      if (length(colors) < length(pathways)) {
        colors <- rep(colors, length.out = length(pathways))
      }

      shiny::incProgress(0.1, detail = "提取GSEA任务...")

      plot_list <- list()

      # 🔧 排列模式：严格按用户选择顺序（支持A_vs_B和B_vs_A同时存在）
      for (i in seq_along(contrasts)) {
        contrast_id <- contrasts[i]

        task_obj <- tryCatch({
          extract_gsea_task(gsea_res, contrast_id, "ALL")
        }, error = function(e) {
          message(sprintf("提取 %s 失败: %s", contrast_id, e$message))
          NULL
        })

        if (is.null(task_obj)) {
          p <- ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0.5, y = 0.5,
                              label = sprintf("未找到:\n%s", contrast_id),
                              size = 3, color = "red") +
            ggplot2::theme_void()
          plot_list[[i]] <- p
          next
        }

        # 🔧 美学统一：与multi_plot完全一致
        main_title <- sprintf(
          "%s [%d 通路]",
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
                              label = sprintf("错误:\n%s", substr(e$message, 1, 80)),
                              size = 3, color = "red") +
            ggplot2::theme_void()
        })

        plot_list[[i]] <- p

        if (i %% 2 == 0 || i == length(contrasts)) {
          shiny::incProgress(0.6 * i / length(contrasts),
                             detail = sprintf("绘制 %d/%d ...", i, length(contrasts)))
        }
      }

      if (length(plot_list) == 0) {
        shiny::showNotification("无有效图表可生成", type = "error")
        return()
      }

      n_plots <- length(plot_list)
      # 🔧 删除最大行数限制：自动计算
      actual_nrow <- ceiling(n_plots / ncol_val)

      # 🔧 美学优化：与multi_plot一致
      combined_plot <- patchwork::wrap_plots(
        plot_list,
        ncol = ncol_val,
        nrow = actual_nrow,
        byrow = TRUE,
        guides = "collect"
      ) + patchwork::plot_annotation(
        title = sprintf("联合GSEA画布: %d 对比组 × %d 通路", n_plots, length(pathways)),
        subtitle = sprintf("排列: %s", paste(contrasts, collapse = " → ")),
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

      shiny::incProgress(1.0, detail = "完成!")
    })

    # 渲染画布
    output$canvas_plot <- shiny::renderPlot({
      shiny::req(canvas_result())
      canvas_result()$plot
    }, height = function() {
      if (is.null(canvas_result())) return(900)
      nrow <- canvas_result()$nrow
      # 🔧 每行400px，无上限
      return(max(900, nrow * 400))
    }, width = 1600, res = 72)

    # 画布信息
    output$canvas_info <- shiny::renderUI({
      shiny::req(canvas_result())
      info <- canvas_result()

      shiny::tags$div(
        style = "margin-bottom: 15px; padding: 10px; background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; border-radius: 8px;",
        shiny::HTML(sprintf(
          "<strong>📊 画布信息:</strong> %d个对比组 | 布局: %d列 × %d行 | 共 %d 条通路<br>
           <small>排列顺序: %s</small>",
          info$n_plots, info$ncol, info$nrow, length(info$pathways),
          paste(info$contrasts, collapse = " → ")
        ))
      )
    })

    # 🔧 删除：保存画布功能已移除（根据用户要求）

  })
}
