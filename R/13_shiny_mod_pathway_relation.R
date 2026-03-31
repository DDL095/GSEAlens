#' @title Pathway Relationship Exploration Module UI
#' @description Network visualization module with dual-mode pathway selection
#' @keywords internal

mod_pathway_relation_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      # ===== 左侧控制面板 =====
      shiny::column(
        3,
        shiny::div(
          class = "well",
          style = "padding: 15px;",

          # ---- Mode 选择器 ----
          shiny::h4("Select Analysis Mode"),
          shiny::radioButtons(
            ns("network_mode"),
            label = NULL,
            choices = c(
              "Top N from Current Set" = "mode_topN",
              "Selected from Main Table" = "mode_select"
            ),
            selected = "mode_topN",
            width = "100%"
          ),

          shiny::hr(),

          # ---- 共享参数面板 ----
          shiny::h4("Shared Parameters"),

          # FDR 阈值
          shiny::numericInput(
            ns("fdr_threshold"),
            label = "FDR Threshold:",
            min = 0,
            max = 1.0,
            value = 0.25,
            step = 0.01
          ),

          # 最小共享基因数
          shiny::numericInput(
            ns("min_shared"),
            label = "Min Shared Core Genes:",
            value = 3,
            min = 1,
            max = 9999,
            step = 1
          ),

          # 最大节点数
          shiny::sliderInput(
            ns("max_nodes"),
            label = "Max Nodes:",
            min = 5,
            max = 100,
            value = 30,
            step = 5
          ),

          # 布局算法
          shiny::selectInput(
            ns("network_layout"),
            label = "Layout Algorithm:",
            choices = c(
              "Fruchterman-Reingold" = "fr",
              "Kamada-Kawai" = "kk",
              "Circle" = "circle"
            ),
            selected = "fr"
          ),

          # Seed 参数
          shiny::numericInput(
            ns("seed"),
            label = "Seed (for reproducibility):",
            value = 42,
            min = 1,
            max = 9999,
            step = 1
          ),

          shiny::helpText(
            style = "color: #666; font-size: 11px;",
            "Seed ensures reproducible layout when parameters change"
          ),

          shiny::hr(),

          # ---- TopN 专属参数 ----
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] == 'mode_topN'", ns("network_mode")),
            shiny::div(
              style = "background: #e3f2fd; padding: 10px; border-radius: 5px;",
              shiny::h5("Top N Configuration"),
              shiny::numericInput(
                ns("topN_count"),
                label = "Top N Count:",
                value = 20,
                min = 1,
                max = 9999,
                step = 1
              )
            )
          ),

          # ---- Select 模式提示 ----
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] == 'mode_select'", ns("network_mode")),
            shiny::div(
              style = "background: #f3e5f5; padding: 10px; border-radius: 5px;",
              shiny::h5("Main Table Selection"),
              shiny::helpText(
                style = "color: #6a1b9a; font-size: 11px;",
                "Check pathways in the main table 'Joint Plot' column"
              )
            )
          ),

          shiny::hr(),

          # ---- 通路预览列表 ----
          shiny::h4("Pathways to Plot"),
          shiny::uiOutput(ns("pathway_preview_list"))
        )
      ),

      # ===== 右侧面板 =====
      shiny::column(
        9,
        shiny::tabsetPanel(
          id = ns("active_tab"),
          type = "tabs",

          # ---- DotPlot Tab ----
          shiny::tabPanel(
            title = "DotPlot",
            value = "dotplot",
            shiny::div(
              class = "white-box",
              style = "padding: 15px; margin-top: 15px;",
              shiny::uiOutput(ns("dotplot_status")),
              shiny::selectInput(
                ns("dotplot_color_mode"),
                label = "Color by:",
                choices = c(
                  "-log10(FDR)" = "padj",
                  "-log10(P-value)" = "pval",
                  "NES" = "nes"
                ),
                selected = "padj"
              ),
              shiny::selectInput(
                ns("dotplot_size_mode"),
                label = "Size by:",
                choices = c(
                  "Core Genes Count" = "core_size",
                  "Set Size" = "setsize"
                ),
                selected = "core_size"
              ),
              plotly::plotlyOutput(ns("plot_dotplot"), height = "600px") %>%
                shinycssloaders::withSpinner(type = 6, color = "#28a745")
            )
          ),

          # ---- Network Tab ----
          shiny::tabPanel(
            title = "Network",
            value = "network",
            shiny::div(
              class = "white-box",
              style = "padding: 15px; margin-top: 15px;",
              shiny::uiOutput(ns("network_status")),
              plotly::plotlyOutput(ns("plot_network"), height = "700px") %>%
                shinycssloaders::withSpinner(type = 6, color = "#28a745")
            )
          )
        )
      )
    )
  )
}


#' @title Pathway Relationship Exploration Module Server
#' @description Dual-mode pathway network visualization with defensive rendering
#' @keywords internal

mod_pathway_relation_server <- function(id, data_prep_list, gsea_res, table_result = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ============================================================
    # 1. 模式状态管理
    # ============================================================

    network_mode <- shiny::reactiveVal("mode_topN")

    shiny::observeEvent(input$network_mode, {
      new_mode <- input$network_mode
      if (new_mode != network_mode()) {
        network_mode(new_mode)
        message(sprintf("[Network] Mode switched to: %s", new_mode))
      }
    })

    # ============================================================
    # 2. 数据源 reactive
    # ============================================================

    # 2a. TopN 数据源
    topN_candidates <- shiny::reactive({
      if (network_mode() != "mode_topN") return(character(0))

      data_list <- data_prep_list$data()
      shiny::req(data_list)

      df <- data_list$df
      shiny::req(nrow(df) > 0)

      top_n <- input$topN_count
      if (is.null(top_n)) top_n <- 50
      top_n <- max(1, min(top_n, nrow(df)))

      df[1:top_n, "ID"]
    })

    # 2b. Select 数据源
    select_candidates <- shiny::reactive({
      if (network_mode() != "mode_select") return(character(0))

      if (is.null(table_result) || is.null(table_result$selected_pathways)) {
        return(character(0))
      }

      selected <- table_result$selected_pathways()
      if (is.null(selected)) return(character(0))
      return(selected)
    })

    # 2c. 原始候选通路
    candidate_raw <- shiny::reactive({
      switch(network_mode(),
             "mode_topN" = topN_candidates(),
             "mode_select" = select_candidates(),
             character(0)
      )
    })

    # 2d. 应用 FDR 阈值筛选后的候选通路
    candidate_filtered <- shiny::reactive({
      pathways <- candidate_raw()

      if (length(pathways) == 0) {
        return(character(0))
      }

      fdr_thresh <- input$fdr_threshold
      if (is.null(fdr_thresh)) fdr_thresh <- 0.25

      data_list <- data_prep_list$data()
      shiny::req(data_list)

      df <- data_list$df
      fdr_vec <- df$p.adjust[match(pathways, df$ID)]
      names(fdr_vec) <- pathways

      filtered <- pathways[!is.na(fdr_vec) & fdr_vec < fdr_thresh]

      message(sprintf("[Network] FDR filter: %d/%d pathways pass (threshold=%.2f)",
                      length(filtered), length(pathways), fdr_thresh))

      return(filtered)
    })

    # ============================================================
    # 3. 最终通路列表管理
    # ============================================================

    final_pathways <- shiny::reactiveVal(character(0))

    shiny::observeEvent(candidate_filtered(), {
      new_candidates <- candidate_filtered()

      if (length(new_candidates) > 0) {
        max_n <- input$max_nodes
        if (is.null(max_n)) max_n <- 999

        final <- new_candidates[1:min(length(new_candidates), max_n)]
        final_pathways(final)

        message(sprintf("[Network] Updated pathways: %d (max_nodes=%d)",
                        length(final), max_n))
      } else {
        final_pathways(character(0))
      }
    })

    # ============================================================
    # 4. 通路预览列表 UI
    # ============================================================

    output$pathway_preview_list <- shiny::renderUI({
      pathways <- final_pathways()

      if (length(pathways) == 0) {
        return(shiny::div(
          style = "background: #fff3cd; padding: 10px; border-radius: 5px; color: #856404;",
          shiny::strong("No pathways available"),
          shiny::br(),
          shiny::small("Adjust parameters or select pathways in Main Table")
        ))
      }

      data_list <- data_prep_list$data()
      df <- data_list$df

      pathway_info <- lapply(pathways, function(pid) {
        row_idx <- which(df$ID == pid)
        if (length(row_idx) == 0) return(NULL)

        row <- df[row_idx[1], ]
        fdr <- row$p.adjust
        nes <- row$NES
        fdr_str <- if (!is.na(fdr) && is.numeric(fdr)) sprintf("%.2e", fdr) else "N/A"
        nes_str <- if (!is.na(nes) && is.numeric(nes)) sprintf("%.2f", nes) else "N/A"

        shiny::div(
          style = "background: #f8f9fa; padding: 8px; margin-bottom: 5px; border-radius: 4px; border-left: 3px solid #007bff;",
          shiny::div(
            style = "font-size: 12px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap;",
            shiny::strong(pid)
          ),
          shiny::div(
            style = "font-size: 11px; color: #666;",
            sprintf("FDR: %s | NES: %s", fdr_str, nes_str)
          )
        )
      })

      shiny::tagList(pathway_info)
    })

    # ============================================================
    # 5. DotPlot 状态提示
    # ============================================================

    output$dotplot_status <- shiny::renderUI({
      pathways <- final_pathways()

      status_parts <- c(
        sprintf("Mode: %s", ifelse(network_mode() == "mode_topN", "TopN", "Select")),
        sprintf("| Pathways: %d", length(pathways)),
        sprintf("| FDR < %.2f", input$fdr_threshold)
      )

      if (length(pathways) == 0) {
        div_style <- "background: #fff3cd; padding: 10px; border-radius: 5px; margin-bottom: 10px; color: #856404;"
      } else {
        div_style <- "background: #d4edda; padding: 10px; border-radius: 5px; margin-bottom: 10px; color: #155724;"
      }

      shiny::div(style = div_style, shiny::HTML(paste(status_parts, collapse = " ")))
    })

    # ============================================================
    # 6. Network 状态提示
    # ============================================================

    output$network_status <- shiny::renderUI({
      pathways <- final_pathways()

      status_parts <- c(
        sprintf("Mode: %s", ifelse(network_mode() == "mode_topN", "TopN", "Select")),
        sprintf("| Pathways: %d", length(pathways)),
        sprintf("| FDR < %.2f", input$fdr_threshold),
        sprintf("| min_shared: %d", input$min_shared),
        sprintf("| seed: %d", input$seed)
      )

      if (length(pathways) == 0) {
        div_style <- "background: #fff3cd; padding: 10px; border-radius: 5px; margin-bottom: 10px; color: #856404;"
      } else {
        div_style <- "background: #d4edda; padding: 10px; border-radius: 5px; margin-bottom: 10px; color: #155724;"
      }

      shiny::div(style = div_style, shiny::HTML(paste(status_parts, collapse = " ")))
    })

    # ============================================================
    # 7. DotPlot 绘图函数
    # ============================================================

    output$plot_dotplot <- plotly::renderPlotly({

      pathways <- final_pathways()

      if (length(pathways) == 0) {
        return(plotly::plot_ly() %>%
                 plotly::layout(
                   title = list(
                     text = "No pathways to display<br><sup>Adjust parameters or select pathways</sup>",
                     font = list(size = 14),
                     x = 0.5
                   ),
                   xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                   yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE)
                 ))
      }

      data_list <- data_prep_list$data()
      shiny::req(data_list)

      df <- data_list$df
      plot_df <- df[df$ID %in% pathways, ]

      if (nrow(plot_df) == 0) {
        return(plotly::plot_ly() %>%
                 plotly::layout(
                   title = list(text = "No matching pathways found", font = list(size = 14)),
                   xaxis = list(showgrid = FALSE),
                   yaxis = list(showgrid = FALSE)
                 ))
      }

      # 获取核心基因数量
      task <- list(
        gsea_res = data_list$gsea_res,
        meta = list(
          left_group = data_list$left_group,
          right_group = data_list$right_group,
          contrast_id = data_list$contrast_id
        )
      )
      class(task) <- "GseaTask"

      core_list <- tryCatch({
        get_core_genes_list(task, plot_df$ID)
      }, error = function(e) NULL)

      plot_df$CoreCount <- sapply(plot_df$ID, function(pid) {
        if (is.null(core_list) || is.null(core_list[[pid]])) return(0)
        length(core_list[[pid]])
      })

      # 颜色映射
      color_mode <- input$dotplot_color_mode
      color_vals <- switch(color_mode,
                           "padj" = -log10(plot_df$p.adjust),
                           "pval" = -log10(plot_df$pvalue),
                           "nes" = abs(plot_df$NES),
                           -log10(plot_df$p.adjust))

      color_title <- switch(color_mode,
                           "padj" = "-log10(FDR)",
                           "pval" = "-log10(P-value)",
                           "nes" = "|NES|")

      # 大小映射
      size_mode <- input$dotplot_size_mode
      size_vals <- switch(size_mode,
                          "core_size" = plot_df$CoreCount,
                          "setsize" = plot_df$setSize,
                          plot_df$CoreCount)

      size_range <- c(5, 25)
      if (max(size_vals) > min(size_vals)) {
        size_scaled <- size_range[1] + (size_vals - min(size_vals)) /
          (max(size_vals) - min(size_vals)) * (size_range[2] - size_range[1])
      } else {
        size_scaled <- mean(size_range)
      }

      # 按 NES 排序
      plot_df <- plot_df[order(plot_df$NES, decreasing = TRUE), ]
      size_scaled <- size_scaled[order(plot_df$NES, decreasing = TRUE)]
      color_vals <- color_vals[order(plot_df$NES, decreasing = TRUE)]

      # 构建 hover 文本
      hover_text <- sprintf(
        "<b>%s</b><br>FDR: %.2e<br>P-value: %.2e<br>NES: %.2f<br>Core Genes: %d<br>Set Size: %d",
        plot_df$ID,
        plot_df$p.adjust,
        plot_df$pvalue,
        plot_df$NES,
        plot_df$CoreCount,
        plot_df$setSize
      )

      # 绘图
      plotly::plot_ly(
        data = plot_df,
        x = plot_df$NES,
        y = seq_len(nrow(plot_df)),
        type = "scatter",
        mode = "markers",
        marker = list(
          size = size_scaled,
          color = color_vals,
          colorscale = "RdYlBu_r",
          showscale = TRUE,
          colorbar = list(title = list(text = color_title, font = list(size = 12))),
          line = list(color = "black", width = 1)
        ),
        text = hover_text,
        hoverinfo = "text",
        hovertemplate = "%{text}<extra></extra>"
      ) %>% plotly::layout(
        title = list(
          text = sprintf("Pathway DotPlot: %s vs %s<br><sub>%d pathways | Size: Core Genes | Color: %s</sub>",
                         data_list$left_group, data_list$right_group,
                         nrow(plot_df), color_title),
          font = list(size = 12),
          x = 0.5,
          xanchor = "center"
        ),
        xaxis = list(title = "NES", zeroline = TRUE),
        yaxis = list(
          title = "",
          tickmode = "array",
          tickvals = seq_len(nrow(plot_df)),
          ticktext = plot_df$ID,
          tickfont = list(size = 9)
        ),
        margin = list(l = 250, r = 50, t = 80, b = 50),
        showlegend = FALSE
      ) %>% plotly::config(
        displayModeBar = TRUE,
        displaylogo = FALSE
      )
    })

    # ============================================================
    # 8. Network 绘图函数（5 阶段防御性渲染）
    # ============================================================

    # 创建全局变量用于边点击事件
    edge_list <<- NULL
    node_df <<- NULL
    core_list <<- NULL

    output$plot_network <- plotly::renderPlotly({

      # ---- Phase 1: 数据获取与验证 ----

      pathways <- final_pathways()

      if (length(pathways) == 0) {
        return(plotly::plot_ly() %>%
                 plotly::layout(
                   title = list(
                     text = "No pathways to display<br><sup>Adjust parameters or select pathways</sup>",
                     font = list(size = 14),
                     x = 0.5
                   ),
                   xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                   yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE)
                 ))
      }

      data_list <- data_prep_list$data()
      shiny::req(data_list)

      task <- list(
        gsea_res = data_list$gsea_res,
        meta = list(
          left_group = data_list$left_group,
          right_group = data_list$right_group,
          contrast_id = data_list$contrast_id
        )
      )
      class(task) <- "GseaTask"

      # ---- Phase 2: 获取核心基因列表 ----

      core_list <<- tryCatch({
        get_core_genes_list(task, pathways)
      }, error = function(e) {
        message(sprintf("[Network] Error extracting core genes: %s", e$message))
        NULL
      })

      if (is.null(core_list) || length(core_list) == 0) {
        return(plotly::plot_ly() %>%
                 plotly::layout(
                   title = list(
                     text = "Selected pathways have no core genes<br><sup>Cannot build network</sup>",
                     font = list(size = 14)
                   ),
                   xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                   yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE)
                 ))
      }

      valid_pathways <- names(core_list)[sapply(core_list, function(x) length(x) > 0)]

      if (length(valid_pathways) == 0) {
        return(plotly::plot_ly() %>%
                 plotly::layout(
                   title = list(
                     text = "All selected pathways have empty core genes",
                     font = list(size = 14)
                   ),
                   xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                   yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE)
                 ))
      }

      # ---- Phase 3: 构建边列表 ----

      min_shared <- input$min_shared
      if (is.null(min_shared)) min_shared <- 3

      edge_list <<- tryCatch({
        build_edge_list_safely(core_list[valid_pathways], min_shared_genes = min_shared)
      }, error = function(e) {
        message(sprintf("[Network] Edge building failed: %s", e$message))
        NULL
      })

      if (is.null(edge_list) || nrow(edge_list) == 0) {
        return(plotly::plot_ly() %>%
                 plotly::layout(
                   title = list(
                     text = sprintf("No connected pathways found<br><sup>min_shared=%d. Try reducing the threshold.</sup>", min_shared),
                     font = list(size = 14),
                     x = 0.5
                   ),
                   xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                   yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE)
                 ))
      }

      # ---- Phase 4: 构建 igraph 对象 ----

      g <- tryCatch({
        igraph::graph_from_data_frame(
          edge_list,
          directed = FALSE,
          vertices = valid_pathways
        )
      }, error = function(e) {
        message(sprintf("[Network] Graph construction failed: %s", e$message))
        NULL
      })

      if (is.null(g)) {
        return(plotly::plot_ly() %>%
                 plotly::layout(
                   title = list(
                     text = "Network graph construction failed",
                     font = list(size = 14)
                   ),
                   xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                   yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE)
                 ))
      }

      # ---- Phase 5: 布局计算 ----

      layout_algo <- input$network_layout
      seed_val <- input$seed

      layout_coords <- tryCatch({
        switch(layout_algo,
               "fr" = {
                 set.seed(seed_val)
                 igraph::layout_with_fr(g)
               },
               "kk" = igraph::layout_with_kk(g),
               "circle" = igraph::layout_in_circle(g),
               {
                 set.seed(seed_val)
                 igraph::layout_with_fr(g)
               }
        )
      }, error = function(e) {
        message(sprintf("[Network] Layout calculation failed: %s", e$message))
        NULL
      })

      if (is.null(layout_coords)) {
        return(plotly::plot_ly() %>%
                 plotly::layout(
                   title = list(
                     text = "Layout calculation failed",
                     font = list(size = 14)
                   ),
                   xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                   yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE)
                 ))
      }

      # ---- Phase 6: Plotly 绘制 ----

      node_df <<- data.frame(
        name = igraph::V(g)$name,
        x = layout_coords[, 1],
        y = layout_coords[, 2],
        stringsAsFactors = FALSE
      )

      res_df <- as.data.frame(task$gsea_res@result)
      node_info <- res_df[match(node_df$name, res_df$ID), ]
      node_df$NES <- node_info$NES
      node_df$FDR <- node_info$p.adjust
      node_df$CoreCount <- sapply(node_df$name, function(n) {
        if (is.null(core_list[[n]])) return(0)
        length(core_list[[n]])
      })
      node_df$color_val <- ifelse(is.na(node_df$NES), 0, node_df$NES)

      # 构建节点 hover 文本
      node_hover_text <- sprintf(
        "<b>%s</b><br>NES: %.2f<br>FDR: %.2e<br>Core Genes: %d",
        node_df$name,
        node_df$NES,
        node_df$FDR,
        node_df$CoreCount
      )

      # 构建边 hover 文本
      edge_hover_text_list <- lapply(seq_len(nrow(edge_list)), function(i) {
        sprintf(
          "<b>%s &harr; %s</b><br>Shared: %d genes<br>Jaccard: %.3f<br><i style='color:#666;'>Click edge for details</i>",
          edge_list$from[i], edge_list$to[i],
          edge_list$shared[i], edge_list$weight[i]
        )
      })

      title_text <- sprintf(
        "Pathway Network: %s vs %s<br><sub>%d nodes, %d edges | min_shared=%d | seed=%d</sub>",
        task$meta$left_group, task$meta$right_group,
        igraph::vcount(g), igraph::ecount(g),
        min_shared, seed_val
      )

      # ============================================================
      # 使用 add_trace 构建图形（更好的 hover 控制）
      # ============================================================

      # 初始化 plotly 对象
      p <- plotly::plot_ly(source = ns("network_plot"))

      # ---- 添加边 (trace 0 到 n-1) ----
      for (i in seq_len(nrow(edge_list))) {
        from_node <- node_df[node_df$name == edge_list$from[i], ]
        to_node <- node_df[node_df$name == edge_list$to[i], ]

        p <- p %>% plotly::add_trace(
          type = "scatter",
          mode = "lines",
          x = c(from_node$x, to_node$x, NA),
          y = c(from_node$y, to_node$y, NA),
          line = list(
            color = "rgba(100, 100, 100, 0.6)",
            width = 3 + edge_list$weight[i] * 8
          ),
          hoverinfo = "text",
          text = edge_hover_text_list[[i]],
          showlegend = FALSE,
          inherit = TRUE
        )
      }

      # ---- 添加节点 (trace n) ----
      p <- p %>% plotly::add_trace(
        type = "scatter",
        mode = "markers",
        x = node_df$x,
        y = node_df$y,
        marker = list(
          size = 18,
          color = node_df$color_val,
          colorscale = list(
            list(0, "#377EB8"),
            list(0.5, "white"),
            list(1, "#E41A1C")
          ),
          cauto = FALSE,
          cmin = -3,
          cmax = 3,
          showscale = TRUE,
          colorbar = list(
            title = list(text = "NES", font = list(size = 12)),
            tickfont = list(size = 10),
            len = 0.5,
            y = 0.5
          ),
          line = list(color = "black", width = 2)
        ),
        text = node_df$name,
        hovertemplate = paste(node_hover_text, "<extra></extra>"),
        showlegend = FALSE,
        key = node_df$name,
        inherit = TRUE
      )

      # ---- 添加节点标签 (trace n+1) ----
      p <- p %>% plotly::add_trace(
        type = "scatter",
        mode = "text",
        x = node_df$x,
        y = node_df$y + 0.12,
        text = node_df$name,
        textposition = "top center",
        textfont = list(size = 10, color = "#333"),
        hoverinfo = "skip",
        showlegend = FALSE,
        inherit = TRUE
      )

      # ---- 布局设置 ----
      x_range <- c(min(layout_coords[, 1]) - 0.5, max(layout_coords[, 1]) + 0.5)
      y_range <- c(min(layout_coords[, 2]) - 0.5, max(layout_coords[, 2]) + 0.5)

      p <- p %>% plotly::layout(
        title = list(
          text = title_text,
          font = list(size = 12),
          x = 0.5,
          xanchor = "center"
        ),
        xaxis = list(
          title = "",
          showgrid = FALSE,
          showticklabels = FALSE,
          zeroline = FALSE,
          range = x_range
        ),
        yaxis = list(
          title = "",
          showgrid = FALSE,
          showticklabels = FALSE,
          zeroline = FALSE,
          range = y_range,
          scaleanchor = "x",
          scaleratio = 1
        ),
        hovermode = "closest",
        dragmode = "pan",
        showlegend = FALSE,
        margin = list(l = 50, r = 50, t = 80, b = 50)
      ) %>% plotly::config(
        displayModeBar = TRUE,
        displaylogo = FALSE,
        modeBarButtonsToRemove = c("lasso2d", "select2d")
      )


    })

    # ============================================================
    # 9. 边点击事件监听
    # ============================================================

    shiny::observeEvent(plotly::event_data("plotly_click", source = ns("network_plot")), {
      click_data <- plotly::event_data("plotly_click", source = ns("network_plot"))

      if (is.null(click_data) || is.null(click_data$curveNumber)) {
        return()
      }

      n_edges <- nrow(edge_list)

      # curveNumber: 0 到 n-1 = 边, n = 节点, n+1 = 标签
      if (click_data$curveNumber < n_edges) {
        edge_idx <- click_data$curveNumber + 1

        from_pw <- as.character(edge_list$from[edge_idx])
        to_pw <- as.character(edge_list$to[edge_idx])
        shared_count <- edge_list$shared[edge_idx]
        jaccard <- edge_list$weight[edge_idx]
        shared_genes <- unlist(edge_list$shared_genes[[edge_idx]])

        pathway_a_genes <- core_list[[from_pw]]
        pathway_b_genes <- core_list[[to_pw]]

        session$sendCustomMessage(
          type = "network_edge_clicked",
          message = list(
            from = from_pw,
            to = to_pw,
            shared = shared_count,
            jaccard = jaccard,
            shared_genes = shared_genes,
            pathway_a_genes = pathway_a_genes,
            pathway_b_genes = pathway_b_genes
          )
        )

        message(sprintf("[Network] Edge clicked: %s <-> %s (shared=%d)",
                        from_pw, to_pw, shared_count))
      }
    })

    # ============================================================
    # 10. 返回值
    # ============================================================

    return(list(
      final_pathways = final_pathways,
      network_mode = network_mode
    ))

  })
}
