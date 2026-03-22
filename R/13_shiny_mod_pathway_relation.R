#' @title 通路关系探索模块 UI
#' @description 包含DotPlot、Chord、Network、UpSet四个子Tab
#' @keywords internal
mod_pathway_relation_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      # 左侧控制面板
      shiny::column(3,
                    shiny::div(class = "well",
                               shiny::h4("🎨 可视化控制"),

                               # DotPlot专用控制
                               shiny::conditionalPanel(
                                 condition = sprintf("input['%s'] == 'dotplot'", ns("active_tab")),
                                 shiny::selectInput(
                                   ns("ratio_source"),
                                   "Ratio计算方式:",
                                   choices = c("ORA (交集/通路基因)" = "ora",
                                               "Leading Edge (交集/DE基因)" = "leading"),
                                   selected = "ora"
                                 ),          # DotPlot专用控制 - 优化版
                                 shiny::conditionalPanel(
                                   condition = sprintf("input['%s'] == 'dotplot'", ns("active_tab")),
                                   shiny::selectInput(
                                     ns("ratio_source"),
                                     "Ratio计算方式:",
                                     choices = c("ORA (交集/通路基因)" = "ora",
                                                 "Leading Edge (交集/DE基因)" = "leading"),
                                     selected = "ora"
                                   ),
                                   shiny::selectInput(
                                     ns("stat_color_mode"),
                                     "颜色映射:",
                                     choices = c("-log10(P-value)" = "pval",
                                                 "-log10(FDR)" = "padj",
                                                 "NES" = "nes"),
                                     selected = "padj"
                                   ),
                                   shiny::selectInput(
                                     ns("size_mode"),
                                     "气泡大小:",
                                     choices = c("Core Genes数" = "core_size",
                                                 "Set Size" = "setsize",
                                                 "Ratio值" = "ratio"),
                                     selected = "core_size"
                                   ),
                                   # 新增：大小范围控制
                                   shiny::sliderInput(
                                     ns("size_range"),
                                     "气泡大小范围:",
                                     min = 1, max = 20, value = c(5, 15), step = 1
                                   ),
                                   # 新增：颜色范围控制（截断极端值）
                                   shiny::sliderInput(
                                     ns("color_cap"),
                                     "颜色值上限 (-log10):",
                                     min = 5, max = 50, value = 20, step = 1
                                   ),
                                   shiny::sliderInput(
                                     ns("dot_alpha"),
                                     "透明度:",
                                     min = 0.1, max = 1, value = 0.8, step = 0.1
                                   ),
                                   # 新增：数据来源提示
                                   shiny::helpText(
                                     style = "color: #666; font-size: 11px;",
                                     "数据来源：主工作台当前选中的通路集合"
                                   )
                                 ),
                                 shiny::selectInput(
                                   ns("stat_color_mode"),
                                   "颜色映射:",
                                   choices = c("-log10(P-value)" = "pval",
                                               "-log10(FDR)" = "padj",
                                               "NES" = "nes"),
                                   selected = "padj"
                                 ),
                                 shiny::selectInput(
                                   ns("size_mode"),
                                   "气泡大小:",
                                   choices = c("Core Genes数" = "core_size",
                                               "Set Size" = "setsize",
                                               "Ratio值" = "ratio"),
                                   selected = "core_size"
                                 ),
                                 shiny::sliderInput(
                                   ns("dot_alpha"),
                                   "透明度:",
                                   min = 0.1, max = 1, value = 0.8, step = 0.1
                                 )
                               ),

                               # Network/UpSet/Chord共享控制
                               shiny::conditionalPanel(
                                 condition = sprintf("input['%s'] != 'dotplot'", ns("active_tab")),
                                 shiny::selectInput(
                                   ns("set_definition"),
                                   "基因集合定义:",
                                   choices = c("Core Genes (Leading Edge)" = "core"),
                                   selected = "core"
                                 ),
                                 shiny::sliderInput(
                                   ns("min_shared"),
                                   "最小共享基因数:",
                                   min = 1, max = 20, value = 3, step = 1
                                 )
                               ),

                               # Network专用
                               shiny::conditionalPanel(
                                 condition = sprintf("input['%s'] == 'network'", ns("active_tab")),
                                 shiny::sliderInput(
                                   ns("max_nodes"),
                                   "最大节点数:",
                                   min = 5, max = 50, value = 20, step = 1
                                 ),
                                 shiny::selectInput(
                                   ns("network_layout"),
                                   "布局算法:",
                                   choices = c("Fruchterman-Reingold" = "fr",
                                               "Kamada-Kawai" = "kk",
                                               "Circular" = "circle"),
                                   selected = "fr"
                                 )
                               ),

                               shiny::hr(),
                               shiny::actionButton(
                                 ns("refresh_plot"),
                                 "🔄 更新图表",
                                 class = "btn-primary",
                                 style = "width: 100%;"
                               ),

                               shiny::helpText(
                                 style = "margin-top: 10px; color: #666;",
                                 "基于主工作台选中的通路集合进行绘制"
                               )
                    )
      ),

      # 右侧绘图区
      shiny::column(9,
                    shiny::tabsetPanel(
                      id = ns("active_tab"),
                      type = "tabs",

                      # Sub 1: DotPlot
                      shiny::tabPanel(
                        title = shiny::HTML("🔴 DotPlot"),
                        value = "dotplot",
                        shiny::div(class = "white-box", style = "margin-top: 15px;",
                                   plotly::plotlyOutput(ns("plot_dotplot"), height = "600px")
                        )
                      ),

                      # Sub 3: Network
                      shiny::tabPanel(
                        title = shiny::HTML("🕸️ Network"),
                        value = "network",
                        shiny::div(class = "white-box", style = "margin-top: 15px;",
                                   plotly::plotlyOutput(ns("plot_network"), height = "600px")
                        )
                      )

                    )
      )
    )
  )
}

#' @title 通路关系探索模块 Server
#' @description 实现DotPlot、Chord、Network、UpSet逻辑，仅hover无click
#' @keywords internal
mod_pathway_relation_server <- function(id, data_prep_list, gsea_res) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # 获取选中的通路（来自主工作台）
    selected_pathways <- shiny::reactive({
      data_list <- data_prep_list$data()
      if (is.null(data_list)) return(character(0))

      # 如果有选中的通路（从主表格传递），使用之；否则使用全部显著通路
      # 这里假设data_list中有一个selected字段，或从上游模块传入
      # 实际实现中应通过shared reactive连接主表格的选择
      pathways <- data_list$df$ID[data_list$df$p.adjust < 0.05]
      if (length(pathways) > 50) pathways <- pathways[1:50]  # 安全限制
      return(pathways)
    })

    # 获取当前对比的GseaTask对象
    current_task <- shiny::reactive({
      data_list <- data_prep_list$data()
      shiny::req(data_list)

      # 重建GseaTask对象（简化版）
      task <- list(
        gsea_res = data_list$gsea_res,
        meta = list(
          left_group = data_list$left_group,
          right_group = data_list$right_group,
          contrast_id = data_list$contrast_id
        )
      )
      class(task) <- "GseaTask"
      return(task)
    })

    # 1. DotPlot实现（经典GSEA气泡图）- 带对数比例尺
    output$plot_dotplot <- plotly::renderPlotly({
      shiny::req(selected_pathways(), current_task())
      shiny::req(input$refresh_plot)

      pathways <- selected_pathways()
      task <- current_task()
      gsea_obj <- task$gsea_res
      res_df <- as.data.frame(gsea_obj@result)
      res_df <- res_df[res_df$ID %in% pathways, ]

      if (nrow(res_df) == 0) {
        return(plotly::plot_ly() %>% plotly::layout(title = "无可用通路数据"))
      }

      # 获取DE基因（用于ORA计算）
      de_table <- tryCatch({
        get_de_table(gsea_res, task$meta$contrast_id)
      }, error = function(e) NULL)

      de_genes_sig <- character(0)
      if (!is.null(de_table) && "pvalue" %in% colnames(de_table)) {
        de_genes_sig <- de_table$gene_symbol[de_table$pvalue < 0.05]
      }

      # 计算Ratio和Core Genes
      plot_data <- lapply(1:nrow(res_df), function(i) {
        row <- res_df[i, ]
        pid <- row$ID

        # Core genes
        core_genes <- get_core_genes_for_pathway(task, pid)

        # Term genes（完整基因集）
        term_genes <- get_term_genes(gsea_res, pid)

        # 计算Ratio
        if (input$ratio_source == "ora") {
          # ORA模式: Core Genes / Term Genes (通路基因)
          if (length(term_genes) > 0) {
            ratio <- length(core_genes) / length(term_genes)
          } else {
            ratio <- 0
          }
        } else {
          # Leading Edge模式: (Core Genes ∩ DE Genes) / DE Genes
          if (length(de_genes_sig) > 0) {
            overlap <- intersect(core_genes, toupper(de_genes_sig))
            ratio <- length(overlap) / length(de_genes_sig)
          } else {
            ratio <- 0
          }
        }

        # Size映射 - 使用对数变换
        raw_size <- switch(input$size_mode,
                           "core_size" = length(core_genes),
                           "setsize" = row$setSize,
                           "ratio" = ratio * 100,
                           length(core_genes))

        # 对数变换用于气泡大小（避免极端值导致比例失调）
        log_size <- log10(raw_size + 1)  # +1避免log(0)

        # Color映射
        raw_color <- switch(input$stat_color_mode,
                            "pval" = -log10(row$pvalue),
                            "padj" = -log10(row$p.adjust),
                            "nes" = abs(row$NES),
                            -log10(row$p.adjust))

        # 颜色截断（避免极端值）
        color_cap <- input$color_cap %||% 20
        color_val <- min(raw_color, color_cap)

        data.frame(
          Pathway = pid,
          Description = ifelse(is.na(row$Description), pid, row$Description),
          Ratio = ratio,
          RawSize = raw_size,      # 原始大小用于显示
          LogSize = log_size,      # 对数变换后的大小用于绘图
          Size = raw_size,         # 保持兼容性
          Color = color_val,
          RawColor = raw_color,    # 原始颜色值
          NES = row$NES,
          PValue = row$pvalue,
          FDR = row$p.adjust,
          CoreCount = length(core_genes),
          stringsAsFactors = FALSE
        )
      })

      plot_df <- do.call(rbind, plot_data)

      # 根据NES排序（正向在上）
      plot_df$Pathway <- factor(plot_df$Pathway,
                                levels = plot_df$Pathway[order(plot_df$NES, decreasing = TRUE)])

      # 颜色标题
      color_title <- switch(input$stat_color_mode,
                            "pval" = "-log10(P-value)",
                            "padj" = "-log10(FDR)",
                            "nes" = "|NES|")

      # 计算气泡大小范围（对数尺度）
      size_range <- input$size_range %||% c(5, 20)
      min_log_size <- min(plot_df$LogSize)
      max_log_size <- max(plot_df$LogSize)

      # 线性映射到像素大小
      if (max_log_size > min_log_size) {
        plot_df$MarkerSize <- size_range[1] +
          (plot_df$LogSize - min_log_size) / (max_log_size - min_log_size) *
          (size_range[2] - size_range[1])
      } else {
        plot_df$MarkerSize <- mean(size_range)
      }

      # 构建标题
      contrast_id <- task$meta$contrast_id
      left_group <- task$meta$left_group
      right_group <- task$meta$right_group

      # 计算比例尺刻度（以10的倍数）
      max_raw_size <- max(plot_df$RawSize)
      # 找到比最大值大的最近的10的幂
      scale_max <- 10^ceiling(log10(max_raw_size + 1))
      # 生成比例尺刻度：1, 10, 100, 1000... 或 1, 5, 10, 50...
      if (scale_max > 100) {
        scale_breaks <- c(1, 10, 100, 1000)
        scale_breaks <- scale_breaks[scale_breaks <= scale_max * 10]
      } else {
        scale_breaks <- c(1, 5, 10, 50, 100)
        scale_breaks <- scale_breaks[scale_breaks <= scale_max]
      }

      # 计算每个刻度对应的像素大小
      scale_sizes <- sapply(scale_breaks, function(s) {
        log_s <- log10(s + 1)
        if (max_log_size > min_log_size) {
          size_range[1] + (log_s - min_log_size) / (max_log_size - min_log_size) *
            (size_range[2] - size_range[1])
        } else {
          mean(size_range)
        }
      })

      # 绘制（仅hover，无click事件监听）
      p <- plotly::plot_ly(
        data = plot_df,
        x = ~Ratio,
        y = ~Pathway,
        type = "scatter",
        mode = "markers",
        marker = list(
          size = ~MarkerSize,
          color = ~Color,
          colorscale = "RdBu",
          showscale = TRUE,
          colorbar = list(
            title = list(text = color_title, font = list(size = 12)),
            tickfont = list(size = 10)
          ),
          opacity = input$dot_alpha %||% 0.8,
          line = list(color = "black", width = 1)
        ),
        text = ~sprintf(
          "<b>%s</b><br>Ratio: %.3f<br>Core Genes: %d<br>NES: %.2f<br>FDR: %.2e",
          Description, Ratio, CoreCount, NES, FDR
        ),
        hoverinfo = "text"
      ) %>% plotly::layout(
        xaxis = list(
          title = list(
            text = ifelse(input$ratio_source == "ora",
                          "Overlap Ratio (Core/Term)",
                          "Overlap Ratio (Core/DE)"),
            font = list(size = 12)
          ),
          tickfont = list(size = 10)
        ),
        yaxis = list(
          title = "",
          tickfont = list(size = 9)
        ),
        showlegend = FALSE,
        dragmode = FALSE,  # 禁止拖拽，确保无click联动
        title = list(
          text = sprintf("DotPlot: %s vs %s (%d pathways)", left_group, right_group, nrow(plot_df)),
          font = list(size = 14),
          x = 0.5,
          xanchor = "center"
        ),
        margin = list(l = 200),  # 给y轴标签留空间
        # 添加自定义比例尺注释
        annotations = list(
          list(
            x = 0.98,
            y = 0.02,
            xref = "paper",
            yref = "paper",
            text = "<b>Bubble Size Scale</b><br>(Core Genes)",
            showarrow = FALSE,
            font = list(size = 10),
            align = "right",
            bgcolor = "rgba(255,255,255,0.8)",
            bordercolor = "gray",
            borderwidth = 1
          )
        )
      )

      # 添加比例尺气泡（在右下角）
      for (i in seq_along(scale_breaks)) {
        p <- p %>% plotly::add_trace(
          x = max(plot_df$Ratio) * 1.05,  # 放在图右侧
          y = nrow(plot_df) * (0.1 + i * 0.08),  # 垂直排列
          mode = "markers+text",
          marker = list(
            size = scale_sizes[i],
            color = "lightgray",
            line = list(color = "black", width = 1)
          ),
          text = paste0(" ", scale_breaks[i]),
          textposition = "middle right",
          textfont = list(size = 9),
          hoverinfo = "skip",
          showlegend = FALSE
        )
      }

      return(p)
    }) %>% shiny::bindEvent(input$refresh_plot, ignoreNULL = FALSE)


    # 3. Network（通路关系网络，基于core_genes Jaccard相似度）- Plotly交互版
    output$plot_network <- plotly::renderPlotly({
      shiny::req(selected_pathways(), current_task())
      shiny::req(input$refresh_plot)

      pathways <- selected_pathways()
      task <- current_task()

      # 限制节点数
      max_nodes <- validate_param(input$max_nodes, 20, 2, 50, "max_nodes")
      if (length(pathways) > max_nodes) {
        pathways <- pathways[1:max_nodes]
      }

      # 计算Core Genes的Jaccard相似度矩阵
      core_list <- get_core_genes_list(task, pathways)

      # 构建边列表
      edge_list <- data.frame(
        from = character(0),
        to = character(0),
        weight = numeric(0),
        shared = integer(0),
        stringsAsFactors = FALSE
      )

      n <- length(pathways)
      for (i in 1:(n-1)) {
        for (j in (i+1):n) {
          genes_i <- core_list[[pathways[i]]]
          genes_j <- core_list[[pathways[j]]]

          if (length(genes_i) == 0 || length(genes_j) == 0) next

          intersection <- length(intersect(genes_i, genes_j))
          if (intersection < input$min_shared) next

          union <- length(union(genes_i, genes_j))
          jaccard <- if (union > 0) intersection / union else 0

          edge_list <- rbind(edge_list, data.frame(
            from = pathways[i],
            to = pathways[j],
            weight = jaccard,
            shared = intersection,
            stringsAsFactors = FALSE
          ))
        }
      }

      if (nrow(edge_list) == 0) {
        return(plotly::plot_ly() %>% plotly::layout(
          title = list(
            text = "无足够连接的通路（请降低最小共享基因数）",
            font = list(size = 14)
          )
        ))
      }

      # 使用igraph构建网络
      g <- igraph::graph_from_data_frame(edge_list, directed = FALSE, vertices = pathways)

      # 布局
      layout_func <- switch(input$network_layout,
                            "fr" = igraph::layout_with_fr,
                            "kk" = igraph::layout_with_kk,
                            "circle" = igraph::layout_in_circle,
                            igraph::layout_with_fr)

      layout_coords <- layout_func(g)

      # 构建节点数据框
      node_df <- data.frame(
        name = igraph::V(g)$name,
        x = layout_coords[, 1],
        y = layout_coords[, 2],
        stringsAsFactors = FALSE
      )

      # 获取节点属性（NES方向）
      res_df <- as.data.frame(task$gsea_res@result)
      node_info <- res_df[match(node_df$name, res_df$ID), ]
      node_df$NES <- node_info$NES
      node_df$color_val <- ifelse(is.na(node_df$NES), 0, node_df$NES)

      # 提取对比组信息用于标题
      contrast_id <- task$meta$contrast_id
      left_group <- task$meta$left_group
      right_group <- task$meta$right_group

      # 创建plotly图形 - 使用add_trace而不是scatter
      p <- plotly::plot_ly()

      # 添加边（使用add_segments）
      for (i in 1:nrow(edge_list)) {
        from_node <- node_df[node_df$name == edge_list$from[i], ]
        to_node <- node_df[node_df$name == edge_list$to[i], ]

        p <- p %>% plotly::add_segments(
          x = from_node$x,
          y = from_node$y,
          xend = to_node$x,
          yend = to_node$y,
          line = list(
            color = "rgba(128, 128, 128, 0.6)",
            width = 1 + edge_list$weight[i] * 3
          ),
          hoverinfo = "text",
          text = sprintf("Shared: %d genes<br>Jaccard: %.3f",
                         edge_list$shared[i], edge_list$weight[i]),
          showlegend = FALSE,
          inherit = FALSE
        )
      }

      # 添加节点
      p <- p %>% plotly::add_markers(
        data = node_df,
        x = ~x,
        y = ~y,
        marker = list(
          size = 15,
          color = ~color_val,
          colorscale = list(
            c(0, "#377EB8"),    # 负值 - 蓝色
            c(0.5, "white"),    # 中点 - 白色
            c(1, "#E41A1C")     # 正值 - 红色
          ),
          showscale = TRUE,
          colorbar = list(
            title = list(text = "NES", font = list(size = 12)),
            tickfont = list(size = 10)
          ),
          line = list(color = "black", width = 2)
        ),
        text = ~name,
        hoverinfo = "text",
        hovertext = ~sprintf(
          "<b>%s</b><br>NES: %.2f<br>Core Genes: %d",
          name,
          ifelse(is.na(NES), 0, NES),
          sapply(name, function(n) length(core_list[[n]]))
        ),
        showlegend = FALSE
      )

      # 添加节点标签
      p <- p %>% plotly::add_text(
        data = node_df,
        x = ~x,
        y = ~y,
        text = ~name,
        textposition = "top center",
        textfont = list(size = 9, color = "black"),
        hoverinfo = "skip",
        showlegend = FALSE
      )

      # 构建标题
      title_text <- sprintf(
        "Pathway Network: %s → %s<br><sub>%d nodes, %d edges (min_shared=%d)</sub>",
        left_group, right_group,
        nrow(node_df), nrow(edge_list), input$min_shared
      )

      # 设置布局
      p %>% plotly::layout(
        title = list(
          text = title_text,
          font = list(size = 14),
          x = 0.5,
          xanchor = "center"
        ),
        xaxis = list(
          title = "",
          showgrid = FALSE,
          showticklabels = FALSE,
          zeroline = FALSE
        ),
        yaxis = list(
          title = "",
          showgrid = FALSE,
          showticklabels = FALSE,
          zeroline = FALSE
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


  })
}
