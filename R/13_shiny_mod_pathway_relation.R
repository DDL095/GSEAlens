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

                               # UpSet专用
                               shiny::conditionalPanel(
                                 condition = sprintf("input['%s'] == 'upset'", ns("active_tab")),
                                 shiny::sliderInput(
                                   ns("max_sets"),
                                   "最大集合数:",
                                   min = 2, max = 15, value = 8, step = 1
                                 ),
                                 shiny::sliderInput(
                                   ns("min_intersection"),
                                   "最小交集大小:",
                                   min = 1, max = 50, value = 5, step = 1
                                 )
                               ),

                               # Chord专用
                               shiny::conditionalPanel(
                                 condition = sprintf("input['%s'] == 'chord'", ns("active_tab")),
                                 shiny::sliderInput(
                                   ns("genes_per_pathway"),
                                   "每通路显示基因数:",
                                   min = 5, max = 50, value = 20, step = 1
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

                      # Sub 2: Chord
                      shiny::tabPanel(
                        title = shiny::HTML("🎻 EnrichmentChord"),
                        value = "chord",
                        shiny::div(class = "white-box", style = "margin-top: 15px;",
                                   shiny::plotOutput(ns("plot_chord"), height = "600px")
                        )
                      ),

                      # Sub 3: Network
                      shiny::tabPanel(
                        title = shiny::HTML("🕸️ Network"),
                        value = "network",
                        shiny::div(class = "white-box", style = "margin-top: 15px;",
                                   plotly::plotlyOutput(ns("plot_network"), height = "600px")
                        )
                      ),

                      # Sub 4: UpSet
                      shiny::tabPanel(
                        title = shiny::HTML("📊 UpSet"),
                        value = "upset",
                        shiny::div(class = "white-box", style = "margin-top: 15px;",
                                   # 使用ComplexUpset输出
                                   shiny::plotOutput(ns("plot_upset"), height = "600px")
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

    # 1. DotPlot实现（经典GSEA气泡图）
    output$plot_dotplot <- plotly::renderPlotly({
      shiny::req(selected_pathways(), current_task())
      shiny::req(input$refresh_plot)  # 依赖更新按钮

      pathways <- selected_pathways()
      task <- current_task()
      gsea_obj <- task$gsea_res
      res_df <- as.data.frame(gsea_obj@result)
      res_df <- res_df[res_df$ID %in% pathways, ]

      if (nrow(res_df) == 0) {
        return(plotly::plot_ly() %>% plotly::layout(title = "无可用通路数据"))
      }

      # 获取DE基因（用于ORA计算）
      # 从gsea_res中提取DE表（需要实现辅助函数）
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
        if (input$ratio_source == "ora" && length(term_genes) > 0) {
          ratio <- length(core_genes) / length(term_genes)
        } else if (length(de_genes_sig) > 0) {
          overlap <- intersect(core_genes, toupper(de_genes_sig))
          ratio <- length(overlap) / length(de_genes_sig)
        } else {
          ratio <- length(core_genes) / row$setSize
        }

        # Size映射
        size_val <- switch(input$size_mode,
                           "core_size" = length(core_genes),
                           "setsize" = row$setSize,
                           "ratio" = ratio * 100,
                           length(core_genes))

        # Color映射
        color_val <- switch(input$stat_color_mode,
                            "pval" = -log10(row$pvalue),
                            "padj" = -log10(row$p.adjust),
                            "nes" = abs(row$NES),
                            -log10(row$p.adjust))

        data.frame(
          Pathway = pid,
          Description = ifelse(is.na(row$Description), pid, row$Description),
          Ratio = ratio,
          Size = size_val,
          Color = color_val,
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

      # 绘制（仅hover，无click事件监听）
      p <- plotly::plot_ly(
        data = plot_df,
        x = ~Ratio,
        y = ~Pathway,
        type = "scatter",
        mode = "markers",
        marker = list(
          size = ~Size,
          color = ~Color,
          colorscale = "RdBu",
          showscale = TRUE,
          colorbar = list(title = color_title),
          opacity = input$dot_alpha,
          line = list(color = "black", width = 1)
        ),
        text = ~sprintf(
          "<b>%s</b><br>Ratio: %.3f<br>Core Genes: %d<br>NES: %.2f<br>FDR: %.2e",
          Description, Ratio, CoreCount, NES, FDR
        ),
        hoverinfo = "text"
      ) %>% plotly::layout(
        xaxis = list(title = ifelse(input$ratio_source == "ora",
                                    "Overlap Ratio (Core/Term)",
                                    "Overlap Ratio (Core/DE)")),
        yaxis = list(title = ""),
        showlegend = FALSE,
        dragmode = FALSE  # 禁止拖拽，确保无click联动
      )

      return(p)
    }) %>% shiny::bindEvent(input$refresh_plot, ignoreNULL = FALSE)

    # 2. EnrichmentChord（基因-通路弦图）
    # 2. EnrichmentChord（基因-通路弦图）- 修复空间不足问题
    output$plot_chord <- shiny::renderPlot({
      shiny::req(selected_pathways(), current_task())
      shiny::req(input$refresh_plot)

      pathways <- selected_pathways()
      task <- current_task()

      # 限制规模，避免gap.degree错误
      max_pathways <- min(length(pathways), 8)  # 最多8个通路
      if (length(pathways) > max_pathways) {
        pathways <- pathways[1:max_pathways]
      }

      # 构建基因-通路关联
      genes_list <- list()
      total_genes <- 0

      for (pid in pathways) {
        core_genes <- get_core_genes_for_pathway(task, pid)
        # 严格限制每个通路的基因数
        n_genes <- min(length(core_genes), input$genes_per_pathway %||% 10)
        if (n_genes > 0) {
          genes_list[[pid]] <- core_genes[1:n_genes]
          total_genes <- total_genes + n_genes
        }
      }

      # 安全检查：如果总节点数过多，进一步截断
      if (total_genes > 50) {
        # 重新分配配额
        genes_per_pathway <- floor(50 / length(pathways))
        for (pid in names(genes_list)) {
          if (length(genes_list[[pid]]) > genes_per_pathway) {
            genes_list[[pid]] <- genes_list[[pid]][1:genes_per_pathway]
          }
        }
      }

      # 构建边列表
      edges <- data.frame(
        from = character(0),
        to = character(0),
        stringsAsFactors = FALSE
      )

      for (pid in names(genes_list)) {
        genes <- genes_list[[pid]]
        if (length(genes) > 0) {
          edges <- rbind(edges, data.frame(
            from = genes,
            to = rep(pid, length(genes)),
            stringsAsFactors = FALSE
          ))
        }
      }

      if (nrow(edges) == 0 || length(unique(edges$from)) < 2) {
        graphics::plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        graphics::text(1, 1, "无足够基因数据\n(每个通路至少需要1个core gene)", col = "red")
        return()
      }

      # 构建邻接矩阵（仅包含出现在edges中的基因）
      all_genes <- unique(edges$from)
      all_paths <- unique(edges$to)

      # 如果维度太大，进一步截断
      if (length(all_genes) > 30 || length(all_paths) > 8) {
        graphics::plot(1, type = "n", axes = FALSE)
        graphics::text(1, 1, "数据维度太大\n请减少通路数量或降低每通路基因数", col = "red")
        return()
      }

      mat <- matrix(0, nrow = length(all_genes), ncol = length(all_paths))
      rownames(mat) <- all_genes
      colnames(mat) <- all_paths

      for (i in 1:nrow(edges)) {
        if (edges$from[i] %in% rownames(mat) && edges$to[i] %in% colnames(mat)) {
          mat[edges$from[i], edges$to[i]] <- 1
        }
      }

      # 设置circos参数（关键修复：减小gap.degree）
      circlize::circos.clear()
      circlize::circos.par(
        gap.degree = 2,  # 从默认5降到2，节省空间
        start.degree = 90,
        track.margin = c(-0.1, 0.1),
        points.overflow.warning = FALSE
      )

      # 准备颜色
      path_colors <- RColorBrewer::brewer.pal(max(3, length(all_paths)), "Set2")[1:length(all_paths)]
      names(path_colors) <- all_paths

      gene_color <- "grey70"

      grid_colors <- c(
        rep(gene_color, length(all_genes)),
        path_colors[all_paths]
      )

      # 绘制弦图（带错误处理）
      tryCatch({
        circlize::chordDiagram(
          mat,
          transparency = 0.4,
          annotationTrack = c("grid", "axis"),
          preAllocateTracks = list(track.height = 0.1),
          grid.col = grid_colors,
          link.border = "white",
          link.lwd = 0.5,
          direction = 1,
          diffHeight = 0.05
        )

        # 添加标签（仅通路，基因太多不显示）
        circlize::circos.track(
          track.index = 1,
          panel.fun = function(x, y) {
            sector.name <- get.cell.meta.data("sector.index")
            if (sector.name %in% all_paths) {
              circlize::circos.text(
                CELL_META$xcenter,
                CELL_META$ylim[1],
                sector.name,
                facing = "clockwise",
                niceFacing = TRUE,
                adj = c(0, 0.5),
                cex = 0.8,
                font = 2,
                col = path_colors[sector.name]
              )
            }
          },
          bg.border = NA
        )

        graphics::title(main = sprintf("Gene-Pathway Chord\n(%d genes × %d pathways)",
                                       length(all_genes), length(all_paths)))

      }, error = function(e) {
        graphics::text(0.5, 0.5, sprintf("Chord Diagram Error:\n%s\n\n建议：减少通路数量至5个以下",
                                         e$message), col = "red")
      })

      circlize::circos.clear()

    }, height = 600, width = 800)

    # 3. Network（通路关系网络，基于core_genes Jaccard相似度）
    # 3. Network（通路关系网络）- 修复连线显示
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
          title = "无足够连接的通路（请降低最小共享基因数）"
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

      # 构建边数据框（关键修复：为plotly创建线段数据）
      edge_segments <- data.frame(
        x = numeric(0),
        y = numeric(0),
        xend = numeric(0),
        yend = numeric(0),
        weight = numeric(0),
        shared = integer(0),
        stringsAsFactors = FALSE
      )

      for (i in 1:nrow(edge_list)) {
        from_pos <- node_df[node_df$name == edge_list$from[i], c("x", "y")]
        to_pos <- node_df[node_df$name == edge_list$to[i], c("x", "y")]

        if (nrow(from_pos) == 1 && nrow(to_pos) == 1) {
          edge_segments <- rbind(edge_segments, data.frame(
            x = from_pos$x,
            y = from_pos$y,
            xend = to_pos$x,
            yend = to_pos$y,
            weight = edge_list$weight[i],
            shared = edge_list$shared[i]
          ))
        }
      }

      # 使用ggplot2绘制（更稳定的线段渲染）
      p <- ggplot2::ggplot() +
        # 边（线段）
        ggplot2::geom_segment(
          data = edge_segments,
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
                       alpha = weight, size = weight),
          color = "gray50",
          show.legend = FALSE
        ) +
        # 节点
        ggplot2::geom_point(
          data = node_df,
          ggplot2::aes(x = x, y = y, color = color_val),
          size = 8,
          shape = 21,
          fill = "white",
          stroke = 2
        ) +
        # 节点标签
        ggplot2::geom_text(
          data = node_df,
          ggplot2::aes(x = x, y = y, label = name),
          size = 3,
          vjust = -1.5,
          check_overlap = TRUE
        ) +
        ggplot2::scale_color_gradient2(
          low = "#377EB8", mid = "white", high = "#E41A1C", midpoint = 0,
          name = "NES"
        ) +
        ggplot2::scale_size_continuous(range = c(0.5, 3)) +
        ggplot2::scale_alpha_continuous(range = c(0.3, 0.8)) +
        ggplot2::coord_fixed() +
        ggplot2::theme_void() +
        ggplot2::theme(
          legend.position = "right",
          plot.title = ggplot2::element_text(hjust = 0.5)
        ) +
        ggplot2::labs(
          title = sprintf("Pathway Network (edges=%d, min_shared=%d)",
                          nrow(edge_segments), input$min_shared)
        )

      # 转换为plotly（保留悬停，移除click）
      ply <- plotly::ggplotly(p, tooltip = c("label", "color_val")) %>%
        plotly::layout(
          dragmode = FALSE,
          hovermode = "closest"
        ) %>%
        plotly::config(displayModeBar = FALSE)

      return(ply)
    }) %>% shiny::bindEvent(input$refresh_plot, ignoreNULL = FALSE)

    # 4. UpSet（交集分析）
    # 4. UpSet（交集分析）- 修复参数名
    output$plot_upset <- shiny::renderPlot({
      shiny::req(selected_pathways(), current_task())
      shiny::req(input$refresh_plot)

      # 检查ComplexUpset是否可用
      if (!requireNamespace("ComplexUpset", quietly = TRUE)) {
        graphics::plot(1, type = "n", axes = FALSE)
        graphics::text(1, 1, "请安装ComplexUpset包\ninstall.packages('ComplexUpset')", col = "red")
        return()
      }

      pathways <- selected_pathways()
      task <- current_task()

      # 限制集合数
      max_sets <- validate_param(input$max_sets, 8, 2, 15, "max_sets")
      if (length(pathways) > max_sets) {
        pathways <- pathways[1:max_sets]
      }

      # 构建集合列表
      core_list <- get_core_genes_list(task, pathways)

      # 转换为数据框（ComplexUpset格式）
      all_genes <- unique(unlist(core_list))
      if (length(all_genes) == 0) {
        graphics::plot(1, type = "n", axes = FALSE)
        graphics::text(1, 1, "无Core Genes数据", col = "red")
        return()
      }

      membership_df <- data.frame(
        gene = all_genes,
        stringsAsFactors = FALSE
      )

      for (pid in pathways) {
        membership_df[[pid]] <- membership_df$gene %in% core_list[[pid]]
      }

      # 修复：ComplexUpset 1.0+ 使用 mode = 'distinct' 替代 min_intersection_size
      # 或使用 intersections 参数控制显示哪些交集
      min_size <- input$min_intersection

      # 方法：先计算所有交集，过滤后再绘图
      tryCatch({
        ComplexUpset::upset(
          membership_df,
          intersect = pathways,
          min_size = min_size,  # 使用min_size而非min_intersection_size
          width_ratio = 0.2,
          name = "Pathway Intersections",
          base_annotations = list(
            'Intersection size' = ComplexUpset::intersection_size(
              text = ComplexUpset::aes_(color = 'black'),
              text_mapping = ComplexUpset::aes_(label = ggplot2::aes(label = !!ggplot2::sym('size')))
            )
          ),
          themes = ComplexUpset::upset_default_themes(
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank()
          )
        )
      }, error = function(e) {
        # 降级方案：使用基础 upset
        graphics::plot(1, type = "n", axes = FALSE)
        graphics::text(1, 1, sprintf("UpSet Error:\n%s", e$message), col = "red")
      })

    }, height = 600, width = 800)

  })
}
