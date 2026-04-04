# =============================================================================
# HubGene Network Module (完全模仿 Pathway Relationship Exploration 设计)
# =============================================================================

#' @title HubGene Network Module UI
#' @description User interface for HubGene Network visualization
#' @param id Module ID
#' @return Shiny UI tagList
#' @keywords internal

mod_hubgene_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::fluidRow(
      # ─── 左侧控制面板 ───
      shiny::column(
        width = 3,
        shiny::div(
          class = "well",
          style = "padding: 15px;",

          # ---- Mode 选择器 ----
          shiny::h4("Pathway Source Mode"),
          shiny::radioButtons(
            ns("hubgene_mode"),
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

          shiny::numericInput(
            ns("fdr_threshold"),
            label = "FDR Threshold:",
            min = 0,
            max = 1.0,
            value = 0.25,
            step = 0.01
          ),

          shiny::hr(),

          # ---- Top N 模式专属参数 ----
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] == 'mode_topN'", ns("hubgene_mode")),
            shiny::div(
              style = "background: #e3f2fd; padding: 10px; border-radius: 5px;",
              shiny::h5("Top N Configuration"),
              shiny::numericInput(
                ns("default_n_pathways"),
                label = "Top N Count:",
                value = 5,
                min = 3,
                max = 999,
                step = 1
              ),
              shiny::helpText(
                style = "color: #666; font-size: 11px;",
                "Display top N pathways by |NES|"
              )
            )
          ),

          # ---- Select 模式专属提示 ----
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] == 'mode_select'", ns("hubgene_mode")),
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

          # ---- Hub 基因阈值 ----
          shiny::h4("Hub Gene Filter"),
          shiny::numericInput(
            ns("min_hub_degree"),
            label = "Min. degree:",
            value = 2,
            min = 1,
            max = 20,
            step = 1
          ),
          shiny::helpText(
            style = "color: #666; font-size: 10px;",
            "Show genes appearing in >= N pathways"
          ),

          shiny::hr(),

          # ---- 节点大小 ----
          shiny::h4("Node Size"),
          shiny::selectInput(
            ns("size_mode"),
            label = NULL,
            choices = c(
              "By Hub Degree" = "degree",
              "By -log10(FDR)" = "fdr"
            ),
            selected = "degree"
          ),

          shiny::hr(),

          # ---- 布局算法 ----
          shiny::h4("Layout"),
          shiny::selectInput(
            ns("layout_algo"),
            label = NULL,
            choices = c(
              "Force-Directed (FR)" = "fr",
              "Kamada-Kawai" = "kk",
              "Circle" = "circle"
            ),
            selected = "fr"
          ),
          shiny::numericInput(
            ns("seed"),
            label = "Seed:",
            value = 42,
            min = 1,
            max = 9999,
            step = 1
          ),

          shiny::hr(),

          # ---- 显示选项 ----
          shiny::h4("Display"),
          shiny::checkboxInput(
            ns("show_gene_labels"),
            label = "Show gene labels",
            value = TRUE
          ),
          shiny::checkboxInput(
            ns("show_pathway_labels"),
            label = "Show pathway labels",
            value = TRUE
          ),

          shiny::hr(),

          # ---- 统计信息 ----
          shiny::h4("Statistics"),
          shiny::verbatimTextOutput(ns("network_stats")),

          shiny::hr(),

          # ---- 当前通路列表预览 ----
          shiny::h4("Pathways to Plot"),
          shiny::uiOutput(ns("pathway_preview_list"))
        )
      ),

      # ─── 主绘图区域 ───
      shiny::column(
        width = 9,
        shiny::div(
          class = "white-box",
          style = "min-height: 800px;",

          # 图例
          shiny::div(
            id = ns("legend_area"),
            style = "background: #f8f9fa; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
            shiny::fluidRow(
              shiny::column(
                width = 6,
                shiny::tags$strong("Pathways: "),
                shiny::span(
                  style = "display:inline-block; width:18px; height:18px; background:#E41A1C; margin-left:8px; border:2px solid black; vertical-align:middle;"
                ),
                shiny::span(" Up (NES>0) "),
                shiny::span(
                  style = "display:inline-block; width:18px; height:18px; background:#377EB8; margin-left:8px; border:2px solid black; vertical-align:middle;"
                ),
                shiny::span(" Down (NES<0)")
              ),
              shiny::column(
                width = 6,
                shiny::tags$strong("Genes: "),
                shiny::span(
                  style = "display:inline-block; width:14px; height:14px; background:#E41A1C; border-radius:50%; margin-left:8px; vertical-align:middle;"
                ),
                shiny::span(" Up "),
                shiny::span(
                  style = "display:inline-block; width:14px; height:14px; background:#377EB8; border-radius:50%; margin-left:8px; vertical-align:middle;"
                ),
                shiny::span(" Down")
              )
            )
          ),

          # Leading Edge 图例
          shiny::div(
            id = ns("leading_edge_legend"),
            style = "background: #f8f9fa; padding: 10px; border-radius: 5px; margin-top: 10px;",
            shiny::fluidRow(
              shiny::column(
                width = 6,
                shiny::tags$strong("Leading Edge: "),
                shiny::span(
                  style = "display:inline-block; width:20px; height:4px; background:#E41A1C; margin-left:8px; vertical-align:middle;"
                ),
                shiny::span(" Edge 3px, Opacity 100% "),
                shiny::span(
                  style = "display:inline-block; width:14px; height:14px; background:#E41A1C; border-radius:50%; margin-left:8px; vertical-align:middle;"
                ),
                shiny::span(" Node Full")
              ),
              shiny::column(
                width = 6,
                shiny::tags$strong("Non-Leading Edge: "),
                shiny::span(
                  style = "display:inline-block; width:20px; height:1px; background:#E41A1C; opacity:0.5; margin-left:8px; vertical-align:middle;"
                ),
                shiny::span(" Edge 1px, Opacity 50% "),
                shiny::span(
                  style = "display:inline-block; width:14px; height:14px; background:#E41A1C; border-radius:50%; opacity:0.5; margin-left:8px; vertical-align:middle;"
                ),
                shiny::span(" Node Dimmed")
              )
            )
          ),

          # 状态信息
          shiny::div(
            id = ns("mode_status"),
            style = "background: #d4edda; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
            shiny::uiOutput(ns("selection_summary"))
          ),

          # 绘图
          plotly::plotlyOutput(
            ns("hubgene_network"),
            height = "700px"
          ) %>%
            shinycssloaders::withSpinner(type = 6, color = "#28a745")
        )
      )
    )
  )
}


#' @title HubGene Network Module Server (完全模仿 Pathway Relationship Exploration 设计)
#' @description Server logic for HubGene Network with dual-mode pathway selection
#' @param id Module ID
#' @param data_prep_list Reactive data from data prep module
#' @param table_controller Table controller with selected_pathways
#' @keywords internal

mod_hubgene_server <- function(id, data_prep_list, table_controller) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ============================================================
    # 1. 模式状态管理
    # ============================================================

    hubgene_mode <- shiny::reactiveVal("mode_topN")

    shiny::observeEvent(input$hubgene_mode, {
      new_mode <- input$hubgene_mode
      if (!is.null(new_mode) && new_mode != hubgene_mode()) {
        hubgene_mode(new_mode)
        message(sprintf("[HubGene] Mode changed to: %s", new_mode))
      }
    })

    # ============================================================
    # 2. 数据源 reactive
    # ============================================================

    # Top N 候选通路
    topN_candidates <- shiny::reactive({
      if (hubgene_mode() != "mode_topN") return(character(0))
      data_list <- data_prep_list$data()
      shiny::req(data_list)
      df <- data_list$df
      shiny::req(nrow(df) > 0)
      top_n <- input$default_n_pathways
      if (is.null(top_n)) top_n <- 20
      top_n <- max(3, min(top_n, nrow(df)))
      df[1:top_n, "ID"]
    })

    # 主表选择候选通路
    select_candidates <- shiny::reactive({
      if (hubgene_mode() != "mode_select") return(character(0))
      sel <- table_controller$selected_pathways()
      if (is.null(sel) || length(sel) == 0) {
        return(character(0))
      }
      return(sel)
    })

    # 原始候选通路（根据模式选择）
    candidate_raw <- shiny::reactive({
      switch(hubgene_mode(),
             "mode_topN" = topN_candidates(),
             "mode_select" = select_candidates(),
             character(0))
    })

    # FDR 过滤后的候选通路
    candidate_filtered <- shiny::reactive({
      pathways <- candidate_raw()
      if (length(pathways) == 0) return(character(0))

      fdr_thresh <- input$fdr_threshold
      if (is.null(fdr_thresh)) fdr_thresh <- 0.25

      data_list <- data_prep_list$data()
      shiny::req(data_list)
      df <- data_list$df

      fdr_vec <- df$p.adjust[match(pathways, df$ID)]
      names(fdr_vec) <- pathways
      filtered <- pathways[!is.na(fdr_vec) & fdr_vec < fdr_thresh]
      return(filtered)
    })

    # ============================================================
    # 3. 最终通路列表管理
    # ============================================================

    final_pathways <- shiny::reactiveVal(character(0))

    shiny::observeEvent(candidate_filtered(), {
      new_candidates <- candidate_filtered()
      if (length(new_candidates) > 0) {
        final_pathways(new_candidates)
        message(sprintf("[HubGene] Final pathways updated: %d", length(new_candidates)))
      } else {
        final_pathways(character(0))
      }
    })

    # ============================================================
    # 4. 模式状态显示
    # ============================================================

    output$selection_summary <- shiny::renderUI({
      mode <- hubgene_mode()
      pathways <- final_pathways()

      if (mode == "mode_topN") {
        mode_label <- "Top N Mode"
        mode_style <- "background: #e3f2fd; padding: 10px; border-radius: 5px;"
        mode_color <- "#004085"
        top_n <- input$default_n_pathways %||% 20
        if (length(pathways) == 0) {
          detail <- sprintf("No pathways passed FDR < %.2f filter", input$fdr_threshold %||% 0.25)
        } else {
          detail <- sprintf("Showing top %d pathways by |NES|, filtered by FDR < %.2f", top_n, input$fdr_threshold %||% 0.25)
        }
      } else {
        mode_label <- "Select Mode"
        mode_style <- "background: #f3e5f5; padding: 10px; border-radius: 5px;"
        mode_color <- "#6a1b9a"
        sel <- table_controller$selected_pathways()
        if (is.null(sel) || length(sel) == 0) {
          detail <- "No pathways selected in main table. Please check 'Joint Plot' column."
        } else {
          detail <- sprintf("%d pathways selected from main table, filtered by FDR < %.2f", length(pathways), input$fdr_threshold %||% 0.25)
        }
      }

      shiny::div(
        style = mode_style,
        shiny::HTML(sprintf(
          "<strong style='color: %s;'>[%s]</strong> | %s",
          mode_color, mode_label, detail
        ))
      )
    })

    # ============================================================
    # 5. 通路预览列表 UI
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
        fdr_str <- if (!is.na(fdr)) sprintf("%.2e", fdr) else "N/A"
        nes_str <- if (!is.na(nes)) sprintf("%.2f", nes) else "N/A"

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
    # 6. 网络统计信息
    # ============================================================

    output$network_stats <- shiny::renderPrint({
      pathways <- final_pathways()
      if (length(pathways) == 0) {
        cat("No pathways to display\n")
        return()
      }

      cat("HubGene Network\n")
      cat(paste(rep("-", 25), collapse = ""), "\n")
      cat(sprintf("Pathways: %d\n", length(pathways)))
      cat(sprintf("Mode: %s\n", ifelse(hubgene_mode() == "mode_topN", "Top N", "Select")))
      cat(sprintf("FDR < %.2f\n", input$fdr_threshold %||% 0.25))
      cat(sprintf("Min Hub Degree: %d\n", input$min_hub_degree %||% 2))
    })

    # ============================================================
    # 7. Reactive: 获取通路列表（兼容旧接口）
    # ============================================================

    selected_pathways <- shiny::reactive({
      return(final_pathways())
    })

    # ============================================================
    # 8. Reactive: 构建网络数据
    # ============================================================

    network_data <- shiny::reactive({
      data_list <- data_prep_list$data()
      shiny::req(!is.null(data_list))

      pathway_ids <- final_pathways()
      shiny::req(length(pathway_ids) > 0)

      # 获取 DE 表
      de_df <- tryCatch({
        get_de_table(data_prep_list$gsea_res, data_list$contrast_id)
      }, error = function(e) {
        warning("[HubGene] Failed to get DE table: ", e$message)
        NULL
      })

      res_df <- data_list$df

      # 构建网络
      net <- build_hubgene_network(
        gsea_task = list(
          gsea_res = data_list$gsea_res,
          meta = list(
            left_group = data_list$left_group,
            right_group = data_list$right_group
          )
        ),
        pathway_ids = pathway_ids,
        min_hub_degree = input$min_hub_degree %||% 2,
        de_df = de_df,
        res_df = res_df
      )

      return(net)
    })

    # ============================================================
    # 9. Reactive: 准备节点数据
    # ============================================================

    node_data <- shiny::reactive({
      net <- network_data()
      shiny::req(!is.null(net))

      node_pos <- prepare_hubgene_nodes(
        network_data = net,
        layout = input$layout_algo %||% "fr",
        seed = input$seed %||% 42
      )

      shiny::req(!is.null(node_pos))

      size_mode <- input$size_mode %||% "degree"
      data_list <- data_prep_list$data()
      left_group <- data_list$left_group
      right_group <- data_list$right_group

      pathway_nodes <- node_pos$pathway
      gene_nodes <- node_pos$gene

      # 通路着色（NES必然非零）
      pathway_nodes$color <- ifelse(pathway_nodes$NES > 0, "#E41A1C", "#377EB8")

      # 基因着色（log2FC可能为零）
      gene_nodes$color <- ifelse(
        is.na(gene_nodes$stat), "#999999",
        ifelse(gene_nodes$stat > 0, "#E41A1C",
               ifelse(gene_nodes$stat < 0, "#377EB8", "#999999"))
      )

      # 通路节点大小
      if (size_mode == "fdr") {
        pathway_nodes$size <- pmin(pmax(-log10(pathway_nodes$FDR + 1e-100) * 5 + 18, 18), 45)
      } else {
        pathway_nodes$size <- pmin(pmax(-log10(pathway_nodes$FDR + 1e-100) * 5 + 18, 18), 45)
      }

      # 基因节点大小（按度数）
      if (nrow(gene_nodes) > 0) {
        gene_nodes$size <- pmin(pmax(gene_nodes$degree * 3 + 8, 8), 30)
      }

      return(list(
        pathway = pathway_nodes,
        gene = gene_nodes
      ))
    })

    # ============================================================
    # 10. 输出: 主绘图
    # ============================================================

    output$hubgene_network <- plotly::renderPlotly({
      data_list <- data_prep_list$data()
      shiny::req(!is.null(data_list))

      nodes <- node_data()
      net <- network_data()

      shiny::req(!is.null(nodes), !is.null(net))

      edges <- net$edges
      pathway_nodes <- nodes$pathway
      gene_nodes <- nodes$gene

      if (nrow(edges) == 0) {
        return(plotly::plot_ly() %>%
                 plotly::layout(title = list(
                   text = "No connections to display. Try adjusting parameters.",
                   font = list(size = 14)
                 )))
      }

      left_group <- data_list$left_group
      right_group <- data_list$right_group

      # 颜色定义
      color_up <- "#E41A1C"
      color_down <- "#377EB8"
      color_neutral <- "#999999"

      # 预计算通路颜色映射
      pathway_color_map <- setNames(
        ifelse(pathway_nodes$NES > 0, color_up, color_down),
        pathway_nodes$id
      )

      # 预计算基因颜色
      gene_color_map <- setNames(
        ifelse(gene_nodes$stat > 0, color_up, color_down),
        gene_nodes$id
      )

      p <- plotly::plot_ly(source = ns("hubgene_source"))

      # 添加边
      for (i in seq_len(nrow(edges))) {
        edge <- edges[i, ]

        gene_row <- gene_nodes[gene_nodes$id == edge$source, ]
        if (nrow(gene_row) == 0) next
        x0 <- gene_row$x
        y0 <- gene_row$y

        pw_row <- pathway_nodes[pathway_nodes$id == edge$target, ]
        if (nrow(pw_row) == 0) next
        x1 <- pw_row$x
        y1 <- pw_row$y

        # 边的颜色: 基于基因和通路的方向关系
        gene_dir <- gene_row$stat > 0
        pw_dir <- pw_row$NES > 0

        if (gene_dir == pw_dir) {
          if (gene_row$stat > 0) {
            edge_color <- color_up
          } else {
            edge_color <- color_down
          }
        } else {
          edge_color <- "#6B8E9F"
        }

        edge_width <- ifelse(edge$is_leading_edge, 3, 1)
        edge_opacity <- ifelse(edge$is_leading_edge, 0.9, 0.5)

        edge_hover <- sprintf(
          "<b>Gene:</b> %s<br><b>Pathway:</b> %s<br><b>Leading Edge:</b> %s",
          edge$source,
          edge$target,
          ifelse(edge$is_leading_edge, "YES", "NO")
        )

        p <- p %>% plotly::add_trace(
          type = "scatter",
          mode = "lines",
          x = c(x0, x1, NA),
          y = c(y0, y1, NA),
          line = list(color = edge_color, width = edge_width),
          opacity = edge_opacity,
          hoverinfo = "text",
          text = edge_hover,
          showlegend = FALSE,
          inherit = FALSE
        )
      }

      # 添加通路节点
      pathway_hover <- sapply(seq_len(nrow(pathway_nodes)), function(i) {
        row <- pathway_nodes[i, ]
        nes <- as.numeric(row$NES)
        fdr <- as.numeric(row$FDR)
        n_core <- as.integer(row$n_core)
        n_total <- as.integer(row$n_total)
        core_ratio <- if (n_total > 0) n_core / n_total * 100 else 0
        direction <- if (nes > 0) paste0("Up in ", left_group) else paste0("Up in ", right_group)

        sprintf("<b style='font-size:13px;'>%s</b>
            <b>NES:</b> %.3f
            <b>FDR:</b> %.2e
            <b>Direction:</b> %s
            <b>Core:</b> %d/%d (%.0f%%)",
                row$id, nes, fdr, direction, n_core, n_total, core_ratio)
      })

      p <- p %>% plotly::add_trace(
        type = "scatter",
        mode = "markers",
        x = pathway_nodes$x,
        y = pathway_nodes$y,
        marker = list(
          size = pathway_nodes$size,
          color = pathway_nodes$color,
          symbol = "square",
          line = list(color = "#333333", width = 1.5)
        ),
        text = pathway_hover,
        hovertemplate = "%{text}<extra></extra>",
        name = "Pathways",
        showlegend = TRUE,
        inherit = FALSE
      )

      # 添加基因节点
      if (nrow(gene_nodes) > 0) {
        gene_leading_status <- sapply(gene_nodes$id, function(g) {
          gene_edges <- edges[edges$source == g, ]
          any(gene_edges$is_leading_edge)
        })

        gene_opacity <- ifelse(gene_leading_status, 1.0, 0.5)

        gene_hover <- sapply(seq_len(nrow(gene_nodes)), function(i) {
          row <- gene_nodes[i, ]
          stat_val <- as.numeric(row$stat)
          degree <- as.integer(row$degree)
          pathways_str <- as.character(row$pathways)
          color <- as.character(gene_color_map[row$id])
          is_leading <- gene_leading_status[i]

          direction <- if (stat_val > 0) {
            paste0("Up in ", left_group)
          } else {
            paste0("Up in ", right_group)
          }

          sprintf("<b style='font-size:13px; color:%s;'>%s</b><br><b>Stat:</b> %.3f
              <b>Direction:</b> %s
              <b>Hub Degree:</b> %d
              <b>Leading Edge:</b> %s
              <b>Pathways:</b>
              <span style='font-size:10px; word-break:break-all;'>%s</span>",
                  color, row$id, stat_val, direction, degree,
                  ifelse(is_leading, "YES", "NO"),
                  pathways_str)
        })

        p <- p %>% plotly::add_trace(
          type = "scatter",
          mode = "markers",
          x = gene_nodes$x,
          y = gene_nodes$y,
          marker = list(
            size = gene_nodes$size,
            color = gene_nodes$color,
            symbol = "circle",
            opacity = gene_opacity,
            line = list(color = "white", width = 0.5)
          ),
          text = gene_hover,
          hovertemplate = "%{text}<extra></extra>",
          name = "Genes",
          showlegend = TRUE,
          inherit = FALSE
        )
      }

      # 添加标签
      if (isTRUE(input$show_pathway_labels)) {
        p <- p %>% plotly::add_trace(
          type = "scatter",
          mode = "text",
          x = pathway_nodes$x,
          y = pathway_nodes$y + 0.1,
          text = pathway_nodes$id,
          textposition = "top center",
          textfont = list(size = 7, color = "#333333"),
          hoverinfo = "skip",
          showlegend = FALSE,
          inherit = FALSE
        )
      }

      if (isTRUE(input$show_gene_labels) && nrow(gene_nodes) > 0) {
        p <- p %>% plotly::add_trace(
          type = "scatter",
          mode = "text",
          x = gene_nodes$x,
          y = gene_nodes$y + 0.06,
          text = gene_nodes$id,
          textposition = "top center",
          textfont = list(size = 6, color = "#555555"),
          hoverinfo = "skip",
          showlegend = FALSE,
          inherit = FALSE
        )
      }

      # 布局设置
      all_x <- c(pathway_nodes$x, gene_nodes$x)
      all_y <- c(pathway_nodes$y, gene_nodes$y)

      x_range <- c(min(all_x) - 1.2, max(all_x) + 1.2)
      y_range <- c(min(all_y) - 1.2, max(all_y) + 1.2)

      n_leading_edges <- sum(edges$is_leading_edge)
      n_total_edges <- nrow(edges)

      p <- p %>% plotly::layout(
        title = list(
          text = sprintf("HubGene Network: %d Pathways, %d Hub Genes, %d Edges (LE: %d)",
                         nrow(pathway_nodes), nrow(gene_nodes), n_total_edges, n_leading_edges),
          font = list(size = 13),
          x = 0.5,
          xanchor = "center"
        ),
        xaxis = list(
          showgrid = FALSE,
          showticklabels = FALSE,
          zeroline = FALSE,
          range = x_range,
          scaleanchor = "x",
          scaleratio = 1
        ),
        yaxis = list(
          showgrid = FALSE,
          showticklabels = FALSE,
          zeroline = FALSE,
          range = y_range
        ),
        dragmode = "pan",
        showlegend = TRUE,
        legend = list(
          orientation = "h",
          x = 0.5,
          y = -0.08,
          xanchor = "center"
        ),
        margin = list(l = 20, r = 20, t = 50, b = 60),
        paper_bgcolor = "white",
        plot_bgcolor = "white",
        hoverlabel = list(
          bgcolor = "white",
          bordercolor = "#666666",
          font = list(size = 11)
        )
      ) %>% plotly::config(
        displayModeBar = TRUE,
        displaylogo = FALSE,
        modeBarButtonsToRemove = c("lasso2d", "select2d"),
        responsive = TRUE
      )

      return(p)
    })

    # ============================================================
    # 11. 返回值（供其他模块调用）
    # ============================================================

    return(list(
      final_pathways = final_pathways,
      hubgene_mode = hubgene_mode,
      selected_pathways = selected_pathways
    ))

  })
}


# ==============================================================================
# 工具函数（保持不变，从 utils_hubgene.R 加载）
# ==============================================================================

# 以下函数已在 R/utils_hubgene.R 中定义，此处不再重复：
# - build_hubgene_network()
# - extract_hub_genes()
# - prepare_hubgene_nodes()
# - color_by_direction()
# - generate_hubgene_hover_text()
# - get_hubgene_legend()
