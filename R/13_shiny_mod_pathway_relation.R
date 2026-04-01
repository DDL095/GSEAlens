#' @title Pathway Relationship Exploration Module UI
#' @description Network visualization module with dual-mode pathway selection
#' @keywords internal

mod_pathway_relation_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      # ===== 左侧控制面板 =====
      shiny::column(
        2,
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

          shiny::numericInput(
            ns("fdr_threshold"),
            label = "FDR Threshold:",
            min = 0,
            max = 1.0,
            value = 0.25,
            step = 0.01
          ),

          shiny::numericInput(
            ns("min_shared"),
            label = "Min Shared Core Genes:",
            value = 3,
            min = 1,
            max = 9999,
            step = 1
          ),

          shiny::sliderInput(
            ns("max_nodes"),
            label = "Max Nodes:",
            min = 5,
            max = 100,
            value = 30,
            step = 5
          ),

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

          # ---- Hover 显示参数 ----
          shiny::h4("Hover Display Settings"),

          shiny::sliderInput(
            ns("hover_max_genes"),
            label = "Max Genes in Hover:",
            min = 3,
            max = 30,
            value = 5,
            step = 1
          ),

          shiny::helpText(
            style = "color: #666; font-size: 11px;",
            "Number of shared genes to display in edge hover tooltip"
          ),

          shiny::hr(),

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

          shiny::h4("Pathways to Plot"),
          shiny::uiOutput(ns("pathway_preview_list"))
        )
      ),

      # ===== 右侧面板 =====
      shiny::column(
        10,
        shiny::tabsetPanel(
          id = ns("active_tab"),
          type = "tabs",

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
              plotly::plotlyOutput(ns("plot_dotplot"), height = "800px") %>%
                shinycssloaders::withSpinner(type = 6, color = "#28a745")
            )
          ),

          shiny::tabPanel(
            title = "Network",
            value = "network",
            shiny::div(
              class = "white-box",
              style = "padding: 15px; margin-top: 15px;",

              shiny::uiOutput(ns("network_status")),

              shiny::div(
                id = ns("selection_panel"),
                style = "background: #e8f4fd; padding: 12px; border-radius: 8px; margin-bottom: 15px;",
                shiny::fluidRow(
                  shiny::column(12,
                                shiny::uiOutput(ns("selection_display"))
                  )
                ),
                shiny::fluidRow(
                  shiny::column(6,
                                shiny::actionButton(
                                  ns("show_edge_detail"),
                                  label = "Show Edge Detail",
                                  class = "btn-primary",
                                  style = "width: 100%;",
                                  disabled = NA
                                )
                  ),
                  shiny::column(6,
                                shiny::actionButton(
                                  ns("clear_selection"),
                                  label = "Clear Selection",
                                  class = "btn-warning",
                                  style = "width: 100%;"
                                )
                  )
                ),
                shiny::div(
                  style = "margin-top: 8px; font-size: 11px; color: #666;",
                  shiny::HTML("<b>Instruction:</b> Click two nodes to select, then click 'Show Edge Detail'")
                )
              ),

              plotly::plotlyOutput(ns("plot_network"), height = "1200px") %>%
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
      }
    })

    # ============================================================
    # 2. 数据源 reactive
    # ============================================================

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

    select_candidates <- shiny::reactive({
      if (network_mode() != "mode_select") return(character(0))
      if (is.null(table_result) || is.null(table_result$selected_pathways)) return(character(0))
      selected <- table_result$selected_pathways()
      if (is.null(selected)) return(character(0))
      return(selected)
    })

    candidate_raw <- shiny::reactive({
      switch(network_mode(),
             "mode_topN" = topN_candidates(),
             "mode_select" = select_candidates(),
             character(0)
      )
    })

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
        max_n <- input$max_nodes
        if (is.null(max_n)) max_n <- 999
        final <- new_candidates[1:min(length(new_candidates), max_n)]
        final_pathways(final)
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
        fdr_str <- if (!is.na(fdr)) sprintf("%.2e", fdr) else "N/A"
        nes_str <- if (!is.na(nes)) sprintf("%.2f", nes) else "N/A"
        shiny::div(
          style = "background: #f8f9fa; padding: 8px; margin-bottom: 5px; border-radius: 4px; border-left: 3px solid #007bff;",
          shiny::div(style = "font-size: 12px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap;", shiny::strong(pid)),
          shiny::div(style = "font-size: 11px; color: #666;", sprintf("FDR: %s | NES: %s", fdr_str, nes_str))
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
        sprintf("Pathways: %d", length(pathways)),
        sprintf("FDR < %.2f", input$fdr_threshold)
      )
      if (length(pathways) == 0) {
        div_style <- "background: #fff3cd; padding: 10px; border-radius: 5px; margin-bottom: 10px; color: #856404;"
      } else {
        div_style <- "background: #d4edda; padding: 10px; border-radius: 5px; margin-bottom: 10px; color: #155724;"
      }
      shiny::div(style = div_style, shiny::HTML(paste(status_parts, collapse = " | ")))
    })

    # ============================================================
    # 6. Network 状态提示
    # ============================================================

    output$network_status <- shiny::renderUI({
      pathways <- final_pathways()
      status_parts <- c(
        sprintf("Mode: %s", ifelse(network_mode() == "mode_topN", "TopN", "Select")),
        sprintf("Pathways: %d", length(pathways)),
        sprintf("FDR < %.2f", input$fdr_threshold),
        sprintf("min_shared: %d", input$min_shared),
        sprintf("seed: %d", input$seed)
      )
      if (length(pathways) == 0) {
        div_style <- "background: #fff3cd; padding: 10px; border-radius: 5px; margin-bottom: 10px; color: #856404;"
      } else {
        div_style <- "background: #d4edda; padding: 10px; border-radius: 5px; margin-bottom: 10px; color: #155724;"
      }
      shiny::div(style = div_style, shiny::HTML(paste(status_parts, collapse = " | ")))
    })

    # ============================================================
    # 7. DotPlot 绘图函数
    # ============================================================

    output$plot_dotplot <- plotly::renderPlotly({
      pathways <- final_pathways()
      if (length(pathways) == 0) {
        return(plotly::plot_ly() %>% plotly::layout(
          title = list(text = "No pathways to display", font = list(size = 14), x = 0.5),
          xaxis = list(showgrid = FALSE), yaxis = list(showgrid = FALSE)
        ))
      }
      data_list <- data_prep_list$data()
      shiny::req(data_list)
      df <- data_list$df
      plot_df <- df[df$ID %in% pathways, ]
      if (nrow(plot_df) == 0) {
        return(plotly::plot_ly() %>% plotly::layout(title = list(text = "No matching pathways", font = list(size = 14))))
      }
      task <- list(
        gsea_res = data_list$gsea_res,
        meta = list(
          left_group = data_list$left_group,
          right_group = data_list$right_group,
          contrast_id = data_list$contrast_id
        )
      )
      class(task) <- "GseaTask"
      core_list <- tryCatch({get_core_genes_list(task, plot_df$ID)}, error = function(e) NULL)
      plot_df$CoreCount <- sapply(plot_df$ID, function(pid) {
        if (is.null(core_list) || is.null(core_list[[pid]])) return(0)
        length(core_list[[pid]])
      })
      color_mode <- input$dotplot_color_mode
      color_vals <- switch(color_mode,
                           "padj" = -log10(plot_df$p.adjust),
                           "pval" = -log10(plot_df$pvalue),
                           "nes" = abs(plot_df$NES),
                           -log10(plot_df$p.adjust))
      color_title <- switch(color_mode, "padj" = "-log10(FDR)", "pval" = "-log10(P-value)", "nes" = "|NES|")
      size_mode <- input$dotplot_size_mode
      size_vals <- switch(size_mode, "core_size" = plot_df$CoreCount, "setsize" = plot_df$setSize, plot_df$CoreCount)
      size_range <- c(5, 25)
      if (max(size_vals) > min(size_vals)) {
        size_scaled <- size_range[1] + (size_vals - min(size_vals)) / (max(size_vals) - min(size_vals)) * (size_range[2] - size_range[1])
      } else {
        size_scaled <- mean(size_range)
      }
      plot_df <- plot_df[order(plot_df$NES, decreasing = TRUE), ]
      size_scaled <- size_scaled[order(plot_df$NES, decreasing = TRUE)]
      color_vals <- color_vals[order(plot_df$NES, decreasing = TRUE)]
      hover_text <- sprintf(
        "<b>%s</b><br>FDR: %.2e<br>NES: %.2f<br>Core Genes: %d",
        plot_df$ID, plot_df$p.adjust, plot_df$NES, plot_df$CoreCount
      )
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
          text = sprintf("Pathway DotPlot: %s vs %s", data_list$left_group, data_list$right_group),
          font = list(size = 12), x = 0.5, xanchor = "center"
        ),
        xaxis = list(title = "NES", zeroline = TRUE),
        yaxis = list(title = "", tickmode = "array", tickvals = seq_len(nrow(plot_df)), ticktext = plot_df$ID, tickfont = list(size = 9)),
        margin = list(l = 250, r = 50, t = 80, b = 50),
        showlegend = FALSE
      ) %>% plotly::config(displayModeBar = TRUE, displaylogo = FALSE)
    })

    # ============================================================
    # 8. 节点选择状态管理
    # ============================================================

    selected_nodes <- shiny::reactiveVal(character(0))

    output$selection_display <- shiny::renderUI({
      sel <- selected_nodes()
      if (length(sel) == 0) {
        shiny::div(style = "color: #666; font-style: italic; text-align: center;", "No nodes selected. Click on network nodes.")
      } else if (length(sel) == 1) {
        shiny::div(style = "color: #007bff; font-weight: bold; text-align: center;", sprintf("Selected: %s", sel[1]))
      } else {
        shiny::div(style = "color: #28a745; font-weight: bold; text-align: center;", sprintf("Selected: %s with %s", sel[1], sel[2]))
      }
    })

    shiny::observe({
      sel <- selected_nodes()
      if (length(sel) == 2) {
        shinyjs::enable("show_edge_detail")
      } else {
        shinyjs::disable("show_edge_detail")
      }
    })

    shiny::observeEvent(input$clear_selection, {
      selected_nodes(character(0))
    })

    shiny::observeEvent(input$show_edge_detail, {
      sel <- selected_nodes()
      if (length(sel) != 2) {
        shiny::showNotification("Please select exactly 2 nodes first", type = "warning")
        return()
      }

      from_pw <- sel[1]
      to_pw <- sel[2]

      edge_idx <- which(
        (edge_list$from == from_pw & edge_list$to == to_pw) |
          (edge_list$from == to_pw & edge_list$to == from_pw)
      )

      if (length(edge_idx) == 0) {
        shiny::showNotification(
          "No edge found between selected nodes. Try reducing min_shared threshold.",
          type = "warning",
          duration = 5
        )
        return()
      }

      i <- edge_idx[1]
      shared_count <- edge_list$shared[i]
      jaccard <- edge_list$weight[i]
      overlap_coef <- edge_list$overlap_coef[i]
      dice_coef <- edge_list$dice_coef[i]
      shared_genes <- unlist(edge_list$shared_genes[[i]])

      pathway_a_genes <- core_list[[from_pw]]
      pathway_b_genes <- core_list[[to_pw]]

      # ============ 修改开始：showdetail 中显示所有基因 ============
      # 在 showdetail 中不限制基因数量，显示所有共享基因
      # 使用逗号和空格隔离基因

      # 对所有基因进行格式化显示
      gene_buttons <- sapply(shared_genes, function(g) {
        sprintf(
          '<span style="display:inline-block;background:#e3f2fd;padding:4px 8px;margin:2px;border-radius:4px;font-size:12px;">%s</span>',
          g
        )
      })

      # 不再显示"+N more"提示
      gene_display <- paste(gene_buttons, collapse = ",")

      # Modal UI
      shiny::showModal(shiny::modalDialog(
        title = HTML(paste0("<b>Edge Detail:</b> ", from_pw, " with ", to_pw)),
        size = "l",
        easyClose = TRUE,
        footer = shiny::modalButton("Close"),

        # 相似度指标
        HTML("<div style='background:#f8f9fa;padding:15px;border-radius:8px;margin-bottom:15px;'>"),
        HTML("<h4 style='margin-top:0;'>Similarity Metrics</h4>"),
        HTML("<table style='width:100%;text-align:center;'>"),
        HTML("<tr><td><h3>", sprintf("%.3f", jaccard), "</h3><small>Jaccard</small></td>"),
        HTML("<td><h3>", sprintf("%.3f", overlap_coef), "</h3><small>Overlap</small></td>"),
        HTML("<td><h3>", sprintf("%.3f", dice_coef), "</h3><small>Dice</small></td>"),
        HTML("<td><h3>", shared_count, "</h3><small>Shared</small></td></tr>"),
        HTML("</table></div>"),

        # 共享基因
        HTML("<div style='background:#fff;padding:15px;border-radius:8px;border:1px solid #dee2e6;'>"),
        HTML("<h4 style='margin-top:0;'>Shared Core Genes (", shared_count, ")</h4>"),
        HTML(gene_display),
        HTML("</div>")
      ))
    })

    # ============================================================
    # 9. Network 绘图函数
    # ============================================================

    edge_list <- NULL
    node_df <- NULL
    core_list <- NULL

    output$plot_network <- plotly::renderPlotly({
      pathways <- final_pathways()

      if (length(pathways) == 0) {
        return(plotly::plot_ly() %>% plotly::layout(
          title = list(text = "No pathways to display", font = list(size = 14), x = 0.5),
          xaxis = list(showgrid = FALSE, showticklabels = FALSE),
          yaxis = list(showgrid = FALSE, showticklabels = FALSE)
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

      core_list <<- tryCatch({get_core_genes_list(task, pathways)}, error = function(e) NULL)

      if (is.null(core_list) || length(core_list) == 0) {
        return(plotly::plot_ly() %>% plotly::layout(
          title = list(text = "No core genes found", font = list(size = 14)),
          xaxis = list(showgrid = FALSE), yaxis = list(showgrid = FALSE)
        ))
      }

      valid_pathways <- names(core_list)[sapply(core_list, function(x) length(x) > 0)]
      if (length(valid_pathways) == 0) {
        return(plotly::plot_ly() %>% plotly::layout(
          title = list(text = "All pathways have empty core genes", font = list(size = 14)),
          xaxis = list(showgrid = FALSE), yaxis = list(showgrid = FALSE)
        ))
      }

      min_shared <- input$min_shared
      if (is.null(min_shared)) min_shared <- 3

      # 获取 hover 显示基因数上限
      hover_max_genes <- input$hover_max_genes
      if (is.null(hover_max_genes)) hover_max_genes <- 10

      # 构建边列表（现在包含 jaccard, overlap_coef, dice_coef）
      edge_list <<- tryCatch({
        build_edge_list_safely(core_list[valid_pathways], min_shared_genes = min_shared)
      }, error = function(e) NULL)

      if (is.null(edge_list) || nrow(edge_list) == 0) {
        return(plotly::plot_ly() %>% plotly::layout(
          title = list(text = paste0("No edges found (min_shared=", min_shared, ")"), font = list(size = 14), x = 0.5),
          xaxis = list(showgrid = FALSE, showticklabels = FALSE),
          yaxis = list(showgrid = FALSE, showticklabels = FALSE)
        ))
      }

      # ============================================================
      # 关键改进：基于 Jaccard 排名映射边缘粗细到 1-5 像素
      # ============================================================

      # 按 Jaccard (weight) 降序排列
      edge_list <- edge_list[order(edge_list$weight, decreasing = TRUE), ]

      # 计算粗细等级
      n_edges <- nrow(edge_list)
      edge_width_mapping <- function(rank, n) {
        # 根据排名分配 1-5 像素
        # rank 1 -> 5px, rank n -> 1px
        width <- 5 - (rank - 1) * (4 / max(1, n - 1))
        return(max(1, min(5, round(width))))
      }

      edge_list$width_rank <- sapply(seq_len(n_edges), function(i) {
        edge_width_mapping(i, n_edges)
      })

      # 普通边粗细范围：1-5 像素
      edge_list$edge_width_normal <- edge_list$width_rank

      # 选中边固定粗细：比普通最粗大 7-8 像素 (5 + 8 = 13)
      edge_list$edge_width_selected <- 13

      # 构建 igraph 对象
      g <- tryCatch({
        igraph::graph_from_data_frame(edge_list, directed = FALSE, vertices = valid_pathways)
      }, error = function(e) NULL)

      if (is.null(g)) {
        return(plotly::plot_ly() %>% plotly::layout(
          title = list(text = "Graph construction failed", font = list(size = 14)),
          xaxis = list(showgrid = FALSE), yaxis = list(showgrid = FALSE)
        ))
      }

      layout_algo <- input$network_layout
      seed_val <- input$seed

      layout_coords <- tryCatch({
        switch(layout_algo,
               "fr" = {set.seed(seed_val); igraph::layout_with_fr(g)},
               "kk" = igraph::layout_with_kk(g),
               "circle" = igraph::layout_in_circle(g),
               {set.seed(seed_val); igraph::layout_with_fr(g)}
        )
      }, error = function(e) NULL)

      if (is.null(layout_coords)) {
        return(plotly::plot_ly() %>% plotly::layout(
          title = list(text = "Layout calculation failed", font = list(size = 14)),
          xaxis = list(showgrid = FALSE), yaxis = list(showgrid = FALSE)
        ))
      }

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

      sel_nodes <- selected_nodes()

      node_hover_text <- sprintf(
        "<b>%s</b><br>NES: %.2f<br>FDR: %.2e<br>Core Genes: %d",
        node_df$name, node_df$NES, node_df$FDR, node_df$CoreCount
      )

      title_text <- sprintf(
        "Pathway Network: %s vs %s<br><sub>%d nodes, %d edges | Width: Jaccard rank (1-5px)</sub>",
        task$meta$left_group, task$meta$right_group,
        igraph::vcount(g), igraph::ecount(g)
      )

      p <- plotly::plot_ly(source = ns("network_plot"))

      # ============================================================
      # 添加边（使用 Jaccard 排名映射的粗细）
      # ============================================================

      for (i in seq_len(nrow(edge_list))) {
        from_node <- node_df[node_df$name == edge_list$from[i], ]
        to_node <- node_df[node_df$name == edge_list$to[i], ]
        is_selected <- (edge_list$from[i] %in% sel_nodes) && (edge_list$to[i] %in% sel_nodes)

        if (is_selected) {
          # 选中边：橙色 + 固定粗细 13px
          edge_color <- "#FF6600"
          edge_width <- edge_list$edge_width_selected[i]
        } else {
          # 普通边：灰色 + Jaccard 排名粗细 1-5px
          edge_color <- "rgba(150, 150, 150, 0.5)"
          edge_width <- edge_list$edge_width_normal[i]
        }

        p <- p %>% plotly::add_trace(
          type = "scatter", mode = "lines",
          x = c(from_node$x, to_node$x, NA),
          y = c(from_node$y, to_node$y, NA),
          line = list(color = edge_color, width = edge_width),
          hoverinfo = "skip", showlegend = FALSE, inherit = TRUE
        )
      }

      # ============================================================
      # 添加边 hover（增强版：显示共享基因列表）
      # ============================================================

      for (i in seq_len(nrow(edge_list))) {

        from_node <- node_df[node_df$name == edge_list$from[i], ]
        to_node   <- node_df[node_df$name == edge_list$to[i], ]

        is_selected <- (edge_list$from[i] %in% sel_nodes) && (edge_list$to[i] %in% sel_nodes)

        shared_genes_vec <- edge_list$shared_genes[[i]]
        shared_count <- edge_list$shared[i]

        if (length(shared_genes_vec) > 0) {
          display_genes <- shared_genes_vec[1:min(length(shared_genes_vec), hover_max_genes)]
          genes_display <- paste(display_genes, collapse = ", ")
          if (length(shared_genes_vec) > hover_max_genes) {
            remaining <- length(shared_genes_vec) - hover_max_genes
            genes_display <- paste0(
              genes_display,
              sprintf(" <span style='color:#ccc;'>(+%d more)</span>", remaining)
            )
          }
        } else {
          genes_display <- "(none)"
        }

        # 简化格式，减少换行
        edge_hover <- sprintf(
          "<b>%s with %s</b><br>Shared Genes (%d): %s<br>Jaccard: %.4f | Overlap: %.4f | Dice: %.4f",
          edge_list$from[i], edge_list$to[i], shared_count,
          genes_display,
          edge_list$weight[i], edge_list$overlap_coef[i], edge_list$dice_coef[i]
        )

        hover_bg <- if (is_selected) "#FF8C00" else "#333"

        p <- p %>% plotly::add_trace(
          type = "scatter", mode = "markers",
          x = c(mean(c(from_node$x, to_node$x))),
          y = c(mean(c(from_node$y, to_node$y))),
          marker = list(size = 12, opacity = 0, color = "transparent"),
          text = edge_hover,
          hoverinfo = "text",
          hoverlabel = list(
            bgcolor = hover_bg,
            font = list(color = "white", size = 12),
            align = "left",     # ← 左对齐
            bordercolor = hover_bg
          ),
          showlegend = FALSE,
          inherit = TRUE
        )
      }

      # ============================================================
      # 添加节点
      # ============================================================

      node_colors <- ifelse(node_df$name %in% sel_nodes, "#FFD700", node_df$color_val)
      node_sizes <- ifelse(node_df$name %in% sel_nodes, 25, 15)
      node_border <- ifelse(node_df$name %in% sel_nodes, 4, 1.5)

      p <- p %>% plotly::add_trace(
        type = "scatter", mode = "markers",
        x = node_df$x, y = node_df$y,
        marker = list(
          size = node_sizes,
          color = node_colors,
          colorscale = list(list(0, "#377EB8"), list(0.5, "white"), list(1, "#E41A1C")),
          cauto = FALSE, cmin = -3, cmax = 3, showscale = TRUE,
          colorbar = list(title = list(text = "NES", font = list(size = 12)), len = 0.5, y = 0.5),
          line = list(color = "black", width = node_border)
        ),
        text = node_df$name,
        hovertemplate = paste(node_hover_text, "<extra></extra>"),
        showlegend = FALSE, key = node_df$name, inherit = TRUE
      )

      # ============================================================
      # 添加标签（选中节点更大更粗）
      # ============================================================

      label_size <- ifelse(node_df$name %in% sel_nodes, 16, 9)  # 从 12/9 改为 16/9
      label_color <- ifelse(node_df$name %in% sel_nodes, "#FF6600", "#333")
      label_bold <- ifelse(node_df$name %in% sel_nodes, "bold", "normal")

      p <- p %>% plotly::add_trace(
        type = "scatter", mode = "text",
        x = node_df$x, y = node_df$y + 0.12,
        text = node_df$name,
        textposition = "top center",
        textfont = list(size = label_size, color = label_color, font = label_bold),
        hoverinfo = "skip", showlegend = FALSE, inherit = TRUE
      )

      x_range <- c(min(layout_coords[, 1]) - 0.5, max(layout_coords[, 1]) + 0.5)
      y_range <- c(min(layout_coords[, 2]) - 0.5, max(layout_coords[, 2]) + 0.5)

      p %>% plotly::layout(
        title = list(text = title_text, font = list(size = 12), x = 0.5, xanchor = "center"),
        xaxis = list(
          title = "",
          showgrid = FALSE,
          showticklabels = FALSE,
          zeroline = FALSE,
          range = x_range,
          scaleanchor = "x",      # 保持宽高比
          scaleratio = 1
        ),
        yaxis = list(
          title = "",
          showgrid = FALSE,
          showticklabels = FALSE,
          zeroline = FALSE,
          range = y_range,
          scaleanchor = "x",      # 保持宽高比
          scaleratio = 1
        ),
        hovermode = "closest",
        dragmode = "pan",
        showlegend = FALSE,
        margin = list(l = 50, r = 50, t = 80, b = 50),
        paper_bgcolor = 'rgba(0,0,0,0)',   # 透明背景
        plot_bgcolor = 'rgba(0,0,0,0)',    # 透明画布
        autosize = TRUE                    # ← 关键：自动调整大小
      ) %>% plotly::config(
        displayModeBar = TRUE,
        displaylogo = FALSE,
        modeBarButtonsToRemove = c("lasso2d", "select2d"),
        responsive = TRUE                  # ← 关键：响应式
      )
    })

    # ============================================================
    # 10. 节点点击事件监听
    # ============================================================

    shiny::observeEvent(plotly::event_data("plotly_click", source = ns("network_plot")), {
      click_data <- plotly::event_data("plotly_click", source = ns("network_plot"))
      if (is.null(click_data) || is.null(click_data$pointNumber)) return()

      n_edges <- nrow(edge_list)
      point_idx <- click_data$pointNumber + 1

      if (point_idx > nrow(node_df) || point_idx < 1) return()

      node_name <- as.character(node_df$name[point_idx])
      current_sel <- selected_nodes()

      if (length(current_sel) < 2) {
        if (!(node_name %in% current_sel)) {
          selected_nodes(c(current_sel, node_name))
        }
      } else {
        selected_nodes(node_name)
      }
    })

    # ============================================================
    # 11. 返回值
    # ============================================================

    return(list(
      final_pathways = final_pathways,
      network_mode = network_mode
    ))

  })
}
