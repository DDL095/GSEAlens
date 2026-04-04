# ==============================================================================
# 文件：R/18_shiny_mod_hubgene_vis.R
# 功能：HubGene Network - visNetwork 交互版
# 依赖：visNetwork, 复用 build_hubgene_network() 等现有函数
# ==============================================================================

#' @title HubGene Network (visNetwork) Module UI
#' @keywords internal

mod_hubgene_vis_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::fluidRow(
      # ─── 左侧控制面板 ───
      shiny::column(
        width = 3,
        shiny::div(
          class = "well",
          style = "padding: 15px; max-height: 90vh; overflow-y: auto;",

          shiny::h4("HubGene Network (visNetwork)"),

          shiny::hr(),

          shiny::h5("Pathway Source Mode"),
          shiny::radioButtons(
            ns("vis_mode"),
            label = NULL,
            choices = c(
              "Top N from Current Set" = "mode_topN",
              "Selected from Main Table" = "mode_select"
            ),
            selected = "mode_topN",
            width = "100%"
          ),

          shiny::hr(),

          shiny::conditionalPanel(
            condition = sprintf("input['%s'] == 'mode_topN'", ns("vis_mode")),
            shiny::div(
              style = "background: #e3f2fd; padding: 10px; border-radius: 5px;",
              shiny::numericInput(
                ns("vis_topN"),
                label = "Top N Count:",
                value = 5,
                min = 3,
                max = 50,
                step = 1
              )
            )
          ),

          shiny::conditionalPanel(
            condition = sprintf("input['%s'] == 'mode_select'", ns("vis_mode")),
            shiny::div(
              style = "background: #f3e5f5; padding: 10px; border-radius: 5px;",
              shiny::helpText(
                style = "color: #6a1b9a; font-size: 11px;",
                "Check pathways in main table 'Joint Plot' column"
              )
            )
          ),

          shiny::hr(),

          shiny::h5("Network Parameters"),

          shiny::numericInput(
            ns("vis_fdr"),
            label = "FDR Threshold:",
            value = 0.25,
            min = 0,
            max = 1,
            step = 0.01
          ),

          shiny::numericInput(
            ns("vis_min_hub"),
            label = "Min Hub Degree:",
            value = 2,
            min = 1,
            max = 20,
            step = 1
          ),

          shiny::hr(),

          shiny::h5("Physics & Interaction"),

          shiny::checkboxInput(
            ns("vis_physics"),
            label = "Enable Physics (Real-time Simulation)",
            value = TRUE
          ),

          shiny::conditionalPanel(
            condition = sprintf("input['%s'] == true", ns("vis_physics")),
            shiny::div(
              style = "background: #fff3e0; padding: 10px; border-radius: 5px; margin-top: 10px;",
              shiny::sliderInput(
                ns("vis_physics_speed"),
                label = "Simulation Speed:",
                min = 0.1,
                max = 2,
                value = 0.5,
                step = 0.1
              ),
              shiny::sliderInput(
                ns("vis_physics_repul"),
                label = "Repulsion Force:",
                min = -500,
                max = -50,
                value = -200,
                step = 10
              )
            )
          ),

          shiny::hr(),

          shiny::checkboxInput(
            ns("vis_nodes_draggable"),
            label = "Allow Node Dragging",
            value = TRUE
          ),

          shiny::checkboxInput(
            ns("vis_smooth_edges"),
            label = "Smooth Curved Edges",
            value = TRUE
          ),

          shiny::hr(),

          shiny::h5("Node Sizing"),

          shiny::selectInput(
            ns("vis_size_mode"),
            label = NULL,
            choices = c(
              "By Hub Degree" = "degree",
              "By -log10(FDR)" = "fdr"
            ),
            selected = "degree"
          ),

          shiny::hr(),

          shiny::h5("Layout"),

          shiny::selectInput(
            ns("vis_layout"),
            label = "Initial Layout:",
            choices = c(
              "Force-Directed (FR)" = "force",
              "Hierarchical" = "hierarchical",
              "LBHU (Large Graph)" = "lbhu"
            ),
            selected = "force"
          ),

          shiny::numericInput(
            ns("vis_seed"),
            label = "Seed:",
            value = 42,
            min = 1,
            max = 9999,
            step = 1
          ),

          shiny::hr(),

          shiny::actionButton(
            ns("vis_refresh"),
            label = "Refresh Network",
            icon = shiny::icon("refresh"),
            class = "btn-primary",
            style = "width: 100%;"
          ),

          shiny::actionButton(
            ns("vis_fix_layout"),
            label = "Fix Current Layout",
            icon = shiny::icon("lock"),
            class = "btn-warning",
            style = "width: 100%; margin-top: 5px;"
          ),

          shiny::actionButton(
            ns("vis_freeze"),
            label = "Stop Physics",
            icon = shiny::icon("pause"),
            class = "btn-danger",
            style = "width: 100%; margin-top: 5px;"
          ),

          shiny::hr(),

          shiny::h5("Export"),

          shiny::downloadButton(
            ns("vis_export_png"),
            label = "Export PNG",
            class = "btn-success",
            style = "width: 100%;"
          ),

          shiny::downloadButton(
            ns("vis_export_json"),
            label = "Export JSON",
            class = "btn-info",
            style = "width: 100%; margin-top: 5px;"
          ),

          shiny::hr(),

          shiny::h5("Network Statistics"),

          shiny::verbatimTextOutput(ns("vis_stats")) %>%
            shiny::tagAppendAttributes(style = "font-size: 11px; max-height: 150px; overflow-y: auto;"),

          shiny::hr(),

          shiny::h5("Pathways Preview"),

          shiny::uiOutput(ns("vis_pathway_list")) %>%
            shiny::tagAppendAttributes(style = "max-height: 200px; overflow-y: auto;")

        )
      ),

      # ─── 主网络绘图区域 ───
      shiny::column(
        width = 9,
        shiny::div(
          class = "white-box",
          style = "min-height: 800px;",

          shiny::div(
            id = ns("vis_status_bar"),
            style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); padding: 10px; border-radius: 5px; margin-bottom: 15px; color: white;",
            shiny::uiOutput(ns("vis_status_text"))
          ),

          visNetwork::visNetworkOutput(
            ns("vis_network"),
            height = "750px",
            width = "100%"
          ) %>%
            shinycssloaders::withSpinner(type = 6, color = "#28a745"),

          shiny::hr(),

          # ─── 调试信息窗口 ───
          shiny::div(
            class = "white-box",
            style = "background: #1e1e1e; color: #d4d4d4; padding: 15px; border-radius: 8px;",
            shiny::h4(
              style = "color: #569cd6; margin-top: 0;",
              "Debug Console"
            ),
            shiny::div(
              id = ns("vis_debug_container"),
              style = "max-height: 200px; overflow-y: auto; font-family: 'Courier New', monospace; font-size: 12px;",
              shiny::verbatimTextOutput(ns("vis_debug")) %>%
                shiny::tagAppendAttributes(style = "background: transparent; border: none; color: #9cdcfe;")
            )
          ),

          # ─── 节点/边详情面板 ───
          shiny::fluidRow(
            shiny::column(
              6,
              shiny::div(
                class = "white-box",
                shiny::h4("Node Details"),
                DT::dataTableOutput(ns("vis_node_detail")) %>%
                  shinycssloaders::withSpinner(type = 4, color = "#28a745")
              )
            ),
            shiny::column(
              6,
              shiny::div(
                class = "white-box",
                shiny::h4("Edge Details"),
                DT::dataTableOutput(ns("vis_edge_detail")) %>%
                  shinycssloaders::withSpinner(type = 4, color = "#28a745")
              )
            )
          )

        )
      )
    )
  )
}


#' @title HubGene Network (visNetwork) Module Server
#' @keywords internal

#' @title HubGene Network (visNetwork) Module Server
#' @keywords internal

mod_hubgene_vis_server <- function(id, data_prep_list, table_controller) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ──────────────────────────────────────────────────────────────
    # 调试日志
    # ──────────────────────────────────────────────────────────────
    .debug_log <<- character(0)

    add_debug <- function(msg, level = "INFO") {
      timestamp <- format(Sys.time(), "%H:%M:%S")
      new_entry <- sprintf("[%s] [%s] %s", timestamp, level, msg)
      .debug_log <<- c(new_entry, .debug_log)[1:100]
      message(sprintf("[HubGene-vis] %s", new_entry))
    }

    output$vis_debug <- shiny::renderText({
      paste(.debug_log, collapse = "\n")
    })

    add_debug("Module initialized")

    # ──────────────────────────────────────────────────────────────
    # 1. 模式状态
    # ──────────────────────────────────────────────────────────────

    vis_mode <- shiny::reactiveVal("mode_topN")

    shiny::observeEvent(input$vis_mode, {
      vis_mode(input$vis_mode)
      add_debug(sprintf("Mode: %s", input$vis_mode))
    })

    # ──────────────────────────────────────────────────────────────
    # 2. 数据源
    # ──────────────────────────────────────────────────────────────

    topN_candidates <- shiny::reactive({
      if (vis_mode() != "mode_topN") return(character(0))
      data_list <- data_prep_list$data()
      shiny::req(data_list)
      df <- data_list$df
      top_n <- input$vis_topN %||% 5
      top_n <- max(3, min(top_n, nrow(df)))
      df[1:top_n, "ID"]
    })

    select_candidates <- shiny::reactive({
      if (vis_mode() != "mode_select") return(character(0))
      sel <- table_controller$selected_pathways()
      if (is.null(sel) || length(sel) == 0) return(character(0))
      return(sel)
    })

    candidate_raw <- shiny::reactive({
      switch(vis_mode(),
             "mode_topN" = topN_candidates(),
             "mode_select" = select_candidates(),
             character(0))
    })

    candidate_filtered <- shiny::reactive({
      pathways <- candidate_raw()
      if (length(pathways) == 0) return(character(0))
      fdr_thresh <- input$vis_fdr %||% 0.25
      data_list <- data_prep_list$data()
      shiny::req(data_list)
      df <- data_list$df
      fdr_vec <- df$p.adjust[match(pathways, df$ID)]
      names(fdr_vec) <- pathways
      pathways[!is.na(fdr_vec) & fdr_vec < fdr_thresh]
    })

    final_pathways <- shiny::reactiveVal(character(0))

    shiny::observeEvent(candidate_filtered(), {
      final_pathways(candidate_filtered())
    })

    # ──────────────────────────────────────────────────────────────
    # 3. 构建网络数据
    # ──────────────────────────────────────────────────────────────

    network_data <- shiny::reactive({
      data_list <- data_prep_list$data()
      shiny::req(!is.null(data_list))
      pathway_ids <- final_pathways()
      shiny::req(length(pathway_ids) > 0)

      add_debug(sprintf("Building network: %d pathways", length(pathway_ids)))

      de_df <- tryCatch({
        get_de_table(data_prep_list$gsea_res, data_list$contrast_id)
      }, error = function(e) NULL)

      net <- build_hubgene_network(
        gsea_task = list(
          gsea_res = data_list$gsea_res,
          meta = list(
            left_group = data_list$left_group,
            right_group = data_list$right_group
          )
        ),
        pathway_ids = pathway_ids,
        min_hub_degree = input$vis_min_hub %||% 2,
        de_df = de_df,
        res_df = data_list$df
      )

      add_debug(sprintf("Network built: nodes=%s, edges=%s",
                        ifelse(is.null(net$nodes), "NULL",
                               paste(nrow(net$nodes$pathway), "+", nrow(net$nodes$gene))),
                        ifelse(is.null(net$edges), "NULL", nrow(net$edges))))
      return(net)
    })

    # ──────────────────────────────────────────────────────────────
    # 4. 准备 visNetwork 数据（简化版）
    # ──────────────────────────────────────────────────────────────

    vis_data <- shiny::reactive({
      net <- network_data()
      shiny::req(!is.null(net), !is.null(net$nodes))

      add_debug("Preparing visNetwork nodes and edges...")

      # ── 节点 ──
      nodes_list <- list()

      # 通路节点
      if (!is.null(net$nodes$pathway) && nrow(net$nodes$pathway) > 0) {
        for (i in seq_len(nrow(net$nodes$pathway))) {
          row <- net$nodes$pathway[i, ]
          nodes_list[[length(nodes_list) + 1]] <- list(
            id = paste0("pw_", row$id),
            label = gsub("HALLMARK_", "", row$id),
            group = "pathway",
            shape = "diamond",
            color = ifelse(row$NES > 0, "#E41A1C", "#377EB8"),
            size = 30,
            title = sprintf("NES: %.2f<br>FDR: %.2e", row$NES, row$FDR)
          )
        }
      }

      # 基因节点
      if (!is.null(net$nodes$gene) && nrow(net$nodes$gene) > 0) {
        for (i in seq_len(nrow(net$nodes$gene))) {
          row <- net$nodes$gene[i, ]
          nodes_list[[length(nodes_list) + 1]] <- list(
            id = paste0("gene_", row$id),
            label = row$id,
            group = "gene",
            shape = "dot",
            color = ifelse(row$stat > 0, "#E41A1C", "#377EB8"),
            size = 15 + row$degree * 2,
            title = sprintf("Degree: %d", row$degree)
          )
        }
      }

      nodes_df <- shiny::tagList(nodes_list)

      # ── 边 ──
      edges_list <- list()

      if (!is.null(net$edges) && nrow(net$edges) > 0) {
        for (i in seq_len(nrow(net$edges))) {
          row <- net$edges[i, ]
          edges_list[[length(edges_list) + 1]] <- list(
            from = paste0("gene_", row$source),
            to = paste0("pw_", row$target),
            color = ifelse(row$is_leading_edge, "#333333", "#AAAAAA"),
            width = ifelse(row$is_leading_edge, 3, 1),
            dashes = !row$is_leading_edge
          )
        }
      }

      edges_df <- shiny::tagList(edges_list)

      add_debug(sprintf("visData ready: %d nodes, %d edges",
                        length(nodes_list), length(edges_list)))

      list(nodes = nodes_list, edges = edges_list)
    })

    # ──────────────────────────────────────────────────────────────
    # 5. 渲染 visNetwork
    # ──────────────────────────────────────────────────────────────

    output$vis_network <- visNetwork::renderVisNetwork({
      vd <- vis_data()
      shiny::req(length(vd$nodes) > 0)

      add_debug("Rendering...")

      # 转换列表为数据框
      nodes_df <- as.data.frame(do.call(rbind, lapply(vd$nodes, function(x) {
        as.data.frame(x, stringsAsFactors = FALSE)
      })), stringsAsFactors = FALSE)

      edges_df <- as.data.frame(do.call(rbind, lapply(vd$edges, function(x) {
        as.data.frame(x, stringsAsFactors = FALSE)
      })), stringsAsFactors = FALSE)

      add_debug(sprintf("DataFrames: %d nodes, %d edges", nrow(nodes_df), nrow(edges_df)))

      # 构建网络
      vis <- visNetwork::visNetwork(nodes_df, edges_df, width = "100%", height = "700px")

      # 物理引擎
      if (input$vis_physics %||% TRUE) {
        vis <- vis %>% visNetwork::visPhysics(
          enabled = TRUE,
          barnesHut = list(
            gravitationalConstant = input$vis_physics_repul %||% -200,
            centralGravity = 0.01,
            springLength = 150,
            springConstant = 0.01
          ),
          stabilization = list(enabled = TRUE, iterations = 200),
          solver = "barnesHut"
        )
      } else {
        vis <- vis %>% visNetwork::visPhysics(enabled = FALSE)
      }

      # 布局
      vis <- vis %>% visNetwork::visLayout(randomSeed = input$vis_seed %||% 42)

      # 交互
      vis <- vis %>% visNetwork::visInteraction(
        dragNodes = TRUE,
        dragView = TRUE,
        zoomView = TRUE,
        hover = TRUE
      )

      # 图例
      vis <- vis %>% visNetwork::visLegend(
        enabled = TRUE,
        position = "right",
        addNodes = list(
          list(label = "Pathway", shape = "diamond", color = "#E41A1C"),
          list(label = "Gene", shape = "dot", color = "#377EB8")
        )
      )

      add_debug("Render complete")
      vis
    })

    # ──────────────────────────────────────────────────────────────
    # 6. 事件
    # ──────────────────────────────────────────────────────────────

    shiny::observeEvent(input$vis_click, {
      click_data <- input$vis_click
      shiny::req(click_data$nodeId)
      add_debug(sprintf("Click: %s", click_data$nodeId), "ACTION")
    })

    shiny::observeEvent(input$vis_refresh, {
      add_debug("Refresh", "ACTION")
    })

    # ──────────────────────────────────────────────────────────────
    # 7. 统计
    # ──────────────────────────────────────────────────────────────

    output$vis_stats <- shiny::renderPrint({
      vd <- tryCatch(vis_data(), error = function(e) NULL)
      if (is.null(vd)) {
        cat("No data\n")
        return()
      }
      pw_count <- sum(sapply(vd$nodes, function(x) x$group == "pathway"))
      gene_count <- sum(sapply(vd$nodes, function(x) x$group == "gene"))
      cat(sprintf("Pathways: %d\n", pw_count))
      cat(sprintf("Genes: %d\n", gene_count))
      cat(sprintf("Edges: %d\n", length(vd$edges)))
    })

    # ──────────────────────────────────────────────────────────────
    # 8. 状态栏
    # ──────────────────────────────────────────────────────────────

    output$vis_status_text <- shiny::renderUI({
      vd <- tryCatch(vis_data(), error = function(e) NULL)
      if (is.null(vd) || length(vd$nodes) == 0) {
        return(shiny::span("No network data", style = "color: #ffcccc;"))
      }
      pw_count <- sum(sapply(vd$nodes, function(x) x$group == "pathway"))
      gene_count <- sum(sapply(vd$nodes, function(x) x$group == "gene"))
      shiny::tagList(
        shiny::strong(sprintf("HubGene Network | %d Pathways + %d Genes",
                              pw_count, gene_count)),
        htmltools::tags$br(),
        htmltools::tags$small(sprintf("Physics: %s | FDR < %.2f",
                                      ifelse(input$vis_physics, "ON", "OFF"),
                                      input$vis_fdr %||% 0.25))
      )
    })

    # ──────────────────────────────────────────────────────────────
    # 9. 通路预览
    # ──────────────────────────────────────────────────────────────

    output$vis_pathway_list <- shiny::renderUI({
      pathways <- final_pathways()
      if (length(pathways) == 0) {
        return(shiny::div("No pathways", style = "color: #856404;"))
      }
      data_list <- data_prep_list$data()
      df <- data_list$df

      lapply(pathways, function(pid) {
        row_idx <- which(df$ID == pid)
        if (length(row_idx) == 0) return(NULL)
        row <- df[row_idx[1], ]
        shiny::div(
          style = "background: #f8f9fa; padding: 5px; margin-bottom: 3px; border-left: 3px solid #007bff; font-size: 11px;",
          shiny::strong(gsub("HALLMARK_", "", pid)),
          shiny::br(),
          sprintf("NES: %.2f | FDR: %.2e", as.numeric(row$NES), as.numeric(row$p.adjust))
        )
      })
    })

    # ──────────────────────────────────────────────────────────────
    # 10. 节点/边详情
    # ──────────────────────────────────────────────────────────────

    output$vis_node_detail <- DT::renderDataTable({
      vd <- tryCatch(vis_data(), error = function(e) NULL)
      if (is.null(vd) || length(vd$nodes) == 0) {
        return(DT::datatable(data.frame(Message = "No nodes"), options = list(dom = "t")))
      }

      display_df <- do.call(rbind, lapply(vd$nodes, function(x) {
        data.frame(
          ID = x$label,
          Type = x$group,
          Shape = x$shape,
          Color = x$color,
          Size = x$size,
          stringsAsFactors = FALSE
        )
      }))

      DT::datatable(display_df, rownames = FALSE, options = list(pageLength = 10, dom = "t"))
    })

    output$vis_edge_detail <- DT::renderDataTable({
      vd <- tryCatch(vis_data(), error = function(e) NULL)
      if (is.null(vd) || length(vd$edges) == 0) {
        return(DT::datatable(data.frame(Message = "No edges"), options = list(dom = "t")))
      }

      display_df <- do.call(rbind, lapply(vd$edges, function(x) {
        data.frame(
          From = gsub("gene_", "", x$from),
          To = gsub("pw_", "", x$to),
          Leading = ifelse(x$dashes, "No", "Yes"),
          Width = x$width,
          stringsAsFactors = FALSE
        )
      }))

      DT::datatable(display_df, rownames = FALSE, options = list(pageLength = 10, dom = "t"))
    })

    # ──────────────────────────────────────────────────────────────
    # 11. 导出
    # ──────────────────────────────────────────────────────────────

    output$vis_export_png <- shiny::downloadHandler(
      filename = function() paste0("hubgene_", Sys.time(), ".png"),
      content = function(file) {
        vis <- visNetwork::visNetwork(
          as.data.frame(do.call(rbind, lapply(vis_data()$nodes, as.data.frame,list))),
          as.data.frame(do.call(rbind, lapply(vis_data()$edges, as.data.frame,list)))
        )
        visNetwork::visSave(vis, file = file, type = "png", background = "white")
      }
    )

    add_debug("Server ready")
    return(list(final_pathways = final_pathways))

  })
}


# ==============================================================================
# END OF FILE
# ==============================================================================
