# ==============================================================================
# 文件：R/18_shiny_mod_hubgene_vis.R
# 功能：HubGene Network - visNetwork 交互版（精简版）
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

          # ── 模式选择 ──
          shiny::h5("Pathway Source Mode"),
          shiny::radioButtons(
            ns("vis_mode"),
            label = NULL,
            choices = c(
              "Top N" = "mode_topN",
              "Selected from Table" = "mode_select"
            ),
            selected = "mode_topN",
            width = "100%"
          ),

          shiny::hr(),

          # ── Top N 配置 ──
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] == 'mode_topN'", ns("vis_mode")),
            shiny::numericInput(
              ns("vis_topN"),
              label = "Top N Count:",
              value = 5,
              min = 3,
              max = 50,
              step = 1
            )
          ),

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

          # ── 物理引擎控制 ──
          shiny::h5("Physics & Interaction"),

          shiny::checkboxInput(
            ns("vis_physics"),
            label = "Enable Physics",
            value = TRUE
          ),

          shiny::conditionalPanel(
            condition = sprintf("input['%s'] == true", ns("vis_physics")),

            shiny::div(
              style = "background: #f8f9fa; padding: 10px; border-radius: 5px; margin-top: 10px;",

              shiny::selectInput(
                ns("vis_solver"),
                label = "Physics Solver:",
                choices = c(
                  "Barnes-Hut (Fast)" = "barnesHut",
                  "Repulsion" = "repulsion",
                  "Hierarchical Repulsion" = "hierarchicalRepulsion"
                ),
                selected = "barnesHut"
              ),

              shiny::sliderInput(
                ns("vis_gravitational"),
                label = "Gravitational Constant:",
                min = -2000,
                max = -50,
                value = -400,
                step = 10,
                post = ""
              ),
              shiny::helpText(
                style = "font-size: 10px; color: #666;",
                "越负值，节点间排斥越强"
              ),

              shiny::sliderInput(
                ns("vis_central_gravity"),
                label = "Central Gravity:",
                min = 0.01,
                max = 0.5,
                value = 0.3,
                step = 0.01
              ),
              shiny::helpText(
                style = "font-size: 10px; color: #666;",
                "将节点拉向中心，防止分散"
              ),

              shiny::sliderInput(
                ns("vis_spring_length"),
                label = "Spring Length:",
                min = 25,
                max = 500,
                value = 150,
                step = 5
              ),
              shiny::helpText(
                style = "font-size: 10px; color: #666;",
                "边的理想长度（像素）"
              ),

              shiny::sliderInput(
                ns("vis_spring_constant"),
                label = "Spring Constant:",
                min = 0.001,
                max = 1.0,
                value = 0.01,
                step = 0.001
              ),
              shiny::helpText(
                style = "font-size: 10px; color: #666;",
                "边的弹性系数"
              ),

              shiny::sliderInput(
                ns("vis_damping"),
                label = "Damping:",
                min = 0.01,
                max = 0.9,
                value = 0.09,
                step = 0.01
              ),
              shiny::helpText(
                style = "font-size: 10px; color: #666;",
                "运动阻尼，防止震荡"
              )
            )
          ),

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

          # ── 节点大小 ──
          shiny::h5("Node Sizing"),

          shiny::sliderInput(
            ns("vis_pw_size"),
            label = "Pathway Node Size:",
            min = 15,
            max = 60,
            value = 30,
            step = 1
          ),
          shiny::helpText(
            style = "font-size: 10px; color: #666;",
            "节点越大，排斥越强"
          ),

          shiny::sliderInput(
            ns("vis_gene_size"),
            label = "Gene Node Size:",
            min = 5,
            max = 30,
            value = 12,
            step = 1
          ),

          shiny::hr(),

          # ── 随机种子 ──
          shiny::numericInput(
            ns("vis_seed"),
            label = "Random Seed:",
            value = 42,
            min = 1,
            max = 9999,
            step = 1
          ),

          shiny::hr(),

          # ── 统计 ──
          shiny::h5("Network Statistics"),
          shiny::verbatimTextOutput(ns("vis_stats")) %>%
            shiny::tagAppendAttributes(style = "font-size: 10px; max-height: 100px; overflow-y: auto;"),

          shiny::hr(),

          # ── 通路预览 ──
          shiny::h5("Pathways Preview"),
          shiny::uiOutput(ns("vis_pathway_list")) %>%
            shiny::tagAppendAttributes(style = "max-height: 150px; overflow-y: auto;")

        )
      ),

      # ─── 主绘图区域 ───
      shiny::column(
        width = 9,
        shiny::div(
          class = "white-box",
          style = "min-height: 850px;",

          # 状态栏
          shiny::div(
            style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); padding: 12px; border-radius: 5px; margin-bottom: 15px; color: white;",
            shiny::uiOutput(ns("vis_status_text"))
          ),

          # visNetwork 输出
          visNetwork::visNetworkOutput(
            ns("vis_network"),
            height = "650px",
            width = "100%"
          ) %>%
            shinycssloaders::withSpinner(type = 6, color = "#28a745"),

          shiny::hr(),

          # ── Debug Console ──
          shiny::div(
            style = "background: #1e1e1e; color: #d4d4d4; padding: 15px; border-radius: 8px; margin-top: 15px;",
            shiny::h4(
              style = "color: #569cd6; margin-top: 0; font-size: 14px;",
              "Debug Console"
            ),
            shiny::div(
              style = "max-height: 150px; overflow-y: auto; font-family: 'Courier New', monospace; font-size: 11px;",
              shiny::verbatimTextOutput(ns("vis_debug")) %>%
                shiny::tagAppendAttributes(style = "background: transparent; border: none; color: #9cdcfe;")
            )
          ),

          # ── 节点/边详情 ──
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


#' @title HubGene Network (visNetwork) Module Server (Simplified)
#' @keywords internal

mod_hubgene_vis_server <- function(id, data_prep_list, table_controller) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ──────────────────────────────────────────────────────────────
    # 调试日志（使用 session$userData）
    # ──────────────────────────────────────────────────────────────
    session$userData$debug_log <- character(0)

    add_debug <- function(msg, level = "INFO") {
      timestamp <- format(Sys.time(), "%H:%M:%S")
      new_entry <- sprintf("[%s] [%s] %s", timestamp, level, msg)
      session$userData$debug_log <- c(new_entry, session$userData$debug_log)[1:50]
      message(sprintf("[HubGene-vis] %s", new_entry))
    }

    output$vis_debug <- shiny::renderText({
      paste(session$userData$debug_log, collapse = "\n")
    })

    add_debug("Module started")

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
      top_n <- input$vis_topN
      if (is.null(top_n)) top_n <- 5
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
      fdr_thresh <- input$vis_fdr
      if (is.null(fdr_thresh)) fdr_thresh <- 0.25
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
      add_debug(sprintf("Pathways: %d", length(candidate_filtered())))
    })

    # ──────────────────────────────────────────────────────────────
    # 3. 网络数据
    # ──────────────────────────────────────────────────────────────

    net_data <- shiny::reactive({
      data_list <- data_prep_list$data()
      shiny::req(!is.null(data_list))
      pathway_ids <- final_pathways()
      shiny::req(length(pathway_ids) > 0)

      add_debug(sprintf("Building: %d pathways", length(pathway_ids)))

      net <- build_hubgene_network(
        gsea_task = list(
          gsea_res = data_list$gsea_res,
          meta = list(
            left_group = data_list$left_group,
            right_group = data_list$right_group
          )
        ),
        pathway_ids = pathway_ids,
        min_hub_degree = input$vis_min_hub,
        de_df = NULL,
        res_df = data_list$df
      )

      if (is.null(net)) {
        add_debug("Build failed", "ERROR")
        return(NULL)
      }

      add_debug(sprintf("Built: %d pw + %d gene + %d edges",
                        ifelse(is.null(net$nodes$pathway), 0, nrow(net$nodes$pathway)),
                        ifelse(is.null(net$nodes$gene), 0, nrow(net$nodes$gene)),
                        ifelse(is.null(net$edges), 0, nrow(net$edges))))
      return(net)
    })

    # ──────────────────────────────────────────────────────────────
    # 4. 渲染网络
    # ──────────────────────────────────────────────────────────────

    output$vis_network <- visNetwork::renderVisNetwork({
      net <- net_data()
      shiny::req(!is.null(net))

      add_debug("Preparing visNetwork...")

      # 获取分组信息用于 tooltip
      data_list <- data_prep_list$data()
      left_group <- data_list$left_group
      right_group <- data_list$right_group

      # 节点（添加 title 列用于悬停）
      nodes <- data.frame(
        id = character(),
        label = character(),
        group = character(),
        shape = character(),
        color = character(),
        size = numeric(),
        title = character(),
        stringsAsFactors = FALSE
      )

      pw_size <- ifelse(is.null(input$vis_pw_size), 30, input$vis_pw_size)
      gene_size <- ifelse(is.null(input$vis_gene_size), 12, input$vis_gene_size)

      # 通路节点 tooltip（添加 leading edge 基因数计算）
      if (!is.null(net$nodes$pathway) && nrow(net$nodes$pathway) > 0) {
        for (i in seq_len(nrow(net$nodes$pathway))) {
          row <- net$nodes$pathway[i, ]
          direction <- ifelse(row$NES > 0, left_group, right_group)

          # 计算该通路的 leading edge 连接数
          if (!is.null(net$edges) && nrow(net$edges) > 0) {
            le_count <- sum(net$edges$target == row$id & net$edges$is_leading_edge, na.rm = TRUE)
          } else {
            le_count <- 0
          }

          nodes <- rbind(nodes, data.frame(
            id = paste0("pw_", row$id),
            label = gsub("HALLMARK_", "", row$id),
            group = "pathway",
            shape = "diamond",
            color = ifelse(row$NES > 0, "#E41A1C", "#377EB8"),
            size = pw_size,
            title = sprintf(
              "Pathway: %s\nDirection: %s\nNES: %.3f\nFDR: %.2e\nCore (LE): %d",
              row$id, direction, row$NES, row$FDR, le_count
            ),
            stringsAsFactors = FALSE
          ))
        }
      }

      # 基因节点 tooltip（添加 LE 状态和连接的通路）
      if (!is.null(net$nodes$gene) && nrow(net$nodes$gene) > 0) {
        for (i in seq_len(nrow(net$nodes$gene))) {
          row <- net$nodes$gene[i, ]
          gene_label <- row$id  # 保持原始大小写

          if (row$stat > 0) {
            direction <- sprintf("Up in %s", left_group)
          } else if (row$stat < 0) {
            direction <- sprintf("Up in %s", right_group)
          } else {
            direction <- "Neutral"
          }

          # 判断是否为 Leading Edge（至少有一条 leading edge 边）
          if (!is.null(net$edges) && nrow(net$edges) > 0) {
            is_le <- any(net$edges$source == row$id & net$edges$is_leading_edge, na.rm = TRUE)
            # 获取连接的通路列表
            connected_pws <- unique(net$edges$target[net$edges$source == row$id])
            if (length(connected_pws) > 0) {
              pw_labels <- gsub("HALLMARK_", "", connected_pws)
              if (length(pw_labels) <= 5) {
                pw_str <- paste(pw_labels, collapse = ", ")
              } else {
                pw_str <- paste(pw_labels[1:5], collapse = ", ")
                pw_str <- paste0(pw_str, sprintf(" (+%d)", length(pw_labels) - 5))
              }
            } else {
              pw_str <- "None"
            }
          } else {
            is_le <- FALSE
            pw_str <- "None"
          }

          nodes <- rbind(nodes, data.frame(
            id = paste0("gene_", gene_label),
            label = gene_label,
            group = "gene",
            shape = "dot",
            color = ifelse(row$stat > 0, "#E41A1C", "#377EB8"),
            size = gene_size + row$degree * 2,
            title = sprintf(
              "Gene: %s\nDirection: %s\nStat: %.3f\nHub Degree: %d\nLE: %s\nConnected: %s",
              gene_label, direction, row$stat, row$degree,
              ifelse(is_le, "YES", "NO"), pw_str
            ),
            stringsAsFactors = FALSE
          ))
        }
      }

      # 边
      edges <- data.frame(
        from = character(),
        to = character(),
        color = character(),
        width = numeric(),
        dashes = logical(),
        title = character(),
        stringsAsFactors = FALSE
      )

      if (!is.null(net$edges) && nrow(net$edges) > 0) {
        for (i in seq_len(nrow(net$edges))) {
          row <- net$edges[i, ]
          # row$source 和 row$target 保持原始大小写
          edges <- rbind(edges, data.frame(
            from = paste0("gene_", row$source),  # 保持原始大小写
            to = paste0("pw_", row$target),
            color = ifelse(row$is_leading_edge, "#333333", "#AAAAAA"),
            width = ifelse(row$is_leading_edge, 3, 1),
            dashes = !row$is_leading_edge,
            title = sprintf("Gene: %s -> Pathway: %s\nLeading Edge: %s",
                            row$source, row$target,
                            ifelse(row$is_leading_edge, "YES", "NO")),
            stringsAsFactors = FALSE
          ))
        }
      }

      add_debug(sprintf("Nodes: %d, Edges: %d", nrow(nodes), nrow(edges)))

      # 物理引擎
      physics_enabled <- isTRUE(input$vis_physics)
      solver <- ifelse(is.null(input$vis_solver), "barnesHut", input$vis_solver)
      grav <- ifelse(is.null(input$vis_gravitational), -400, input$vis_gravitational)

      vis <- visNetwork::visNetwork(nodes, edges, width = "100%", height = "650px")

      if (physics_enabled) {
        vis <- vis %>% visNetwork::visPhysics(
          enabled = TRUE,
          solver = solver,
          barnesHut = list(
            gravitationalConstant = grav,
            centralGravity = 0.3,
            springLength = 150,
            springConstant = 0.01,
            damping = 0.09
          )
        )
        add_debug(sprintf("Physics: ON, solver=%s, grav=%.0f", solver, grav))
      } else {
        vis <- vis %>% visNetwork::visPhysics(enabled = FALSE)
        add_debug("Physics: OFF")
      }

      # 布局
      seed <- ifelse(is.null(input$vis_seed), 42, input$vis_seed)
      vis <- vis %>% visNetwork::visLayout(randomSeed = seed)

      # 交互
      vis <- vis %>% visNetwork::visInteraction(
        dragNodes = TRUE,
        dragView = TRUE,
        zoomView = TRUE,
        hover = TRUE,
        tooltipDelay = 200
      )

      # 图例
      vis <- vis %>% visNetwork::visLegend(
        enabled = TRUE,
        position = "right",
        addNodes = list(
          list(label = "Pathway Up", shape = "diamond", color = "#E41A1C"),
          list(label = "Pathway Down", shape = "diamond", color = "#377EB8"),
          list(label = "Gene Up", shape = "dot", color = "#E41A1C"),
          list(label = "Gene Down", shape = "dot", color = "#377EB8")
        )
      )

      # 点击事件
      vis <- vis %>% visNetwork::visEvents(
        click = sprintf("function(props) {
      var n = props.nodes[0];
      if(n) Shiny.setInputValue('%s', {id: n, ts: Date.now()}, {priority:'event'});
    }", ns("vis_click"))
      )

      add_debug("Render complete")
      vis
    })

    # ──────────────────────────────────────────────────────────────
    # 5. 事件
    # ──────────────────────────────────────────────────────────────

    shiny::observeEvent(input$vis_click, {
      add_debug(sprintf("Click: %s", input$vis_click$id), "CLICK")
    })

    # ──────────────────────────────────────────────────────────────
    # 6. 统计
    # ──────────────────────────────────────────────────────────────

    output$vis_stats <- shiny::renderPrint({
      net <- tryCatch(net_data(), error = function(e) NULL)
      if (is.null(net)) {
        cat("No data\n")
        return()
      }
      pw <- ifelse(is.null(net$nodes$pathway), 0, nrow(net$nodes$pathway))
      gene <- ifelse(is.null(net$nodes$gene), 0, nrow(net$nodes$gene))
      edge <- ifelse(is.null(net$edges), 0, nrow(net$edges))
      cat(sprintf("Pathways: %d\n", pw))
      cat(sprintf("Hub Genes: %d\n", gene))
      cat(sprintf("Edges: %d\n", edge))
    })

    # ──────────────────────────────────────────────────────────────
    # 7. 状态栏
    # ──────────────────────────────────────────────────────────────

    output$vis_status_text <- shiny::renderUI({
      net <- tryCatch(net_data(), error = function(e) NULL)
      if (is.null(net)) {
        return(shiny::span("No data", style = "color: #ffcccc;"))
      }
      pw <- ifelse(is.null(net$nodes$pathway), 0, nrow(net$nodes$pathway))
      gene <- ifelse(is.null(net$nodes$gene), 0, nrow(net$nodes$gene))
      edge <- ifelse(is.null(net$edges), 0, nrow(net$edges))
      physics <- if (isTRUE(input$vis_physics)) "ON" else "OFF"

      shiny::tagList(
        shiny::strong(sprintf("HubGene Network | %d Pathways + %d Hub Genes + %d Edges", pw, gene, edge)),
        htmltools::tags$br(),
        htmltools::tags$small(sprintf("Physics: %s | FDR < %.2f", physics, input$vis_fdr %||% 0.25))
      )
    })

    # ──────────────────────────────────────────────────────────────
    # 8. 通路预览
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
          style = "background: #f8f9fa; padding: 5px; margin-bottom: 3px; border-left: 3px solid #007bff; font-size: 10px;",
          shiny::strong(gsub("HALLMARK_", "", pid)),
          htmltools::tags$br(),
          sprintf("NES: %.2f", as.numeric(row$NES))
        )
      })
    })

    # ──────────────────────────────────────────────────────────────
    # 9. 详情表
    # ──────────────────────────────────────────────────────────────

    output$vis_node_detail <- DT::renderDataTable({
      net <- tryCatch(net_data(), error = function(e) NULL)
      if (is.null(net) || is.null(net$nodes$pathway) && is.null(net$nodes$gene)) {
        return(DT::datatable(data.frame(Message = "No nodes"), options = list(dom = "t")))
      }

      rows <- list()
      if (!is.null(net$nodes$pathway) && nrow(net$nodes$pathway) > 0) {
        for (i in seq_len(nrow(net$nodes$pathway))) {
          r <- net$nodes$pathway[i, ]
          rows[[length(rows) + 1]] <- data.frame(
            ID = r$id, Type = "Pathway",
            NES = round(r$NES, 3),
            Degree = NA,
            stringsAsFactors = FALSE
          )
        }
      }
      if (!is.null(net$nodes$gene) && nrow(net$nodes$gene) > 0) {
        for (i in seq_len(nrow(net$nodes$gene))) {
          r <- net$nodes$gene[i, ]
          rows[[length(rows) + 1]] <- data.frame(
            ID = r$id, Type = "Gene",
            NES = round(r$stat, 3),
            Degree = r$degree,
            stringsAsFactors = FALSE
          )
        }
      }

      df <- do.call(rbind, rows)
      DT::datatable(df, rownames = FALSE, options = list(pageLength = 10, dom = "t"))
    })

    output$vis_edge_detail <- DT::renderDataTable({
      net <- tryCatch(net_data(), error = function(e) NULL)
      if (is.null(net) || is.null(net$edges) || nrow(net$edges) == 0) {
        return(DT::datatable(data.frame(Message = "No edges"), options = list(dom = "t")))
      }

      rows <- lapply(seq_len(nrow(net$edges)), function(i) {
        r <- net$edges[i, ]
        data.frame(
          Gene = r$source,
          Pathway = r$target,
          LeadingEdge = ifelse(r$is_leading_edge, "Yes", "No"),
          stringsAsFactors = FALSE
        )
      })

      df <- do.call(rbind, rows)
      DT::datatable(df, rownames = FALSE, options = list(pageLength = 10, dom = "t"))
    })

    add_debug("Ready")
    return(list(final_pathways = final_pathways))

  })
}
