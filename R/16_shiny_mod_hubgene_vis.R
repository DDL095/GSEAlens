# ==============================================================================
# 文件：R/18_shiny_mod_hubgene_vis.R
# 功能：HubGene Network - visNetwork 交互版
# 美学调整：Gene透明度、Edge颜色跟随基因stat、highlightNearest
# 修复：添加 symbol_map 定义，恢复基因名原始大小写
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

              shiny::sliderInput(
                ns("vis_gravitational"),
                label = "Gravitational Constant:",
                min = -2000,
                max = -50,
                value = -400,
                step = 10
              ),

              shiny::sliderInput(
                ns("vis_spring_length"),
                label = "Spring Length:",
                min = 25,
                max = 500,
                value = 150,
                step = 5
              ),

              shiny::sliderInput(
                ns("vis_spring_constant"),
                label = "Spring Constant:",
                min = 0,
                max = 1.0,
                value = 0.001,
                step = 0.001
              ),
              shiny::sliderInput(
                ns("vis_central_gravity"),
                label = "Central Gravity:",
                min = 0,
                max = 1,
                value = 0.3,
                step = 0.05
              ),

              shiny::sliderInput(
                ns("vis_damping"),
                label = "Damping:",
                min = 0.01,
                max = 0.9,
                value = 0.09,
                step = 0.01
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

          shiny::sliderInput(
            ns("vis_gene_size"),
            label = "Gene Node Base Size:",
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
            height = "1200px",
            width = "100%"
          ) %>%
            shinycssloaders::withSpinner(type = 6, color = "#28a745"),

          shiny::hr()

        )
      )
    )
  )
}


#' @title HubGene Network (visNetwork) Module Server
#' @keywords internal

mod_hubgene_vis_server <- function(id, data_prep_list, table_controller, gsea_res) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ──────────────────────────────────────────────────────────────
    # 调试日志
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

      # 获取分组信息
      data_list <- data_prep_list$data()
      left_group <- data_list$left_group
      right_group <- data_list$right_group

      # ──────────────────────────────────────────────────────────────
      # 修复：重建基因名大小写映射
      # ──────────────────────────────────────────────────────────────
      symbol_map <- NULL
      if (!is.null(gsea_res) && !is.null(data_list$contrast_id)) {
        symbol_map <- tryCatch({
          .rebuild_symbol_map(gsea_res, data_list$contrast_id)
        }, error = function(e) {
          message("[HubGene-vis] symbol_map build failed: ", e$message)
          NULL
        })
        add_debug(sprintf("symbol_map: %d genes mapped", length(symbol_map)))
      } else {
        add_debug("symbol_map skipped: gsea_res or contrast_id is NULL", "WARN")
      }

      # 节点大小参数
      pw_size <- ifelse(is.null(input$vis_pw_size), 30, input$vis_pw_size)
      gene_size <- ifelse(is.null(input$vis_gene_size), 12, input$vis_gene_size)

      # ──────────────────────────────────────────────────────────────
      # 节点构建（添加opacity列）
      # ──────────────────────────────────────────────────────────────

      nodes <- data.frame(
        id = character(),
        label = character(),
        group = character(),
        shape = character(),
        color = character(),
        size = numeric(),
        opacity = numeric(),
        title = character(),
        stringsAsFactors = FALSE
      )

      # 通路节点
      if (!is.null(net$nodes$pathway) && nrow(net$nodes$pathway) > 0) {
        for (i in seq_len(nrow(net$nodes$pathway))) {
          row <- net$nodes$pathway[i, ]
          direction <- ifelse(row$NES > 0, left_group, right_group)

          # 计算LE连接数
          le_count <- 0
          if (!is.null(net$edges) && nrow(net$edges) > 0) {
            le_count <- sum(net$edges$target == row$id & net$edges$is_leading_edge, na.rm = TRUE)
          }

          nodes <- rbind(nodes, data.frame(
            id = paste0("pw_", row$id),
            label = gsub("^[^_]*_", "", row$id),
            group = "pathway",
            shape = "diamond",
            color = ifelse(row$NES > 0, "#E41A1C", "#377EB8"),
            size = pw_size,
            opacity = 1.0,
            title = sprintf(
              "<b>%s</b><br>Direction: %s<br>NES: %.3f<br>FDR: %.2e<br>LE Genes: %d",
              row$id, direction, row$NES, row$FDR, le_count
            ),
            stringsAsFactors = FALSE
          ))
        }
      }

      # 基因节点（添加opacity区分）
      if (!is.null(net$nodes$gene) && nrow(net$nodes$gene) > 0) {
        for (i in seq_len(nrow(net$nodes$gene))) {
          row <- net$nodes$gene[i, ]
          # 恢复原始大小写
          gene_label <- .get_display_symbol(row$id, symbol_map)

          # 判断是否为Leading Edge
          is_le <- FALSE
          pw_str <- "None"
          if (!is.null(net$edges) && nrow(net$edges) > 0) {
            is_le <- any(net$edges$source == row$id & net$edges$is_leading_edge, na.rm = TRUE)
            connected_pws <- unique(net$edges$target[net$edges$source == row$id])
            if (length(connected_pws) > 0) {
              pw_labels <- gsub("^[^_]*_", "", connected_pws)
              if (length(pw_labels) <= 5) {
                pw_str <- paste(pw_labels, collapse = ", ")
              } else {
                pw_str <- paste(pw_labels[1:5], collapse = ", ")
                pw_str <- paste0(pw_str, sprintf(" (+%d)", length(pw_labels) - 5))
              }
            }
          }

          # 方向
          if (row$stat > 0) {
            direction <- sprintf("Up in %s", left_group)
          } else if (row$stat < 0) {
            direction <- sprintf("Up in %s", right_group)
          } else {
            direction <- "Neutral"
          }

          # 透明度：LE=1.0, non-LE=0.5
          gene_opacity <- ifelse(is_le, 1.0, 0.5)

          nodes <- rbind(nodes, data.frame(
            id = paste0("gene_", gene_label),
            label = gene_label,
            group = "gene",
            shape = "dot",
            color = ifelse(row$stat > 0, "#E41A1C", ifelse(row$stat < 0, "#377EB8", "#999999")),
            size = gene_size + row$degree * 3,
            opacity = gene_opacity,
            title = sprintf(
              "<b>%s</b><br>Direction: %s<br>Stat: %.3f<br>Hub Degree: %d<br>LE: %s<br>Connected: %s",
              gene_label, direction, row$stat, row$degree,
              ifelse(is_le, "YES", "NO"), pw_str
            ),
            stringsAsFactors = FALSE
          ))
        }
      }

      # ──────────────────────────────────────────────────────────────
      # 边构建（颜色跟随基因stat）
      # ──────────────────────────────────────────────────────────────

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

          # 获取源基因的stat值来决定颜色
          gene_stat <- 0
          if (!is.null(net$nodes$gene) && nrow(net$nodes$gene) > 0) {
            gene_row <- net$nodes$gene[net$nodes$gene$id == row$source, ]
            if (nrow(gene_row) > 0) {
              gene_stat <- gene_row$stat[1]
            }
          }

          # 颜色跟随基因stat：浅色（带透明度）
          if (gene_stat > 0) {
            edge_color <- "#E41A1C80"
          } else if (gene_stat < 0) {
            edge_color <- "#377EB880"
          } else {
            edge_color <- "#99999980"
          }

          # 恢复原始大小写
          gene_display <- .get_display_symbol(row$source, symbol_map)

          edges <- rbind(edges, data.frame(
            from = paste0("gene_", gene_display),
            to = paste0("pw_", row$target),
            color = edge_color,
            width = ifelse(row$is_leading_edge, 3, 1.5),
            dashes = !row$is_leading_edge,
            title = sprintf(
              "<b>Gene:</b> %s<br><b>Pathway:</b> %s<br><b>Leading Edge:</b> %s",
              gene_display, row$target,
              ifelse(row$is_leading_edge, "YES", "NO")
            ),
            stringsAsFactors = FALSE
          ))
        }
      }

      add_debug(sprintf("Nodes: %d, Edges: %d", nrow(nodes), nrow(edges)))

      # ──────────────────────────────────────────────────────────────
      # 构建visNetwork对象
      # ──────────────────────────────────────────────────────────────

      vis <- visNetwork::visNetwork(nodes, edges, width = "100%", height = "650px")

      # 物理引擎
      physics_enabled <- isTRUE(input$vis_physics)
      grav <- ifelse(is.null(input$vis_gravitational), -400, input$vis_gravitational)

      spring_length   <- ifelse(is.null(input$vis_spring_length), 150,   input$vis_spring_length)
      spring_constant <- ifelse(is.null(input$vis_spring_constant), 0.01, input$vis_spring_constant)
      damping         <- ifelse(is.null(input$vis_damping), 0.09,          input$vis_damping)
      central_gravity <- ifelse(is.null(input$vis_central_gravity), 0.3, input$vis_central_gravity)

      if (physics_enabled) {
        vis <- vis %>% visNetwork::visPhysics(
          enabled = TRUE,
          solver = "barnesHut",
          barnesHut = list(
            gravitationalConstant = grav,
            centralGravity = central_gravity,
            springLength = spring_length,
            springConstant = spring_constant,
            damping = damping
          )
        )
        add_debug(sprintf("Physics: ON, solver=barnesHut, grav=%.0f", grav))
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
        tooltipDelay = 300,
        navigationButtons = TRUE
      )

      # 点击高亮相关节点和边
      vis <- vis %>% visNetwork::visOptions(
        highlightNearest = TRUE,
        nodesIdSelection = TRUE
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
      le_edges <- ifelse(is.null(net$edges), 0, sum(net$edges$is_leading_edge, na.rm = TRUE))
      cat(sprintf("Pathways: %d\n", pw))
      cat(sprintf("Hub Genes: %d\n", gene))
      cat(sprintf("Total Edges: %d\n", edge))
      cat(sprintf("LE Edges: %d\n", le_edges))
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
      le_edges <- ifelse(is.null(net$edges), 0, sum(net$edges$is_leading_edge, na.rm = TRUE))
      physics <- if (isTRUE(input$vis_physics)) "ON" else "OFF"

      shiny::tagList(
        shiny::strong(sprintf("HubGene Network | %d Pathways + %d Hub Genes + %d Edges (LE: %d)",
                              pw, gene, edge, le_edges)),
        htmltools::tags$br(),
        htmltools::tags$small(sprintf("Physics: %s (Barnes-Hut) | FDR < %.2f | tooltipDelay: 300ms ",
                                      physics, input$vis_fdr %||% 0.25))
      )
    })

    # ──────────────────────────────────────────────────────────────
    # 8. 通路预览
    # ──────────────────────────────────────────────────────────────

    output$vis_pathway_list <- shiny::renderUI({
      pathways <- final_pathways()
      net <- tryCatch(net_data(), error = function(e) NULL)

      if (length(pathways) == 0) {
        return(shiny::div("No pathways", style = "color: #856404;"))
      }
      data_list <- data_prep_list$data()
      df <- data_list$df

      lapply(pathways, function(pid) {
        row_idx <- which(df$ID == pid)
        if (length(row_idx) == 0) return(NULL)
        row <- df[row_idx[1], ]
        le_count <- 0
        if (!is.null(net) && !is.null(net$edges)) {
          le_count <- sum(net$edges$target == pid & net$edges$is_leading_edge, na.rm = TRUE)
        }
        shiny::div(
          style = "background: #f8f9fa; padding: 5px; margin-bottom: 3px; border-left: 3px solid #007bff; font-size: 10px;",
          shiny::strong(gsub("^[^_]*_", "", pid)),
          htmltools::tags$br(),
          sprintf("NES: %.2f | LE: %d", as.numeric(row$NES), le_count)
        )
      })
    })

    # ──────────────────────────────────────────────────────────────
    # 9. 详情表
    # ──────────────────────────────────────────────────────────────

    output$vis_node_detail <- DT::renderDataTable({
      net <- tryCatch(net_data(), error = function(e) NULL)
      if (is.null(net) || (is.null(net$nodes$pathway) && is.null(net$nodes$gene))) {
        return(DT::datatable(data.frame(Message = "No nodes"), options = list(dom = "t")))
      }

      # 获取 symbol_map 用于显示原始大小写
      symbol_map <- NULL
      data_list <- tryCatch(data_prep_list$data(), error = function(e) NULL)
      if (!is.null(data_list) && !is.null(data_prep_list$gsea_res) && !is.null(data_list$contrast_id)) {
        symbol_map <- tryCatch({
          .rebuild_symbol_map(data_prep_list$gsea_res, data_list$contrast_id)
        }, error = function(e) NULL)
      }

      rows <- list()

      # Pathway节点
      if (!is.null(net$nodes$pathway) && nrow(net$nodes$pathway) > 0) {
        for (i in seq_len(nrow(net$nodes$pathway))) {
          r <- net$nodes$pathway[i, ]
          le_count <- sum(net$edges$target == r$id & net$edges$is_leading_edge, na.rm = TRUE)
          total_count <- sum(net$edges$target == r$id, na.rm = TRUE)
          rows[[length(rows) + 1]] <- data.frame(
            ID = r$id,
            Type = "Pathway",
            NES = round(r$NES, 3),
            Direction = ifelse(r$NES > 0, "Up", "Down"),
            Degree = NA_integer_,
            LE = le_count,
            Total = total_count,
            stringsAsFactors = FALSE
          )
        }
      }

      # Gene节点
      if (!is.null(net$nodes$gene) && nrow(net$nodes$gene) > 0) {
        for (i in seq_len(nrow(net$nodes$gene))) {
          r <- net$nodes$gene[i, ]
          is_le <- any(net$edges$source == r$id & net$edges$is_leading_edge, na.rm = TRUE)
          gene_display <- .get_display_symbol(r$id, symbol_map)
          rows[[length(rows) + 1]] <- data.frame(
            ID = gene_display,
            Type = "Gene",
            NES = round(r$stat, 3),
            Direction = ifelse(r$stat > 0, "Up", ifelse(r$stat < 0, "Down", "Neutral")),
            Degree = r$degree,
            LE = ifelse(is_le, "YES", "NO"),
            Total = NA_integer_,
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

      # 获取 symbol_map 用于显示原始大小写
      symbol_map <- NULL
      data_list <- tryCatch(data_prep_list$data(), error = function(e) NULL)
      if (!is.null(data_list) && !is.null(data_prep_list$gsea_res) && !is.null(data_list$contrast_id)) {
        symbol_map <- tryCatch({
          .rebuild_symbol_map(data_prep_list$gsea_res, data_list$contrast_id)
        }, error = function(e) NULL)
      }

      rows <- lapply(seq_len(nrow(net$edges)), function(i) {
        r <- net$edges[i, ]
        gene_display <- .get_display_symbol(r$source, symbol_map)
        data.frame(
          Gene = gene_display,
          Pathway = r$target,
          LeadingEdge = ifelse(r$is_leading_edge, "YES", "NO"),
          Width = ifelse(r$is_leading_edge, 3, 1.5),
          Style = ifelse(r$is_leading_edge, "Solid", "Dashed"),
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
