# =============================================================================
# HubGene Network Module
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

          # 通路数量设置
          shiny::h4("Pathway Settings"),
          shiny::numericInput(
            ns("default_n_pathways"),
            label = "Number of pathways:",
            value = 5,
            min = 3,
            max = 999,
            step = 1
          ),
          shiny::helpText(
            style = "color: #666; font-size: 11px;",
            "Top N pathways by |NES| to display"
          ),
          shiny::div(
            style = "background: #e8f4fd; padding: 10px; border-radius: 5px; margin-top: 10px;",
            shiny::textOutput(ns("selection_summary"), inline = TRUE)
          ),

          shiny::hr(),

          # Hub 基因阈值
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

          # 节点大小
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

          # 布局算法
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

          # 显示选项
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

          # 统计信息
          shiny::h4("Statistics"),
          shiny::verbatimTextOutput(ns("network_stats"))
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


#' @title HubGene Network Module Server
#' @description Server logic for HubGene Network
#' @param id Module ID
#' @param data_prep_list Reactive data from data prep module
#' @param table_controller Table controller with selected_pathways
#' @keywords internal

mod_hubgene_server <- function(id, data_prep_list, table_controller) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ─── Reactive: 获取通路列表 ───
    selected_pathways <- shiny::reactive({
      sel <- table_controller$selected_pathways()

      if (!is.null(sel) && length(sel) > 0) {
        return(sel)
      }

      data_list <- data_prep_list$data()
      shiny::req(!is.null(data_list))

      n <- input$default_n_pathways %||% 12
      n <- max(3, min(n, nrow(data_list$df)))

      return(data_list$df$ID[1:n])
    })

    # ─── Reactive: 构建网络数据 ───
    network_data <- shiny::reactive({
      data_list <- data_prep_list$data()
      shiny::req(!is.null(data_list))

      pathway_ids <- selected_pathways()
      shiny::req(length(pathway_ids) > 0)

      de_df <- tryCatch({
        get_de_table(data_list$gsea_res, data_list$contrast_id)
      }, error = function(e) NULL)

      res_df <- data_list$df

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

    # ─── Reactive: 准备节点数据 ───
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

      # ─── 通路着色（NES必然非零）───
      pathway_nodes$color <- ifelse(pathway_nodes$NES > 0, "#E41A1C", "#377EB8")

      # ─── 基因着色（log2FC可能为零）───
      gene_nodes$color <- ifelse(
        is.na(gene_nodes$stat), "#999999",
        ifelse(gene_nodes$stat > 0, "#E41A1C",
               ifelse(gene_nodes$stat < 0, "#377EB8", "#999999"))
      )

      # ─── 通路节点大小 ───
      if (size_mode == "fdr") {
        pathway_nodes$size <- pmin(pmax(-log10(pathway_nodes$FDR + 1e-100) * 5 + 18, 18), 45)
      } else {
        pathway_nodes$size <- pmin(pmax(-log10(pathway_nodes$FDR + 1e-100) * 5 + 18, 18), 45)
      }

      # ─── 基因节点大小（按度数）───
      if (nrow(gene_nodes) > 0) {
        gene_nodes$size <- pmin(pmax(gene_nodes$degree * 3 + 8, 8), 30)
      }

      return(list(
        pathway = pathway_nodes,
        gene = gene_nodes
      ))
    })

    # ─── 输出: 选区摘要 ───
    output$selection_summary <- shiny::renderText({
      n <- length(selected_pathways())
      sel <- table_controller$selected_pathways()
      if (!is.null(sel) && length(sel) > 0) {
        return(sprintf("%d from Main Table", n))
      }
      return(sprintf("Top %d by |NES|", n))
    })

    # ─── 输出: 网络统计 ───
    output$network_stats <- shiny::renderPrint({
      net <- network_data()

      if (is.null(net) || is.null(net$edges)) {
        cat("No data\n")
        return()
      }

      cat("HubGene Network\n")
      cat(paste(rep("-", 25), collapse = ""), "\n")
      cat(sprintf("Pathways: %d\n", nrow(net$nodes$pathway)))
      cat(sprintf("Genes: %d\n", nrow(net$nodes$gene)))
      cat(sprintf("Edges: %d\n", nrow(net$edges)))

      if (!is.null(net$hub_df) && nrow(net$hub_df) > 0) {
        hub_only <- net$hub_df[net$hub_df$is_hub, ]
        if (nrow(hub_only) > 0) {
          cat(sprintf("Hub degree: %d-%d\n", min(hub_only$degree), max(hub_only$degree)))
        }
      }
    })

    # ─── 输出: 主绘图 ───
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
                 plotly::layout(title = "No connections"))
      }

      left_group <- data_list$left_group
      right_group <- data_list$right_group

      p <- plotly::plot_ly(source = ns("hubgene_source"))

      # ─── 添加边 ───
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

        p <- p %>% plotly::add_trace(
          type = "scatter",
          mode = "lines",
          x = c(x0, x1, NA),
          y = c(y0, y1, NA),
          line = list(color = "rgba(150,150,150,0.35)", width = 1.2),
          hoverinfo = "skip",
          showlegend = FALSE,
          inherit = FALSE
        )
      }

      # ─── 通路 hover text ───
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

      # ─── 添加通路节点 ───
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

      # ─── 基因 hover text ───
      if (nrow(gene_nodes) > 0) {
        gene_hover <- sapply(seq_len(nrow(gene_nodes)), function(i) {
          row <- gene_nodes[i, ]
          stat_val <- as.numeric(row$stat)
          degree <- as.integer(row$degree)
          pathways <- as.character(row$pathways)
          color <- as.character(row$color)
          direction <- if (is.na(stat_val) || stat_val == 0) {
            "No direction"
          } else if (stat_val > 0) {
            paste0("Up in ", left_group)
          } else {
            paste0("Up in ", right_group)
          }

          sprintf("<b style='font-size:13px; color:%s;'>%s</b><br><b>Stat:</b> %.3f
                  <b>Direction:</b> %s
                  <b>Degree:</b> %d
                  <b>Pathways:</b>
                  <span style='font-size:10px; word-break:break-all;'>%s</span>",
                  color, row$id, stat_val, direction, degree, pathways)
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
            line = list(color = "white", width = 0.5)
          ),
          text = gene_hover,
          hovertemplate = "%{text}<extra></extra>",
          name = "Genes",
          showlegend = TRUE,
          inherit = FALSE
        )
      }

      # ─── 添加标签 ───
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

      # ─── 布局 ───
      all_x <- c(pathway_nodes$x, gene_nodes$x)
      all_y <- c(pathway_nodes$y, gene_nodes$y)

      x_range <- c(min(all_x) - 1.2, max(all_x) + 1.2)
      y_range <- c(min(all_y) - 1.2, max(all_y) + 1.2)

      p <- p %>% plotly::layout(
        title = list(
          text = sprintf("HubGene Network: %d Pathways, %d Genes, %d Edges",
                         nrow(pathway_nodes), nrow(gene_nodes), nrow(edges)),
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

  })
}
