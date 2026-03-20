#' @title 四重联动 UI
#' @keywords internal
mod_quadrant_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      shiny::column(6, shiny::div(class = "white-box",
                                  shiny::h4(.tr("quadrant.title_pathway_volcano")),
                                  plotly::plotlyOutput(ns("volcano_pathway"), height = "450px"))),
      shiny::column(6, shiny::div(class = "white-box",
                                  shiny::h4(.tr("quadrant.title_rank_dist")),
                                  plotly::plotlyOutput(ns("volcano_gene"), height = "450px")))
    ),
    shiny::fluidRow(
      shiny::column(6, shiny::div(class = "white-box",
                                  shiny::h4(.tr("quadrant.title_deg_volcano")),
                                  plotly::plotlyOutput(ns("de_volcano"), height = "450px"))),
      shiny::column(6, shiny::div(class = "white-box",
                                  shiny::h4(.tr("quadrant.title_expr_box")),
                                  plotly::plotlyOutput(ns("gene_expr_box"), height = "450px")))
    )
  )
}


#' @title 四重联动 Server（修复版）
#' @description 修复：DESeq2列名适配 + 宏观微观联动高亮
#' @keywords internal
mod_quadrant_server <- function(id, data_prep, gsea_res) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    COLOR_LEFT <- "#E41A1C"
    COLOR_RIGHT <- "#377EB8"
    COLOR_NS <- "#C0C0C0"
    COLOR_HIGHLIGHT <- "#FF9800"

    selected_labels <- shiny::reactiveVal(character(0))
    # 🆕 新增：存储选中通路的基因列表，用于差异火山图高亮
    selected_pathway_genes <- shiny::reactiveVal(character(0))

    # 1. 通路火山图
    output$volcano_pathway <- plotly::renderPlotly({
      data_list <- data_prep()
      shiny::req(data_list)

      df <- data_list$df
      left_grp <- data_list$left_group
      right_grp <- data_list$right_group

      df$color <- ifelse(df$NES > 0, COLOR_LEFT, COLOR_RIGHT)

      p <- plotly::plot_ly(
        data = df,
        x = ~NES,
        y = ~-log10(p.adjust),
        type = "scatter",
        mode = "markers",
        marker = list(color = ~color, size = 10, opacity = 0.8, line = list(color = "white", width = 1)),
        text = ~sprintf("%s<br>NES: %.2f<br>FDR: %.2e", ID, NES, p.adjust),
        hoverinfo = "text",
        key = ~ID,
        source = ns("pathway_volcano")
      ) %>%
        plotly::event_register('plotly_click') %>%
        plotly::layout(
          xaxis = list(title = "NES", zeroline = FALSE),
          yaxis = list(title = "-log10 (FDR)", zeroline = FALSE),
          dragmode = "select",
          showlegend = FALSE
        )

      labs <- selected_labels()
      if (length(labs) > 0) {
        lab_df <- df[df$ID %in% labs, ]
        if (nrow(lab_df) > 0) {
          p <- p %>% plotly::add_annotations(
            x = lab_df$NES, y = -log10(lab_df$p.adjust),
            text = lab_df$ID, showarrow = TRUE, arrowhead = 2,
            ax = 20, ay = -30, font = list(size = 10, color = "black"),
            bgcolor = "rgba(255,255,255,0.9)", bordercolor = "#333", borderwidth = 1
          )
        }
      }
      p
    })

    # 🔧 修改：点击事件同时更新标签和存储通路基因
    shiny::observeEvent(plotly::event_data("plotly_click", source = ns("pathway_volcano")), {
      click <- plotly::event_data("plotly_click", source = ns("pathway_volcano"))
      if (!is.null(click$key)) {
        # 更新标签（原有功能）
        current <- selected_labels()
        if (click$key %in% current) {
          selected_labels(setdiff(current, click$key))
        } else {
          selected_labels(c(current, click$key))
        }

        # 🆕 新增：存储该通路的基因，用于差异火山图高亮
        data_list <- data_prep()
        if (!is.null(data_list)) {
          pathway_genes <- data_list$gsea_res@geneSets[[click$key]]
          selected_pathway_genes(toupper(pathway_genes))
          message(sprintf("🎯 已选中通路: %s，包含 %d 个基因", click$key, length(pathway_genes)))
        }
      }
    })

    # 2. 基因 Rank 分布
    output$volcano_gene <- plotly::renderPlotly({
      data_list <- data_prep()
      shiny::req(data_list)

      genelist <- data_list$gsea_res@geneList

      rank_df <- data.frame(
        Rank = seq_along(genelist),
        Metric = as.numeric(genelist),
        Gene = names(genelist),
        stringsAsFactors = FALSE
      )

      click <- plotly::event_data("plotly_click", source = ns("pathway_volcano"))
      rank_df$Color <- COLOR_NS
      rank_df$Size <- 4

      if (!is.null(click) && !is.null(click$key)) {
        pathway_genes <- data_list$gsea_res@geneSets[[click$key]]
        match_idx <- which(toupper(rank_df$Gene) %in% toupper(pathway_genes))
        if (length(match_idx) > 0) {
          rank_df$Color[match_idx] <- COLOR_HIGHLIGHT
          rank_df$Size[match_idx] <- 12
        }
      }

      plotly::plot_ly(
        data = rank_df, x = ~Rank, y = ~Metric,
        type = "scattergl", mode = "markers",
        marker = list(color = ~Color, size = ~Size, opacity = 0.8, line = list(width = 0)),
        text = ~Gene, hoverinfo = "text"
      ) %>% plotly::layout(
        xaxis = list(title = "Gene Rank"),
        yaxis = list(title = "Ranking Metric (Stat)"),
        showlegend = FALSE
      )
    })

    # 🔧 修复4：差异表达火山图（DESeq2适配 + NA处理+高亮基因置顶）
    output$de_volcano <- plotly::renderPlotly({
      data_list <- data_prep()
      shiny::req(data_list)

      contrast_id <- data_list$contrast_id

      de_df <- tryCatch({
        get_de_table(gsea_res, contrast_id)
      }, error = function(e) {
        message(sprintf("DE表获取失败: %s", e$message))
        return(NULL)
      })

      shiny::req(de_df)

      # 列名标准化
      if (!"logFC" %in% colnames(de_df) && "log2FoldChange" %in% colnames(de_df)) {
        de_df$logFC <- de_df$log2FoldChange
      }
      if (!"pvalue" %in% colnames(de_df) && "p.value" %in% colnames(de_df)) {
        de_df$pvalue <- de_df$p.value
      }
      if (!"padj" %in% colnames(de_df) && "p.adjust" %in% colnames(de_df)) {
        de_df$padj <- de_df$p.adjust
      }

      # 移除NA
      de_df <- de_df[!is.na(de_df$logFC) & !is.na(de_df$pvalue) & !is.na(de_df$padj), ]

      if (nrow(de_df) == 0) return(NULL)

      de_df$x_axis <- de_df$logFC
      de_df$y_axis <- -log10(de_df$pvalue)

      # 处理Inf
      inf_y <- is.infinite(de_df$y_axis)
      if (any(inf_y)) {
        max_y <- max(de_df$y_axis[!inf_y], na.rm = TRUE)
        de_df$y_axis[inf_y] <- max_y * 1.1
      }

      # 标记显著性
      de_df$significant <- abs(de_df$logFC) > 0.5 & de_df$padj < 0.05

      # 基础设置
      de_df$color <- COLOR_NS
      de_df$size <- 6

      # 高亮通路基因
      highlight_genes <- selected_pathway_genes()
      is_highlight <- toupper(de_df$gene_symbol) %in% highlight_genes

      if (any(is_highlight)) {
        de_df$color[is_highlight] <- COLOR_HIGHLIGHT
        de_df$size[is_highlight] <- 12
      }

      # 显著性颜色（非高亮基因）
      sig_not_highlight <- de_df$significant & de_df$color == COLOR_NS
      if (any(sig_not_highlight)) {
        de_df$color[sig_not_highlight] <- ifelse(de_df$logFC[sig_not_highlight] > 0, COLOR_LEFT, COLOR_RIGHT)
      }

      # 🔧 关键修复：按color列排序，确保高亮基因（橙色）最后绘制，置顶显示
      # 排序优先级: NS(灰) -> 显著(红/蓝) -> 高亮(橙)
      color_order <- ifelse(de_df$color == COLOR_NS, 1,
                            ifelse(de_df$color == COLOR_HIGHLIGHT, 3, 2))
      de_df <- de_df[order(color_order), ]

      message(sprintf("🎨 火山图绘制顺序: %d灰 -> %d显著 -> %d高亮",
                      sum(color_order==1), sum(color_order==2), sum(color_order==3)))

      plotly::plot_ly(
        data = de_df,
        x = ~x_axis,
        y = ~y_axis,
        type = "scattergl",
        mode = "markers",
        marker = list(
          color = ~color,
          size = ~size,
          opacity = 0.9,  # 略微增加不透明度
          line = list(color = "white", width = 1)
        ),
        text = ~sprintf("%s<br>logFC: %.2f<br>FDR: %.2e%s",
                        gene_symbol, logFC, padj,
                        ifelse(color == COLOR_HIGHLIGHT, "<br>⭐ PATHWAY GENE", "")),
        hoverinfo = "text",
        key = ~toupper(gene_symbol),
        source = ns("deg_volcano")  # 🔧 关键添加：确保图4能接收点击事件
      ) %>%
        plotly::layout(
          xaxis = list(title = "log2 Fold Change", zeroline = FALSE),
          yaxis = list(title = "-log10 P-value", zeroline = FALSE),
          showlegend = FALSE
        )
    })

    # 4. 表达箱线图
    output$gene_expr_box <- plotly::renderPlotly({
      gene_click <- plotly::event_data("plotly_click", source = ns("deg_volcano"))

      if (is.null(gene_click) || is.null(gene_click$key)) {
        return(plotly::plot_ly() %>% plotly::layout(
          title = list(text = .tr("quadrant.hint_click_gene"), font = list(size = 14))
        ))
      }

      data_list <- data_prep()
      shiny::req(data_list)

      expr_mat <- get_expr_matrix(gsea_res, type = data_list$expression_type)
      sample_meta <- get_sample_meta(gsea_res)

      if (is.null(expr_mat) || is.null(sample_meta)) {
        return(plotly::plot_ly() %>% plotly::layout(title = "表达矩阵不可用"))
      }

      # 查找基因
      target_gene_upper <- gene_click$key
      gene_names_upper <- toupper(rownames(expr_mat))
      match_idx <- which(gene_names_upper == target_gene_upper)

      if (length(match_idx) == 0) {
        return(plotly::plot_ly() %>% plotly::layout(
          title = sprintf("基因 '%s' 未找到", target_gene_upper)
        ))
      }

      actual_gene <- rownames(expr_mat)[match_idx[1]]
      expr_values <- expr_mat[actual_gene, ]

      # 构建绘图数据
      sample_names <- names(expr_values)
      group_info <- sample_meta$group[match(sample_names, rownames(sample_meta))]

      plot_data <- data.frame(
        Sample = sample_names,
        Expression = as.numeric(expr_values),
        Group = group_info,
        stringsAsFactors = FALSE
      )

      # 移除NA
      plot_data <- plot_data[!is.na(plot_data$Group), ]

      if (nrow(plot_data) == 0) {
        return(plotly::plot_ly() %>% plotly::layout(title = "无匹配的样本分组数据"))
      }

      left_grp <- data_list$left_group
      right_grp <- data_list$right_group

      group_colors <- c(COLOR_LEFT, COLOR_RIGHT)
      names(group_colors) <- c(left_grp, right_grp)

      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Group, y = Expression, fill = Group)) +
        ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        ggplot2::geom_jitter(width = 0.2, size = 3, ggplot2::aes(text = Sample)) +
        ggplot2::scale_fill_manual(values = group_colors) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::labs(title = actual_gene, y = data_list$expression_type, x = "") +
        ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

      plotly::ggplotly(p, tooltip = c("text", "y"))
    })
  })
}
