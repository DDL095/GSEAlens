#' @title 四重联动 UI（修复版）
#' @keywords internal
mod_quadrant_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      shiny::column(6, shiny::div(class = "white-box",
                                  shiny::h4("1. 宏观: 通路火山图"),
                                  shiny::div(
                                    style = "position: relative;",
                                    plotly::plotlyOutput(ns("volcano_pathway"), height = "450px"),
                                    shiny::div(
                                      style = "position: absolute; top: 10px; right: 10px; background: rgba(255,255,255,0.9); padding: 5px 10px; border-radius: 4px; font-size: 11px; color: #666;",
                                      "点击通路查看详情"
                                    )
                                  ))),
      shiny::column(6, shiny::div(class = "white-box",
                                  shiny::h4("2. 微观: 基因 Rank 分布"),
                                  plotly::plotlyOutput(ns("volcano_gene"), height = "450px")))
    ),
    shiny::fluidRow(
      shiny::column(6, shiny::div(class = "white-box",
                                  shiny::h4("3. 差异表达火山图"),
                                  plotly::plotlyOutput(ns("de_volcano"), height = "450px"))),
      shiny::column(6, shiny::div(class = "white-box",
                                  shiny::h4("4. 全量表达箱线图"),
                                  shiny::div(
                                    style = "position: relative;",
                                    plotly::plotlyOutput(ns("gene_expr_box"), height = "450px"),
                                    shiny::uiOutput(ns("boxplot_order_status"))
                                  )))
    )
  )
}

#' @title 四重联动 Server（连续点击多选 + 排序根治版）
#' @description 支持连续点击标记多个通路，Boxplot排序强制锁定
#' @keywords internal
mod_quadrant_server <- function(id, data_prep_list, gsea_res) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    COLOR_LEFT <- "#E41A1C"
    COLOR_RIGHT <- "#377EB8"
    COLOR_NS <- "#C0C0C0"
    COLOR_HIGHLIGHT <- "#FF9800"

    selected_pathway_ids <- shiny::reactiveVal(character(0))
    selected_pathway_genes <- shiny::reactiveVal(character(0))

    data_prep_data <- data_prep_list$data
    highlight_genes_reactive <- data_prep_list$highlight_genes

    # 🔧 修复：直接持有reactiveVal引用，确保实时同步
    boxplot_order_ref <- data_prep_list$boxplot_order

    shiny::observeEvent(data_prep_data(), {
      selected_pathway_ids(character(0))
      selected_pathway_genes(character(0))

      # 🔴 关键修复：当对比组切换时，重置排序为默认
      boxplot_order_ref("default")
      message("🔄 [联动] 对比组已切换，排序已重置为默认")
    })

    # 1. 通路火山图（连续点击多选）
    output$volcano_pathway <- plotly::renderPlotly({
      data_list <- data_prep_data()
      shiny::req(data_list)

      df <- data_list$df
      left_grp <- data_list$left_group
      right_grp <- data_list$right_group

      current_selections <- selected_pathway_ids()

      df$color <- ifelse(df$ID %in% current_selections, COLOR_HIGHLIGHT,
                         ifelse(df$NES > 0, COLOR_LEFT, COLOR_RIGHT))
      df$size <- ifelse(df$ID %in% current_selections, 18, 10)
      df$opacity <- ifelse(length(current_selections) == 0, 0.8,
                           ifelse(df$ID %in% current_selections, 1.0, 0.35))
      df$linewidth <- ifelse(df$ID %in% current_selections, 3, 1)

      annotations_list <- list()
      if (length(current_selections) > 0) {
        selected_df <- df[df$ID %in% current_selections, ]
        for (i in 1:nrow(selected_df)) {
          row <- selected_df[i, ]
          ay_offset <- -40 - ((i-1) %% 3) * 25
          ax_offset <- ifelse(i %% 2 == 1, 0, 60)

          annotations_list[[i]] <- list(
            x = row$NES,
            y = -log10(row$p.adjust),
            text = row$ID,
            showarrow = TRUE,
            arrowhead = 2,
            arrowsize = 1,
            arrowwidth = 2,
            arrowcolor = COLOR_HIGHLIGHT,
            ax = ax_offset,
            ay = ay_offset,
            font = list(size = 11, color = COLOR_HIGHLIGHT, family = "Arial"),
            bgcolor = "rgba(255,255,255,0.95)",
            bordercolor = COLOR_HIGHLIGHT,
            borderwidth = 2,
            borderpad = 4
          )
        }
      }

      plotly::plot_ly(
        data = df,
        x = ~NES,
        y = ~-log10(p.adjust),
        type = "scatter",
        mode = "markers",
        marker = list(
          color = ~color,
          size = ~size,
          opacity = ~opacity,
          line = list(color = "white", width = ~linewidth)
        ),
        text = ~sprintf("%s<br>NES: %.2f<br>FDR: %.2e", ID, NES, p.adjust),
        hoverinfo = "text",
        key = ~ID,
        source = ns("pathway_volcano")
      ) %>%
        plotly::layout(
          xaxis = list(title = "NES", zeroline = FALSE),
          yaxis = list(title = "-log10 (FDR)", zeroline = FALSE),
          showlegend = FALSE,
          dragmode = "pan",
          annotations = annotations_list
        )
    })

    shiny::observeEvent(plotly::event_data("plotly_click", source = ns("pathway_volcano")), {
      click <- plotly::event_data("plotly_click", source = ns("pathway_volcano"))

      if (is.null(click) || is.null(click$key)) return()

      clicked_id <- click$key
      current <- selected_pathway_ids()

      if (clicked_id %in% current) {
        new_selection <- setdiff(current, clicked_id)
        selected_pathway_ids(new_selection)
        message(sprintf("❌ 移除: %s (剩余%d个)", clicked_id, length(new_selection)))
      } else {
        new_selection <- c(current, clicked_id)
        selected_pathway_ids(new_selection)
        message(sprintf("✅ 添加: %s (共%d个)", clicked_id, length(new_selection)))

        data_list <- data_prep_data()
        if (!is.null(data_list)) {
          pathway_genes <- data_list$gsea_res@geneSets[[clicked_id]]
          current_genes <- selected_pathway_genes()
          selected_pathway_genes(unique(c(current_genes, toupper(pathway_genes))))
        }
      }
    })

    # 2. 基因 Rank 分布
    output$volcano_gene <- plotly::renderPlotly({
      data_list <- data_prep_data()
      shiny::req(data_list)

      genelist <- data_list$gsea_res@geneList

      rank_df <- data.frame(
        Rank = seq_along(genelist),
        Metric = as.numeric(genelist),
        Gene = names(genelist),
        stringsAsFactors = FALSE
      )

      current_pws <- selected_pathway_ids()
      rank_df$Color <- COLOR_NS
      rank_df$Size <- 4

      if (length(current_pws) > 0) {
        pathway_genes <- selected_pathway_genes()
        match_idx <- which(toupper(rank_df$Gene) %in% pathway_genes)
        if (length(match_idx) > 0) {
          rank_df$Color[match_idx] <- COLOR_HIGHLIGHT
          rank_df$Size[match_idx] <- 12
        }
      }

      plotly::plot_ly(
        data = rank_df, x = ~Rank, y = ~Metric,
        type = "scattergl",
        mode = "markers",
        marker = list(color = ~Color, size = ~Size, opacity = 0.8, line = list(width = 0)),
        text = ~Gene, hoverinfo = "text"
      ) %>% plotly::layout(
        xaxis = list(title = "Gene Rank"),
        yaxis = list(title = "Ranking Metric (Stat)"),
        showlegend = FALSE,
        title = list(
          text = ifelse(length(current_pws) > 0,
                        sprintf("已选%d个通路 | 共%d个基因", length(current_pws), length(selected_pathway_genes())),
                        "点击上方火山图标记通路"),
          font = list(size = 12)
        )
      )
    })

    # 3. 差异表达火山图
    output$de_volcano <- plotly::renderPlotly({
      data_list <- data_prep_data()
      shiny::req(data_list)

      contrast_id <- data_list$contrast_id

      de_df <- tryCatch({
        get_de_table(gsea_res, contrast_id)
      }, error = function(e) {
        message(sprintf("DE表获取失败: %s", e$message))
        return(NULL)
      })

      shiny::req(de_df)

      col_mapping <- list(
        logFC = c("logFC", "log2FoldChange"),
        pvalue = c("pvalue", "p.value", "P.Value"),
        padj = c("padj", "p.adjust", "adj.P.Val")
      )

      for (std_name in names(col_mapping)) {
        if (!(std_name %in% colnames(de_df))) {
          for (alt_name in col_mapping[[std_name]]) {
            if (alt_name %in% colnames(de_df)) {
              de_df[[std_name]] <- de_df[[alt_name]]
              break
            }
          }
        }
      }

      de_df <- de_df[!is.na(de_df$logFC) & !is.na(de_df$pvalue) & !is.na(de_df$padj), ]
      if (nrow(de_df) == 0) return(NULL)

      de_df$x_axis <- de_df$logFC
      de_df$y_axis <- -log10(de_df$pvalue)

      inf_y <- is.infinite(de_df$y_axis)
      if (any(inf_y)) {
        max_y <- max(de_df$y_axis[!inf_y], na.rm = TRUE)
        de_df$y_axis[inf_y] <- max_y * 1.1
      }

      user_genes <- toupper(highlight_genes_reactive())
      pathway_genes <- selected_pathway_genes()
      de_df$gene_upper <- toupper(de_df$gene_symbol)

      de_df$category <- "normal"
      de_df$category[de_df$gene_upper %in% pathway_genes] <- "pathway"
      de_df$category[de_df$gene_upper %in% user_genes] <- "user"

      de_df$color <- COLOR_NS
      de_df$size <- 6
      de_df$significant <- abs(de_df$logFC) > 0.5 & de_df$padj < 0.05

      user_idx <- which(de_df$category == "user")
      if (length(user_idx) > 0) {
        de_df$color[user_idx] <- ifelse(de_df$logFC[user_idx] > 0, COLOR_LEFT, COLOR_RIGHT)
        de_df$size[user_idx] <- 15
      }

      pathway_idx <- which(de_df$category == "pathway" & de_df$color == COLOR_NS)
      if (length(pathway_idx) > 0) {
        de_df$color[pathway_idx] <- COLOR_HIGHLIGHT
        de_df$size[pathway_idx] <- 12
      }

      sig_idx <- which(de_df$significant & de_df$color == COLOR_NS)
      if (length(sig_idx) > 0) {
        de_df$color[sig_idx] <- ifelse(de_df$logFC[sig_idx] > 0, COLOR_LEFT, COLOR_RIGHT)
        de_df$size[sig_idx] <- 8
      }

      de_df$priority <- ifelse(de_df$category == "user", 4,
                               ifelse(de_df$category == "pathway", 3,
                                      ifelse(de_df$significant, 2, 1)))
      de_df <- de_df[order(de_df$priority), ]

      user_genes_df <- de_df[de_df$category == "user", ]

      p <- plotly::plot_ly(
        data = de_df,
        x = ~x_axis,
        y = ~y_axis,
        type = "scattergl",
        mode = "markers",
        marker = list(
          color = ~color,
          size = ~size,
          opacity = 0.9,
          line = list(color = "white", width = 1)
        ),
        text = ~sprintf("%s<br>logFC: %.2f<br>FDR: %.2e%s",
                        gene_symbol, logFC, padj,
                        ifelse(category == "user", " ⭐ USER",
                               ifelse(category == "pathway", " 🔥 PATHWAY", ""))),
        hoverinfo = "text",
        key = ~gene_upper,
        source = ns("deg_volcano")
      ) %>%
        plotly::layout(
          xaxis = list(title = "log2 Fold Change", zeroline = FALSE),
          yaxis = list(title = "-log10 P-value", zeroline = FALSE),
          showlegend = FALSE
        )

      if (nrow(user_genes_df) > 0) {
        annotations <- lapply(1:nrow(user_genes_df), function(i) {
          gene <- user_genes_df[i, ]
          base_ay <- -40
          stagger <- (i %% 3) * 15
          ax <- ifelse(gene$logFC > 0, 30 + (i %% 2) * 20, -30 - (i %% 2) * 20)

          list(
            x = gene$x_axis,
            y = gene$y_axis,
            text = gene$gene_symbol,
            showarrow = TRUE,
            arrowhead = 0,
            arrowsize = 1,
            arrowwidth = 1,
            arrowcolor = ifelse(gene$logFC > 0, COLOR_LEFT, COLOR_RIGHT),
            ax = ax,
            ay = base_ay - stagger,
            bgcolor = "rgba(255,255,255,0.9)",
            bordercolor = ifelse(gene$logFC > 0, COLOR_LEFT, COLOR_RIGHT),
            borderwidth = 1,
            font = list(size = 10, color = ifelse(gene$logFC > 0, COLOR_LEFT, COLOR_RIGHT))
          )
        })
        p <- p %>% plotly::layout(annotations = annotations)
      }

      p
    })

    # 4. 全量表达箱线图（完全根治版）
    output$gene_expr_box <- plotly::renderPlotly({
      # ✅ 第一步：最前端立即读取排序选择
      current_pending_order <- input$boxplot_order_pending  # 从UI读取
      current_confirmed_order <- boxplot_order_ref()       # 从reactiveVal读取

      # 使用已确认的排序（优先级更高）
      final_order_to_use <- if (!is.null(current_confirmed_order) &&
                                current_confirmed_order != "default" &&
                                current_confirmed_order != "") {
        current_confirmed_order
      } else {
        "default"
      }

      message(sprintf("📋 [Boxplot Render] 读取排序: pending='%s' | confirmed='%s' | 最终使用='%s'",
                      current_pending_order, current_confirmed_order, final_order_to_use))

      # 获取点击事件
      gene_click <- plotly::event_data("plotly_click", source = ns("deg_volcano"))

      if (is.null(gene_click) || is.null(gene_click$key)) {
        return(plotly::plot_ly() %>% plotly::layout(
          title = list(text = "👈 请在左侧火山图点击基因", font = list(size = 14))
        ))
      }

      data_list <- data_prep_data()
      shiny::req(data_list)

      tryCatch({
        expr_mat <- get_expr_matrix(gsea_res, type = data_list$expression_type)
        sample_meta <- get_sample_meta(gsea_res)

        if (is.null(expr_mat) || is.null(sample_meta)) {
          return(plotly::plot_ly() %>% plotly::layout(title = "表达矩阵不可用"))
        }

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

        sample_names <- names(expr_values)
        group_info <- sample_meta$group[match(sample_names, rownames(sample_meta))]

        plot_data <- data.frame(
          Sample = sample_names,
          Expression = as.numeric(expr_values),
          Group = group_info,
          stringsAsFactors = FALSE
        )

        plot_data <- plot_data[!is.na(plot_data$Group), ]

        if (nrow(plot_data) == 0) {
          return(plotly::plot_ly() %>% plotly::layout(title = "无有效分组数据"))
        }

        # ✅ 关键第二步：解析排序设置，转换为x_categories排列顺序
        actual_groups <- unique(as.character(plot_data$Group))
        x_categories <- NULL

        if (final_order_to_use != "default" && final_order_to_use != "" && !is.na(final_order_to_use)) {
          # 将逗号分隔的字符串转换为组名向量
          order_parts <- strsplit(final_order_to_use, ",")[[1]]
          order_parts <- trimws(order_parts)  # 移除空格

          message(sprintf("   🔍 解析排序字符串: [%s]", paste(order_parts, collapse=" | ")))

          # 验证解析结果
          valid_parts <- order_parts[order_parts %in% actual_groups]
          missing_parts <- setdiff(order_parts, actual_groups)

          if (length(missing_parts) > 0) {
            message(sprintf("   ⚠️ 警告：排序包含不存在的组[%s]，自动修正",
                            paste(missing_parts, collapse=",")))
          }

          # 严格按照指定顺序排列
          x_categories <- valid_parts
          # 追加未指定的组
          x_categories <- c(x_categories, setdiff(actual_groups, x_categories))

          message(sprintf("   ✅ 最终排序: [%s]", paste(x_categories, collapse=" → ")))
        } else {
          x_categories <- actual_groups
          message(sprintf("   ✅ 使用默认排序: [%s]", paste(x_categories, collapse=" → ")))
        }

        # ✅ 关键第三步：确保Group因子级别与x_categories完全一致
        plot_data <- plot_data[plot_data$Group %in% x_categories, ]
        plot_data$Group <- factor(plot_data$Group, levels = x_categories, ordered = TRUE)

        # 颜色设置
        unique_groups <- levels(plot_data$Group)
        if (length(unique_groups) == 2) {
          group_colors <- c(COLOR_LEFT, COLOR_RIGHT)
          names(group_colors) <- unique_groups
        } else {
          group_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                            "#FF7F00", "#A65628", "#F781BF", "#999999")[1:length(unique_groups)]
          names(group_colors) <- unique_groups
        }

        # ✅ 关键第四步：ggplot2层面强制顺序
        p <- ggplot2::ggplot(plot_data,
                             ggplot2::aes(x = Group, y = Expression, fill = Group)) +
          ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA) +
          ggplot2::geom_jitter(width = 0.2, size = 3, alpha = 0.6) +
          ggplot2::scale_fill_manual(values = group_colors) +
          # 🔴 最关键：使用 limits 参数固定顺序
          ggplot2::scale_x_discrete(limits = x_categories, drop = FALSE) +
          ggplot2::theme_bw(base_size = 12) +
          ggplot2::labs(title = sprintf("%s (排序: %s)",
                                        actual_gene,
                                        if(final_order_to_use == "default") "默认" else "自定义"),
                        y = data_list$expression_type,
                        x = NULL) +
          ggplot2::theme(
            legend.position = "none",
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
          )

        # 转换为plotly，并在plotly层面也强制顺序
        ply <- plotly::ggplotly(p)

        ply <- ply %>% plotly::layout(
          xaxis = list(
            categoryorder = "array",
            categoryarray = x_categories,
            title = ""
          )
        )

        return(ply)

      }, error = function(e) {
        message(sprintf("❌ Boxplot错误: %s", e$message))
        return(plotly::plot_ly() %>% plotly::layout(
          title = sprintf("错误: %s", e$message)
        ))
      })
    })

    # 显示排序状态指示器
    output$boxplot_order_status <- shiny::renderUI({
      order <- boxplot_order_ref()
      if (is.null(order) || order == "default" || order == "" || is.na(order)) {
        shiny::div(
          style = "position: absolute; top: 10px; right: 10px; background: rgba(255,255,255,0.9); padding: 5px 10px; border-radius: 4px; font-size: 11px; color: #666;",
          "排序: 默认"
        )
      } else {
        shiny::div(
          style = "position: absolute; top: 10px; right: 10px; background: rgba(76,175,80,0.2); padding: 5px 10px; border-radius: 4px; font-size: 11px; color: #2e7d32; font-weight: bold;",
          sprintf("排序: %s ✓", gsub(",", " → ", order))
        )
      }
    })

    # ✅ 修复 6：显示当前排序状态（实时更新）
    output$boxplot_order_status <- shiny::renderUI({
      order <- boxplot_order_ref()  # 实时读取
      if (is.null(order) || is.na(order) ||
          order == "default" || order == "" || order == "limited") {
        shiny::div(
          style = "position: absolute; top: 10px; right: 10px; background: rgba(255,255,255,0.9); padding: 5px 10px; border-radius: 4px; font-size: 11px; color: #666;",
          "排序: 默认"
        )
      } else {
        shiny::div(
          style = "position: absolute; top: 10px; right: 10px; background: rgba(212,237,218,0.9); padding: 5px 10px; border-radius: 4px; font-size: 11px; color: #155724;",
          sprintf("排序: %s", gsub(",", " → ", order))
        )
      }
    })
  })
}
