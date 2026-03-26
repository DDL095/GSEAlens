#' @title Quadrant Linkage UI
#' @keywords internal

mod_quadrant_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(


      shiny::column(6, shiny::div(class = "white-box",
                                  shiny::h4("1. Macro: Pathway Volcano"),
                                  shiny::div(
                                    style = "position: relative;",
                                    plotly::plotlyOutput(ns("volcano_pathway"), height = "450px"),
                                    shiny::div(
                                      style = "position: absolute; top: 10px; left: 10px; background: rgba(255,255,255,0.9); padding: 5px 10px; border-radius: 4px; font-size: 11px; color: #666;",
                                      "Click pathway to view details"
                                    )
                                  ))),


      shiny::column(6, shiny::div(class = "white-box",
                                  shiny::h4("2. Micro: Gene Rank Distribution"),
                                  plotly::plotlyOutput(ns("volcano_gene"), height = "450px")))
    ),


    shiny::fluidRow(
      shiny::column(6, shiny::div(class = "white-box",
                                  shiny::h4("3. Differential Expression Volcano"),
                                  shiny::div(
                                    style = "position: absolute; top: 10px; right: 50%; z-index: 100;",
                                    shiny::actionButton(
                                      ns("toggle_volcano_settings"),
                                      label = "⚙️",
                                      style = "padding: 2px 8px; font-size: 12px;"
                                    )
                                  ),
                                  shiny::conditionalPanel(
                                    condition = sprintf("input['%s'] %% 2 == 1", ns("toggle_volcano_settings")),
                                    shiny::div(
                                      style = "background: #f8f9fa; padding: 10px; border-radius: 5px; margin-bottom: 10px;",
                                      shiny::numericInput(
                                        ns("volcano_logfc_thresh"),
                                        "logFC Threshold:",
                                        value = 1,
                                        min = 0,
                                        max = 22,
                                        step = 0.5
                                      ),
                                      shiny::numericInput(
                                        ns("volcano_pval_thresh"),
                                        "P-value Threshold:",
                                        value = 0.05,
                                        min = 0.001,
                                        max = 1,
                                        step = 0.01
                                      ),
                                      shiny::helpText("Click volcano area to refresh after setting")
                                    )
                                  ),
                                  plotly::plotlyOutput(ns("de_volcano"), height = "450px"))),


      shiny::column(6, shiny::div(class = "white-box",
                                  shiny::h4("4. Full Expression Distribution"),
                                  shiny::div(
                                    style = "position: relative;",
                                    # 添加切换控件
                                    shiny::div(
                                      style = "position: absolute; top: -40px; left: 30%; z-index: 100;",
                                      # 添加 0 基准点切换
                                      shiny::checkboxInput(
                                        ns("zero_baseline"),
                                        label = "Use 0 as baseline",
                                        value = FALSE
                                      )
                                    ),
                                    plotly::plotlyOutput(ns("gene_expr_box"), height = "450px"),
                                    shiny::uiOutput(ns("boxplot_order_status"))
                                  )))
    )
  )
}

#' @title Quadrant Linkage Server (Multi-select + Sort Fix)
#' @description Supports continuous clicking to mark multiple pathways, Boxplot sorting forced lock
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
      message("🔄 [Linkage] Comparison group switched, sorting reset to default")
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
          title = list(
            text = sprintf("Pathway Volcano: %s vs %s<br><sub>%d pathways | %d significant (FDR<0.25)</sub>",
                           data_list$left_group, data_list$right_group,
                           nrow(df), sum(df$p.adjust < 0.25, na.rm = TRUE)),
            font = list(size = 14),
            x = 0.5,
            xanchor = "center"
          ),
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
        message(sprintf("❌ Removed: %s (remaining %d)", clicked_id, length(new_selection)))
      } else {
        new_selection <- c(current, clicked_id)
        selected_pathway_ids(new_selection)
        message(sprintf("✅ Added: %s (total %d)", clicked_id, length(new_selection)))

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
                        sprintf("Selected %d pathways | %d genes total", length(current_pws), length(selected_pathway_genes())),
                        "Click pathway in volcano above to mark"),
          font = list(size = 12)
        )
      )
    })

    # 3. 差异表达火山图（恢复美学 + 可调阈值 + 修复点击联动）
    output$de_volcano <- plotly::renderPlotly({
      data_list <- data_prep_data()
      shiny::req(data_list)

      contrast_id <- data_list$contrast_id
      left_group <- data_list$left_group
      right_group <- data_list$right_group

      de_df <- tryCatch({
        get_de_table(gsea_res, contrast_id)
      }, error = function(e) {
        message(sprintf("DE table retrieval failed: %s", e$message))
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

      # 获取阈值设置（使用输入值或默认值）
      logfc_thresh <- if (!is.null(input$volcano_logfc_thresh)) input$volcano_logfc_thresh else 0
      pval_thresh <- if (!is.null(input$volcano_pval_thresh)) input$volcano_pval_thresh else 0.05

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

      # 分类基因 - 优先级：user > pathway > normal
      de_df$category <- "normal"
      de_df$category[de_df$gene_upper %in% pathway_genes] <- "pathway"
      de_df$category[de_df$gene_upper %in% user_genes] <- "user"

      # 根据阈值判断差异基因
      de_df$is_significant <- abs(de_df$logFC) > logfc_thresh & de_df$pvalue < pval_thresh

      # 统计差异基因数量
      n_up <- sum(de_df$is_significant & de_df$logFC > 0, na.rm = TRUE)
      n_down <- sum(de_df$is_significant & de_df$logFC < 0, na.rm = TRUE)
      n_not_sig <- sum(!de_df$is_significant, na.rm = TRUE)

      # 颜色设置
      COLOR_LEFT <- "#E41A1C"    # 红色 - 上调
      COLOR_RIGHT <- "#377EB8"   # 蓝色 - 下调
      COLOR_NS <- "#C0C0C0"      # 灰色 - 不显著
      COLOR_USER <- "#4DAF4A"    # 绿色- 用户基因
      COLOR_PATHWAY <- "#FF9800" # 橙色 - 通路基因

      # 设置基础颜色
      de_df$color <- dplyr::case_when(
        de_df$category == "user" ~ COLOR_USER,
        de_df$category == "pathway" ~ COLOR_PATHWAY,
        de_df$is_significant & de_df$logFC > 0 ~ COLOR_LEFT,
        de_df$is_significant & de_df$logFC < 0 ~ COLOR_RIGHT,
        TRUE ~ COLOR_NS
      )

      # 大小设置
      de_df$size <- dplyr::case_when(
        de_df$category == "user" ~ 17,
        de_df$category == "pathway" ~ 13,
        de_df$is_significant ~ 9,
        TRUE ~ 4
      )

      # 透明度设置
      de_df$opacity <- dplyr::case_when(
        de_df$category %in% c("user", "pathway") ~ 1.0,
        de_df$is_significant ~ 0.7,
        TRUE ~ 0.5
      )

      # 边框宽度
      de_df$linewidth <- dplyr::case_when(
        de_df$category %in% c("user", "pathway") ~ 0.8,
        TRUE ~ 0
      )

      # 关键修复：按plot_order排序，确保正确的绘制顺序
      de_df$plot_order <- dplyr::case_when(
        de_df$category == "user" ~ 3,
        de_df$category == "pathway" ~ 2,
        TRUE ~ 1
      )
      de_df <- de_df[order(de_df$plot_order), ]

      # 构建标题
      title_text <- sprintf(
        "%s vs %s<br><sup>↑ %d | ↓ %d | NS %d (logFC>|%.1f|, p<%.3f)</sup>",
        left_group, right_group,
        n_up, n_down, n_not_sig,
        logfc_thresh, pval_thresh
      )

      # 准备注释
      annotations_list <- list()

      if (logfc_thresh > 0) {
        annotations_list[[length(annotations_list) + 1]] <- list(
          x = logfc_thresh,
          y = max(de_df$y_axis) * 0.95,
          text = sprintf("logFC=%.1f", logfc_thresh),
          showarrow = FALSE,
          font = list(size = 10, color = "gray")
        )
        annotations_list[[length(annotations_list) + 1]] <- list(
          x = -logfc_thresh,
          y = max(de_df$y_axis) * 0.95,
          text = sprintf("-logFC=%.1f", logfc_thresh),
          showarrow = FALSE,
          font = list(size = 10, color = "gray")
        )
      }

      # 用户基因标签
      user_genes_df <- de_df[de_df$category == "user", ]
      if (nrow(user_genes_df) > 0) {
        annotations_list <- c(annotations_list, lapply(1:nrow(user_genes_df), function(i) {
          gene <- user_genes_df[i, ]
          list(
            x = gene$x_axis,
            y = gene$y_axis,
            text = gene$gene_symbol,
            showarrow = TRUE,
            arrowhead = 0,
            arrowsize = 1,
            arrowwidth = 2,
            arrowcolor = COLOR_USER,
            ax = ifelse(gene$logFC > 0, 50, -50),
            ay = -35,
            bgcolor = "rgba(255,255,255,0.85)",
            bordercolor = COLOR_USER,
            borderwidth = 0.5,
            font = list(size = 13, color = COLOR_USER)
          )
        }))
      }

      # 关键修复：使用plot_ly初始化，确保key正确绑定
      # 将数据框按plot_order排序后一次性传入
      p <- plotly::plot_ly(
        data = de_df,
        x = ~x_axis,
        y = ~y_axis,
        type = "scatter",
        mode = "markers",
        marker = list(
          color = ~color,
          size = ~size,
          opacity = ~opacity,
          line = list(color = "white", width = ~linewidth)
        ),
        text = ~sprintf("%s<br>logFC: %.2f<br>-log10(p): %.2f<br>FDR: %.2e<br>Category: %s",
                        gene_symbol, logFC, y_axis, padj, category),
        hoverinfo = "text",
        key = ~gene_upper,  # 关键：确保key正确绑定到gene_upper列
        source = ns("deg_volcano"),
        showlegend = FALSE
      ) %>% plotly::layout(
        title = list(
          text = title_text,
          font = list(size = 14),
          x = 0.5,
          xanchor = "center"
        ),
        xaxis = list(
          title = "logFC",
          zeroline = FALSE,
          showgrid = TRUE,
          gridcolor = "lightgray"
        ),
        yaxis = list(
          title = "-log10 (P-value)",
          zeroline = FALSE,
          showgrid = TRUE,
          gridcolor = "lightgray"
        ),
        showlegend = FALSE,
        dragmode = "pan",
        annotations = annotations_list,
        shapes = list(
          list(
            type = "line",
            x0 = logfc_thresh,
            x1 = logfc_thresh,
            y0 = 0,
            y1 = max(de_df$y_axis) * 1.05,
            line = list(color = "gray", dash = "dash", width = 1)
          ),
          list(
            type = "line",
            x0 = -logfc_thresh,
            x1 = -logfc_thresh,
            y0 = 0,
            y1 = max(de_df$y_axis) * 1.05,
            line = list(color = "gray", dash = "dash", width = 1)
          ),
          list(
            type = "line",
            x0 = min(de_df$x_axis) * 1.1,
            x1 = max(de_df$x_axis) * 1.1,
            y0 = -log10(pval_thresh),
            y1 = -log10(pval_thresh),
            line = list(color = "gray", dash = "dash", width = 1)
          )
        )
      )

      p
    })

    # 4. 全量表达箱线图（带0基准线选项）
    output$gene_expr_box <- plotly::renderPlotly({
      # 读取排序设置（保持与原有逻辑一致）
      current_confirmed_order <- boxplot_order_ref()
      final_order_to_use <- if (!is.null(current_confirmed_order) &&
                                current_confirmed_order != "default" &&
                                current_confirmed_order != "") {
        current_confirmed_order
      } else {
        "default"
      }

      # 获取点击事件（来自DE火山图）
      gene_click <- plotly::event_data("plotly_click", source = ns("deg_volcano"))

      if (is.null(gene_click) || is.null(gene_click$key)) {
        return(plotly::plot_ly() %>% plotly::layout(
          title = list(text = "👈 Please click gene in left volcano plot", font = list(size = 14))
        ))
      }

      data_list <- data_prep_data()
      shiny::req(data_list)

      tryCatch({
        expr_mat <- get_expr_matrix(gsea_res, type = data_list$expression_type)
        sample_meta <- get_sample_meta(gsea_res)

        if (is.null(expr_mat) || is.null(sample_meta)) {
          return(plotly::plot_ly() %>% plotly::layout(title = "Expression matrix not available"))
        }





        target_gene_upper <- gene_click$key


        # 🧬 多策略基因匹配逻辑（SYMBOL ↔ Ensembl ID 双向匹配）

        match_idx <- integer(0)

        # 策略1: 直接匹配表达矩阵行名（可能是 SYMBOL 或 Ensembl ID）
        gene_names_upper <- toupper(rownames(expr_mat))
        match_idx <- which(gene_names_upper == target_gene_upper)

        # 策略2: 如果失败，通过 gene_meta 进行 SYMBOL -> Ensembl ID 映射
        if (length(match_idx) == 0) {
          gene_meta <- gsea_res$expr_bundle$gene_meta

          if (!is.null(gene_meta) && nrow(gene_meta) > 0) {
            # 确保 gene_meta 有行名（与 expr_mat 对应）
            if (is.null(rownames(gene_meta))) {
              rownames(gene_meta) <- rownames(expr_mat)
            }

            # 查找可能的 SYMBOL 列（多候选兼容不同注释来源）
            symbol_candidates <- c("SYMBOL", "symbol", "Gene", "gene_name",
                                   "gene_symbol", "Gene.Symbol", "GeneName",
                                   "SYMBOL_char", "symbol_char")
            symbol_col <- intersect(symbol_candidates, colnames(gene_meta))[1]

            if (!is.na(symbol_col)) {
              # 在 gene_meta 中查找匹配的 SYMBOL（不区分大小写）
              meta_symbols_upper <- toupper(as.character(gene_meta[[symbol_col]]))
              meta_matches <- which(meta_symbols_upper == target_gene_upper)

              if (length(meta_matches) > 0) {
                # 获取对应的 Ensembl ID（行名）
                ensembl_ids <- rownames(gene_meta)[meta_matches]

                # 在表达矩阵中查找这些 Ensembl ID（优先找第一个存在的）
                for (eid in ensembl_ids) {
                  idx_in_expr <- which(rownames(expr_mat) == eid)
                  if (length(idx_in_expr) > 0) {
                    match_idx <- idx_in_expr[1]
                    break
                  }
                }

                if (length(match_idx) > 0) {
                  message(sprintf("✅ Mapped via gene_meta$%s: %s -> %s",
                                  symbol_col, target_gene_upper, ensembl_ids[1]))
                }
              }
            }
          }
        }

        # 策略3: 尝试部分匹配（容错处理，适用于不完整输入）
        if (length(match_idx) == 0) {
          # 在行名中部分匹配
          match_idx <- grep(target_gene_upper, gene_names_upper, ignore.case = TRUE)[1]

          # 如果仍失败，在 gene_meta 的所有字符列中搜索
          if (is.na(match_idx) && !is.null(gene_meta)) {
            for (col in colnames(gene_meta)) {
              if (is.character(gene_meta[[col]]) || is.factor(gene_meta[[col]])) {
                col_values <- toupper(as.character(gene_meta[[col]]))
                found_idx <- which(col_values == target_gene_upper)[1]
                if (!is.na(found_idx)) {
                  ensembl_id <- rownames(gene_meta)[found_idx]
                  match_idx <- which(rownames(expr_mat) == ensembl_id)[1]
                  if (!is.na(match_idx)) break
                }
              }
            }
          }
        }

        if (length(match_idx) == 0 || is.na(match_idx)) {
          return(plotly::plot_ly() %>% plotly::layout(
            title = sprintf("Gene '%s' not found (tried matching ID and SYMBOL)", target_gene_upper)
          ))
        }

        actual_gene <- rownames(expr_mat)[match_idx]
        # 🎯 关键：保留原始点击的 SYMBOL 作为显示名称（而非 Ensembl ID）
        display_gene_name <- gene_click$key

        message(sprintf("📊 Boxplot data: %s (matched to row: %s)", display_gene_name, actual_gene))




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
          return(plotly::plot_ly() %>% plotly::layout(title = "No valid group data"))
        }

        # 解析排序（保持与原有逻辑完全一致）
        actual_groups <- unique(as.character(plot_data$Group))
        x_categories <- NULL

        if (final_order_to_use != "default" && final_order_to_use != "" && !is.na(final_order_to_use)) {
          sep <- if (grepl("→", final_order_to_use, fixed = TRUE)) "→" else ","
          order_parts <- strsplit(final_order_to_use, sep)[[1]]
          order_parts <- trimws(order_parts)
          valid_parts <- order_parts[order_parts %in% actual_groups]
          x_categories <- c(valid_parts, setdiff(actual_groups, valid_parts))
        } else {
          x_categories <- actual_groups
        }

        # 强制factor levels（关键：保持排序稳定性）
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

        # 获取0基准线设置
        use_zero_baseline <- input$zero_baseline %||% FALSE

        # 计算y轴范围
        y_min <- min(plot_data$Expression, na.rm = TRUE)
        y_max <- max(plot_data$Expression, na.rm = TRUE)

        # 如果启用0基准线，确保包含0
        if (use_zero_baseline) {
          y_min <- min(y_min, 0)
          y_max <- max(y_max, 0)
        }

        # 添加一些边距
        y_range <- y_max - y_min
        y_min <- y_min - y_range * 0.1
        y_max <- y_max + y_range * 0.1

        # 绘制箱线图
        p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Group, y = Expression, fill = Group)) +
          ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA) +
          ggplot2::geom_jitter(width = 0.2, size = 3, alpha = 0.6,
                               ggplot2::aes(
                                 text = sprintf(
                                   "<b>Sample:</b> %s<br><b>Group:</b> %s<br><b>Expression:</b> %.3f",
                                   Sample, Group, Expression
                                 )
                               )
          )+
          ggplot2::scale_fill_manual(values = group_colors) +
          ggplot2::scale_x_discrete(limits = x_categories, drop = FALSE) +
          ggplot2::coord_cartesian(ylim = c(y_min, y_max)) +
          ggplot2::theme_bw(base_size = 12) +
          ggplot2::labs(
            title = sprintf("%s",
                            display_gene_name),
            y = data_list$expression_type,
            x = NULL
          ) +
          ggplot2::theme(
            legend.position = "none",
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
          )

        # 如果启用0基准线，添加虚线
        if (use_zero_baseline) {
          p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                                       color = "red", alpha = 0.7, size = 0.8)
        }

        # 转换为plotly并强制X轴顺序
        ply <- plotly::ggplotly(p, tooltip = "text")
        #ply <- plotly::ggplotly(p, tooltip = c("x", "y", "Sample"))

        ply <- ply %>% plotly::layout(
          xaxis = list(
            categoryorder = "array",
            categoryarray = x_categories,
            title = ""
          ),
          dragmode = FALSE  # 禁止拖拽，避免误触改变顺序
        )

        return(ply)

      }, error = function(e) {
        message(sprintf("❌ Boxplot error: %s", e$message))
        return(plotly::plot_ly() %>% plotly::layout(
          title = sprintf("Error: %s", e$message)
        ))
      })
    })



  })
}
