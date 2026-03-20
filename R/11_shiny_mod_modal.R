#' @title 详情弹窗 UI
#' @keywords internal
mod_pathway_modal_ui <- function(id) {
  ns <- shiny::NS(id)
  NULL  # 弹窗通过 server 动态生成
}

#' @title 详情弹窗 Server
#' @description ComplexHeatmap重构，完全复刻pheatmap美学：Leading Edge排序、Z-score截断、细胞数值、行名样式
#' @param id 模块 ID
#' @param data_prep 数据流
#' @param trigger_event 来自表格的点击事件 (reactive)
#' @param gsea_res GseaRes 对象 (用于获取表达矩阵)
#' @keywords internal
mod_pathway_modal_server <- function(id, data_prep, trigger_event, gsea_res) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    current_data <- shiny::reactiveVal(NULL)

    # 监听弹窗触发
    shiny::observeEvent(trigger_event(), {
      pathway_id <- trigger_event()
      shiny::req(pathway_id)

      data_list <- data_prep()
      shiny::req(data_list)

      # 提取通路基因和核心基因
      gsea_obj <- data_list$gsea_res
      pathway_genes <- gsea_obj@geneSets[[pathway_id]]

      res_df <- as.data.frame(gsea_obj@result)
      core_str <- res_df$core_enrichment[res_df$ID == pathway_id]
      core_genes <- if (length(core_str) > 0 && !is.na(core_str[1])) {
        unlist(strsplit(as.character(core_str[1]), "/"))
      } else {
        character(0)
      }

      # 获取NES值判断富集方向
      nes_val <- res_df$NES[res_df$ID == pathway_id]

      current_data(list(
        pathway_id = pathway_id,
        pathway_genes = pathway_genes,
        core_genes = core_genes,
        nes = nes_val,
        data_list = data_list
      ))

      # 显示弹窗
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML(sprintf("<b>%s</b><br><small>%s vs %s | NES: %.2f</small>",
                                    pathway_id,
                                    data_list$left_group,
                                    data_list$right_group,
                                    nes_val)),
        size = "l",
        easyClose = TRUE,
        shiny::fluidRow(
          shiny::column(5, shiny::div(class = "white-box",
                                      shiny::h4(.tr("modal.title_gsea_profile")),
                                      shiny::plotOutput(ns("modal_gsea_plot"), height = "400px"))),
          shiny::column(7, shiny::div(class = "white-box",
                                      shiny::h4(.tr("modal.title_heatmap")),
                                      # 关键：外层div限制最大高度650px，支持内部滚动
                                      shiny::div(
                                        style = "height: 650px; overflow-y: auto; overflow-x: auto; border: 1px solid #ddd;",
                                        shiny::plotOutput(ns("modal_heatmap"), height = "auto", width = "100%")
                                      )))
        ),
        shiny::hr(),
        shiny::div(class = "white-box",
                   shiny::h4(.tr("modal.title_gene_table")),
                   DT::dataTableOutput(ns("modal_gene_table"))),
        footer = shiny::modalButton(.tr("modal.btn_close"))
      ))
    })

    # 渲染 GSEA 图
    output$modal_gsea_plot <- shiny::renderPlot({
      pdata <- current_data()
      shiny::req(pdata)

      tryCatch({
        print(plot_directional_gsea(
          directional_gsea_obj = list(
            gsea_res = pdata$data_list$gsea_res,
            meta = list(
              left_group = pdata$data_list$left_group,
              right_group = pdata$data_list$right_group
            )
          ),
          target_pathways = pdata$pathway_id,
          subPlot = 3,
          add_pval = TRUE
        ))
      }, error = function(e) {
        graphics::plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        graphics::text(1, 1, sprintf("绘图失败: %s", e$message), col = "red")
      })
    })

    # 渲染热图 - ComplexHeatmap完全复刻pheatmap美学
    output$modal_heatmap <- shiny::renderPlot({
      pdata <- current_data()
      shiny::req(pdata)

      # 使用访问器获取表达矩阵（Phase 2修复后的版本）
      expr_mat <- get_expr_matrix(gsea_res, type = pdata$data_list$expression_type)
      sample_meta <- get_sample_meta(gsea_res)

      shiny::req(expr_mat, sample_meta)

      # 筛选样本（仅当前对比组的两组）
      left_grp <- pdata$data_list$left_group
      right_grp <- pdata$data_list$right_group
      groups <- c(left_grp, right_grp)

      # 获取样本索引（使用访问器统一处理后的group列）
      sample_idx <- which(sample_meta$group %in% groups)
      if (length(sample_idx) == 0) return()

      target_samples <- rownames(sample_meta)[sample_idx]
      expr_mat <- expr_mat[, target_samples, drop = FALSE]

      # 筛选通路基因
      pathway_genes <- pdata$pathway_genes
      gene_idx <- which(toupper(rownames(expr_mat)) %in% toupper(pathway_genes))
      if (length(gene_idx) < 2) return()

      plot_genes <- rownames(expr_mat)[gene_idx]
      plot_mat <- expr_mat[plot_genes, , drop = FALSE]

      # 移除方差为0的基因（避免scale报错）
      plot_mat <- plot_mat[apply(plot_mat, 1, var) > 1e-6, , drop = FALSE]
      if (nrow(plot_mat) < 2) return()

      # 获取GSEA的geneList用于排序（Leading Edge逻辑）
      gene_list <- pdata$data_list$gsea_res@geneList

      # 计算每个基因在geneList中的统计量（用于排序）
      gene_metrics <- sapply(rownames(plot_mat), function(g) {
        idx <- match(toupper(g), toupper(names(gene_list)))
        if (is.na(idx)) return(0)
        return(gene_list[idx])
      })

      # 判断是否为Leading Edge基因
      is_leading <- toupper(rownames(plot_mat)) %in% toupper(pdata$core_genes)

      # 关键排序逻辑（复刻原始美学）：
      # - NES > 0（富集在左组）：Leading Edge基因排在上方（高统计量在前）
      # - NES < 0（富集在右组）：Leading Edge基因排在下方（低统计量在前，即翻转）
      if (pdata$nes > 0) {
        # 正向富集：按统计量降序，Leading Edge自然在上方
        sort_order <- order(is_leading, gene_metrics, decreasing = c(FALSE, TRUE))
      } else {
        # 负向富集：按统计量升序，Leading Edge在下方
        sort_order <- order(is_leading, gene_metrics, decreasing = c(FALSE, FALSE))
      }

      plot_mat <- plot_mat[sort_order, , drop = FALSE]
      is_leading <- is_leading[sort_order]
      gene_metrics <- gene_metrics[sort_order]

      # Z-score标准化（行方向）
      z_mat <- t(scale(t(plot_mat)))
      z_mat[is.na(z_mat)] <- 0

      # 关键：严格截断到-1到1（复刻原始pheatmap的Z-score范围）
      z_mat[z_mat > 1] <- 1
      z_mat[z_mat < -1] <- -1

      # 获取CPM数值用于显示（复刻原始美学：显示整数CPM）
      cpm_mat <- tryCatch({
        if (pdata$data_list$backend == "limma_voom" &&
            !is.null(gsea_res$expr_bundle$dge_list)) {
          edgeR::cpm(gsea_res$expr_bundle$dge_list, log = FALSE)[rownames(plot_mat), target_samples, drop = FALSE]
        } else {
          # DESeq2或其他：手动计算CPM
          raw_counts <- gsea_res$expr_bundle$raw_counts[rownames(plot_mat), target_samples, drop = FALSE]
          t(t(raw_counts) / colSums(raw_counts)) * 1e6
        }
      }, error = function(e) {
        # 如果CPM计算失败，返回原始值
        plot_mat
      })

      display_numbers <- round(cpm_mat)

      # 构建ComplexHeatmap（完全复刻pheatmap美学）

      # 1. 颜色映射：淡蓝-白-亮粉橙，对应-1, 0, 1
      col_fun <- circlize::colorRamp2(
        c(-1, 0, 1),
        c("#67a9cf", "#f7f7f7", "#ef8a62")
      )

      # 2. 列注释（分组条）：红蓝对比
      grp_col <- c("#E41A1C", "#377EB8")
      names(grp_col) <- c(left_grp, right_grp)

      # 构建分组因子（确保顺序：左组在前，右组在后）
      group_factor <- factor(
        sample_meta$group[sample_idx],
        levels = c(left_grp, right_grp)
      )

      top_ann <- ComplexHeatmap::HeatmapAnnotation(
        Group = group_factor,
        col = list(Group = grp_col),
        annotation_name_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        simple_anno_size = grid::unit(0.6, "cm"),
        show_legend = TRUE
      )

      # 3. 行注释：Leading Edge标记（右侧）
      leading_status <- ifelse(is_leading, "YES", "NO")
      leading_colors <- c("YES" = "#FF9800", "NO" = "#E0E0E0")

      right_ann <- ComplexHeatmap::rowAnnotation(
        LeadingEdge = leading_status,
        col = list(LeadingEdge = leading_colors),
        annotation_name_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        simple_anno_size = grid::unit(0.4, "cm"),
        show_legend = TRUE
      )

      # 4. 细胞渲染函数：显示CPM整数，黑色粗体（复刻原始美学）
      cell_fun <- function(j, i, x, y, width, height, fill) {
        val <- display_numbers[i, j]
        # 固定黑色字体，字号13，粗体（完全复刻原始设置）
        grid::grid.text(
          val,
          x,
          y,
          gp = grid::gpar(fontsize = 13, col = "black", fontface = "bold")
        )
      }

      # 5. 构建Heatmap对象（复刻所有原始参数）
      ht <- ComplexHeatmap::Heatmap(
        z_mat,
        name = "Z-Score",
        col = col_fun,
        cluster_rows = FALSE,           # 禁用行聚类（保持GSEA排序）
        cluster_columns = FALSE,        # 禁用列聚类
        column_split = group_factor,    # 按组分割列
        cluster_column_slices = FALSE,  # 不聚类分组
        top_annotation = top_ann,
        right_annotation = right_ann,
        cell_fun = cell_fun,            # 自定义细胞显示
        # 行名样式：Leading Edge加粗+橙色，其他普通
        row_names_gp = grid::gpar(
          fontsize = 15,
          fontface = ifelse(is_leading, "bold", "plain"),
          col = ifelse(is_leading, "#FF9800", "black")
        ),
        # 列名样式
        column_names_gp = grid::gpar(fontsize = 15, fontface = "bold"),
        # 白色网格线（复刻原始）
        rect_gp = grid::gpar(col = "white", lwd = 1),
        show_heatmap_legend = TRUE,
        # 图例设置
        heatmap_legend_param = list(
          title = "Z-Score",
          at = c(-1, 0, 1),
          labels = c("-1", "0", "1")
        ),
        # 尺寸自适应
        width = NULL,
        height = NULL
      )

      # 绘制（带标题）
      title_text <- sprintf(
        "Row-Scaled Z-Score [-1, 1]\nEnriched in: %s",
        ifelse(pdata$nes > 0, left_grp, right_grp)
      )

      ComplexHeatmap::draw(
        ht,
        merge_legend = TRUE,
        column_title = title_text,
        column_title_gp = grid::gpar(fontsize = 14, fontface = "bold")
      )

    }, height = function() {
      # 动态高度计算：每基因25px + 基础高度
      pdata <- current_data()
      if (is.null(pdata)) return(400)
      n_genes <- length(pdata$pathway_genes)
      # 限制最大高度（避免过大），最小400px
      height_px <- max(400, min(n_genes * 25 + 150, 2000))
      return(height_px)
    })

    # 基因统计表（Leading Edge表）
    output$modal_gene_table <- DT::renderDataTable({
      pdata <- current_data()
      shiny::req(pdata)

      genelist <- pdata$data_list$gsea_res@geneList
      all_genes <- pdata$pathway_genes
      core_genes <- pdata$core_genes

      # 构建完整的基因排名表
      gene_df <- data.frame(
        Gene = names(genelist),
        Rank_Metric = round(as.numeric(genelist), 3),
        Rank_in_List = seq_along(genelist),
        stringsAsFactors = FALSE
      )

      # 筛选通路内基因
      gene_df <- gene_df[toupper(gene_df$Gene) %in% toupper(all_genes), ]

      # 标记Leading Edge
      gene_df$Is_Core <- ifelse(
        toupper(gene_df$Gene) %in% toupper(core_genes),
        "✅ YES",
        "—"
      )

      # 排序：Leading Edge在前，其次按统计量（绝对值）降序
      gene_df <- gene_df[order(
        gene_df$Is_Core == "✅ YES",
        abs(gene_df$Rank_Metric),
        decreasing = TRUE
      ), ]

      # 重命名列（国际化）
      colnames(gene_df) <- c(
        .tr("modal.col_gene"),
        .tr("modal.col_rank_metric"),
        .tr("modal.col_rank_in_list"),
        .tr("modal.col_is_core")
      )

      DT::datatable(
        gene_df,
        rownames = FALSE,
        escape = FALSE,
        extensions = c('Scroller'),
        options = list(
          pageLength = -1,
          scrollY = "40vh",
          scroller = TRUE,
          dom = 'frtip'
        )
      ) %>%
        DT::formatStyle(
          columns = .tr("modal.col_is_core"),
          backgroundColor = DT::styleEqual(c("✅ YES", "—"), c("#FF9800", "transparent")),
          fontWeight = DT::styleEqual("✅ YES", "bold"),
          color = DT::styleEqual("✅ YES", "white")
        )
    })
  })
}
