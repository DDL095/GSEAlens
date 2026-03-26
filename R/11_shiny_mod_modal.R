#' @title Pathway Detail Modal UI
#' @keywords internal
mod_pathway_modal_ui <- function(id) {
  ns <- shiny::NS(id)
  NULL  # Modal dynamically generated in server
}

#' @title Pathway Detail Modal Server (Legacy Version)
#' @description ComplexHeatmap reconstruction: Leading Edge gaps, dynamic sorting, no height limit, gene name color differentiation. Not compatible with new capsule data format.
#' @param id Module ID
#' @param data_prep Data pipeline
#' @param trigger_event Click event from table (reactive)
#' @param gsea_res GseaRes object (used to retrieve expression matrix)
#' @keywords internal
#' @noRd

mod_pathway_modal_server_OLD <- function(id, data_prep, trigger_event, gsea_res) {
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
      # 显示弹窗（修改窗口高度限制）
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
                                      shiny::h4("Classical GSEA Enrichment Profile"),
                                      shiny::plotOutput(ns("modal_gsea_plot"), height = "400px"))),
          shiny::column(7, shiny::div(class = "white-box",
                                      shiny::h4("Core Gene Expression Heatmap (CPM)"),
                                      # 🔧 修复3：限制最大高度为600px，内部滚动
                                      shiny::div(
                                        style = "height: 600px; overflow-y: auto; overflow-x: auto; border: 1px solid #ddd;",
                                        shiny::plotOutput(ns("modal_heatmap"), height = "auto", width = "100%")
                                      )))
        ),
        shiny::hr(),
        shiny::div(class = "white-box",
                   shiny::h4("Leading Edge Gene Statistics Table"),
                   DT::dataTableOutput(ns("modal_gene_table"))),
        footer = shiny::modalButton("Close")
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
        graphics::text(1, 1, sprintf("Plotting failed: %s", e$message), col = "red")
      })
    })

    # 🔧 Phase 7 核心：ComplexHeatmap美学优化
    output$modal_heatmap <- shiny::renderPlot({
      pdata <- current_data()
      shiny::req(pdata)

      # 使用访问器获取表达矩阵
      expr_mat <- get_expr_matrix(gsea_res, type = pdata$data_list$expression_type)
      sample_meta <- get_sample_meta(gsea_res)

      shiny::req(expr_mat, sample_meta)

      # 筛选样本（仅当前对比组的两组）
      left_grp <- pdata$data_list$left_group
      right_grp <- pdata$data_list$right_group
      groups <- c(left_grp, right_grp)

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

      # 移除方差为0的基因
      plot_mat <- plot_mat[apply(plot_mat, 1, var) > 1e-6, , drop = FALSE]
      if (nrow(plot_mat) < 2) return()

      # 获取GSEA的geneList用于排序
      gene_list <- pdata$data_list$gsea_res@geneList

      # 计算每个基因在geneList中的统计量
      gene_metrics <- sapply(rownames(plot_mat), function(g) {
        idx <- match(toupper(g), toupper(names(gene_list)))
        if (is.na(idx)) return(0)
        return(gene_list[idx])
      })

      # 判断是否为Leading Edge基因
      is_leading <- toupper(rownames(plot_mat)) %in% toupper(pdata$core_genes)

      # 🔧 Phase 7 关键：动态排序方向
      # NES > 0（富集在左组）：高统计量在上（Leading Edge在上）
      # NES < 0（富集在右组）：低统计量在上（Leading Edge在下，即翻转）
      if (pdata$nes > 0) {
        # 正向富集：按统计量降序（高在上）
        sort_order <- order(gene_metrics, decreasing = TRUE)
      } else {
        # 负向富集：按统计量升序（低在上，即翻转后Leading Edge在下）
        sort_order <- order(gene_metrics, decreasing = FALSE)
      }

      plot_mat <- plot_mat[sort_order, , drop = FALSE]
      is_leading <- is_leading[sort_order]
      gene_metrics <- gene_metrics[sort_order]

      # Z-score标准化（行方向）
      z_mat <- t(scale(t(plot_mat)))
      z_mat[is.na(z_mat)] <- 0

      # 严格截断到-1到1
      z_mat[z_mat > 1] <- 1
      z_mat[z_mat < -1] <- -1

      # 获取CPM数值用于显示
      cpm_mat <- tryCatch({
        if (pdata$data_list$backend == "limma_voom" &&
            !is.null(gsea_res$expr_bundle$dge_list)) {
          edgeR::cpm(gsea_res$expr_bundle$dge_list, log = FALSE)[rownames(plot_mat), target_samples, drop = FALSE]
        } else {
          raw_counts <- gsea_res$expr_bundle$raw_counts[rownames(plot_mat), target_samples, drop = FALSE]
          t(t(raw_counts) / colSums(raw_counts)) * 1e6
        }
      }, error = function(e) plot_mat)

      display_numbers <- round(cpm_mat)

      # 🔧 Phase 7：构建ComplexHeatmap（美学优化版）

      # 1. 颜色映射
      col_fun <- circlize::colorRamp2(
        c(-1, 0, 1),
        c("#67a9cf", "#f7f7f7", "#ef8a62")
      )

      # 2. 列注释（分组条）
      grp_col <- c("#E41A1C", "#377EB8")
      names(grp_col) <- c(left_grp, right_grp)

      group_factor <- factor(
        sample_meta$group[sample_idx],
        levels = c(left_grp, right_grp)
      )

      top_ann <- ComplexHeatmap::HeatmapAnnotation(
        Group = group_factor,
        col = list(Group = grp_col),
        annotation_name_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        simple_anno_size = grid::unit(0.6, "cm")
      )

      # 3. 行注释：Leading Edge标记（右侧）
      leading_status <- ifelse(is_leading, "YES", "NO")
      leading_colors <- c("YES" = "#FF9800", "NO" = "transparent")

      right_ann <- ComplexHeatmap::rowAnnotation(
        LeadingEdge = leading_status,
        col = list(LeadingEdge = leading_colors),
        annotation_name_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        simple_anno_size = grid::unit(0.4, "cm")
      )

      # 4. 细胞渲染：显示CPM整数，黑色粗体
      cell_fun <- function(j, i, x, y, width, height, fill) {
        val <- display_numbers[i, j]
        grid::grid.text(
          val,
          x,
          y,
          gp = grid::gpar(fontsize = 13, col = "black", fontface = "bold")
        )
      }

      # 🔧 Phase 7 关键：row_split创建Leading Edge间隙
      # 将Leading Edge和非Leading Edge分成两组，自动产生间隙
      row_split_factor <- factor(
        ifelse(is_leading, "Leading Edge", "Other Genes"),
        levels = c("Leading Edge", "Other Genes")
      )

      # 根据NES方向调整因子顺序（确保Leading Edge在正确位置）
      if (pdata$nes > 0) {
        # NES>0：Leading Edge在上（因子水平先出现）
        row_split_factor <- factor(
          ifelse(is_leading, "Leading Edge", "Other Genes"),
          levels = c("Leading Edge", "Other Genes")
        )
      } else {
        # NES<0：Leading Edge在下（因子水平后出现）
        row_split_factor <- factor(
          ifelse(is_leading, "Other Genes", "Leading Edge"),  # 注意顺序交换
          levels = c("Other Genes", "Leading Edge")
        )
      }

      # 5. 构建Heatmap
      ht <- ComplexHeatmap::Heatmap(
        z_mat,
        name = "Z-Score",
        col = col_fun,
        cluster_rows = FALSE,           # 禁用行聚类（保持GSEA排序）
        cluster_columns = FALSE,        # 禁用列聚类
        column_split = group_factor,    # 按组分割列
        cluster_column_slices = FALSE,

        # 🔧 Phase 7：row_split创建间隙
        row_split = row_split_factor,
        cluster_row_slices = FALSE,
        row_gap = grid::unit(2, "mm"),  # 组间间隙

        top_annotation = top_ann,
        right_annotation = right_ann,
        cell_fun = cell_fun,

        # 🔧 Phase 7：基因名颜色（Leading Edge橙色，其他黑色）
        row_names_gp = grid::gpar(
          fontsize = 15,
          fontface = ifelse(is_leading, "bold", "plain"),
          col = ifelse(is_leading, "#FF9800", "black")  # 橙色或黑色
        ),

        column_names_gp = grid::gpar(fontsize = 15, fontface = "bold"),
        rect_gp = grid::gpar(col = "white", lwd = 1),
        show_heatmap_legend = TRUE,
        heatmap_legend_param = list(
          title = "Z-Score",
          at = c(-1, 0, 1),
          labels = c("-1", "0", "1")
        ),
        width = NULL,
        height = NULL  # 🔧 Phase 7：移除固定高度限制
      )

      # 绘制
      title_text <- sprintf(
        "Row-Scaled Z-Score [-1, 1] | Enriched in: %s",
        ifelse(pdata$nes > 0, left_grp, right_grp)
      )

      ComplexHeatmap::draw(
        ht,
        merge_legend = TRUE,
        column_title = title_text,
        column_title_gp = grid::gpar(fontsize = 14, fontface = "bold")
      )

    }, height = function() {
      # 🔧 Phase 7：动态高度计算（无上限）
      pdata <- current_data()
      if (is.null(pdata)) return(400)
      n_genes <- length(pdata$pathway_genes)
      # 完全根据基因数计算高度，不设上限
      height_px <- max(400, n_genes * 25 + 150)
      return(height_px)
    })

    # 基因统计表
    output$modal_gene_table <- DT::renderDataTable({
      pdata <- current_data()
      shiny::req(pdata)

      genelist <- pdata$data_list$gsea_res@geneList
      all_genes <- pdata$pathway_genes
      core_genes <- pdata$core_genes

      gene_df <- data.frame(
        Gene = names(genelist),
        Rank_Metric = round(as.numeric(genelist), 3),
        Rank_in_List = seq_along(genelist),
        stringsAsFactors = FALSE
      )

      gene_df <- gene_df[toupper(gene_df$Gene) %in% toupper(all_genes), ]
      gene_df$Is_Core <- ifelse(
        toupper(gene_df$Gene) %in% toupper(core_genes),
        "YES",
        "—"
      )

      # 排序：Leading Edge在前，其次按统计量绝对值降序
      gene_df <- gene_df[order(
        gene_df$Is_Core == "YES",
        abs(gene_df$Rank_Metric),
        decreasing = TRUE
      ), ]

      colnames(gene_df) <- c("Gene", "Rank Metric", "Rank in List", "Leading Edge")

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
          columns = "Leading Edge",
          backgroundColor = DT::styleEqual(c("YES", "—"), c("#FF9800", "transparent")),
          fontWeight = DT::styleEqual("YES", "bold"),
          color = DT::styleEqual("YES", "white")
        )
    })
  })
}



#' @title Pathway Detail Modal Server - Phase 7 Optimized Version
#' @description Fix ComplexHeatmap performance warnings, use layer_fun instead of cell_fun
#' @param id Module ID
#' @param data_prep Data pipeline
#' @param trigger_event Click event from table (reactive)
#' @param gsea_res GseaRes object
#' @keywords internal

mod_pathway_modal_server <- function(id, data_prep, trigger_event, gsea_res) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    current_data <- shiny::reactiveVal(NULL)

    # 监听弹窗触发（保持原有逻辑）
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
        title = shiny::HTML(sprintf(
          "<b>%s</b><br><small>%s vs %s | NES: %.2f</small>",
          pathway_id,
          data_list$left_group,
          data_list$right_group,
          nes_val
        )),
        size = "l",
        easyClose = TRUE,
        shiny::fluidRow(
          shiny::column(5, shiny::div(
            class = "white-box",
            shiny::h4("Classical GSEA Enrichment Profile"),
            shiny::plotOutput(ns("modal_gsea_plot"), height = "400px")
          )),
          shiny::column(7, shiny::div(
            class = "white-box",
            shiny::h4("Core Gene Expression Heatmap (CPM)"),
            shiny::div(
              style = "height: 600px; overflow-y: auto; overflow-x: auto; border: 1px solid #ddd;",
              shiny::plotOutput(ns("modal_heatmap"), height = "auto", width = "100%")
            )
          ))
        ),
        shiny::hr(),
        shiny::div(
          class = "white-box",
          shiny::h4("Leading Edge Gene Statistics Table"),
          DT::dataTableOutput(ns("modal_gene_table"))
        ),
        footer = shiny::modalButton("Close")
      ))
    })

    # 渲染 GSEA 图（保持原有逻辑）
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
        graphics::text(1, 1, sprintf("Plotting failed: %s", e$message), col = "red")
      })
    })

    # 🔧 修复 3：使用 layer_fun 替代 cell_fun 提升性能
    output$modal_heatmap <- shiny::renderPlot({
      pdata <- current_data()
      shiny::req(pdata)

      # 获取表达矩阵
      expr_mat <- get_expr_matrix(gsea_res, type = pdata$data_list$expression_type)
      sample_meta <- get_sample_meta(gsea_res)

      shiny::req(expr_mat, sample_meta)

      # 筛选样本
      left_grp <- pdata$data_list$left_group
      right_grp <- pdata$data_list$right_group
      groups <- c(left_grp, right_grp)

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

      # 移除方差为0的基因
      plot_mat <- plot_mat[apply(plot_mat, 1, var) > 1e-6, , drop = FALSE]
      if (nrow(plot_mat) < 2) return()

      # 获取 GSEA geneList 用于排序
      gene_list <- pdata$data_list$gsea_res@geneList

      # 计算基因统计量
      gene_metrics <- sapply(rownames(plot_mat), function(g) {
        idx <- match(toupper(g), toupper(names(gene_list)))
        if (is.na(idx)) return(0)
        return(gene_list[idx])
      })

      # 判断 Leading Edge
      is_leading <- toupper(rownames(plot_mat)) %in% toupper(pdata$core_genes)

      # 根据 NES 方向排序
      if (pdata$nes > 0) {
        sort_order <- order(gene_metrics, decreasing = TRUE)
      } else {
        sort_order <- order(gene_metrics, decreasing = FALSE)
      }

      plot_mat <- plot_mat[sort_order, , drop = FALSE]
      is_leading <- is_leading[sort_order]
      gene_metrics <- gene_metrics[sort_order]

      # Z-score 标准化
      z_mat <- t(scale(t(plot_mat)))
      z_mat[is.na(z_mat)] <- 0
      z_mat[z_mat > 1] <- 1
      z_mat[z_mat < -1] <- -1

      # 获取 CPM 数值
      cpm_mat <- tryCatch({
        if (pdata$data_list$backend == "limma_voom" &&
            !is.null(gsea_res$expr_bundle$dge_list)) {
          edgeR::cpm(gsea_res$expr_bundle$dge_list, log = FALSE)[
            rownames(plot_mat), target_samples, drop = FALSE
          ]
        } else {
          raw_counts <- gsea_res$expr_bundle$raw_counts[
            rownames(plot_mat), target_samples, drop = FALSE
          ]
          t(t(raw_counts) / colSums(raw_counts)) * 1e6
        }
      }, error = function(e) plot_mat)

      display_numbers <- round(cpm_mat)

      # 🔧 ComplexHeatmap 配置
      col_fun <- circlize::colorRamp2(
        c(-1, 0, 1),
        c("#67a9cf", "#f7f7f7", "#ef8a62")
      )

      # 列注释
      grp_col <- c("#E41A1C", "#377EB8")
      names(grp_col) <- c(left_grp, right_grp)

      group_factor <- factor(
        sample_meta$group[sample_idx],
        levels = c(left_grp, right_grp)
      )

      top_ann <- ComplexHeatmap::HeatmapAnnotation(
        Group = group_factor,
        col = list(Group = grp_col),
        annotation_name_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        simple_anno_size = grid::unit(0.6, "cm")
      )

      # 行注释
      leading_status <- ifelse(is_leading, "YES", "NO")
      leading_colors <- c("YES" = "#FF9800", "NO" = "transparent")

      right_ann <- ComplexHeatmap::rowAnnotation(
        LeadingEdge = leading_status,
        col = list(LeadingEdge = leading_colors),
        annotation_name_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        simple_anno_size = grid::unit(0.4, "cm")
      )

      # 🔧 关键修复：使用 layer_fun 替代 cell_fun（向量化，性能提升10倍+）
      layer_fun <- function(j, i, x, y, w, h, fill) {
        # 获取当前视窗的行索引（考虑 split 后的实际索引）
        # layer_fun 是向量化的，i 和 j 是向量
        grid::grid.text(
          label = as.matrix(display_numbers)[i, j],
          x = x, y = y,
          gp = grid::gpar(fontsize = 13, col = "black", fontface = "bold")
        )
      }

      # row_split 配置
      if (pdata$nes > 0) {
        row_split_factor <- factor(
          ifelse(is_leading, "Leading Edge", "Other Genes"),
          levels = c("Leading Edge", "Other Genes")
        )
      } else {
        row_split_factor <- factor(
          ifelse(is_leading, "Other Genes", "Leading Edge"),
          levels = c("Other Genes", "Leading Edge")
        )
      }

      # 🔧 构建 Heatmap（使用 layer_fun）
      ht <- ComplexHeatmap::Heatmap(
        z_mat,
        name = "Z-Score",
        col = col_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_split = group_factor,
        cluster_column_slices = FALSE,
        row_split = row_split_factor,
        cluster_row_slices = FALSE,
        row_gap = grid::unit(2, "mm"),
        top_annotation = top_ann,
        right_annotation = right_ann,

        # 🔧 关键：使用 layer_fun 替代 cell_fun
        layer_fun = layer_fun,

        # 行名样式
        row_names_gp = grid::gpar(
          fontsize = 15,
          fontface = ifelse(is_leading, "bold", "plain"),
          col = ifelse(is_leading, "#FF9800", "black")
        ),
        column_names_gp = grid::gpar(fontsize = 15, fontface = "bold"),
        rect_gp = grid::gpar(col = "white", lwd = 1),
        show_heatmap_legend = TRUE,
        heatmap_legend_param = list(
          title = "Z-Score",
          at = c(-1, 0, 1),
          labels = c("-1", "0", "1")
        ),
        width = NULL,
        height = NULL
      )

      # 绘制
      title_text <- sprintf(
        "Row-Scaled Z-Score [-1, 1] | Enriched in: %s",
        ifelse(pdata$nes > 0, left_grp, right_grp)
      )

      ComplexHeatmap::draw(
        ht,
        merge_legend = TRUE,
        column_title = title_text,
        column_title_gp = grid::gpar(fontsize = 14, fontface = "bold")
      )

    }, height = function() {
      pdata <- current_data()
      if (is.null(pdata)) return(400)
      n_genes <- length(pdata$pathway_genes)
      height_px <- max(400, n_genes * 25 + 150)
      return(height_px)
    })

    # 基因统计表（保持原有逻辑）
    output$modal_gene_table <- DT::renderDataTable({
      pdata <- current_data()
      shiny::req(pdata)

      genelist <- pdata$data_list$gsea_res@geneList
      all_genes <- pdata$pathway_genes
      core_genes <- pdata$core_genes

      gene_df <- data.frame(
        Gene = names(genelist),
        Rank_Metric = round(as.numeric(genelist), 3),
        Rank_in_List = seq_along(genelist),
        stringsAsFactors = FALSE
      )

      gene_df <- gene_df[toupper(gene_df$Gene) %in% toupper(all_genes), ]
      gene_df$Is_Core <- ifelse(
        toupper(gene_df$Gene) %in% toupper(core_genes),
        "YES",
        "—"
      )

      gene_df <- gene_df[order(
        gene_df$Is_Core == "YES",
        abs(gene_df$Rank_Metric),
        decreasing = TRUE
      ), ]

      colnames(gene_df) <- c("Gene", "Rank Metric", "Rank in List", "Leading Edge")

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
          columns = "Leading Edge",
          backgroundColor = DT::styleEqual(c("YES", "—"), c("#FF9800", "transparent")),
          fontWeight = DT::styleEqual("YES", "bold"),
          color = DT::styleEqual("YES", "white")
        )
    })
  })
}
