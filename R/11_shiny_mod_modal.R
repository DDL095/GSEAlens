#' @title Pathway Detail Modal UI
#' @keywords internal
mod_pathway_modal_ui <- function(id) {
  ns <- shiny::NS(id)
  NULL  # Modal dynamically generated in server
}



#' @title Pathway Detail Modal Server - Fixed for DESeq2 DFrame compatibility
#' @description Forces raw CPM usage, handles DFrame/DESeq2 colData properly,
#'   ensures sample name matching between expression matrix and metadata
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
            shiny::h4("Core Gene Expression Heatmap"),
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

    # 通路核心基因热图绘制
    output$modal_heatmap <- shiny::renderPlot({
      pdata <- current_data()
      shiny::req(pdata)

      # ==================== 1. 样本元数据获取与标准化 ====================
      sample_meta_raw <- gsea_res$expr_bundle$sample_meta

      # 处理 DFrame（DESeq2 的 colData 返回类型）
      if (inherits(sample_meta_raw, "DFrame") || inherits(sample_meta_raw, "DataFrame")) {
        sample_names <- rownames(sample_meta_raw)
        sample_meta <- as.data.frame(sample_meta_raw, stringsAsFactors = FALSE)
        if (!is.null(sample_names) && (is.null(rownames(sample_meta)) || all(rownames(sample_meta) == ""))) {
          rownames(sample_meta) <- sample_names
        }
      } else {
        sample_meta <- as.data.frame(sample_meta_raw, stringsAsFactors = FALSE)
      }

      group_col <- NULL
      if (!is.null(gsea_res$backend_info$target_factor)) {
        # 优先使用构建时记录的目标因子（如"分组"）
        target <- gsea_res$backend_info$target_factor
        if (target %in% colnames(sample_meta)) {
          group_col <- target
        }
      }

      # 回退：检测常见分组列名
      if (is.null(group_col)) {
        candidates <- c("group", "Group", "condition", "Condition",
                        "treatment", "Treatment", "type", "Type")
        for (col in candidates) {
          if (col %in% colnames(sample_meta)) {
            group_col <- col
            break
          }
        }
      }

      if (is.null(group_col)) {
        stop(sprintf("Cannot find group column. Available columns: %s",
                     paste(colnames(sample_meta), collapse = ", ")))
      }

      # 统一转换为字符型（处理 Factor）
      sample_meta$group <- as.character(sample_meta[[group_col]])

      # ==================== 2. 样本筛选 ====================
      left_grp <- pdata$data_list$left_group
      right_grp <- pdata$data_list$right_group

      sample_idx <- which(sample_meta$group %in% c(left_grp, right_grp))
      if (length(sample_idx) == 0) {
        stop(sprintf("No samples found for groups '%s' or '%s'. Available groups: %s",
                     left_grp, right_grp,
                     paste(unique(sample_meta$group), collapse = ", ")))
      }

      target_samples <- rownames(sample_meta)[sample_idx]

      # ==================== 3. CPM 矩阵计算（保持原有逻辑） ====================
      expr_bundle <- gsea_res$expr_bundle
      cpm_mat <- NULL

      # 策略 A：从 dds_obj 计算（DESeq2 流程）
      if (!is.null(expr_bundle$dds_obj)) {
        tryCatch({
          raw_counts <- DESeq2::counts(expr_bundle$dds_obj, normalized = FALSE)
          # 确保样本匹配
          common_samples <- intersect(target_samples, colnames(raw_counts))
          if (length(common_samples) > 0) {
            raw_counts <- raw_counts[, common_samples, drop = FALSE]
            lib_sizes <- colSums(raw_counts)
            cpm_mat <- t(t(raw_counts) / lib_sizes) * 1e6
            target_samples <- common_samples  # 更新为实际可用的样本
            target_samples <- common_samples
          }
        }, error = function(e) {
          message("DDS CPM calculation failed: ", e$message)
        })
      }

      # 策略 B：从 raw_counts 计算（通用流程）
      if (is.null(cpm_mat) && !is.null(expr_bundle$raw_counts)) {
        common_samples <- intersect(target_samples, colnames(expr_bundle$raw_counts))
        if (length(common_samples) > 0) {
          rc <- expr_bundle$raw_counts[, common_samples, drop = FALSE]
          lib_sizes <- colSums(rc)
          cpm_mat <- t(t(rc) / lib_sizes) * 1e6
          target_samples <- common_samples
        }
      }

      if (is.null(cpm_mat) || ncol(cpm_mat) == 0) {
        stop("Cannot compute CPM matrix. Please check expression data")
      }

      pathway_genes <- pdata$pathway_genes

      # 获取基因标识符：优先使用 gene_meta 行名，否则使用 cpm_mat 行名
      gene_identifiers <- rownames(cpm_mat)
      if (!is.null(expr_bundle$gene_meta) && !is.null(rownames(expr_bundle$gene_meta))) {
        gene_identifiers <- rownames(expr_bundle$gene_meta)
      }

      # 🔧 防御性检查：检测 Ensembl ID（禁止直接使用，必须转换为 SYMBOL）
      if (any(grepl("^ENS(MUS)?G[0-9]+", gene_identifiers, ignore.case = TRUE))) {
        stop("Detected Ensembl IDs (e.g., ENSMUSG...) as row names.\n",
             "Heatmap requires gene symbols. Please ensure row names are gene symbols.")
      }

      # 🔧 跨物种匹配：统一转为大写进行匹配（处理人类/小鼠大小写差异，如 Cyp2e1 <-> CYP2E1）
      expr_genes_upper <- toupper(gene_identifiers)
      pathway_genes_upper <- toupper(pathway_genes)

      matched_idx <- which(expr_genes_upper %in% pathway_genes_upper)

      if (length(matched_idx) < 2) {
        stop(sprintf("Gene matching failed: only %d/%d pathway genes matched.\n",
                     length(matched_idx), length(pathway_genes)),
             "Possible reasons:\n",
             "1. Species mismatch\n",
             "2. Row names are not gene symbols\n",
             "3. Gene name case mismatch")
      }

      message(sprintf("Successfully matched %d/%d pathway genes", length(matched_idx), length(pathway_genes)))

      # 提取表达矩阵（保留原始 CPM 值用于显示）
      plot_mat <- cpm_mat[matched_idx, , drop = FALSE]
      # 确保行名与通路基因一致（使用原始大小写中第一个匹配的）
      rownames(plot_mat) <- gene_identifiers[matched_idx]

      # ==================== 5. 数据标准化与可视化（保持原有 ComplexHeatmap 逻辑） ====================
      # 移除低方差基因
      gene_vars <- apply(plot_mat, 1, var, na.rm = TRUE)
      valid_genes <- gene_vars > 1e-6 & !is.na(gene_vars)
      if (sum(valid_genes) < 2) {
        stop("Insufficient genes with variance after filtering")
      }
      plot_mat <- plot_mat[valid_genes, , drop = FALSE]

      # 保存原始 CPM 用于热图细胞格显示（取整）
      display_numbers <- round(plot_mat)

      # Z-score 标准化（行方向，用于颜色映射）
      z_mat <- t(scale(t(plot_mat)))
      z_mat[is.na(z_mat)] <- 0
      z_mat[z_mat > 1] <- 1
      z_mat[z_mat < -1] <- -1

      # 获取 GSEA geneList 用于基因排序
      gene_list <- pdata$data_list$gsea_res@geneList

      # 计算每个基因的 Ranking Metric（用于排序）
      gene_metrics <- sapply(rownames(plot_mat), function(g) {
        idx <- match(toupper(g), toupper(names(gene_list)))
        if (is.na(idx)) return(0)
        return(gene_list[idx])
      })

      # 确定 Leading Edge 基因
      res_df <- as.data.frame(pdata$data_list$gsea_res@result)
      core_str <- res_df$core_enrichment[res_df$ID == pdata$pathway_id]
      core_genes <- character(0)
      if (length(core_str) > 0 && !is.na(core_str[1])) {
        core_genes <- unlist(strsplit(as.character(core_str[1]), "/"))
      }
      is_leading <- toupper(rownames(plot_mat)) %in% toupper(core_genes)

      # 根据 NES 方向排序基因
      if (pdata$nes > 0) {
        sort_order <- order(gene_metrics, decreasing = TRUE)
      } else {
        sort_order <- order(gene_metrics, decreasing = FALSE)
      }

      plot_mat <- plot_mat[sort_order, , drop = FALSE]
      z_mat <- z_mat[sort_order, , drop = FALSE]
      display_numbers <- display_numbers[sort_order, , drop = FALSE]
      is_leading <- is_leading[sort_order]
      gene_metrics <- gene_metrics[sort_order]

      # ComplexHeatmap 配置（保持原有美学设置）
      col_fun <- circlize::colorRamp2(
        c(-1, 0, 1),
        c("#67a9cf", "#f7f7f7", "#ef8a62")
      )

      # 分组颜色（保持红蓝配色）
      grp_col <- c("#E41A1C", "#377EB8")
      names(grp_col) <- c(left_grp, right_grp)

      # 确保分组因子水平顺序正确（左组在前）
      group_factor <- factor(
        sample_meta$group[match(colnames(plot_mat), rownames(sample_meta))],
        levels = c(left_grp, right_grp)
      )

      top_ann <- ComplexHeatmap::HeatmapAnnotation(
        Group = group_factor,
        col = list(Group = grp_col),
        annotation_name_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        simple_anno_size = grid::unit(0.6, "cm")
      )

      # Leading Edge 注释
      leading_status <- ifelse(is_leading, "YES", "NO")
      leading_colors <- c("YES" = "#FF9800", "NO" = "transparent")

      right_ann <- ComplexHeatmap::rowAnnotation(
        LeadingEdge = leading_status,
        col = list(LeadingEdge = leading_colors),
        annotation_name_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        simple_anno_size = grid::unit(0.4, "cm")
      )

      # 单元格数值显示函数（显示原始 CPM）
      layer_fun <- function(j, i, x, y, w, h, fill) {
        vals <- ComplexHeatmap::pindex(display_numbers, i, j)
        grid::grid.text(
          label = vals,
          x = x, y = y,
          gp = grid::gpar(fontsize = 13, col = "black", fontface = "bold")
        )
      }

      # 行分割（Leading Edge 置顶）
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

      # 基因名样式（Leading Edge 加粗橙色）
      row_names_gp <- grid::gpar(
        fontsize = 15,
        fontface = ifelse(is_leading, "bold", "plain"),
        col = ifelse(is_leading, "#FF9800", "black")
      )

      # 构建 Heatmap
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
        layer_fun = layer_fun,
        row_names_gp = row_names_gp,
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

      # 标题
      enriched_group <- ifelse(pdata$nes > 0, left_grp, right_grp)
      title_text <- sprintf(
        "%s\nRow-Scaled Z-Score [-1, 1] | Raw CPM Values Shown | Enriched in: %s",
        pdata$pathway_id,
        enriched_group
      )

      ComplexHeatmap::draw(
        ht,
        merge_legend = TRUE,
        column_title = title_text,
        column_title_gp = grid::gpar(fontsize = 14, fontface = "bold")
      )

    }, height = function() {
      # 动态高度计算
      pdata <- current_data()
      if (is.null(pdata)) return(400)
      n_genes <- length(pdata$pathway_genes)
      # 每个基因 25px + 基础高度 150px
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
        "-"
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
          backgroundColor = DT::styleEqual(c("YES", "-"), c("#FF9800", "transparent")),
          fontWeight = DT::styleEqual("YES", "bold"),
          color = DT::styleEqual("YES", "white")
        )
    })
  })
}
