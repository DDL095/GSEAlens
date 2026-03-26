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

    # 🔧 修复：强制使用原始CPM，标题显示通路名称
    output$modal_heatmap <- shiny::renderPlot({
      pdata <- current_data()
      shiny::req(pdata)

      # 获取样本元数据（用于分组）
      sample_meta <- get_sample_meta(gsea_res)
      shiny::req(sample_meta)

      # 筛选样本（仅当前对比组的两组）
      left_grp <- pdata$data_list$left_group
      right_grp <- pdata$data_list$right_group
      groups <- c(left_grp, right_grp)

      sample_idx <- which(sample_meta$group %in% groups)
      if (length(sample_idx) == 0) return()

      target_samples <- rownames(sample_meta)[sample_idx]

      # ==================== 强制使用原始CPM（不受expression_type影响） ====================
      # 优先从dge_list获取CPM，否则从raw_counts计算
      expr_bundle <- gsea_res$expr_bundle

      if (!is.null(expr_bundle$dge_list)) {
        # 使用edgeR::cpm计算原始CPM（log=FALSE确保非log2）
        cpm_mat <- edgeR::cpm(expr_bundle$dge_list, log = FALSE)
        message("📊 热图使用 edgeR::cpm(dge_list, log=FALSE) 计算原始CPM")
      } else if (!is.null(expr_bundle$raw_counts)) {
        # 从raw_counts手动计算CPM
        raw_counts <- expr_bundle$raw_counts[, target_samples, drop = FALSE]
        lib_sizes <- colSums(raw_counts)
        cpm_mat <- t(t(raw_counts) / lib_sizes) * 1e6
        message("📊 热图从 raw_counts 手动计算原始CPM")
      } else {
        # 备用方案：尝试从expression_type获取（但会警告）
        expr_mat <- get_expr_matrix(gsea_res, type = "cpm")
        if (is.null(expr_mat)) {
          stop("无法获取CPM数据：dge_list和raw_counts均不可用")
        }
        cpm_mat <- expr_mat
        warning("使用备用CPM数据源")
      }

      # 确保只保留目标样本
      cpm_mat <- cpm_mat[, target_samples, drop = FALSE]
      # ==================== CPM数据准备完成 ====================

      # 获取基因注释
      gene_meta <- expr_bundle$gene_meta

      # ==================== 基因映射逻辑（基于箱形图策略） ====================
      pathway_genes <- pdata$pathway_genes

      # 检测 pathway_genes 的 ID 类型
      is_ensembl <- all(grepl("^ENS(MUS)?G", pathway_genes))

      if (is_ensembl) {
        # 🧬 策略A: GSEA使用Ensembl ID
        matched_mask <- rownames(cpm_mat) %in% pathway_genes
        if (sum(matched_mask) == 0) return()

        plot_mat <- cpm_mat[matched_mask, , drop = FALSE]

        # 尝试通过 gene_meta 获取SYMBOL（保留原始大小写）
        if (!is.null(gene_meta) && nrow(gene_meta) > 0) {
          symbol_candidates <- c("SYMBOL", "symbol", "Gene", "gene_name",
                                 "gene_symbol", "Gene.Symbol", "GeneName")
          symbol_col <- intersect(symbol_candidates, colnames(gene_meta))[1]

          if (!is.na(symbol_col)) {
            ens_to_symbol <- setNames(
              as.character(gene_meta[[symbol_col]]),
              rownames(gene_meta)
            )

            current_ens_ids <- rownames(plot_mat)
            symbol_names <- ens_to_symbol[current_ens_ids]

            # 处理未匹配的（保留ENS ID）
            na_idx <- is.na(symbol_names)
            symbol_names[na_idx] <- current_ens_ids[na_idx]

            # 处理重复SYMBOL（多个ENS对应一个SYMBOL）：取均值
            if (any(duplicated(symbol_names[!na_idx]))) {
              unique_syms <- unique(symbol_names)
              merged_mat <- do.call(rbind, lapply(unique_syms, function(sym) {
                rows <- plot_mat[symbol_names == sym, , drop = FALSE]
                if (nrow(rows) > 1) {
                  colMeans(rows, na.rm = TRUE)
                } else {
                  as.numeric(rows[1, ])
                }
              }))
              plot_mat <- merged_mat
              rownames(plot_mat) <- unique_syms
              colnames(plot_mat) <- target_samples
            } else {
              rownames(plot_mat) <- symbol_names
            }
          }
        }

      } else {
        # 🧬 策略B: GSEA使用SYMBOL，需通过gene_meta映射到ENS ID
        if (is.null(gene_meta) || nrow(gene_meta) == 0) return()

        symbol_candidates <- c("SYMBOL", "symbol", "Gene", "gene_name",
                               "gene_symbol", "Gene.Symbol", "GeneName")
        symbol_col <- intersect(symbol_candidates, colnames(gene_meta))[1]
        if (is.na(symbol_col)) return()

        meta_symbols <- as.character(gene_meta[[symbol_col]])

        # 为每个 pathway_gene（保留原始大小写）提取表达数据
        gene_expr_list <- lapply(pathway_genes, function(sym) {
          match_idx <- which(toupper(meta_symbols) == toupper(sym))
          if (length(match_idx) == 0) return(NULL)

          ens_ids <- rownames(gene_meta)[match_idx]
          valid_ens <- ens_ids[ens_ids %in% rownames(cpm_mat)]
          if (length(valid_ens) == 0) return(NULL)

          sub_mat <- cpm_mat[valid_ens, , drop = FALSE]

          # 多个ENS对应一个SYMBOL：取均值（na.rm = TRUE）
          if (nrow(sub_mat) > 1) {
            expr_vals <- colMeans(sub_mat, na.rm = TRUE)
          } else {
            expr_vals <- as.numeric(sub_mat[1, ])
          }

          return(list(
            symbol = sym,  # 保留原始大小写
            values = expr_vals
          ))
        })

        # 过滤未匹配的基因
        gene_expr_list <- gene_expr_list[!sapply(gene_expr_list, is.null)]
        if (length(gene_expr_list) == 0) return()

        # 构建表达矩阵（行名为原始SYMBOL）
        plot_mat <- do.call(rbind, lapply(gene_expr_list, `[[`, "values"))
        rownames(plot_mat) <- sapply(gene_expr_list, `[[`, "symbol")
        colnames(plot_mat) <- target_samples
      }

      # 移除方差为0或NA的基因
      gene_vars <- apply(plot_mat, 1, var, na.rm = TRUE)
      plot_mat <- plot_mat[gene_vars > 1e-6 & !is.na(gene_vars), , drop = FALSE]
      if (nrow(plot_mat) < 2) return()
      # ==================== 基因映射逻辑结束 ====================

      # 保存原始CPM用于显示（取整）
      display_numbers <- round(plot_mat)

      # Z-score标准化（用于颜色映射，行方向）
      z_mat <- t(scale(t(plot_mat)))
      z_mat[is.na(z_mat)] <- 0
      z_mat[z_mat > 1] <- 1
      z_mat[z_mat < -1] <- -1

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

      # 根据NES方向排序
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

      # ComplexHeatmap配置
      col_fun <- circlize::colorRamp2(
        c(-1, 0, 1),
        c("#67a9cf", "#f7f7f7", "#ef8a62")
      )

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

      leading_status <- ifelse(is_leading, "YES", "NO")
      leading_colors <- c("YES" = "#FF9800", "NO" = "transparent")

      right_ann <- ComplexHeatmap::rowAnnotation(
        LeadingEdge = leading_status,
        col = list(LeadingEdge = leading_colors),
        annotation_name_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        simple_anno_size = grid::unit(0.4, "cm")
      )

      # 使用pindex正确处理row_split后的索引
      layer_fun <- function(j, i, x, y, w, h, fill) {
        vals <- ComplexHeatmap::pindex(display_numbers, i, j)
        grid::grid.text(
          label = vals,
          x = x, y = y,
          gp = grid::gpar(fontsize = 13, col = "black", fontface = "bold")
        )
      }

      # row_split配置
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

      # 构建Heatmap
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

      # 🔧 修复：动态标题显示通路名称和实际数据类型
      enriched_group <- ifelse(pdata$nes > 0, left_grp, right_grp)
      title_text <- sprintf(
        "%s\nRow-Scaled Z-Score [-1, 1] | Raw CPM Values Shown | Enriched in: %s",
        pdata$pathway_id,  # 加入通路名称
        enriched_group
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



#' @title Pathway Detail Modal Server - Fixed CPM Logic
#' @description Forces raw CPM usage regardless of expression_type selection,
#'   detects and reverses log2 transformation if necessary, displays pathway name in title
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

      # 显示弹窗（修改标题，移除固定的CPM字样）
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
            # 🔧 修改：移除固定的"(CPM)"，改为通用描述
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

    # 🔧 核心修复：强制使用原始CPM，检测并反转换log2数据
    output$modal_heatmap <- shiny::renderPlot({
      pdata <- current_data()
      shiny::req(pdata)

      # 获取样本元数据
      sample_meta <- get_sample_meta(gsea_res)
      shiny::req(sample_meta)

      # 筛选样本（仅当前对比组的两组）
      left_grp <- pdata$data_list$left_group
      right_grp <- pdata$data_list$right_group
      groups <- c(left_grp, right_grp)

      sample_idx <- which(sample_meta$group %in% groups)
      if (length(sample_idx) == 0) return()

      target_samples <- rownames(sample_meta)[sample_idx]

      # ==================== 强制原始CPM获取逻辑 ====================
      expr_bundle <- gsea_res$expr_bundle

      # 策略1：优先使用DGEList计算CPM（最可靠，不受expression_type影响）
      if (!is.null(expr_bundle$dge_list)) {
        # 使用edgeR::cpm计算原始CPM（确保log=FALSE）
        cpm_mat <- edgeR::cpm(expr_bundle$dge_list, log = FALSE)
        message(sprintf("📊 [%s] Using edgeR::cpm(dge_list, log=FALSE)", pdata$pathway_id))

        # 防御性检查：如果看起来像是log2转换过的数据（有负值或值域过小），进行反转换
        if (max(cpm_mat, na.rm = TRUE) < 50 || any(cpm_mat < 0, na.rm = TRUE)) {
          message("⚠️  Detected possible log2-transformed CPM, reversing transformation...")
          cpm_mat <- 2^cpm_mat - 1  # 反转换：假设是log2(x+1)
          # 再次计算CPM确保数值正确
          cpm_mat <- edgeR::cpm(cpm_mat, log = FALSE)
        }

      } else if (!is.null(expr_bundle$raw_counts)) {
        # 策略2：从raw_counts手动计算CPM
        raw_counts <- expr_bundle$raw_counts[, target_samples, drop = FALSE]

        # 检查raw_counts是否可能是log2转换过的（防御性编程）
        if (max(raw_counts, na.rm = TRUE) < 50 || any(raw_counts < 0, na.rm = TRUE)) {
          message("⚠️  raw_counts appears log2-transformed, reversing...")
          raw_counts <- 2^raw_counts - 1
        }

        lib_sizes <- colSums(raw_counts)
        cpm_mat <- t(t(raw_counts) / lib_sizes) * 1e6
        message(sprintf("📊 [%s] Calculated CPM from raw_counts", pdata$pathway_id))

      } else {
        # 备用方案：尝试从display_expr反转换（不推荐，但兼容）
        stop("Cannot generate heatmap: dge_list and raw_counts both unavailable")
      }

      # 确保只保留目标样本
      cpm_mat <- cpm_mat[, target_samples, drop = FALSE]

      # 最终验证：确保不是log2数据（CPM应该都是正数，且通常>1）
      if (any(cpm_mat < 0, na.rm = TRUE) || max(cpm_mat, na.rm = TRUE) < 20) {
        warning(sprintf("Expression values may still be log2-transformed (max=%.2f)",
                        max(cpm_mat, na.rm = TRUE)))
      }
      # ==================== CPM数据准备完成 ====================

      # 获取基因注释
      gene_meta <- expr_bundle$gene_meta
      pathway_genes <- pdata$pathway_genes

      # ==================== 基因映射逻辑（复用箱形图策略） ====================
      is_ensembl <- all(grepl("^ENS(MUS)?G", pathway_genes))

      if (is_ensembl) {
        # GSEA使用Ensembl ID
        matched_mask <- rownames(cpm_mat) %in% pathway_genes
        if (sum(matched_mask) == 0) return()

        plot_mat <- cpm_mat[matched_mask, , drop = FALSE]

        # 尝试映射到SYMBOL
        if (!is.null(gene_meta) && nrow(gene_meta) > 0) {
          symbol_candidates <- c("SYMBOL", "symbol", "Gene", "gene_name",
                                 "gene_symbol", "Gene.Symbol", "GeneName")
          symbol_col <- intersect(symbol_candidates, colnames(gene_meta))[1]

          if (!is.na(symbol_col)) {
            ens_to_symbol <- setNames(
              as.character(gene_meta[[symbol_col]]),
              rownames(gene_meta)
            )

            current_ens_ids <- rownames(plot_mat)
            symbol_names <- ens_to_symbol[current_ens_ids]
            na_idx <- is.na(symbol_names)
            symbol_names[na_idx] <- current_ens_ids[na_idx]

            # 处理重复SYMBOL：取均值
            if (any(duplicated(symbol_names[!na_idx]))) {
              unique_syms <- unique(symbol_names)
              merged_mat <- do.call(rbind, lapply(unique_syms, function(sym) {
                rows <- plot_mat[symbol_names == sym, , drop = FALSE]
                if (nrow(rows) > 1) {
                  colMeans(rows, na.rm = TRUE)
                } else {
                  as.numeric(rows[1, ])
                }
              }))
              plot_mat <- merged_mat
              rownames(plot_mat) <- unique_syms
              colnames(plot_mat) <- target_samples
            } else {
              rownames(plot_mat) <- symbol_names
            }
          }
        }

      } else {
        # GSEA使用SYMBOL，映射到ENS ID
        if (is.null(gene_meta) || nrow(gene_meta) == 0) return()

        symbol_candidates <- c("SYMBOL", "symbol", "Gene", "gene_name",
                               "gene_symbol", "Gene.Symbol", "GeneName")
        symbol_col <- intersect(symbol_candidates, colnames(gene_meta))[1]
        if (is.na(symbol_col)) return()

        meta_symbols <- as.character(gene_meta[[symbol_col]])

        gene_expr_list <- lapply(pathway_genes, function(sym) {
          match_idx <- which(toupper(meta_symbols) == toupper(sym))
          if (length(match_idx) == 0) return(NULL)

          ens_ids <- rownames(gene_meta)[match_idx]
          valid_ens <- ens_ids[ens_ids %in% rownames(cpm_mat)]
          if (length(valid_ens) == 0) return(NULL)

          sub_mat <- cpm_mat[valid_ens, , drop = FALSE]

          if (nrow(sub_mat) > 1) {
            expr_vals <- colMeans(sub_mat, na.rm = TRUE)
          } else {
            expr_vals <- as.numeric(sub_mat[1, ])
          }

          return(list(symbol = sym, values = expr_vals))
        })

        gene_expr_list <- gene_expr_list[!sapply(gene_expr_list, is.null)]
        if (length(gene_expr_list) == 0) return()

        plot_mat <- do.call(rbind, lapply(gene_expr_list, `[[`, "values"))
        rownames(plot_mat) <- sapply(gene_expr_list, `[[`, "symbol")
        colnames(plot_mat) <- target_samples
      }

      # 移除低方差基因
      gene_vars <- apply(plot_mat, 1, var, na.rm = TRUE)
      plot_mat <- plot_mat[gene_vars > 1e-6 & !is.na(gene_vars), , drop = FALSE]
      if (nrow(plot_mat) < 2) return()

      # 保存原始CPM用于显示（取整）
      display_numbers <- round(plot_mat)

      # Z-score标准化（用于颜色）
      z_mat <- t(scale(t(plot_mat)))
      z_mat[is.na(z_mat)] <- 0
      z_mat[z_mat > 1] <- 1
      z_mat[z_mat < -1] <- -1

      # 获取GSEA geneList用于排序
      gene_list <- pdata$data_list$gsea_res@geneList

      gene_metrics <- sapply(rownames(plot_mat), function(g) {
        idx <- match(toupper(g), toupper(names(gene_list)))
        if (is.na(idx)) return(0)
        return(gene_list[idx])
      })

      is_leading <- toupper(rownames(plot_mat)) %in% toupper(pdata$core_genes)

      # 根据NES方向排序
      if (pdata$nes > 0) {
        sort_order <- order(gene_metrics, decreasing = TRUE)
      } else {
        sort_order <- order(gene_metrics, decreasing = FALSE)
      }

      plot_mat <- plot_mat[sort_order, , drop = FALSE]
      z_mat <- z_mat[sort_order, , drop = FALSE]
      display_numbers <- display_numbers[sort_order, , drop = FALSE]
      is_leading <- is_leading[sort_order]

      # ComplexHeatmap配置
      col_fun <- circlize::colorRamp2(
        c(-1, 0, 1),
        c("#67a9cf", "#f7f7f7", "#ef8a62")
      )

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

      leading_status <- ifelse(is_leading, "YES", "NO")
      leading_colors <- c("YES" = "#FF9800", "NO" = "transparent")

      right_ann <- ComplexHeatmap::rowAnnotation(
        LeadingEdge = leading_status,
        col = list(LeadingEdge = leading_colors),
        annotation_name_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        simple_anno_size = grid::unit(0.4, "cm")
      )

      # 使用pindex正确处理row_split后的索引
      layer_fun <- function(j, i, x, y, w, h, fill) {
        vals <- ComplexHeatmap::pindex(display_numbers, i, j)
        grid::grid.text(
          label = vals,
          x = x, y = y,
          gp = grid::gpar(fontsize = 13, col = "black", fontface = "bold")
        )
      }

      # row_split配置
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

      # 构建Heatmap
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

      # 🔧 修复：动态标题显示通路名称和实际数据类型
      enriched_group <- ifelse(pdata$nes > 0, left_grp, right_grp)
      title_text <- sprintf(
        "%s\nRow-Scaled Z-Score [-1, 1] | Raw CPM Values Shown | Enriched in: %s",
        pdata$pathway_id,  # 加入通路名称
        enriched_group
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
