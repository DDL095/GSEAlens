#' @title Pathway Detail Modal UI
#' @keywords internal
mod_pathway_modal_ui <- function(id) {
  ns <- shiny::NS(id)
  NULL  # Modal dynamically generated in server
}



#' @title Pathway Detail Modal Server - Enhanced with Gene Selection
#' @description Adds checkbox selection to Leading Edge Gene table with confirm button.
#'   Selected genes are passed to sidebar pending list (not directly applied).
#'   Buttons are fixed at the top of the table section using sticky positioning.
#' @param id Module ID
#' @param data_prep Reactive data from data prep module
#' @param trigger_event Event that triggers modal opening
#' @param gsea_res GseaRes object
#' @param pending_genes_reactive Reactive expression returning current pending genes
#' @param update_pending_genes Function to update pending genes list
#' @keywords internal

mod_pathway_modal_server <- function(id, data_prep, trigger_event, gsea_res,
                                     pending_genes_reactive, update_pending_genes) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    current_data <- shiny::reactiveVal(NULL)

    # Modal内部基因选择状态
    modal_selected_genes <- shiny::reactiveVal(character(0))

    # 监听弹窗触发
    shiny::observeEvent(trigger_event(), {
      pathway_id <- trigger_event()
      shiny::req(pathway_id)

      # 清空modal内部状态
      modal_selected_genes(character(0))

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

      # =========================================
      # Modal UI - 按钮固定在表格上方
      # =========================================
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

        # 第一行：GSEA图和热图
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

        # 基因选择区域
        shiny::div(
          class = "white-box",
          shiny::h4("Leading Edge Gene Statistics Table"),

          # CSS样式 - 按钮固定定位
          shiny::tags$head(shiny::tags$style(shiny::HTML(sprintf("
            #%s-modal-gene-section {
              position: relative;
            }
            #%s-modal-button-row {
              position: sticky;
              top: 0;
              background-color: #fafafa;
              padding: 15px 0;
              z-index: 10;
              border-bottom: 2px solid #dee2e6;
              margin-bottom: 10px;
            }
            .modal-gene-checkbox {
              width: 18px;
              height: 18px;
              cursor: pointer;
            }
          ", ns(""), ns(""))))),

          # 按钮行（sticky）
          shiny::div(
            id = ns("modal-button-row"),
            shiny::fluidRow(
              shiny::column(6,
                            shiny::div(
                              style = "color: #666; font-size: 13px;",
                              shiny::HTML("<b>Instruction:</b> Check genes in table below, then click <b>Confirm</b> to add to sidebar list")
                            )
              ),
              shiny::column(6,
                            shiny::div(
                              style = "display: flex; gap: 10px; justify-content: flex-end;",
                              shiny::actionButton(
                                ns("select_all_leading"),
                                label = "Select All Leading Edge",
                                class = "btn-info"
                              ),
                              shiny::actionButton(
                                ns("clear_modal_selection"),
                                label = "Clear",
                                class = "btn-warning"
                              ),
                              shiny::actionButton(
                                ns("confirm_gene_selection"),
                                label = "Confirm Selection & Add to Sidebar",
                                class = "btn-success btn-lg",
                                style = "font-weight: bold;"
                              )
                            )
              )
            )
          ),

          # 选择状态显示
          shiny::uiOutput(ns("modal_selection_status")),

          shiny::hr(),

          # 表格区域（可滚动）
          shiny::div(
            style = "max-height: 50vh; overflow-y: auto; border: 1px solid #ddd; border-radius: 5px; background: white;",
            DT::dataTableOutput(ns("modal_gene_table"))
          )
        ),

        footer = shiny::modalButton("Close")
      ))
    })

    # =========================================
    # Modal基因选择逻辑
    # =========================================

    # 监听checkbox点击
    shiny::observeEvent(input$modal_gene_toggle, {
      toggle <- input$modal_gene_toggle
      if (is.null(toggle) || !is.list(toggle)) return()

      # JavaScript传递的是字符串，需要转换
      is_checked <- identical(toggle$checked, TRUE) || identical(toggle$checked, "true")
      gene_name <- toggle$id

      current <- modal_selected_genes()

      if (is_checked) {
        if (!(gene_name %in% current)) {
          modal_selected_genes(c(current, gene_name))
          message(sprintf("[Modal] Gene selected: %s", gene_name))
        }
      } else {
        modal_selected_genes(setdiff(current, gene_name))
        message(sprintf("[Modal] Gene deselected: %s", gene_name))
      }
    })

    # Select All Leading Edge按钮
    shiny::observeEvent(input$select_all_leading, {
      pdata <- current_data()
      if (is.null(pdata)) return()

      leading_genes <- pdata$core_genes
      if (length(leading_genes) == 0) {
        shiny::showNotification("No leading edge genes found for this pathway", type = "warning")
        return()
      }

      # 更新modal内部选择状态
      modal_selected_genes(leading_genes)
      message(sprintf("[Modal] Selected all %d leading edge genes", length(leading_genes)))

      # 强制刷新表格UI（通过触发input更新）
      shiny::updateTextInput(session, "modal_gene_table_refresh", value = as.character(Sys.time()))
    })

    # Clear Modal Selection按钮
    shiny::observeEvent(input$clear_modal_selection, {
      modal_selected_genes(character(0))
      message("[Modal] Cleared gene selection")
    })

    # Confirm Selection按钮 - 关键逻辑：更新pending_genes
    shiny::observeEvent(input$confirm_gene_selection, {
      selected <- modal_selected_genes()

      if (length(selected) == 0) {
        shiny::showNotification("No genes selected! Please check genes in the table first.",
                                type = "warning", duration = 3)
        return()
      }

      # 获取当前pending genes并执行并集
      update_pending_genes(unique(c(selected, pending_genes_reactive())))

      message(sprintf("[Modal] Confirmed %d genes", length(selected)))

      # 清空modal内部状态
      modal_selected_genes(character(0))

      # 关闭弹窗
      shiny::removeModal()

      # 提示用户
      shiny::showNotification(
        sprintf("%d genes added to 'Select Genes of Interest' list! Click 'Confirm and Apply Gene Markers' in sidebar to apply.",
                length(selected)),
        type = "message",
        duration = 5
      )
    })

    # 显示当前选择状态
    output$modal_selection_status <- shiny::renderUI({
      # 响应 modal_selected_genes 变化
      modal_sel <- modal_selected_genes()
      # 响应 pending_genes 变化
      current_pending <- pending_genes_reactive()

      union_count <- length(unique(c(toupper(modal_sel), toupper(current_pending))))

      if (length(modal_sel) == 0) {
        shiny::div(
          style = "margin: 10px 0; padding: 10px; background: #f8f9fa; border-radius: 5px; color: #666;",
          shiny::HTML(sprintf(
            "<b>Selected:</b> 0 genes | <b>In sidebar pending:</b> %d genes | <b>After confirm:</b> %d genes",
            length(current_pending), union_count
          ))
        )
      } else {
        genes_preview <- paste(head(modal_sel, 6), collapse = ", ")
        if (length(modal_sel) > 6) {
          genes_preview <- paste0(genes_preview, sprintf(" (+%d more)", length(modal_sel) - 6))
        }
        shiny::div(
          style = "margin: 10px 0; padding: 10px; background: #d4edda; border-radius: 5px; color: #155724;",
          shiny::HTML(sprintf(
            "<b>Selected:</b> %d genes | <b>After confirm:</b> %d genes total<br><small>%s</small>",
            length(modal_sel), union_count, genes_preview
          ))
        )
      }
    })

    # =========================================
    # Modal内容渲染
    # =========================================

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

      sample_meta_raw <- gsea_res$expr_bundle$sample_meta

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
        target <- gsea_res$backend_info$target_factor
        if (target %in% colnames(sample_meta)) {
          group_col <- target
        }
      }

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

      sample_meta$group <- as.character(sample_meta[[group_col]])

      left_grp <- pdata$data_list$left_group
      right_grp <- pdata$data_list$right_group

      sample_idx <- which(sample_meta$group %in% c(left_grp, right_grp))
      if (length(sample_idx) == 0) {
        stop(sprintf("No samples found for groups '%s' or '%s'.",
                     left_grp, right_grp))
      }

      target_samples <- rownames(sample_meta)[sample_idx]

      expr_bundle <- gsea_res$expr_bundle
      cpm_mat <- NULL

      if (!is.null(expr_bundle$dds_obj)) {
        tryCatch({
          raw_counts <- DESeq2::counts(expr_bundle$dds_obj, normalized = FALSE)
          common_samples <- intersect(target_samples, colnames(raw_counts))
          if (length(common_samples) > 0) {
            raw_counts <- raw_counts[, common_samples, drop = FALSE]
            lib_sizes <- colSums(raw_counts)
            cpm_mat <- t(t(raw_counts) / lib_sizes) * 1e6
            target_samples <- common_samples
          }
        }, error = function(e) {
          message("DDS CPM calculation failed: ", e$message)
        })
      }

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
        stop("Cannot compute CPM matrix.")
      }

      pathway_genes <- pdata$pathway_genes

      gene_identifiers <- rownames(cpm_mat)
      if (!is.null(expr_bundle$gene_meta) && !is.null(rownames(expr_bundle$gene_meta))) {
        gene_identifiers <- rownames(expr_bundle$gene_meta)
      }

      if (any(grepl("^ENS(MUS)?G[0-9]+", gene_identifiers, ignore.case = TRUE))) {
        stop("Detected Ensembl IDs as row names. Heatmap requires gene symbols.")
      }

      expr_genes_upper <- toupper(gene_identifiers)
      pathway_genes_upper <- toupper(pathway_genes)

      matched_idx <- which(expr_genes_upper %in% pathway_genes_upper)

      if (length(matched_idx) < 2) {
        stop(sprintf("Gene matching failed: only %d/%d pathway genes matched.",
                     length(matched_idx), length(pathway_genes)))
      }

      message(sprintf("[Modal] Matched %d/%d pathway genes", length(matched_idx), length(pathway_genes)))

      plot_mat <- cpm_mat[matched_idx, , drop = FALSE]
      rownames(plot_mat) <- gene_identifiers[matched_idx]

      gene_vars <- apply(plot_mat, 1, var, na.rm = TRUE)
      valid_genes <- gene_vars > 1e-6 & !is.na(gene_vars)
      if (sum(valid_genes) < 2) {
        stop("Insufficient genes with variance after filtering")
      }
      plot_mat <- plot_mat[valid_genes, , drop = FALSE]

      display_numbers <- round(plot_mat)

      z_mat <- t(scale(t(plot_mat)))
      z_mat[is.na(z_mat)] <- 0
      z_mat[z_mat > 1] <- 1
      z_mat[z_mat < -1] <- -1

      gene_list <- pdata$data_list$gsea_res@geneList

      gene_metrics <- sapply(rownames(plot_mat), function(g) {
        idx <- match(toupper(g), toupper(names(gene_list)))
        if (is.na(idx)) return(0)
        return(gene_list[idx])
      })

      res_df <- as.data.frame(pdata$data_list$gsea_res@result)
      core_str <- res_df$core_enrichment[res_df$ID == pdata$pathway_id]
      core_genes <- character(0)
      if (length(core_str) > 0 && !is.na(core_str[1])) {
        core_genes <- unlist(strsplit(as.character(core_str[1]), "/"))
      }
      is_leading <- toupper(rownames(plot_mat)) %in% toupper(core_genes)

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

      col_fun <- circlize::colorRamp2(
        c(-1, 0, 1),
        c("#67a9cf", "#f7f7f7", "#ef8a62")
      )

      grp_col <- c("#E41A1C", "#377EB8")
      names(grp_col) <- c(left_grp, right_grp)

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

      leading_status <- ifelse(is_leading, "YES", "NO")
      leading_colors <- c("YES" = "#FF9800", "NO" = "transparent")

      right_ann <- ComplexHeatmap::rowAnnotation(
        LeadingEdge = leading_status,
        col = list(LeadingEdge = leading_colors),
        annotation_name_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        simple_anno_size = grid::unit(0.4, "cm")
      )

      layer_fun <- function(j, i, x, y, w, h, fill) {
        vals <- ComplexHeatmap::pindex(display_numbers, i, j)
        grid::grid.text(
          label = vals,
          x = x, y = y,
          gp = grid::gpar(fontsize = 13, col = "black", fontface = "bold")
        )
      }

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

      row_names_gp <- grid::gpar(
        fontsize = 15,
        fontface = ifelse(is_leading, "bold", "plain"),
        col = ifelse(is_leading, "#FF9800", "black")
      )

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

      enriched_group <- ifelse(pdata$nes > 0, left_grp, right_grp)
      title_text <- sprintf(
        "%s\nRow-Scaled Z-Score [-1, 1] | Enriched in: %s",
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
      pdata <- current_data()
      if (is.null(pdata)) return(400)
      n_genes <- length(pdata$pathway_genes)
      height_px <- max(400, n_genes * 25 + 150)
      return(height_px)
    })

    # =========================================
    # 基因统计表 - checkbox随modal_selected_genes响应式更新
    # =========================================
    output$modal_gene_table <- DT::renderDataTable({
      # 关键：依赖 modal_selected_genes() 使得checkbox状态能响应式更新
      pdata <- current_data()
      shiny::req(pdata)

      # 获取当前modal选择状态（不再用isolate）
      current_selection <- modal_selected_genes()

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

      # 生成checkbox列（根据current_selection显示勾选状态）
      gene_df$Select <- sapply(gene_df$Gene, function(g) {
        is_checked <- ifelse(g %in% current_selection, 'checked="checked"', '')
        sprintf(
          '<input type="checkbox" class="modal-gene-checkbox" data-gene="%s" %s onclick="Shiny.setInputValue(\'%s\', {id: \'%s\', checked: this.checked}, {priority: \'event\'});"/>',
          g, is_checked, ns("modal_gene_toggle"), g
        )
      })

      # 排序列
      gene_df <- gene_df[, c("Select", "Gene", "Rank_Metric", "Rank_in_List", "Is_Core")]
      colnames(gene_df) <- c("Select", "Gene", "Rank Metric", "Rank in List", "Leading Edge")

      DT::datatable(
        gene_df,
        rownames = FALSE,
        escape = FALSE,
        extensions = c('Scroller'),
        options = list(
          pageLength = -1,
          scrollY = FALSE,  # 禁用DT内置滚动，使用外层div滚动
          scroller = FALSE,
          dom = 'rtip',
          columnDefs = list(
            list(width = '50px', targets = 0, orderable = FALSE, className = 'dt-center'),
            list(width = '100px', targets = 1, className = 'dt-center'),
            list(width = '80px', targets = 2, className = 'dt-center'),
            list(width = '80px', targets = 3, className = 'dt-center'),
            list(width = '80px', targets = 4, className = 'dt-center')
          )
        )
      ) %>%
        DT::formatStyle(
          columns = "Leading Edge",
          backgroundColor = DT::styleEqual(c("YES", "-"), c("#FF9800", "transparent")),
          fontWeight = DT::styleEqual("YES", "bold"),
          color = DT::styleEqual("YES", "white")
        ) %>%
        DT::formatStyle(
          columns = "Select",
          backgroundColor = "#f0f8ff"
        )
    })
  })
}
