#' @title 启动 GSEAlens PRO 2.0 终极全息探索空间
#' @description 彻底终结 [object Object] 报错、修复 ALL 选项、加入 Leading Edge 追踪及宽屏 Dashboard。
#' @param res_capsule 您环境中的 GseaResPro 对象
#' @importFrom magrittr %>%
#' @export
launch_gsea_app <- function(res_capsule) {

  if (!inherits(res_capsule, "GseaResPro")) stop("❌ 严重错误: 传入的对象不是标准的 GseaResPro 计算胶囊！")

  ui <- shiny::fluidPage(
    shiny::tags$head(
      shiny::tags$style(shiny::HTML("
        /* 超宽屏 Dashboard 模态窗配置，占据 80% 宽度，完美 16:9/4:3 比例 */
        .modal-dialog { max-width: 80vw !important; width: 80vw !important; margin: 3vh auto; }
        .modal-content { border-radius: 12px; box-shadow: 0 5px 15px rgba(0,0,0,0.3); }
        .modal-body { min-height: 75vh; overflow-y: auto; background-color: #f8f9fa; }

        .white-box { background-color: white; padding: 15px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-bottom: 15px;}

        /* 强力保障 DataTable 的横向滚动条 */
        .master-table-container { width: 100%; overflow-x: auto !important; margin-bottom: 20px;}
        .multi-plot-container { min-height: 50vh; margin-top: 10px; }
      "))
    ),
    shiny::titlePanel("🧬 GSEAlens PRO 2.0: Holographic Knowledge Space"),

    shiny::sidebarLayout(
      shiny::sidebarPanel(
        width = 3,
        shiny::selectInput("selected_contrast", "⚖️ 1. 选择对比组 (Contrast):", choices = names(res_capsule$results)),
        shiny::hr(),
        shiny::h4("🎯 2. 数据切片与排序"),
        shiny::selectizeInput("selected_collections", "选择基因集亚组 (支持多选):", choices = NULL, multiple = TRUE),
        shiny::selectInput("sort_by", "全局排序策略:", choices = c("按 NES (降序)" = "nes_desc", "按 NES (升序)" = "nes_asc", "按 P-value (升序)" = "pval_asc", "按 FDR (升序)" = "fdr_asc")),
        shiny::helpText("提示：缩小背景集合将自动过滤结果，并重算挽救 FDR (P.adjust)！"),
        shiny::br(),
        shiny::actionButton("run_btn", "🚀 确认配置 / 更新工作台", class = "btn-success", style = "width: 100%; font-weight: bold; font-size: 16px; margin-bottom: 15px;"),
        shiny::hr(),
        shiny::h4("🎨 3. 联合绘图控制台"),
        shiny::selectInput("plot_subtype", "GSEAvis 样式 (Subtype):", choices = c("1: 仅经典富集" = 1, "2: 富集+热图带" = 2, "3: 完整带Rank" = 3), selected = 3),
        shiny::textInput("custom_colors", "多通路自定义颜色:", value = "#E41A1C, #377EB8, #4DAF4A, #984EA3")
      ),
      shiny::mainPanel(
        width = 9,
        shiny::tabsetPanel(
          shiny::tabPanel("📊 主工作台 (Master Table)",
                          shiny::br(),
                          shiny::div(class = "white-box master-table-container", DT::dataTableOutput("master_table")),
                          shiny::conditionalPanel(condition = "input.master_table_rows_selected.length > 0",
                                                  shiny::div(class = "white-box multi-plot-container",
                                                             shiny::h4("🖼️ 联合绘图区 - 勾选上方表格行即可渲染"),
                                                             shiny::plotOutput("multi_gsea_plot", height = "500px")))
          ),
          shiny::tabPanel("🌋 全息四重联动",
                          shiny::br(),
                          shiny::fluidRow(
                            shiny::column(6, shiny::div(class="white-box", shiny::h4("1. 宏观: 通路火山图 👉"), plotly::plotlyOutput("volcano_pathway", height="400px"))),
                            shiny::column(6, shiny::div(class="white-box", shiny::h4("2. 微观: 基因 Rank 分布"), plotly::plotlyOutput("volcano_gene", height="400px")))
                          ),
                          shiny::fluidRow(
                            shiny::column(6, shiny::div(class="white-box", shiny::h4("3. 真实基因表达火山图 👉"), plotly::plotlyOutput("limma_volcano", height="400px"))),
                            shiny::column(6, shiny::div(class="white-box", shiny::h4("4. 全景表达散点图 (全量组别)"), plotly::plotlyOutput("gene_expr_scatter", height="400px")))
                          )
          )
        )
      )
    )
  )

  server <- function(input, output, session) {

    # 🌟 核心修复1：直接提取 GSEA 数据自带的 Collection 进行选项构建，修复 ALL 失败问题
    shiny::observe({
      shiny::req(input$selected_contrast)
      task_info <- res_capsule$results[[input$selected_contrast]]
      if (is.null(task_info$data)) return()

      df <- as.data.frame(task_info$data)
      avail_cols <- c()
      if ("Collection" %in% colnames(df)) avail_cols <- c(avail_cols, unique(df$Collection))
      if ("Combo_Name" %in% colnames(df)) avail_cols <- c(avail_cols, unique(df$Combo_Name))
      if ("Subcollection" %in% colnames(df)) avail_cols <- c(avail_cols, unique(df$Subcollection))

      avail_cols <- unique(avail_cols[!is.na(avail_cols) & avail_cols != ""])

      # 智能默认寻找 Hallmark (H 集合)
      hm_idx <- grep("^(H|H:.*|Hallmark)$", avail_cols, ignore.case = TRUE)
      default_sel <- if (length(hm_idx) > 0) avail_cols[hm_idx[1]] else "ALL"

      choices <- c("ALL", sort(avail_cols))
      shiny::updateSelectizeInput(session, "selected_collections", choices = choices, selected = default_sel)
    })

    current_task_wrapper <- shiny::eventReactive(input$run_btn, {
      shiny::req(input$selected_contrast)
      task_info <- res_capsule$results[[input$selected_contrast]]
      if (is.null(task_info$data)) return(NULL)

      parts <- strsplit(input$selected_contrast, "_vs_")[[1]]
      left_grp <- if(length(parts) >= 1) parts[1] else input$selected_contrast
      right_grp <- if(length(parts) >= 2) parts[2] else "Background"

      list(gsea_res = task_info$data, meta = list(left_group = left_grp, right_group = right_grp, expr_data = res_capsule$expr_data))
    }, ignoreNULL = FALSE)

    display_df <- shiny::reactive({
      task <- current_task_wrapper()
      shiny::req(task, task$gsea_res)
      df <- as.data.frame(task$gsea_res) # 直接 as.data.frame，防止 merge 导致列重复

      df$Display_Collection <- if("Combo_Name" %in% colnames(df)) df$Combo_Name else if("Collection" %in% colnames(df)) df$Collection else "Unknown"

      sel_cols <- shiny::isolate(input$selected_collections)
      if (!is.null(sel_cols) && length(sel_cols) > 0 && !("ALL" %in% sel_cols)) {
        # 兼容匹配大类、小类、混合类
        df <- df %>% dplyr::filter(
          Display_Collection %in% sel_cols |
            (if("Collection" %in% colnames(df)) Collection %in% sel_cols else FALSE) |
            (if("Subcollection" %in% colnames(df)) Subcollection %in% sel_cols else FALSE)
        )
      }
      if(nrow(df) == 0) return(NULL)

      df$p.adjust <- stats::p.adjust(df$pvalue, method = "BH")
      df$Safe_ID <- gsub("'", "\\\\'", df$ID)
      df$Pathway_Link <- if("URL" %in% colnames(df)) {
        ifelse(is.na(df$URL) | df$URL == "", sprintf("<b>%s</b>", df$ID), sprintf('<a href="%s" target="_blank" style="color: #0056b3; text-decoration: none;">%s</a>', df$URL, df$ID))
      } else { sprintf("<b>%s</b>", df$ID) }
      df$Description <- if("Description" %in% colnames(df)) df$Description else df$ID

      df %>%
        dplyr::mutate(
          Enriched_In = factor(ifelse(NES > 0, task$meta$left_group, task$meta$right_group), levels = c(task$meta$left_group, task$meta$right_group)),
          Detail_Page = sprintf('<button class="btn btn-sm btn-success" style="padding: 2px 10px;" onClick="Shiny.setInputValue(\'show_detail\', \'%s\', {priority: \'event\'})">🔍 Dashboard</button>', Safe_ID)
        ) %>%
        dplyr::arrange(
          dplyr::case_when(
            shiny::isolate(input$sort_by) == "nes_desc" ~ dplyr::desc(NES),
            shiny::isolate(input$sort_by) == "nes_asc" ~ NES,
            shiny::isolate(input$sort_by) == "pval_asc" ~ pvalue,
            shiny::isolate(input$sort_by) == "fdr_asc" ~ p.adjust
          )
        ) %>%
        dplyr::mutate(Rank = dplyr::row_number()) %>%
        dplyr::select(Rank, Detail_Page, Pathway = Pathway_Link, Collection = Display_Collection, Enriched_In, Size = setSize, NES = NES, pvalue = pvalue, p.adjust = p.adjust, Description, ID)
    })

    # 横向滚轴表格
    output$master_table <- DT::renderDataTable({
      df_show <- display_df()
      if(is.null(df_show) || nrow(df_show) == 0) return(DT::datatable(data.frame("系统提示" = "🚨 该子集过滤后无数据，请尝试选择其他子集或回到 ALL"), options = list(dom = 't')))

      DT::datatable(
        df_show %>% dplyr::select(-ID) %>% dplyr::mutate(NES=round(NES,3), pvalue=signif(pvalue,3), p.adjust=signif(p.adjust,3)),
        escape = FALSE, selection = "multiple", rownames = FALSE,
        options = list(
          scrollX = TRUE, scrollY = "40vh", scroller = TRUE, paging = FALSE, dom = 'Bfrtip', autoWidth = TRUE,
          columnDefs = list(list(width = '400px', targets = 9))
        )
      ) %>%
        DT::formatStyle('Enriched_In', backgroundColor = DT::styleEqual(c(current_task_wrapper()$meta$left_group, current_task_wrapper()$meta$right_group), c('#fee0d2', '#deebf7'))) %>%
        DT::formatStyle('NES', color = DT::styleInterval(0, c('blue', 'red')), fontWeight = 'bold')
    })

    # 🌟 核心修复2：终结 [object Object] 报错的绝对防御（invisible(NULL)）
    output$multi_gsea_plot <- shiny::renderPlot({
      shiny::req(input$master_table_rows_selected)
      df <- display_df()
      shiny::req(nrow(df) > 0)

      valid_idx <- intersect(input$master_table_rows_selected, seq_len(nrow(df)))
      shiny::req(length(valid_idx) > 0)

      selected_ids <- df$ID[valid_idx]
      colors <- trimws(unlist(strsplit(input$custom_colors, ",")))

      tryCatch({
        p <- plot_directional_gsea(directional_gsea_obj = current_task_wrapper(), target_pathways = selected_ids, subPlot = as.numeric(input$plot_subtype), curveCol = colors, main_title = paste("联合展示:", length(selected_ids), "条通路"))
        print(p)
      }, error = function(e) {
        par(mar=c(1,1,1,1))
        plot(1, type="n", axes=FALSE, xlab="", ylab="")
        text(1, 1, labels=paste("❌ 绘图引擎冲突:\n", e$message), col="red", cex=1.2)
      })
      invisible(NULL) # 必须返回 NULL，阻断 Shiny JS 层尝试解析图层对象导致的 [object Object] 报错！
    })

    # 🌟 超宽屏 Dashboard & Leading Edge 热图追踪
    shiny::observeEvent(input$show_detail, {
      pw_id <- input$show_detail
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML(sprintf("<b style='color:#0056b3; font-size:22px;'><i class='fas fa-microscope'></i> Pathway Dashboard: %s</b>", pw_id)),
        shiny::fluidRow(
          shiny::column(5, shiny::div(class = "white-box", shiny::h4("经典 GSEA 轮廓 (防崩溃版)"), shiny::plotOutput("modal_gsea_plot", height = "600px"))),
          shiny::column(7, shiny::div(class = "white-box", shiny::h4("核心表达模式 & Leading Edge 标记追踪"), shiny::div(style = "height: 600px; overflow-y: auto;", shiny::plotOutput("modal_heatmap", height = "auto"))))
        ), easyClose = TRUE
      ))

      output$modal_gsea_plot <- shiny::renderPlot({
        tryCatch({
          print(plot_directional_gsea(current_task_wrapper(), target_pathways = pw_id, subPlot = 3))
        }, error = function(e) {
          par(mar=c(1,1,1,1)); plot(1, type="n", axes=FALSE, xlab="", ylab=""); text(1, 1, labels=paste("❌ 制图失败:", e$message), col="red")
        })
        invisible(NULL) # 阻断报错
      })

      output$modal_heatmap <- shiny::renderPlot({
        task <- current_task_wrapper()
        dge_list <- task$meta$expr_data; if(is.null(dge_list)) return()
        sample_meta <- dge_list$samples
        target_samples <- rownames(sample_meta)[sample_meta$group %in% c(task$meta$left_group, task$meta$right_group)]

        expr_mat <- edgeR::cpm(dge_list, log = TRUE)[, target_samples, drop = FALSE]
        cpm_mat <- edgeR::cpm(dge_list, log = FALSE)[, target_samples, drop = FALSE]

        all_genes <- task$gsea_res@geneSets[[pw_id]]
        expr_genes <- rownames(expr_mat)

        # 必须转大写进行匹配，因为表达矩阵基因名可能是 Cyp2e1，而字典里是 CYP2E1
        plot_genes_idx <- which(toupper(expr_genes) %in% toupper(all_genes))
        if(length(plot_genes_idx) < 2) return()

        # 精确提取 Leading Edge 基因
        res_df <- as.data.frame(task$gsea_res)
        core_str <- res_df$core_enrichment[res_df$ID == pw_id]
        core_genes <- if(length(core_str) > 0 && !is.na(core_str)) unlist(strsplit(core_str[1], "/")) else character(0)

        plot_genes <- expr_genes[plot_genes_idx]
        gene_metrics <- sapply(plot_genes, function(x) {
          idx <- match(toupper(x), toupper(names(task$gsea_res@geneList)))
          if(is.na(idx)) 0 else task$gsea_res@geneList[idx]
        })
        is_le <- toupper(plot_genes) %in% toupper(core_genes)

        # Leading Edge 基因强势置顶，其他随后，按重要性排序
        sort_order <- order(is_le, gene_metrics, decreasing = TRUE)
        plot_genes <- plot_genes[sort_order]

        plot_mat <- expr_mat[plot_genes, , drop = FALSE]
        plot_mat <- plot_mat[apply(plot_mat, 1, stats::var) > 1e-6, , drop = FALSE]
        if(nrow(plot_mat) < 2) return()

        z_mat <- t(scale(t(plot_mat)))
        z_mat[is.na(z_mat)] <- 0; z_mat[z_mat > 1.5] <- 1.5; z_mat[z_mat < -1.5] <- -1.5

        # 增加双重注释: 分组(Col) + LeadingEdge(Row)
        ann_col <- data.frame(Group = sample_meta[target_samples, "group"], row.names = target_samples)
        ann_row <- data.frame(LeadingEdge = ifelse(toupper(rownames(z_mat)) %in% toupper(core_genes), "YES", "NO"), row.names = rownames(z_mat))

        grp_col <- c("#E41A1C", "#377EB8"); names(grp_col) <- c(task$meta$left_group, task$meta$right_group)
        le_col <- c("YES" = "#FF9800", "NO" = "#FAFAFA") # YES 用醒目的橙色

        pheatmap::pheatmap(
          z_mat, scale = "none", cluster_cols = FALSE, cluster_rows = FALSE,
          gaps_col = sum(sample_meta[target_samples, "group"] == task$meta$left_group),
          color = grDevices::colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
          breaks = seq(-1.5, 1.5, length.out = 101),
          annotation_col = ann_col, annotation_row = ann_row,
          annotation_colors = list(Group = grp_col, LeadingEdge = le_col),
          show_rownames = TRUE, fontsize_row = 10,
          display_numbers = round(cpm_mat[rownames(z_mat), ], 1),
          number_color = "black", fontsize_number = 8
        )
      }, height = function() { max(600, length(which(toupper(rownames(res_capsule$expr_data)) %in% toupper(current_task_wrapper()$gsea_res@geneSets[[pw_id]]))) * 20 + 100) })
    })

    # 联动绘图 - 宏观、微观、真实火山
    output$volcano_pathway <- plotly::renderPlotly({
      df <- display_df(); shiny::req(df)
      plotly::plot_ly(data = df, x = ~NES, y = ~-log10(p.adjust), type = "scatter", mode = "markers", text = ~ID, hoverinfo = "text", key = ~ID, marker = list(color = ifelse(df$NES > 0, "#E41A1C", "#377EB8"), size = 8, opacity = 0.6), source = "pathway_volcano") %>% plotly::layout(xaxis = list(title = "NES"), yaxis = list(title = "-log10 (FDR)"), dragmode = "select")
    })

    output$volcano_gene <- plotly::renderPlotly({
      task <- current_task_wrapper(); gList <- task$gsea_res@geneList
      plot_df <- data.frame(Rank = 1:length(gList), Metric = as.numeric(gList), Gene = names(gList), Color = "Background", Size = 4, Alpha = 0.2, stringsAsFactors = FALSE)
      click_data <- plotly::event_data("plotly_click", source = "pathway_volcano")
      if (!is.null(click_data)) {
        hl_idx <- toupper(plot_df$Gene) %in% toupper(task$gsea_res@geneSets[[click_data$key]])
        plot_df$Color[hl_idx] <- "Target"; plot_df$Size[hl_idx] <- 10; plot_df$Alpha[hl_idx] <- 1
        plot_df <- plot_df[order(plot_df$Color == "Target"), ]
      }
      plotly::plot_ly(data = plot_df, x = ~Rank, y = ~Metric, type = "scattergl", mode = "markers", color = ~Color, colors = c("Background" = "#CFD8DC", "Target" = "#FF9800"), size = ~Size, sizes = c(4, 10), marker = list(opacity = ~Alpha), text = ~Gene, hoverinfo = "text", key = ~Gene, source = "volcano_gene") %>% plotly::layout(title = list(text = "全局基因 Rank", font=list(size=12)), showlegend = FALSE)
    })

    output$limma_volcano <- plotly::renderPlotly({
      shiny::req(res_capsule$limma_fit); fit <- res_capsule$limma_fit; task_name <- shiny::isolate(input$selected_contrast)
      avail_coefs <- if(!is.null(fit$contrasts)) colnames(fit$contrasts) else colnames(fit$coefficients)
      coef_to_use <- avail_coefs[1]; is_reverse <- FALSE

      c_mat <- res_capsule$contrast_matrix
      if (!is.null(c_mat) && "Contrast_Name" %in% colnames(c_mat)) {
        idx <- which(paste0(c_mat$Num, "_vs_", c_mat$Den) == task_name)
        if (length(idx) > 0) { coef_to_use <- c_mat$Contrast_Name[idx[1]] } else {
          idx_rev <- which(paste0(c_mat$Den, "_vs_", c_mat$Num) == task_name)
          if (length(idx_rev) > 0) { coef_to_use <- c_mat$Contrast_Name[idx_rev[1]]; is_reverse <- TRUE }
        }
      } else {
        parts <- strsplit(task_name, "_vs_")[[1]]; if (length(parts) == 2) { c1 <- paste0(parts[1], " - ", parts[2]); c3 <- paste0(parts[2], " - ", parts[1]); if (c1 %in% avail_coefs) coef_to_use <- c1 else if (c3 %in% avail_coefs) { coef_to_use <- c3; is_reverse <- TRUE } }
      }

      tt <- tryCatch({ as.data.frame(limma::topTable(fit, coef = coef_to_use, number = Inf, sort.by = "none")) }, error = function(e) { as.data.frame(limma::topTable(fit, coef = 1, number = Inf, sort.by = "none")) })
      if (is_reverse && "logFC" %in% colnames(tt)) tt$logFC <- -tt$logFC
      if (!"SYMBOL" %in% colnames(tt)) tt$SYMBOL <- rownames(tt)

      tt$Gene <- tt$SYMBOL; pval_col <- if("adj.P.Val" %in% colnames(tt)) "adj.P.Val" else "P.Value"
      tt$Status <- "NotSig"
      if("logFC" %in% colnames(tt) && pval_col %in% colnames(tt)) { tt$Status[tt$logFC > 0.5 & tt[[pval_col]] < 0.05] <- "Up"; tt$Status[tt$logFC < -0.5 & tt[[pval_col]] < 0.05] <- "Down" }
      tt$Size <- 4; tt$Alpha <- 0.4

      click_data <- plotly::event_data("plotly_click", source = "pathway_volcano")
      if (!is.null(click_data)) {
        hl_idx <- toupper(tt$Gene) %in% toupper(current_task_wrapper()$gsea_res@geneSets[[click_data$key]])
        tt$Status[hl_idx] <- "PathwayTarget"; tt$Size[hl_idx] <- 12; tt$Alpha[hl_idx] <- 1
        tt <- tt[order(tt$Status == "PathwayTarget"), ]
      }
      pal <- c("NotSig"="#CFD8DC", "Up"="#F8BBD0", "Down"="#BBDEFB", "PathwayTarget"="#FF1744")
      plotly::plot_ly(data = tt, x = ~logFC, y = ~-log10(P.Value), type = "scattergl", mode = "markers", color = ~Status, colors = pal, size = ~Size, sizes = c(4, 12), marker = list(opacity = ~Alpha, line = list(color="white", width=0.5)), text = ~Gene, hoverinfo = "text", key = ~Gene, source = "gene_volcano") %>% plotly::layout(title = list(text = sprintf("Limma Volcano (coef: %s)", coef_to_use), font=list(size=12)), showlegend = FALSE)
    })

    # 🌟 修复: 解除过滤，全量散点图回归展示所有的 Group 分组！
    output$gene_expr_scatter <- plotly::renderPlotly({
      shiny::req(res_capsule$expr_data)
      clk_gene <- plotly::event_data("plotly_click", source = "volcano_gene")
      clk_limma <- plotly::event_data("plotly_click", source = "gene_volcano")
      target_gene <- if(!is.null(clk_limma)) clk_limma$key else if(!is.null(clk_gene)) clk_gene$key else NULL
      if (is.null(target_gene)) return(plotly::plot_ly() %>% plotly::layout(title = list(text="👈 请在左上方或左侧图中点击一个基因", font=list(size=14))))

      all_cpm <- edgeR::cpm(res_capsule$expr_data, log = TRUE)
      match_idx <- match(toupper(target_gene), toupper(rownames(all_cpm)))
      if (is.na(match_idx)) return(plotly::plot_ly() %>% plotly::layout(title = list(text=sprintf("基因 '%s' 不在表达矩阵中", target_gene), font=list(size=14))))

      actual_gene <- rownames(all_cpm)[match_idx]
      gene_expr <- all_cpm[actual_gene, ]

      samples <- res_capsule$expr_data$samples
      plot_df <- data.frame(Sample = names(gene_expr), Expression = as.numeric(gene_expr), Group = if("group" %in% colnames(samples)) as.character(samples$group) else "Unknown", stringsAsFactors = FALSE)

      p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Group, y = Expression, fill = Group)) + ggplot2::geom_boxplot(alpha = 0.5, outlier.shape = NA) + ggplot2::geom_jitter(width = 0.2, size = 2, ggplot2::aes(text = Sample)) + ggplot2::theme_bw() + ggplot2::labs(title = sprintf("全量组别表达 (logCPM): %s", actual_gene), y = "log2(CPM)") + ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      plotly::ggplotly(p, tooltip = c("text", "y"))
    })
  }

  message("🚀 GSEAlens PRO 2.0 数据解析级重构完毕！(防崩溃已开启 / 全景图 / 标记追踪)")
  shiny::shinyApp(ui, server)
}
