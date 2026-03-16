#' @title 启动 GSEAlens PRO 2.0 终极全息探索空间
#' @description 动态聚合子集合、手动挡刷新防止卡顿、重算 FDR。完美修复了 Limma 下标越界、表达矩阵基因大小写匹配问题，并将表格完全镜像 HTML 报表。
#' @param res_capsule 您环境中的 GseaResPro 对象
#' @importFrom magrittr %>%
#' @export
launch_gsea_app <- function(res_capsule) {

  if (!inherits(res_capsule, "GseaResPro")) {
    stop("❌ 严重错误: 传入的对象不是标准的 GseaResPro 计算胶囊！")
  }

  ui <- shiny::fluidPage(
    shiny::tags$head(
      shiny::tags$style(shiny::HTML("
        .modal-dialog { max-width: 90vw; margin: 3vh auto;}
        .modal-body { max-height: 80vh; overflow-y: auto; background-color: #f8f9fa; }
        .white-box { background-color: white; padding: 15px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-bottom: 15px;}
        .master-table-container { height: 45vh; overflow-y: hidden; }
        .multi-plot-container { height: 50vh; margin-top: 10px; }
      "))
    ),
    shiny::titlePanel("🧬 GSEAlens PRO 2.0: Holographic Knowledge Space"),

    shiny::sidebarLayout(
      shiny::sidebarPanel(
        width = 3,
        shiny::selectInput("selected_contrast", "⚖️ 1. 选择对比组 (Contrast):", choices = names(res_capsule$results)),
        shiny::hr(),
        shiny::h4("🎯 2. 数据切片与排序"),
        shiny::selectizeInput("selected_collections", "选择基因集亚组 (支持多选):", choices = NULL, multiple = TRUE, options = list(placeholder = '留空代表全选 (ALL)...')),
        shiny::selectInput("sort_by", "全局排序策略:", choices = c("按 NES (降序)" = "nes_desc", "按 NES (升序)" = "nes_asc", "按 P-value (升序)" = "pval_asc", "按 FDR (升序)" = "fdr_asc")),
        shiny::helpText("提示：缩小背景集合将自动过滤结果，并重算挽救 FDR (P.adjust)！"),
        shiny::br(),
        # 🌟 核心改进：取消自动更新，加入手动执行按钮
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
                                                             shiny::plotOutput("multi_gsea_plot", height = "400px")))
          ),
          shiny::tabPanel("🌋 全息四重联动",
                          shiny::br(),
                          shiny::fluidRow(
                            shiny::column(6, shiny::div(class="white-box", shiny::h4("1. 宏观: 通路火山图 👉"), plotly::plotlyOutput("volcano_pathway", height="400px"))),
                            shiny::column(6, shiny::div(class="white-box", shiny::h4("2. 微观: 基因 Rank 分布"), plotly::plotlyOutput("volcano_gene", height="400px")))
                          ),
                          shiny::fluidRow(
                            shiny::column(6, shiny::div(class="white-box", shiny::h4("3. 真实基因表达火山图 👉"), plotly::plotlyOutput("limma_volcano", height="400px"))),
                            shiny::column(6, shiny::div(class="white-box", shiny::h4("4. 基因表达散点图 (大小写容错)"), plotly::plotlyOutput("gene_expr_scatter", height="400px")))
                          )
          )
        )
      )
    )
  )

  server <- function(input, output, session) {

    # 根据组别动态更新子集选项，不触发自动执行
    shiny::observeEvent(input$selected_contrast, {
      meta_dict <- as.data.frame(res_capsule$geneset_info$meta_dict)
      if ("Combo_Name" %in% colnames(meta_dict)) {
        avail_cols <- unique(meta_dict$Combo_Name)
      } else if ("Collection" %in% colnames(meta_dict)) {
        avail_cols <- unique(meta_dict$Collection)
      } else {
        avail_cols <- "ALL"
      }
      shiny::updateSelectizeInput(session, "selected_collections", choices = avail_cols, selected = NULL)
    })

    # 🌟 核心改进：只有按下按钮才提取任务
    current_task_wrapper <- shiny::eventReactive(input$run_btn, {
      shiny::req(input$selected_contrast)
      task_info <- res_capsule$results[[input$selected_contrast]]
      if (is.null(task_info$data)) return(NULL)
      gsea_obj <- task_info$data

      c_mat <- res_capsule$contrast_matrix
      left_grp <- "Group1"; right_grp <- "Group2"
      if (!is.null(c_mat) && "Contrast_Name" %in% colnames(c_mat)) {
        idx <- which(paste0(c_mat$Num, "_vs_", c_mat$Den) == input$selected_contrast)
        if (length(idx) > 0) {
          left_grp <- c_mat$Num[idx[1]]; right_grp <- c_mat$Den[idx[1]]
        } else {
          idx_rev <- which(paste0(c_mat$Den, "_vs_", c_mat$Num) == input$selected_contrast)
          if (length(idx_rev) > 0) {
            left_grp <- c_mat$Den[idx_rev[1]]; right_grp <- c_mat$Num[idx_rev[1]]
          }
        }
      }
      list(gsea_res = gsea_obj, meta = list(left_group = left_grp, right_group = right_grp, expr_data = res_capsule$expr_data))
    }, ignoreNULL = FALSE)

    display_df <- shiny::reactive({
      task <- current_task_wrapper()
      shiny::req(task, task$gsea_res)
      df <- as.data.frame(task$gsea_res@result)
      meta_dict <- as.data.frame(res_capsule$geneset_info$meta_dict)
      df <- merge(df, meta_dict, by="ID", all.x=TRUE)

      df$Display_Collection <- if("Combo_Name" %in% colnames(df)) df$Combo_Name else if("Collection" %in% colnames(df)) df$Collection else "Unknown"
      df$Pathway_Link <- if("URL" %in% colnames(df)) {
        ifelse(is.na(df$URL) | df$URL == "", sprintf("<b>%s</b>", df$ID), sprintf('<a href="%s" target="_blank" style="color: #0056b3; text-decoration: none;">%s</a>', df$URL, df$ID))
      } else { sprintf("<b>%s</b>", df$ID) }
      df$Description <- if("long_description_for_html" %in% colnames(df)) df$long_description_for_html else if("Description" %in% colnames(df)) df$Description else df$ID

      # 读取按钮按下时的配置隔离状态
      sel_cols <- shiny::isolate(input$selected_collections)
      if (!is.null(sel_cols) && length(sel_cols) > 0 && !("ALL" %in% sel_cols)) {
        df <- df %>% dplyr::filter(Display_Collection %in% sel_cols)
      }

      df$p.adjust <- stats::p.adjust(df$pvalue, method = "BH")

      # 🌟 完美复刻 HTML 报表字段结构
      df <- df %>%
        dplyr::mutate(
          Enriched_In = factor(ifelse(NES > 0, task$meta$left_group, task$meta$right_group), levels = c(task$meta$left_group, task$meta$right_group)),
          NES_Round = round(NES, 3),
          Pval = signif(pvalue, 4),
          FDR = signif(p.adjust, 4),
          Detail_Page = sprintf('<button class="btn btn-sm btn-success" style="padding: 2px 10px;" onClick="Shiny.setInputValue(\'show_detail\', \'%s\', {priority: \'event\'})">🔍 Dashboard</button>', ID)
        ) %>%
        dplyr::arrange(
          dplyr::case_when(
            shiny::isolate(input$sort_by) == "nes_desc" ~ dplyr::desc(NES),
            shiny::isolate(input$sort_by) == "nes_asc" ~ NES,
            shiny::isolate(input$sort_by) == "pval_asc" ~ pvalue,
            shiny::isolate(input$sort_by) == "fdr_asc" ~ p.adjust
          )
        ) %>%
        dplyr::mutate(Rank = dplyr::row_number())

      df %>% dplyr::select(Rank, Detail_Page, Pathway = Pathway_Link, Collection = Display_Collection, Enriched_In, Size = setSize, NES = NES_Round, pvalue = Pval, p.adjust = FDR, Description, ID)
    })

    output$master_table <- DT::renderDataTable({
      df_show <- display_df()
      shiny::req(nrow(df_show) > 0)

      left_grp <- current_task_wrapper()$meta$left_group
      right_grp <- current_task_wrapper()$meta$right_group

      DT::datatable(
        df_show %>% dplyr::select(-ID),
        escape = FALSE, selection = "multiple", rownames = FALSE,
        colnames = c("Rank", "子网页/操作", "Pathway", "Collection", "Enriched In", "Size", "NES", "P-value", "FDR", "Description"),
        options = list(scrollY = "45vh", scroller = TRUE, paging = FALSE, dom = 'Bfrtip', autoWidth = TRUE)
      ) %>%
        DT::formatStyle('Enriched_In', backgroundColor = DT::styleEqual(c(left_grp, right_grp), c('#fee0d2', '#deebf7'))) %>%
        DT::formatStyle('NES', color = DT::styleInterval(0, c('blue', 'red')), fontWeight = 'bold')
    })

    output$multi_gsea_plot <- shiny::renderPlot({
      shiny::req(input$master_table_rows_selected)
      df <- display_df()
      selected_ids <- df$ID[input$master_table_rows_selected]
      colors <- trimws(unlist(strsplit(input$custom_colors, ",")))
      plot_directional_gsea(directional_gsea_obj = current_task_wrapper(), target_pathways = selected_ids, subPlot = as.numeric(input$plot_subtype), curveCol = colors, main_title = paste("联合分析:", length(selected_ids), "条通路"))
    })

    shiny::observeEvent(input$show_detail, {
      pw_id <- input$show_detail
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML(sprintf("<b style='color:#0056b3;'><i class='fas fa-microscope'></i> Pathway Dashboard:</b> %s", pw_id)),
        shiny::fluidRow(
          shiny::column(5, shiny::div(class = "white-box", shiny::h4("经典 GSEA 轮廓"), shiny::plotOutput("modal_gsea_plot", height = "500px"))),
          shiny::column(7, shiny::div(class = "white-box", shiny::h4("核心表达模式 (CPM 数值)"), shiny::div(style = "height: 500px; overflow-y: auto;", shiny::plotOutput("modal_heatmap", height = "auto"))))
        ), easyClose = TRUE, size = "l"
      ))

      output$modal_gsea_plot <- shiny::renderPlot({ plot_directional_gsea(current_task_wrapper(), target_pathways = pw_id, subPlot = 3) })

      output$modal_heatmap <- shiny::renderPlot({
        task <- current_task_wrapper()
        dge_list <- task$meta$expr_data
        sample_meta <- dge_list$samples
        target_samples <- rownames(sample_meta)[sample_meta$group %in% c(task$meta$left_group, task$meta$right_group)]

        expr_mat <- edgeR::cpm(dge_list, log = TRUE)[, target_samples, drop = FALSE]
        cpm_mat <- edgeR::cpm(dge_list, log = FALSE)[, target_samples, drop = FALSE]

        # 🌟 修复弹窗热图大小写容错匹配问题
        expr_genes <- rownames(expr_mat)
        all_genes <- task$gsea_res@geneSets[[pw_id]]
        plot_genes_idx <- which(toupper(expr_genes) %in% toupper(all_genes))
        if(length(plot_genes_idx) < 2) return()
        plot_genes <- expr_genes[plot_genes_idx]

        gene_metrics <- sapply(plot_genes, function(x) {
          idx <- match(toupper(x), toupper(names(task$gsea_res@geneList)))
          if(is.na(idx)) 0 else task$gsea_res@geneList[idx]
        })

        plot_mat <- expr_mat[plot_genes[order(gene_metrics, decreasing = TRUE)], , drop = FALSE]
        plot_mat <- plot_mat[apply(plot_mat, 1, var) > 1e-6, , drop = FALSE]
        if(nrow(plot_mat) < 2) return()

        z_mat <- t(scale(t(plot_mat)))
        z_mat[is.na(z_mat)] <- 0; z_mat[z_mat > 1] <- 1; z_mat[z_mat < -1] <- -1

        ann_col <- data.frame(Group = sample_meta[target_samples, "group"], row.names = target_samples)
        grp_col <- c("#E41A1C", "#377EB8"); names(grp_col) <- c(task$meta$left_group, task$meta$right_group)

        pheatmap::pheatmap(
          z_mat, scale = "none", cluster_cols = FALSE, cluster_rows = FALSE,
          gaps_col = sum(sample_meta[target_samples, "group"] == task$meta$left_group),
          color = grDevices::colorRampPalette(c("#67a9cf", "#f7f7f7", "#ef8a62"))(100),
          breaks = seq(-1, 1, length.out = 101), annotation_col = ann_col,
          annotation_colors = list(Group = grp_col), show_rownames = TRUE,
          fontsize_row = 10, display_numbers = round(cpm_mat[rownames(z_mat), ], 1),
          number_color = "black", fontsize_number = 8
        )
      }, height = function() { max(500, length(which(toupper(rownames(res_capsule$expr_data)) %in% toupper(current_task_wrapper()$gsea_res@geneSets[[pw_id]]))) * 20 + 100) })
    })

    output$volcano_pathway <- plotly::renderPlotly({
      df <- display_df()
      plotly::plot_ly(data = df, x = ~NES, y = ~-log10(p.adjust), type = "scatter", mode = "markers", text = ~ID, hoverinfo = "text", key = ~ID, marker = list(color = ifelse(df$NES > 0, "#E41A1C", "#377EB8"), size = 8, opacity = 0.6), source = "pathway_volcano") %>%
        plotly::layout(xaxis = list(title = "NES"), yaxis = list(title = "-log10 (FDR)"), dragmode = "select")
    })

    output$volcano_gene <- plotly::renderPlotly({
      task <- current_task_wrapper()
      gList <- task$gsea_res@geneList
      plot_df <- data.frame(Rank = 1:length(gList), Metric = as.numeric(gList), Gene = names(gList), Color = "Background", Size = 4, Alpha = 0.2, stringsAsFactors = FALSE)

      click_data <- plotly::event_data("plotly_click", source = "pathway_volcano")
      if (!is.null(click_data)) {
        hl_idx <- toupper(plot_df$Gene) %in% toupper(task$gsea_res@geneSets[[click_data$key]])
        plot_df$Color[hl_idx] <- "Target"; plot_df$Size[hl_idx] <- 10; plot_df$Alpha[hl_idx] <- 1
        plot_df <- plot_df[order(plot_df$Color == "Target"), ]
      }
      plotly::plot_ly(data = plot_df, x = ~Rank, y = ~Metric, type = "scattergl", mode = "markers", color = ~Color, colors = c("Background" = "#CFD8DC", "Target" = "#FF9800"), size = ~Size, sizes = c(4, 10), marker = list(opacity = ~Alpha), text = ~Gene, hoverinfo = "text") %>%
        plotly::layout(title = list(text = "全局基因 Rank", font=list(size=12)), showlegend = FALSE)
    })

    # 🌟 修复点1: Limma火山图越界防护，还原真实的对比 coef 和 反向翻转逻辑
    output$limma_volcano <- plotly::renderPlotly({
      shiny::req(res_capsule$limma_fit)
      fit <- res_capsule$limma_fit
      task_name <- shiny::isolate(input$selected_contrast) # 读取按钮隔离时的状态

      c_mat <- res_capsule$contrast_matrix
      is_reverse <- FALSE

      if (!is.null(c_mat) && "Contrast_Name" %in% colnames(c_mat)) {
        idx <- which(paste0(c_mat$Num, "_vs_", c_mat$Den) == task_name)
        if (length(idx) > 0) {
          coef_name <- c_mat$Contrast_Name[idx[1]]
        } else {
          idx_rev <- which(paste0(c_mat$Den, "_vs_", c_mat$Num) == task_name)
          if (length(idx_rev) > 0) {
            coef_name <- c_mat$Contrast_Name[idx_rev[1]]
            is_reverse <- TRUE # 逆向计算时需要翻转 FC
          } else { coef_name <- colnames(fit)[1] }
        }
      } else { coef_name <- colnames(fit)[1] }

      # 加一层保险，应对极为罕见的无法匹配
      tt <- tryCatch({ limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
      }, error = function(e) { limma::topTable(fit, coef = 1, number = Inf, sort.by = "none") })

      if (is_reverse && "logFC" %in% colnames(tt)) { tt$logFC <- -tt$logFC }

      tt$Gene <- rownames(tt)
      tt$Status <- "NotSig"
      pval_col <- if("adj.P.Val" %in% colnames(tt)) "adj.P.Val" else "P.Value"
      tt$Status[tt$logFC > 0.5 & tt[[pval_col]] < 0.05] <- "Up"
      tt$Status[tt$logFC < -0.5 & tt[[pval_col]] < 0.05] <- "Down"
      tt$Size <- 4; tt$Alpha <- 0.4

      click_data <- plotly::event_data("plotly_click", source = "pathway_volcano")
      if (!is.null(click_data)) {
        hl_idx <- toupper(tt$Gene) %in% toupper(current_task_wrapper()$gsea_res@geneSets[[click_data$key]])
        tt$Status[hl_idx] <- "PathwayTarget"; tt$Size[hl_idx] <- 12; tt$Alpha[hl_idx] <- 1
        tt <- tt[order(tt$Status == "PathwayTarget"), ]
      }
      pal <- c("NotSig" = "#CFD8DC", "Up" = "#F8BBD0", "Down" = "#BBDEFB", "PathwayTarget" = "#FF1744")
      plotly::plot_ly(data = tt, x = ~logFC, y = ~-log10(P.Value), type = "scattergl", mode = "markers", color = ~Status, colors = pal, size = ~Size, sizes = c(4, 12), marker = list(opacity = ~Alpha, line = list(color="white", width=0.5)), text = ~Gene, hoverinfo = "text", key = ~Gene, source = "gene_volcano") %>%
        plotly::layout(title = list(text = "真实基因火山图 (Limma)", font=list(size=12)), showlegend = FALSE)
    })

    # 🌟 修复点2: 表达散点图忽略大小写强制挂载
    output$gene_expr_scatter <- plotly::renderPlotly({
      shiny::req(res_capsule$expr_data)
      click_gene <- plotly::event_data("plotly_click", source = "gene_volcano")
      if (is.null(click_gene)) return(plotly::plot_ly() %>% plotly::layout(title = list(text="👈 请在左侧真实火山图中点击特定基因查看", font=list(size=14))))

      target_gene <- click_gene$key
      all_cpm <- edgeR::cpm(res_capsule$expr_data, log = TRUE)

      # 极简核心修复：全部转大写寻找坐标！
      match_idx <- match(toupper(target_gene), toupper(rownames(all_cpm)))
      if (is.na(match_idx)) return(plotly::plot_ly() %>% plotly::layout(title = list(text=sprintf("基因 '%s' 不在表达矩阵中", target_gene), font=list(size=14))))

      actual_gene <- rownames(all_cpm)[match_idx]
      gene_expr <- all_cpm[actual_gene, ]

      plot_df <- data.frame(Sample = names(gene_expr), Expression = as.numeric(gene_expr), Group = res_capsule$expr_data$samples$group)
      task <- current_task_wrapper()
      plot_df <- plot_df[plot_df$Group %in% c(task$meta$left_group, task$meta$right_group), ]

      p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Group, y = Expression, fill = Group)) +
        ggplot2::geom_boxplot(alpha = 0.5, outlier.shape = NA) +
        ggplot2::geom_jitter(width = 0.2, size = 2, ggplot2::aes(text = Sample)) +
        ggplot2::theme_bw() +
        ggplot2::labs(title = sprintf("基因表达 (logCPM): %s", actual_gene), y = "log2(CPM)") +
        ggplot2::theme(legend.position = "none")
      plotly::ggplotly(p, tooltip = c("text", "y"))
    })
  }

  message("🚀 GSEAlens PRO 2.0 终极版启动！(手动挡控制已激活)")
  shiny::shinyApp(ui, server)
}
