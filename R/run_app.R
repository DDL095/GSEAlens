#' @title 启动 GSEAlens PRO 交互式探索空间
#' @description 接收一个 GseaResPro 计算胶囊对象，启动本地 Shiny 交互式数据分析网页。
#' @param res_capsule 必须是使用 \code{setup_gsea_env_pro} 等流程生成的 \code{GseaResPro} 对象。
#' @importFrom shiny fluidPage tags HTML titlePanel sidebarLayout sidebarPanel mainPanel tabsetPanel tabPanel br div hr h3 h4 selectInput textInput helpText conditionalPanel plotOutput observeEvent showModal modalDialog fluidRow column modalButton reactive req shinyApp
#' @importFrom DT dataTableOutput renderDataTable datatable formatStyle styleInterval formatRound styleEqual
#' @importFrom plotly plotlyOutput renderPlotly plot_ly layout event_data
#' @importFrom dplyr select mutate everything
#' @importFrom pheatmap pheatmap
#' @importFrom edgeR cpm
#' @export
launch_gsea_app <- function(res_capsule) {

  if (!inherits(res_capsule, "GseaResPro")) {
    stop("❌ 严重错误: 传入的对象不是标准的 GseaResPro 计算胶囊！请检查您的对象 class。")
  }

  # ==========================================
  # 前端 UI 界面架构
  # ==========================================
  ui <- shiny::fluidPage(
    shiny::tags$head(
      shiny::tags$style(shiny::HTML("
        .modal-dialog { max-width: 95vw; margin: 2vh auto;}
        .modal-body { max-height: 85vh; overflow-y: auto; background-color: #f8f9fa; }
        .scrollable-heatmap { max-height: 70vh; overflow-y: auto; overflow-x: hidden; border: 1px solid #ccc; background-color: white;}
        .white-box { background-color: white; padding: 15px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-bottom: 15px;}
      "))
    ),

    shiny::titlePanel("🧬 GSEAlens PRO: Interactive Knowledge Space"),

    shiny::sidebarLayout(
      shiny::sidebarPanel(
        width = 3,
        shiny::selectInput("selected_contrast", "⚖️ 1. 选择对比组 (Contrast):", choices = names(res_capsule$results)),
        shiny::hr(),
        shiny::h4("🎨 动线 A：联合绘图控制台"),
        shiny::selectInput("plot_subtype", "GSEAvis 样式 (Subtype):",
                           choices = c("1: 仅经典富集" = 1, "2: 富集+热图带" = 2, "3: 完整带Rank" = 3), selected = 3),
        shiny::textInput("custom_colors", "多通路自定义颜色 (逗号分隔):", value = "#E41A1C, #377EB8, #4DAF4A, #984EA3"),
        shiny::helpText("操作指南：在右侧表格中左键勾选多个通路，下方将自动渲染出完美的组合图。")
      ),

      shiny::mainPanel(
        width = 9,
        shiny::tabsetPanel(
          shiny::tabPanel("📊 主工作台 (Master Table)",
                          shiny::br(),
                          shiny::div(class = "white-box",
                                     DT::dataTableOutput("master_table")
                          ),
                          shiny::conditionalPanel(
                            condition = "input.master_table_rows_selected.length > 0",
                            shiny::div(class = "white-box",
                                       shiny::h3("🖼️ 联合绘图区 (Multi-Pathway Canvas)"),
                                       shiny::plotOutput("multi_gsea_plot", height = "550px")
                            )
                          )
          ),
          shiny::tabPanel("🌋 全息双重联动 (Volcano & Rank)",
                          shiny::br(),
                          shiny::fluidRow(
                            shiny::column(6,
                                          shiny::div(class = "white-box",
                                                     shiny::h4("宏观：通路级火山图 (点击任意点 👉)"),
                                                     plotly::plotlyOutput("volcano_pathway", height = "600px")
                                          )
                            ),
                            shiny::column(6,
                                          shiny::div(class = "white-box",
                                                     shiny::h4("微观：全局基因 Rank 分布 (WebGL 加速)"),
                                                     plotly::plotlyOutput("volcano_gene", height = "600px")
                                          )
                            )
                          )
          )
        )
      )
    )
  )

  # ==========================================
  # 后端 Server 逻辑核心
  # ==========================================
  server <- function(input, output, session) {

    # 提取核心计算对象 (调用包内同级函数 extract_gsea_task_pro)
    current_task <- shiny::reactive({
      shiny::req(input$selected_contrast)
      extract_gsea_task_pro(gsea_capsule = res_capsule, task_name = input$selected_contrast, target_collection = "ALL")
    })

    # 构建展示表格
    display_df <- shiny::reactive({
      task <- current_task()
      df <- as.data.frame(task$gsea_res)
      df$NES_Round <- round(df$NES, 3)
      df$Pval <- signif(df$pvalue, 3)
      df$FDR <- signif(df$p.adjust, 3)

      # 植入 JS 触发器
      df$Action <- sprintf('<button class="btn btn-sm btn-primary" onClick="Shiny.setInputValue(\'show_detail\', \'%s\', {priority: \'event\'})">🔍 深度解析</button>', df$ID)

      df[, c("Action", "ID", "Description", "setSize", "NES_Round", "Pval", "FDR", "core_enrichment")]
    })

    output$master_table <- DT::renderDataTable({
      df_show <- display_df()
      df_show$core_enrichment <- NULL # 表格中不显示冗长的 core genes

      DT::datatable(
        df_show,
        escape = FALSE, selection = "multiple",
        colnames = c("操作", "ID", "通路描述", "基因数", "NES", "pvalue", "FDR"),
        options = list(pageLength = 10, scrollX = TRUE, dom = 'Bfrtip')
      ) %>%
        DT::formatStyle('NES_Round', color = DT::styleInterval(0, c('blue', 'red')), fontWeight = 'bold')
    })

    # 联合绘图 (调用包内同级函数 plot_directional_gsea)
    output$multi_gsea_plot <- shiny::renderPlot({
      shiny::req(input$master_table_rows_selected)
      task <- current_task()
      df <- display_df()
      selected_ids <- df$ID[input$master_table_rows_selected]
      colors <- trimws(unlist(strsplit(input$custom_colors, ",")))

      plot_directional_gsea(
        directional_gsea_obj = task,
        target_pathways = selected_ids,
        subPlot = as.numeric(input$plot_subtype),
        curveCol = colors
      )
    })

    # 动线 B: 模态弹窗深度解析
    shiny::observeEvent(input$show_detail, {
      pw_id <- input$show_detail

      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML(sprintf("<b style='color:#0056b3;'>Pathway Deep Dive:</b> %s", pw_id)),
        shiny::fluidRow(
          shiny::column(5,
                        shiny::div(class = "white-box",
                                   shiny::h4("经典 GSEA 轮廓"),
                                   shiny::plotOutput("modal_gsea_plot", height = "500px")
                        )
          ),
          shiny::column(7,
                        shiny::div(class = "white-box",
                                   shiny::h4("核心表达模式热图 (Design 完美复原)"),
                                   shiny::div(class = "scrollable-heatmap",
                                              shiny::plotOutput("modal_heatmap", height = "auto")
                                   )
                        )
          )
        ),
        shiny::fluidRow(
          shiny::column(12,
                        shiny::div(class = "white-box",
                                   shiny::h4("Leading Edge Gene Statistics (核心基因雷达)"),
                                   DT::dataTableOutput("modal_gene_table")
                        )
          )
        ),
        easyClose = TRUE,
        footer = shiny::modalButton("关闭 (Close)"),
        size = "l"
      ))

      output$modal_gsea_plot <- shiny::renderPlot({
        plot_directional_gsea(current_task(), target_pathways = pw_id, subPlot = 3)
      })

      output$modal_heatmap <- shiny::renderPlot({
        task <- current_task()
        dge_list <- task$meta$expr_data
        meta_info <- task$meta

        if (is.null(dge_list)) {
          plot.new(); text(0.5, 0.5, "❌ 胶囊内未包含 expr_data，无法绘制热图", cex = 1.2, col="red"); return()
        }

        sample_meta <- dge_list$samples
        left_samples <- rownames(sample_meta)[sample_meta$group == meta_info$left_group]
        right_samples <- rownames(sample_meta)[sample_meta$group == meta_info$right_group]
        target_samples <- c(left_samples, right_samples)

        expr_mat <- edgeR::cpm(dge_list, log = TRUE)[, target_samples, drop = FALSE]
        sample_meta_sub <- sample_meta[target_samples, , drop = FALSE]

        all_genes <- task$gsea_res@geneSets[[pw_id]]
        expr_genes <- rownames(expr_mat)
        plot_genes <- expr_genes[which(toupper(expr_genes) %in% toupper(all_genes))]

        if (length(plot_genes) < 2) {
          plot.new(); text(0.5, 0.5, "⚠️ 在表达矩阵中匹配到的核心基因不足 2 个", cex = 1.2); return()
        }

        # 按照 Rank 排序
        gene_metrics <- sapply(toupper(plot_genes), function(x) {
          idx <- match(x, toupper(names(task$gsea_res@geneList)))
          if (is.na(idx)) 0 else task$gsea_res@geneList[idx]
        })
        plot_genes <- plot_genes[order(gene_metrics, decreasing = TRUE)]

        plot_mat <- expr_mat[plot_genes, , drop = FALSE]
        plot_mat <- plot_mat[apply(plot_mat, 1, var) > 1e-6, , drop = FALSE]

        z_mat <- t(scale(t(plot_mat)))
        z_mat[is.na(z_mat)] <- 0; z_mat[z_mat > 1] <- 1; z_mat[z_mat < -1] <- -1

        ann_col <- data.frame(Group = sample_meta_sub$group, row.names = rownames(sample_meta_sub))
        grp_col <- c("#E41A1C", "#377EB8"); names(grp_col) <- c(meta_info$left_group, meta_info$right_group)

        pheatmap::pheatmap(
          z_mat, scale = "none", cluster_cols = FALSE, cluster_rows = FALSE, gaps_col = length(left_samples),
          color = grDevices::colorRampPalette(c("#67a9cf", "#f7f7f7", "#ef8a62"))(100), breaks = seq(-1, 1, length.out = 101),
          annotation_col = ann_col, annotation_colors = list(Group = grp_col),
          show_rownames = TRUE, fontsize_row = 10
        )
      }, height = function() {
        task <- current_task()
        match_g <- intersect(toupper(rownames(task$meta$expr_data)), toupper(task$gsea_res@geneSets[[pw_id]]))
        max(400, length(match_g) * 16 + 100)
      })

      output$modal_gene_table <- DT::renderDataTable({
        task <- current_task()
        df <- display_df()

        all_genes <- task$gsea_res@geneSets[[pw_id]]
        match_idx <- which(toupper(names(task$gsea_res@geneList)) %in% toupper(all_genes))
        valid_genes <- names(task$gsea_res@geneList)[match_idx]

        core_str <- df$core_enrichment[df$ID == pw_id]
        core_genes <- unlist(strsplit(as.character(core_str), "/"))

        gene_table <- data.frame(
          Gene = valid_genes,
          Rank_Metric = unname(task$gsea_res@geneList[valid_genes]),
          Status = ifelse(toupper(valid_genes) %in% toupper(core_genes), "✅ Core Edge", "Background"),
          stringsAsFactors = FALSE
        )
        gene_table <- gene_table[order(gene_table$Rank_Metric, decreasing = TRUE), ]

        DT::datatable(gene_table, rownames = FALSE, options = list(pageLength = 5, dom = 'ftip')) %>%
          DT::formatRound('Rank_Metric', 3) %>%
          DT::formatStyle('Status', color = DT::styleEqual(c('✅ Core Edge', 'Background'), c('red', 'grey')))
      })
    })

    # ==========================================
    # 双重视界联动机制
    # ==========================================
    output$volcano_pathway <- plotly::renderPlotly({
      df <- display_df()
      plotly::plot_ly(
        data = df, x = ~NES_Round, y = ~-log10(FDR),
        type = "scatter", mode = "markers",
        text = ~ID, hoverinfo = "text",
        key = ~ID,
        marker = list(color = ifelse(df$NES_Round > 0, "#E41A1C", "#377EB8"), size = 10, opacity = 0.7),
        source = "pathway_volcano"
      ) %>% plotly::layout(xaxis = list(title = "NES"), yaxis = list(title = "-log10(FDR)"), dragmode = "select")
    })

    output$volcano_gene <- plotly::renderPlotly({
      task <- current_task()
      gList <- task$gsea_res@geneList

      plot_df <- data.frame(
        Rank = 1:length(gList), Metric = as.numeric(gList), Gene = names(gList),
        Color = "Background", Size = 4, Alpha = 0.2, stringsAsFactors = FALSE
      )

      click_data <- plotly::event_data("plotly_click", source = "pathway_volcano")
      title_text <- "全局基因 Rank 分布 (请在左侧点击通路 👉)"

      if (!is.null(click_data)) {
        pw_id <- click_data$key
        title_text <- sprintf("当前高亮基因集: %s", pw_id)
        pw_genes <- task$gsea_res@geneSets[[pw_id]]
        hl_idx <- toupper(plot_df$Gene) %in% toupper(pw_genes)

        plot_df$Color[hl_idx] <- "Target"
        plot_df$Size[hl_idx] <- 10
        plot_df$Alpha[hl_idx] <- 1
        plot_df <- plot_df[order(plot_df$Color == "Target"), ] # 置顶高亮层
      }

      plotly::plot_ly(
        data = plot_df, x = ~Rank, y = ~Metric, type = "scattergl", mode = "markers",
        color = ~Color, colors = c("Background" = "#CFD8DC", "Target" = "#FF9800"),
        size = ~Size, sizes = c(4, 10), marker = list(opacity = ~Alpha), text = ~Gene, hoverinfo = "text"
      ) %>%
        plotly::layout(title = list(text = title_text, font = list(size = 14)), showlegend = FALSE)
    })
  }

  message("🚀 正在通过 GSEAlens 包唤醒交互空间...")
  shiny::shinyApp(ui, server)
}
