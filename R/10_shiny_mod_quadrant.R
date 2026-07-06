# Section ----



#' @title Quadrant Linkage Module

#' @description Four-quadrant interactive visualization module for pathway and gene analysis.

#'   Provides synchronized pathway selection, gene ranking, differential expression volcano,

#'   and expression boxplot visualizations with bidirectional table-volcano synchronization.

#' @keywords internal



#' @title Quadrant Linkage UI

#' @description UI components for the four-quadrant interactive visualization

#' @param id Module ID

#' @return Shiny UI tagList

#' @keywords internal



# Subsection ----



mod_quadrant_ui <- function(id) {

  ns <- shiny::NS(id)

  shiny::tagList(

    shiny::fluidRow(

      shiny::column(6, shiny::div(

        class = "white-box",

        shiny::h4("1. Pathway Volcano"),

        shiny::div(

          style = "position: relative;",

          plotly::plotlyOutput(ns("volcano_pathway"), height = "400px"),

          shiny::div(

            style = "position: absolute; top: 10px; left: 10px; background: rgba(255,255,255,0.9); padding: 5px 10px; border-radius: 4px; font-size: 11px; color: #666;",

            "Click pathway to highlight | Check 'Joint Plot' in table"

          )

        ),

        shiny::uiOutput(ns("selected_pathways_display")),

        shiny::div(

          style = "margin-top: 10px;",

          shiny::actionButton(

            ns("open_volcano_export"),

            label = "Export Publication Plot",

            icon = shiny::icon("file-image"),

            class = "btn-success btn-sm",

            style = "width: 100%;"

          )

        )

      )),

      shiny::column(6, shiny::div(

        class = "white-box",

        shiny::h4("2. Gene Rank Distribution"),

        plotly::plotlyOutput(ns("volcano_gene"), height = "400px")

      ))

    ),

    shiny::fluidRow(

      shiny::column(6, shiny::div(

        class = "white-box",

        shiny::h4("3. Differential Expression Volcano"),

        shiny::div(

          style = "margin-bottom: 10px;",

          shiny::checkboxInput(ns("toggle_volcano_settings"), label = "Settings", value = FALSE)

        ),

        shiny::conditionalPanel(

          condition = sprintf("input['%s']", ns("toggle_volcano_settings")),

          shiny::div(

            style = "background: #f8f9fa; padding: 10px; border-radius: 5px; margin-bottom: 10px;",

            shiny::numericInput(ns("volcano_logfc_thresh"), "logFC Threshold:", value = 1, min = 0, max = 22, step = 0.5),

            shiny::numericInput(ns("volcano_pval_thresh"), "P-value Threshold:", value = 0.05, min = 0.001, max = 1, step = 0.01)

          )

        ),

        plotly::plotlyOutput(ns("de_volcano"), height = "400px"),

        shiny::div(

          style = "margin-top: 10px;",

          shiny::actionButton(

            ns("open_de_volcano_export"),

            label = "Export Publication Plot",

            icon = shiny::icon("file-image"),

            class = "btn-success",

            style = "width: 100%;"

          )

        )

      )),

      shiny::column(6, shiny::div(

        class = "white-box",

        shiny::h4("4. Full Expression Distribution"),

        shiny::div(

          style = "margin-bottom: 10px;",

          shiny::checkboxInput(ns("zero_baseline"), label = "Use 0 as baseline", value = TRUE)

        ),

        plotly::plotlyOutput(ns("gene_expr_box"), height = "400px"),

        shiny::uiOutput(ns("boxplot_order_status")),

        shiny::div(

          style = "margin-top: 10px;",

          shiny::actionButton(

            ns("export_boxplot_data_btn"),

            label = "Export Boxplot Data",

            icon = shiny::icon("table"),

            class = "btn-info btn-sm",

            style = "width: 100%;"

          ),

          shiny::helpText(

            "Opens a table of the per-sample expression values (long and wide",

            "formats) for the currently selected gene \u2014 copy to Excel /",

            "GraphPad Prism or download as CSV."

          ),

          shiny::hr(),

          shiny::actionButton(

            ns("open_boxplot_image_export"),

            label = "Export Publication Plot",

            icon = shiny::icon("file-image"),

            class = "btn-success btn-sm",

            style = "width: 100%;"

          ),

          shiny::helpText(

            "Opens a publication-grade ggplot2 export modal (PDF/PNG/SVG +",

            "Live Preview + Copy R Code) for the currently selected gene's",

            "boxplot. Matches the interactive plotly styling."

          )

        )

      ))

    ),

    shiny::hr(),

    shiny::fluidRow(

      shiny::column(12, shiny::div(

        class = "white-box",

        shiny::h4("Gene Expression Table (Click View to display boxplot | Click X to remove)"),

        shiny::div(

          style = "margin-bottom: 10px;",

          shiny::actionButton(ns("clear_all_genes_btn_quadrant"), label = "Clear Table", class = "btn-warning btn-sm"),

          shiny::helpText("Click 'Clear Table' to remove all genes from the table")

        ),

        shiny::div(

          style = "margin-bottom: 10px; font-size: 12px; color: #666;",

          shiny::HTML("<b>Gene Sources:</b> Interest genes from sidebar | DE volcano clicks | Pathway genes (orange)")

        ),

        DT::dataTableOutput(ns("gene_expr_table"))

      ))

    ),

  )

}





#' @title Quadrant Linkage Server

#' @description Server logic with unified gene pool and bidirectional table-volcano sync.

#'   Now properly responds to table selection changes.

#' @param id Module ID

#' @param data_prep_list List from data prep module

#' @param gsea_res GseaRes object

#' @param table_controller Table controller with update_selection method

#' @keywords internal



# Subsection ----



mod_quadrant_server <- function(id, data_prep_list, gsea_res, table_controller) {

  shiny::moduleServer(id, function(input, output, session) {

    ns <- session$ns



    # Color Definitions

    COLOR_LEFT <- "#E41A1C"

    COLOR_RIGHT <- "#377EB8"

    COLOR_NS <- "#C0C0C0"

    COLOR_USER <- "#4DAF4A"

    COLOR_PATHWAY <- "#FF9800"

    COLOR_BOTH <- "#9C27B0"



    # Reactive Values

    # 监听 table_controller$selected_pathways 的变化

    selected_pathway_ids <- shiny::reactiveVal(character(0))

    selected_pathway_genes <- shiny::reactiveVal(character(0))

    current_boxplot_gene <- shiny::reactiveVal(NULL)



    # Get Reactive References

    data_prep_data <- data_prep_list$data

    highlight_genes_reactive <- data_prep_list$highlight_genes

    boxplot_order_ref <- data_prep_list$boxplot_order



    # =========================================

    # 关键：响应主Table勾选状态的变化

    # =========================================

    shiny::observe({

      # 监听 table_controller$selected_pathways 的变化

      table_selected <- table_controller$selected_pathways()



      if (!is.null(table_selected) && length(table_selected) > 0) {

        current <- selected_pathway_ids()

        # 只有当与当前状态不同时才更新

        if (!identical(sort(current), sort(table_selected))) {

          selected_pathway_ids(table_selected)

          message(sprintf("[Quad] Table selection changed: %s", paste(table_selected, collapse = ", ")))



          # 更新pathway基因标记（使用最后一个选中的）

          data_list <- data_prep_data()

          if (!is.null(data_list)) {

            last_id <- tail(table_selected, 1)

            if (!is.null(data_list$gsea_res@geneSets[[last_id]])) {

              pathway_genes <- data_list$gsea_res@geneSets[[last_id]]

              selected_pathway_genes(toupper(pathway_genes))

              message(sprintf("[Quad] Updated pathway genes from: %s (%d genes)", last_id, length(pathway_genes)))

            }

          }

        }

      } else {

        # Table清空时同步清空

        if (length(selected_pathway_ids()) > 0) {

          selected_pathway_ids(character(0))

          selected_pathway_genes(character(0))

          message("[Quad] Table selection cleared")

        }

      }

    })



    # Reset on Contrast Switch

    shiny::observeEvent(data_prep_data(), {

      selected_pathway_ids(character(0))

      selected_pathway_genes(character(0))

      current_boxplot_gene(NULL)

      boxplot_order_ref("default")

    })



    # Selected Pathways Display

    output$selected_pathways_display <- shiny::renderUI({

      sel_ids <- selected_pathway_ids()



      if (length(sel_ids) == 0) {

        shiny::div(

          style = "margin-top: 10px; padding: 10px; background: #f8f9fa; border-radius: 5px;",

          shiny::div(style = "color: #999; font-style: italic;", "No pathways selected")

        )

      } else {

        tag_list <- lapply(sel_ids, function(pid) {

          shiny::tags$span(

            style = "display: inline-block; margin-right: 5px; margin-bottom: 5px; padding: 5px 10px; background: #FF9800; color: white; border-radius: 3px; font-size: 12px;",

            pid,

            shiny::tags$button(

              class = "btn btn-xs",

              style = "margin-left: 5px; padding: 0 5px; background: rgba(255,255,255,0.3); border: none; color: white;",

              onclick = sprintf("Shiny.setInputValue('%s', '%s', {priority: 'event'})", ns("remove_pathway"), pid),

              "x"

            )

          )

        })



        tag_list <- c(tag_list, list(

          shiny::tags$button(

            class = "btn btn-sm btn-warning",

            style = "margin-left: 10px;",

            onclick = sprintf("Shiny.setInputValue('%s', 'CLEAR_ALL', {priority: 'event'})", ns("remove_pathway")),

            "Clear All"

          )

        ))



        shiny::div(

          style = "margin-top: 10px; padding: 10px; background: #f8f9fa; border-radius: 5px;",

          shiny::tags$strong("Selected Pathways: "),

          tag_list

        )

      }

    })



    # Remove Pathway Event

    shiny::observeEvent(input$remove_pathway, {

      pathway_to_remove <- input$remove_pathway

      data_list <- data_prep_data()



      if (pathway_to_remove == "CLEAR_ALL") {

        selected_pathway_ids(character(0))

        selected_pathway_genes(character(0))

        # 同步清空主table勾选

        if (!is.null(table_controller$update_selection)) {

          table_controller$update_selection("clear")

        }

        message("[Pathway] Cleared all selected pathways")

      } else {

        current <- selected_pathway_ids()

        new_selection <- setdiff(current, pathway_to_remove)

        selected_pathway_ids(new_selection)



        # 同步更新主table勾选

        if (!is.null(table_controller$update_selection)) {

          table_controller$update_selection("remove", pathway_to_remove)

        }



        if (length(new_selection) > 0) {

          last_id <- tail(new_selection, 1)

          if (!is.null(data_list) && !is.null(data_list$gsea_res@geneSets[[last_id]])) {

            pathway_genes <- data_list$gsea_res@geneSets[[last_id]]

            selected_pathway_genes(toupper(pathway_genes))

          }

        } else {

          selected_pathway_genes(character(0))

        }

        message(sprintf("[Pathway] Removed: %s", pathway_to_remove))

      }

    })



    # Clear All Genes Button

    shiny::observeEvent(input$clear_all_genes_btn_quadrant, {

      highlight_genes_reactive(character(0))

      message("[Gene] Cleared all gene markers from quadrant")

    })



    # 1. Pathway Volcano Plot

    output$volcano_pathway <- plotly::renderPlotly({

      data_list <- data_prep_data()

      shiny::req(data_list)

      df <- data_list$df

      current_selections <- selected_pathway_ids()



      df$color <- ifelse(df$ID %in% current_selections, COLOR_PATHWAY,

        ifelse(df$NES > 0, COLOR_LEFT, COLOR_RIGHT)

      )

      df$size <- ifelse(df$ID %in% current_selections, 18, 10)

      df$opacity <- ifelse(length(current_selections) == 0, 0.8,

        ifelse(df$ID %in% current_selections, 1.0, 0.35)

      )

      df$linewidth <- ifelse(df$ID %in% current_selections, 3, 1)



      annotations_list <- list()

      if (length(current_selections) > 0) {

        selected_df <- df[df$ID %in% current_selections, ]

        n_sel <- nrow(selected_df)

        # Distribute annotation labels in radial directions to minimise

        # overlap (was: all labels at ay=-30 alternating ax 0/50, which

        # stacked and obscured data points).

        angles <- seq(0, 2 * pi, length.out = 8 + 1)[1:8]

        for (i in seq_len(n_sel)) {

          row <- selected_df[i, ]

          ring  <- (i - 1) %/% 8          # 0 for first 8, 1 for next 8, ...

          radius <- 45 + ring * 25

          ang <- angles[(i - 1) %% 8 + 1]

          ax_offset <- round(radius * cos(ang))

          ay_offset <- round(radius * sin(ang))

          annotations_list[[i]] <- list(

            x = row$NES,

            y = -log10(row$pvalue),

            text = as.character(row$ID),

            showarrow = TRUE,

            arrowhead = 2,

            arrowsize = 0.8,

            arrowwidth = 1.5,

            arrowcolor = COLOR_PATHWAY,

            ax = ax_offset,

            ay = ay_offset,

            font = list(size = 9, color = "#222"),

            bgcolor = "rgba(255,255,255,0.92)",

            bordercolor = COLOR_PATHWAY,

            borderwidth = 1,

            borderpad = 3

          )

        }

      }



      plotly::plot_ly(

        data = df,

        x = ~NES,

        y = ~ -log10(pvalue),

        type = "scatter",

        mode = "markers",

        marker = list(

          color = ~color,

          size = ~size,

          opacity = ~opacity,

          line = list(color = "white", width = ~linewidth)

        ),

        text = ~ sprintf("%s<br>NES: %.2f<br>P-value: %.2e<br>FDR: %.2e", ID, NES, pvalue, p.adjust),

        hoverinfo = "text",

        key = ~ID,

        source = ns("pathway_volcano")

      ) |>

        plotly::layout(

          title = list(

            text = sprintf(

              "Pathway Volcano: %s vs %s<br><sub>%d pathways | %d selected | %d significant (P<0.05)</sub>",

              data_list$left_group, data_list$right_group,

              nrow(df), length(current_selections), sum(df$p.adjust < 0.25, na.rm = TRUE)

            ),

            font = list(size = 14),

            x = 0.5,

            xanchor = "center"

          ),

          xaxis = list(title = "NES", zeroline = FALSE),

          yaxis = list(title = "-log10 (P-value)", zeroline = FALSE),

          showlegend = FALSE,

          dragmode = "pan",

          annotations = annotations_list

        )

    })



    # Pathway Volcano Click Event - 双向同步

    shiny::observeEvent(plotly::event_data("plotly_click", source = ns("pathway_volcano")), {

      click <- plotly::event_data("plotly_click", source = ns("pathway_volcano"))

      if (is.null(click) || is.null(click$key)) {

        return()

      }



      clicked_id <- click$key

      current <- selected_pathway_ids()

      data_list <- data_prep_data()



      if (clicked_id %in% current) {

        # 点击已选中的 → 取消选中

        new_selection <- setdiff(current, clicked_id)

        selected_pathway_ids(new_selection)



        # 同步更新主table勾选

        if (!is.null(table_controller$update_selection)) {

          table_controller$update_selection("remove", clicked_id)

        }



        if (length(new_selection) > 0) {

          last_id <- tail(new_selection, 1)

          if (!is.null(data_list) && !is.null(data_list$gsea_res@geneSets[[last_id]])) {

            pathway_genes <- data_list$gsea_res@geneSets[[last_id]]

            selected_pathway_genes(toupper(pathway_genes))

          }

        } else {

          selected_pathway_genes(character(0))

        }

        message(sprintf("[Pathway] Deselected: %s", clicked_id))

      } else {

        # 点击未选中的 → 添加选中

        new_selection <- c(current, clicked_id)

        selected_pathway_ids(new_selection)



        # 同步更新主table勾选

        if (!is.null(table_controller$update_selection)) {

          table_controller$update_selection("add", clicked_id)

        }



        if (!is.null(data_list) && !is.null(data_list$gsea_res@geneSets[[clicked_id]])) {

          pathway_genes <- data_list$gsea_res@geneSets[[clicked_id]]

          selected_pathway_genes(toupper(pathway_genes))

        }



        message(sprintf("[Pathway] Selected: %s", clicked_id))

      }

    })



    # 2. Gene Rank Distribution

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



      pathway_genes <- selected_pathway_genes()

      rank_df$Color <- COLOR_NS

      rank_df$Size <- 4



      if (length(pathway_genes) > 0) {

        match_idx <- which(toupper(rank_df$Gene) %in% pathway_genes)

        if (length(match_idx) > 0) {

          rank_df$Color[match_idx] <- COLOR_PATHWAY

          rank_df$Size[match_idx] <- 12

        }

      }



      plotly::plot_ly(

        data = rank_df,

        x = ~Rank,

        y = ~Metric,

        type = "scattergl",

        mode = "markers",

        marker = list(color = ~Color, size = ~Size, opacity = 0.8, line = list(width = 0)),

        text = ~Gene,

        hoverinfo = "text"

      ) |>

        plotly::layout(

          xaxis = list(title = "Gene Rank"),

          yaxis = list(title = "Ranking Metric (Stat)"),

          showlegend = FALSE,

          title = list(

            text = ifelse(length(pathway_genes) > 0,

              sprintf("Selected: %d pathway genes highlighted", length(pathway_genes)),

              "Click pathway in volcano above to mark"

            ),

            font = list(size = 12)

          )

        )

    })



    # 3. DE Volcano Plot

    output$de_volcano <- plotly::renderPlotly({

      data_list <- data_prep_data()

      shiny::req(data_list)



      contrast_id <- data_list$contrast_id

      left_group <- data_list$left_group

      right_group <- data_list$right_group



      de_df <- tryCatch(

        {

          get_de_table(gsea_res, contrast_id)

        },

        error = function(e) NULL

      )



      shiny::req(de_df)



      if (!"logFC" %in% colnames(de_df)) de_df$logFC <- de_df$log2FoldChange

      if (!"pvalue" %in% colnames(de_df)) de_df$pvalue <- de_df$p.value

      if (!"padj" %in% colnames(de_df)) de_df$padj <- de_df$p.adjust



      de_df <- de_df[!is.na(de_df$logFC) & !is.na(de_df$pvalue), ]

      if (nrow(de_df) == 0) {

        return(NULL)

      }



      logfc_thresh <- if (!is.null(input$volcano_logfc_thresh)) input$volcano_logfc_thresh else 1

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



      de_df$is_user <- de_df$gene_upper %in% user_genes

      de_df$is_pathway <- de_df$gene_upper %in% pathway_genes

      de_df$is_significant <- abs(de_df$logFC) > logfc_thresh & de_df$pvalue < pval_thresh



      n_up <- sum(de_df$is_significant & de_df$logFC > 0, na.rm = TRUE)

      n_down <- sum(de_df$is_significant & de_df$logFC < 0, na.rm = TRUE)

      n_not_sig <- sum(!de_df$is_significant, na.rm = TRUE)

      n_user <- sum(de_df$is_user, na.rm = TRUE)

      n_pathway <- sum(de_df$is_pathway, na.rm = TRUE)

      n_both <- sum(de_df$is_user & de_df$is_pathway, na.rm = TRUE)



      de_df$color <- dplyr::case_when(

        de_df$is_user & de_df$is_pathway ~ COLOR_BOTH,

        de_df$is_user ~ COLOR_USER,

        de_df$is_pathway ~ COLOR_PATHWAY,

        de_df$is_significant & de_df$logFC > 0 ~ COLOR_LEFT,

        de_df$is_significant & de_df$logFC < 0 ~ COLOR_RIGHT,

        TRUE ~ COLOR_NS

      )



      de_df$size <- dplyr::case_when(

        de_df$is_user | de_df$is_pathway ~ 15,

        de_df$is_significant ~ 9,

        TRUE ~ 4

      )



      de_df$opacity <- dplyr::case_when(

        de_df$is_user | de_df$is_pathway ~ 1.0,

        de_df$is_significant ~ 0.7,

        TRUE ~ 0.5

      )



      de_df$linewidth <- dplyr::case_when(

        de_df$is_user | de_df$is_pathway ~ 1.0,

        TRUE ~ 0

      )



      de_df$plot_order <- dplyr::case_when(

        de_df$is_user | de_df$is_pathway ~ 3,

        de_df$is_significant ~ 2,

        TRUE ~ 1

      )

      de_df <- de_df[order(de_df$plot_order), ]



      current_selected <- selected_pathway_ids()

      pathway_label <- if (length(current_selected) > 0) {

        sprintf(" | From: %s", tail(current_selected, 1))

      } else {

        ""

      }



      title_text <- sprintf(

        "%s vs %s<br><sup>Up %d | Down %d | NS %d | Selected %d | Pway %d | Both %d%s</sup>",

        left_group, right_group, n_up, n_down, n_not_sig, n_user, n_pathway, n_both, pathway_label

      )



      annotations_list <- list()

      max_y_val <- max(de_df$y_axis, na.rm = TRUE)



      annotations_list[[1]] <- list(

        x = 0.99, y = 0.99, xref = "paper", yref = "paper",

        text = paste0("<b style='color:", COLOR_LEFT, ";'>High in ", left_group, "</b>"),

        showarrow = FALSE, font = list(size = 14),

        xanchor = "right", yanchor = "top",

        bgcolor = "rgba(255,255,255,0.9)",

        bordercolor = COLOR_LEFT, borderwidth = 2, borderpad = 6

      )

      annotations_list[[2]] <- list(

        x = 0.01, y = 0.99, xref = "paper", yref = "paper",

        text = paste0("<b style='color:", COLOR_RIGHT, ";'>High in ", right_group, "</b>"),

        showarrow = FALSE, font = list(size = 14),

        xanchor = "left", yanchor = "top",

        bgcolor = "rgba(255,255,255,0.9)",

        bordercolor = COLOR_RIGHT, borderwidth = 2, borderpad = 6

      )



      user_genes_df <- de_df[de_df$is_user, ]

      if (nrow(user_genes_df) > 0) {

        for (i in seq_len(min(nrow(user_genes_df), 20))) {

          gene <- user_genes_df[i, ]

          gene_color <- if (gene$is_user && gene$is_pathway) COLOR_BOTH else COLOR_USER

          annotations_list[[length(annotations_list) + 1]] <- list(

            x = gene$x_axis,

            y = gene$y_axis,

            text = gene$gene_symbol,

            showarrow = TRUE,

            arrowhead = 0,

            arrowsize = 1,

            arrowwidth = 2,

            arrowcolor = gene_color,

            ax = ifelse(gene$logFC > 0, 50, -50),

            ay = -35,

            bgcolor = "rgba(255,255,255,0.85)",

            bordercolor = gene_color,

            borderwidth = 0.5,

            font = list(size = 12, color = gene_color)

          )

        }

      }



      if (logfc_thresh > 0) {

        annotations_list[[length(annotations_list) + 1]] <- list(

          x = logfc_thresh, y = max_y_val * 0.95,

          text = sprintf("|logFC|=%.1f", logfc_thresh),

          showarrow = FALSE, font = list(size = 10, color = "gray")

        )

        annotations_list[[length(annotations_list) + 1]] <- list(

          x = -logfc_thresh, y = max_y_val * 0.95,

          text = sprintf("|logFC|=%.1f", logfc_thresh),

          showarrow = FALSE, font = list(size = 10, color = "gray")

        )

      }





      # ============================================================

      # 分层渲染：拆分数据为多个子集（从底到顶）

      # ============================================================



      df_ns <- de_df[!de_df$is_significant, ]

      df_up <- de_df[de_df$is_significant & de_df$logFC > 0, ]

      df_down <- de_df[de_df$is_significant & de_df$logFC < 0, ]

      df_pathway <- de_df[de_df$is_pathway & !de_df$is_user, ]

      df_user <- de_df[de_df$is_user & !de_df$is_pathway, ]

      df_both <- de_df[de_df$is_user & de_df$is_pathway, ]



      # 初始化 plot_ly 时设置 source（点击事件需要）

      p <- plotly::plot_ly(

        data = de_df,

        source = ns("deg_volcano")

      )



      # Layer 1: NS (灰色，最底层)

      if (nrow(df_ns) > 0) {

        p <- p |> plotly::add_trace(

          data = df_ns,

          x = df_ns$x_axis,

          y = df_ns$y_axis,

          type = "scatter",

          mode = "markers",

          marker = list(color = COLOR_NS, size = 4, opacity = 0.5, line = list(color = "white", width = 0)),

          text = sprintf("%s<br>logFC: %.2f<br>-log10(p): %.2f<br>FDR: %.2e", df_ns$gene_symbol, df_ns$logFC, df_ns$y_axis, df_ns$padj),

          hoverinfo = "text",

          key = df_ns$gene_upper,

          showlegend = FALSE,

          inherit = FALSE

        )

      }



      # Layer 2: UP (红色)

      if (nrow(df_up) > 0) {

        p <- p |> plotly::add_trace(

          data = df_up,

          x = df_up$x_axis,

          y = df_up$y_axis,

          type = "scatter",

          mode = "markers",

          marker = list(color = COLOR_LEFT, size = 9, opacity = 0.7, line = list(color = "white", width = 0)),

          text = sprintf("%s<br>logFC: %.2f<br>-log10(p): %.2f<br>FDR: %.2e", df_up$gene_symbol, df_up$logFC, df_up$y_axis, df_up$padj),

          hoverinfo = "text",

          key = df_up$gene_upper,

          showlegend = FALSE,

          inherit = FALSE

        )

      }



      # Layer 3: DOWN (蓝色)

      if (nrow(df_down) > 0) {

        p <- p |> plotly::add_trace(

          data = df_down,

          x = df_down$x_axis,

          y = df_down$y_axis,

          type = "scatter",

          mode = "markers",

          marker = list(color = COLOR_RIGHT, size = 9, opacity = 0.7, line = list(color = "white", width = 0)),

          text = sprintf("%s<br>logFC: %.2f<br>-log10(p): %.2f<br>FDR: %.2e", df_down$gene_symbol, df_down$logFC, df_down$y_axis, df_down$padj),

          hoverinfo = "text",

          key = df_down$gene_upper,

          showlegend = FALSE,

          inherit = FALSE

        )

      }



      # Layer 4: 仅Pathway (橙色)

      if (nrow(df_pathway) > 0) {

        p <- p |> plotly::add_trace(

          data = df_pathway,

          x = df_pathway$x_axis,

          y = df_pathway$y_axis,

          type = "scatter",

          mode = "markers",

          marker = list(color = COLOR_PATHWAY, size = 15, opacity = 1.0, line = list(color = "white", width = 1.5)),

          text = sprintf("%s<br>logFC: %.2f<br>-log10(p): %.2f<br>FDR: %.2e<br>[Pathway]", df_pathway$gene_symbol, df_pathway$logFC, df_pathway$y_axis, df_pathway$padj),

          hoverinfo = "text",

          key = df_pathway$gene_upper,

          showlegend = FALSE,

          inherit = FALSE

        )

      }



      # Layer 5: 仅User (绿色)

      if (nrow(df_user) > 0) {

        p <- p |> plotly::add_trace(

          data = df_user,

          x = df_user$x_axis,

          y = df_user$y_axis,

          type = "scatter",

          mode = "markers",

          marker = list(color = COLOR_USER, size = 15, opacity = 1.0, line = list(color = "white", width = 1.5)),

          text = sprintf("%s<br>logFC: %.2f<br>-log10(p): %.2f<br>FDR: %.2e<br>[Selected]", df_user$gene_symbol, df_user$logFC, df_user$y_axis, df_user$padj),

          hoverinfo = "text",

          key = df_user$gene_upper,

          showlegend = FALSE,

          inherit = FALSE

        )

      }



      # Layer 6: Both (紫色，最顶层)

      if (nrow(df_both) > 0) {

        p <- p |> plotly::add_trace(

          data = df_both,

          x = df_both$x_axis,

          y = df_both$y_axis,

          type = "scatter",

          mode = "markers",

          marker = list(color = COLOR_BOTH, size = 18, opacity = 1.0, line = list(color = "white", width = 2)),

          text = sprintf("%s<br>logFC: %.2f<br>-log10(p): %.2f<br>FDR: %.2e<br>[Both]", df_both$gene_symbol, df_both$logFC, df_both$y_axis, df_both$padj),

          hoverinfo = "text",

          key = df_both$gene_upper,

          showlegend = FALSE,

          inherit = FALSE

        )

      }



      # 布局设置

      p <- p |>

        plotly::layout(

          title = list(text = title_text, font = list(size = 14), x = 0.5, xanchor = "center"),

          xaxis = list(title = "logFC", zeroline = FALSE, showgrid = TRUE, gridcolor = "lightgray"),

          yaxis = list(title = "-log10 (P-value)", zeroline = FALSE, showgrid = TRUE, gridcolor = "lightgray"),

          showlegend = FALSE,

          dragmode = "pan",

          annotations = annotations_list,

          shapes = list(

            list(

              type = "line", x0 = logfc_thresh, x1 = logfc_thresh,

              y0 = 0, y1 = max_y_val * 1.05,

              line = list(color = "gray", dash = "dash", width = 1)

            ),

            list(

              type = "line", x0 = -logfc_thresh, x1 = -logfc_thresh,

              y0 = 0, y1 = max_y_val * 1.05,

              line = list(color = "gray", dash = "dash", width = 1)

            ),

            list(

              type = "line", x0 = min(de_df$x_axis) * 1.1, x1 = max(de_df$x_axis) * 1.1,

              y0 = -log10(pval_thresh), y1 = -log10(pval_thresh),

              line = list(color = "gray", dash = "dash", width = 1)

            )

          )

        )



      p

    })



    # DE Volcano Click Event

    shiny::observeEvent(plotly::event_data("plotly_click", source = ns("deg_volcano")), {

      click <- plotly::event_data("plotly_click", source = ns("deg_volcano"))

      if (is.null(click) || is.null(click$key)) {

        return()

      }



      clicked_gene <- click$key

      current_boxplot_gene(clicked_gene)



      current_applied <- highlight_genes_reactive()

      highlight_genes_reactive(union(toupper(current_applied), toupper(clicked_gene)))

    })



    # 4. Gene Expression Table

    output$gene_expr_table <- DT::renderDataTable({

      genes <- highlight_genes_reactive()

      shiny::req(length(genes) > 0)



      data_list <- data_prep_data()

      shiny::req(data_list)

      symbol_map <- .rebuild_symbol_map(gsea_res, data_list$contrast_id)



      left_grp <- data_list$left_group

      right_grp <- data_list$right_group



      # 获取 DE volcano 的阈值

      pval_thresh <- input$volcano_pval_thresh %||% 0.05

      logfc_thresh <- input$volcano_logfc_thresh %||% 1



      user_genes <- toupper(highlight_genes_reactive())

      pathway_genes <- selected_pathway_genes()



      de_df <- tryCatch(

        {

          get_de_table(gsea_res, data_list$contrast_id)

        },

        error = function(e) NULL

      )



      table_data <- lapply(genes, function(g) {

        gene_upper <- toupper(g)

        is_user <- gene_upper %in% user_genes

        is_pathway <- gene_upper %in% pathway_genes



        logfc <- NA_real_

        pval <- NA_real_

        padj_val <- NA_real_

        if (!is.null(de_df)) {

          if ("gene_symbol" %in% colnames(de_df)) {

            idx <- which(toupper(de_df$gene_symbol) == gene_upper)

          } else {

            idx <- which(rownames(de_df) == g)

          }

          if (length(idx) > 0) {

            row <- de_df[idx[1], ]

            logfc <- if ("logFC" %in% colnames(row)) row$logFC else row$log2FoldChange

            pval <- row$pvalue

            padj_val <- row$padj

          }

        }



        updown <- ifelse(is.na(logfc), "-", ifelse(logfc > 0, "UP", "DOWN"))

        high_in <- ifelse(is.na(logfc), "-", ifelse(logfc > 0, left_grp, right_grp))



        # 保持数值为数值类型，不转字符串

        pval_display <- ifelse(is.na(pval), NA_real_, pval)

        padj_display <- ifelse(is.na(padj_val), NA_real_, padj_val)



        # 使用 DE volcano 阈值的 Significant 判断

        sig_label <- {

          if (is.na(pval) || is.na(logfc)) {

            "-"

          } else if (pval < pval_thresh && logfc > logfc_thresh) {

            "upreg"

          } else if (pval < pval_thresh && logfc < -logfc_thresh) {

            "downreg"

          } else {

            "ns"

          }

        }



        # ✅ 修复后的Delete和View按钮生成

        delete_btn <- sprintf(

          '<button class="btn btn-xs btn-danger" onclick="Shiny.setInputValue(\'%s\', {gene_to_remove: \'%s\', time: Date.now()}, {priority: \'event\'});" style="margin-right: 5px;">X</button>',

          ns("delete_gene_action"), g

        )



        view_btn <- sprintf(

          '<button class="btn btn-xs btn-primary" onclick="Shiny.setInputValue(\'%s\', {gene_to_view: \'%s\', time: Date.now()}, {priority: \'event\'});">View</button>',

          ns("view_boxplot_action"), g

        )



        data.frame(

          Gene = .get_display_symbol(g, symbol_map),

          Log2FC = ifelse(is.na(logfc), "-", sprintf("%.3f", logfc)),

          Pvalue = pval_display,

          Padjust = padj_display,

          UpDown = updown,

          Significant = sig_label,

          HighIn = high_in,

          IsPathway = ifelse(is_pathway, "Yes", "No"),

          Delete = delete_btn,

          ViewBox = view_btn,

          stringsAsFactors = FALSE,

          check.names = FALSE

        )

      })



      table_df <- do.call(rbind, table_data)



      # 按 P-value 升序排序

      table_df <- table_df[order(table_df[["Pvalue"]], na.last = TRUE), ]



      # JS 渲染函数：科学计数法

      js_pval_render <- "

        function(data, type) {

          if (type === 'sort' || type === 'type') return data;

          if (data === null || isNaN(data)) return '-';

          if (Math.abs(data) < 0.001) return data.toExponential(2);

          return data.toFixed(4);

        }

      "



      dt <- DT::datatable(

        table_df,

        escape = FALSE,

        rownames = FALSE,

        selection = "none",

        extensions = c("Scroller"),

        options = list(

          pageLength = 10,

          scrollY = "50vh",

          scroller = TRUE,

          dom = "frtip",

          search = list(caseInsensitive = TRUE),

          columnDefs = list(

            list(orderable = FALSE, targets = c(7, 8)),

            list(className = "dt-center", targets = c(0, 2, 3, 4, 5, 6)),

            list(

              targets = 2,

              render = htmlwidgets::JS(js_pval_render)

            ),

            list(

              targets = 3,

              render = htmlwidgets::JS(js_pval_render)

            )

          ),

          order = list(list(2, "asc"))

        )

      )



      # Up/Down 列样式

      dt <- dt |> DT::formatStyle(

        columns = "UpDown",

        backgroundColor = DT::styleEqual(c("UP", "DOWN", "-"), c("#FFCDD2", "#BBDEFB", "transparent")),

        color = DT::styleEqual(c("UP", "DOWN", "-"), c(COLOR_LEFT, COLOR_RIGHT, "#666")),

        fontWeight = DT::styleEqual(c("UP", "DOWN"), c("bold", "bold"))

      )



      # Significant 列样式

      dt <- dt |> DT::formatStyle(

        columns = "Significant",

        color = DT::styleEqual(c("upreg", "downreg", "ns"), c("#E41A1C", "#377EB8", "#999999")),

        fontWeight = "bold"

      )



      # High In 列样式

      dt <- dt |> DT::formatStyle(

        columns = "HighIn",

        color = DT::styleEqual(c(left_grp, right_grp, "-"), c("#E41A1C", "#377EB8", "#999999"))

      )



      # In Pathway 列样式

      dt <- dt |> DT::formatStyle(

        columns = "IsPathway",

        backgroundColor = DT::styleEqual(c("Yes", "No"), c("#9C27B0", "transparent")),

        color = DT::styleEqual(c("Yes", "No"), c("#FFFFFF", "#999999")),

        fontWeight = DT::styleEqual(c("Yes", "No"), c("bold", "normal"))

      )



      # P-value 列样式

      dt <- dt |> DT::formatStyle(

        columns = "Pvalue",

        color = DT::styleInterval(0.001, c("red", "black")),

        fontWeight = DT::styleInterval(0.05, c("bold", "normal"))

      )



      # Adj P-value 列样式

      dt <- dt |> DT::formatStyle(

        columns = "Padjust",

        color = DT::styleInterval(0.001, c("red", "black")),

        fontWeight = DT::styleInterval(0.05, c("bold", "normal"))

      )



      dt

    })



    # ✅ 修复Delete基因事件监听

    shiny::observeEvent(input$delete_gene_action,

      {

        action_data <- input$delete_gene_action

        if (is.null(action_data) || is.null(action_data$gene_to_remove)) {

          return()

        }



        gene_to_remove <- action_data$gene_to_remove

        current_applied <- highlight_genes_reactive()

        new_applied <- setdiff(toupper(current_applied), toupper(gene_to_remove))

        highlight_genes_reactive(new_applied)



        # 同时清除对应的箱线图

        if (!is.null(current_boxplot_gene()) && toupper(current_boxplot_gene()) == toupper(gene_to_remove)) {

          current_boxplot_gene(NULL)

        }



        message(sprintf("[Gene] Removed from table: %s", gene_to_remove))

        shiny::showNotification(sprintf("Gene removed: %s", gene_to_remove), type = "message", duration = 2)

      },

      ignoreNULL = TRUE

    )



    # ✅ 修复View箱线图事件监听

    shiny::observeEvent(input$view_boxplot_action,

      {

        action_data <- input$view_boxplot_action

        if (is.null(action_data) || is.null(action_data$gene_to_view)) {

          return()

        }



        gene_name <- action_data$gene_to_view

        current_boxplot_gene(gene_name)



        message(sprintf("[Boxplot] Displaying: %s", gene_name))

        shiny::showNotification(sprintf("Showing boxplot for: %s", gene_name), type = "message", duration = 2)

      },

      ignoreNULL = TRUE

    )



    # 5. Expression Boxplot

    output$gene_expr_box <- plotly::renderPlotly({

      current_gene <- current_boxplot_gene()



      if (is.null(current_gene)) {

        return(plotly::plot_ly() |> plotly::layout(

          title = list(text = "Click gene in DE volcano or table to view expression", font = list(size = 14))

        ))

      }



      data_list <- data_prep_data()

      shiny::req(data_list)



      tryCatch(

        {

          expr_mat <- get_expr_matrix(gsea_res, type = data_list$expression_type)

          sample_meta <- get_sample_meta(gsea_res)



          if (is.null(expr_mat) || is.null(sample_meta)) {

            return(plotly::plot_ly() |> plotly::layout(title = "Expression matrix not available"))

          }



          target_gene_upper <- toupper(current_gene)

          gene_names_upper <- toupper(rownames(expr_mat))

          match_idx <- which(gene_names_upper == target_gene_upper)



          if (length(match_idx) == 0) {

            gene_meta <- gsea_res$expr_bundle$gene_meta

            if (!is.null(gene_meta) && nrow(gene_meta) > 0) {

              rownames(gene_meta) <- if (is.null(rownames(gene_meta))) rownames(expr_mat) else rownames(gene_meta)

              symbol_col <- intersect(c("SYMBOL", "symbol", "Gene", "gene_name", "gene_symbol"), colnames(gene_meta))[1]

              if (!is.na(symbol_col)) {

                meta_symbols_upper <- toupper(as.character(gene_meta[[symbol_col]]))

                meta_matches <- which(meta_symbols_upper == target_gene_upper)

                if (length(meta_matches) > 0) {

                  ensembl_id <- rownames(gene_meta)[meta_matches[1]]

                  match_idx <- which(rownames(expr_mat) == ensembl_id)

                }

              }

            }

          }



          if (length(match_idx) == 0 || is.na(match_idx)) {

            return(plotly::plot_ly() |> plotly::layout(

              title = sprintf("Gene '%s' not found", current_gene)

            ))

          }



          symbol_map <- .rebuild_symbol_map(gsea_res, data_list$contrast_id)

          actual_gene <- rownames(expr_mat)[match_idx[1]]

          display_gene_name <- .get_display_symbol(current_gene, symbol_map)



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

            return(plotly::plot_ly() |> plotly::layout(title = "No valid group data"))

          }



          current_confirmed_order <- boxplot_order_ref()

          final_order_to_use <- if (!is.null(current_confirmed_order) &&

            current_confirmed_order != "default" &&

            current_confirmed_order != "") {

            current_confirmed_order

          } else {

            "default"

          }



          actual_groups <- unique(as.character(plot_data$Group))

          x_categories <- NULL



          if (final_order_to_use != "default" && final_order_to_use != "" && !is.na(final_order_to_use)) {

            sep <- if (grepl("->", final_order_to_use, fixed = TRUE)) "->" else ","

            order_parts <- strsplit(final_order_to_use, sep)[[1]]

            order_parts <- trimws(order_parts)

            valid_parts <- order_parts[order_parts %in% actual_groups]

            x_categories <- c(valid_parts, setdiff(actual_groups, valid_parts))

          } else {

            x_categories <- actual_groups

          }



          plot_data <- plot_data[plot_data$Group %in% x_categories, ]

          plot_data$Group <- factor(plot_data$Group, levels = x_categories, ordered = TRUE)



          unique_groups <- levels(plot_data$Group)

          if (length(unique_groups) == 2) {

            group_colors <- c(COLOR_LEFT, COLOR_RIGHT)

            names(group_colors) <- unique_groups

          } else {

            group_colors <- c(

              "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",

              "#FF7F00", "#A65628", "#F781BF", "#999999"

            )[seq_along(unique_groups)]

            names(group_colors) <- unique_groups

          }



          use_zero_baseline <- input$zero_baseline %||% FALSE



          y_min <- min(plot_data$Expression, na.rm = TRUE)

          y_max <- max(plot_data$Expression, na.rm = TRUE)



          if (use_zero_baseline) {

            y_min <- min(y_min, 0)

            y_max <- max(y_max, 0)

          }



          y_range <- y_max - y_min

          y_min <- y_min - y_range * 0.1

          y_max <- y_max + y_range * 0.1



          p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Group, y = Expression, fill = Group)) +

            ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA) +

            ggplot2::geom_jitter(

              width = 0.2, size = 3, alpha = 0.6,

              ggplot2::aes(text = sprintf(

                "<b>Sample:</b> %s<br><b>Group:</b> %s<br><b>Expression:</b> %.3f",

                Sample, Group, Expression

              ))

            ) +

            ggplot2::scale_fill_manual(values = group_colors) +

            ggplot2::scale_x_discrete(limits = x_categories, drop = FALSE) +

            ggplot2::coord_cartesian(ylim = c(y_min, y_max)) +

            ggplot2::theme_bw(base_size = 12) +

            ggplot2::labs(

              title = sprintf("%s", display_gene_name),

              y = data_list$expression_type, x = NULL

            ) +

            ggplot2::theme(

              legend.position = "none",

              axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)

            )



          if (use_zero_baseline) {

            p <- p + ggplot2::geom_hline(

              yintercept = 0, linetype = "dashed",

              color = "red", alpha = 0.7, size = 0.8

            )

          }



          ply <- plotly::ggplotly(p, tooltip = "text")

          ply |> plotly::layout(

            xaxis = list(categoryorder = "array", categoryarray = x_categories, title = ""),

            dragmode = FALSE

          )

        },

        error = function(e) {

          return(plotly::plot_ly() |> plotly::layout(

            title = sprintf("Error: %s", e$message)

          ))

        }

      )

    })



    # Boxplot Order Status

    output$boxplot_order_status <- shiny::renderUI({

      order_info <- boxplot_order_ref()

      if (is.null(order_info) || order_info == "default") {

        shiny::div(style = "font-size: 12px; color: #666;", "Order: Default")

      } else {

        shiny::div(

          style = "font-size: 12px; color: #666;",

          "Order:", gsub(",", " -> ", order_info)

        )

      }

    })



    # ==========================================================================

    # 7. Export Boxplot Data Button (panel 4, "Full Expression Distribution")

    # --------------------------------------------------------------------------

    # Opens a modal with two TSV tables (long + wide) of the per-sample

    # expression values for the currently selected gene. Each table has its

    # own "Copy" / "Download CSV" controls. The user can paste the wide

    # table directly into GraphPad Prism (Column format) or Excel.

    # ==========================================================================

    current_boxplot_data <- shiny::reactiveVal(NULL)

    current_boxplot_filenames <- shiny::reactiveVal(list(long = "long.csv", wide = "wide.csv"))



    # Helper: turn an arbitrary string into a filename-safe token

    .sanitize_filename <- function(x) {

      x <- gsub("[^A-Za-z0-9_.-]+", "_", x %||% "gene")

      x <- gsub("^_+|_+$", "", x)

      if (!nzchar(x)) x <- "gene"

      x

    }



    # Helper: write a string to the system clipboard (uses clipr when available)

    .write_to_clipboard <- function(text) {

      if (requireNamespace("clipr", quietly = TRUE)) {

        clipr::write_clip(text)

        return(invisible(TRUE))

      }

      con <- file("clipboard", "w")

      on.exit(close(con), add = TRUE)

      writeLines(text, con)

      invisible(TRUE)

    }



    shiny::observeEvent(input$export_boxplot_data_btn, {

      data_list <- data_prep_data()

      shiny::req(data_list)



      current_gene <- current_boxplot_gene()

      if (is.null(current_gene) || !nzchar(current_gene)) {

        shiny::showNotification(

          "No gene selected. Click a gene in the DE volcano or in the gene table first.",

          type = "warning",

          duration = 3

        )

        return()

      }



      # Use the same display symbol that the boxplot title shows

      symbol_map <- .rebuild_symbol_map(gsea_res, data_list$contrast_id)

      display_gene <- .get_display_symbol(current_gene, symbol_map)



      result <- tryCatch(

        extract_boxplot_data(

          GSEAlens_res = gsea_res,

          contrast_id  = data_list$contrast_id,

          gene_symbol  = display_gene,

          expr_type    = data_list$expression_type

        ),

        error = function(e) {

          shiny::showNotification(

            sprintf("Failed to extract boxplot data: %s", e$message),

            type = "error",

            duration = 5

          )

          NULL

        }

      )

      if (is.null(result)) return()



      current_boxplot_data(result)



      # Pre-build the two TSV strings and the safe filenames

      tsv_long <- .serialize_df_tsv(result$long)

      tsv_wide <- .serialize_df_tsv(result$wide)

      safe_gene     <- .sanitize_filename(result$gene)

      safe_contrast <- .sanitize_filename(result$contrast_id)

      current_boxplot_filenames(list(

        long = sprintf("GSEAlens_%s_%s_box_long.csv", safe_gene, safe_contrast),

        wide = sprintf("GSEAlens_%s_%s_box_wide.csv", safe_gene, safe_contrast)

      ))



      # Format the sample-count line for the header

      n_samples <- nrow(result$long)

      group_counts <- table(result$long$Group)

      group_summary <- paste(

        sprintf("%s: %d", names(group_counts), as.integer(group_counts)),

        collapse = ", "

      )



      shiny::showModal(shiny::modalDialog(

        title = sprintf(

          "Expression Data: %s (n = %d samples; %s)",

          result$gene, n_samples, group_summary

        ),

        size = "l",

        easyClose = TRUE,

        footer = shiny::modalButton("Close"),

        shiny::tagList(

          # ---- Section 1: Wide format ------------------------------------

          shiny::div(

            style = "margin-top: 4px;",

            shiny::tags$b("Wide format"),

            shiny::helpText(

              style = "margin: 0 0 6px 0;",

              "First row = sample IDs, second row = expression values. ",

              "Left-most cell = gene name. Paste into GraphPad Prism (Column data) ",

              "or Excel."

            )

          ),

          shiny::div(

            style = "display: flex; gap: 8px; align-items: center; margin-bottom: 6px;",

            shiny::actionButton(

              ns("copy_wide_tsv_btn"),

              label = "Copy Wide",

              icon = shiny::icon("copy"),

              class = "btn-success btn-sm"

            ),

            shiny::downloadButton(

              ns("download_wide_csv_btn"),

              label = "Download Wide CSV",

              class = "btn-info btn-sm"

            )

          ),

          shiny::div(

            style = "background: #f5f5f5; padding: 12px; border-radius: 5px; max-height: 160px; overflow: auto;",

            shiny::tags$pre(

              shiny::code(tsv_wide),

              style = "font-size: 11px; white-space: pre; word-break: keep-all; margin: 0;"

            )

          ),



          shiny::hr(),



          # ---- Section 2: Long format ------------------------------------

          shiny::div(

            shiny::tags$b("Long format"),

            shiny::helpText(

              style = "margin: 0 0 6px 0;",

              "One row per sample. Suitable for GraphPad Prism (Grouped data) ",

              "and most R / Python plotting tools."

            )

          ),

          shiny::div(

            style = "display: flex; gap: 8px; align-items: center; margin-bottom: 6px;",

            shiny::actionButton(

              ns("copy_long_tsv_btn"),

              label = "Copy Long",

              icon = shiny::icon("copy"),

              class = "btn-success btn-sm"

            ),

            shiny::downloadButton(

              ns("download_long_csv_btn"),

              label = "Download Long CSV",

              class = "btn-info btn-sm"

            )

          ),

          shiny::div(

            style = "background: #f5f5f5; padding: 12px; border-radius: 5px; max-height: 220px; overflow: auto;",

            shiny::tags$pre(

              shiny::code(tsv_long),

              style = "font-size: 11px; white-space: pre; word-break: keep-all; margin: 0;"

            )

          )

        )

      ))

    })



    # Copy handlers ---------------------------------------------------------

    shiny::observeEvent(input$copy_wide_tsv_btn, {

      result <- current_boxplot_data()

      shiny::req(result)

      tsv <- .serialize_df_tsv(result$wide)

      .write_to_clipboard(tsv)

      shiny::showNotification("Wide table copied to clipboard.", type = "message", duration = 2)

    })



    shiny::observeEvent(input$copy_long_tsv_btn, {

      result <- current_boxplot_data()

      shiny::req(result)

      tsv <- .serialize_df_tsv(result$long)

      .write_to_clipboard(tsv)

      shiny::showNotification("Long table copied to clipboard.", type = "message", duration = 2)

    })



    # Download handlers -----------------------------------------------------

    output$download_wide_csv_btn <- shiny::downloadHandler(

      filename = function() current_boxplot_filenames()$wide,

      content  = function(file) {

        result <- current_boxplot_data()

        tsv <- if (is.null(result)) "" else .serialize_df_tsv(result$wide)

        # writeLines adds a final "\n", matching the modal preview

        writeLines(if (nzchar(tsv)) tsv else "Gene\t(empty)\n", file, useBytes = TRUE)

      }

    )



    output$download_long_csv_btn <- shiny::downloadHandler(

      filename = function() current_boxplot_filenames()$long,

      content  = function(file) {

        result <- current_boxplot_data()

        tsv <- if (is.null(result)) "" else .serialize_df_tsv(result$long)

        writeLines(if (nzchar(tsv)) tsv else "Sample\tGroup\tExpression\n", file, useBytes = TRUE)

      }

    )



    # ============================================================

    # Boxplot Image Export Center (panel 4, "Full Expression Distribution")

    # ------------------------------------------------------------

    # Publication-grade ggplot2 export for the selected gene's boxplot.

    # Mirrors the volcano export modal pattern + internal/external preview.

    # ============================================================

    .build_boxplot_export_modal <- function() {

      shiny::modalDialog(

        title = "Export Publication Plot - Gene Boxplot",

        size = "l", footer = NULL, easyClose = TRUE,

        shiny::fluidRow(

          shiny::column(5,

            shiny::h5("Dimensions"),

            shiny::fluidRow(

              shiny::column(6,

                shiny::numericInput(ns("box_exp_width"),  "Width (inch)",  value = 5, min = 2, max = 20)

              ),

              shiny::column(6,

                shiny::numericInput(ns("box_exp_height"), "Height (inch)", value = 5, min = 2, max = 20)

              )

            ),

            shiny::numericInput(ns("box_exp_dpi"), "DPI", value = 300, min = 72, max = 600),

            shiny::hr(),

            shiny::h5("Download"),

            shiny::fluidRow(

              shiny::column(6,

                shiny::downloadButton(ns("box_exp_download_pdf"),

                  "Download PDF", class = "btn-danger btn-block",

                  icon = shiny::icon("file-pdf"))

              ),

              shiny::column(6,

                shiny::downloadButton(ns("box_exp_download_png"),

                  "Download PNG", class = "btn-primary btn-block",

                  icon = shiny::icon("file-image"))

              )

            ),

            shiny::hr(),

            shiny::fluidRow(

              shiny::column(6,

                shiny::selectInput(ns("box_exp_format"), "Other format",

                  choices = c("SVG" = "svg", "TIFF" = "tiff"), selected = "svg")

              ),

              shiny::column(6,

                shiny::div(style = "margin-top: 22px;",

                  shiny::downloadButton(ns("box_exp_download_other"),

                    "Download", class = "btn-default btn-block")

                )

              )

            ),

            shiny::hr(),

            shiny::fluidRow(

              shiny::column(6,

                shiny::actionButton(ns("box_exp_copy_code"),

                  "Copy R Code", class = "btn-info btn-block",

                  icon = shiny::icon("clipboard"))

              ),

              shiny::column(6,

                shiny::actionButton(ns("box_exp_dismiss"),

                  "Close", class = "btn-default btn-block")

              )

            )

          ),

          shiny::column(7,

            shiny::h5("Live Preview"),

            shiny::div(style = "background:#f5f5f5; border:1px solid #ddd; border-radius:4px; padding:8px; width:100%; height:520px; display:flex; align-items:center; justify-content:center; overflow:hidden;",

              shiny::plotOutput(ns("box_exp_preview")) |>

                shinycssloaders::withSpinner(type = 6, color = "#28a745")

            ),

            shiny::tags$small(style = "color: #666; display:block; margin-top:6px;",

              "Preview is auto-scaled to fit the window while keeping the export aspect ratio (tall figures appear as vertical thumbnails, wide figures as horizontal ones).")

          )

        )

      )

    }



    shiny::observeEvent(input$open_boxplot_image_export, {

      # Auto-extract data if not already cached (so the user doesn't need

      # to click "Export Boxplot Data" first)

      result <- current_boxplot_data()

      if (is.null(result) || is.null(result$long) || nrow(result$long) == 0) {

        data_list <- data_prep_data()

        current_gene <- current_boxplot_gene()

        if (is.null(data_list) || is.null(current_gene) || !nzchar(current_gene)) {

          shiny::showNotification(

            "No gene selected. Click a gene in the DE volcano or gene table first.",

            type = "warning", duration = 4); return()

        }

        symbol_map  <- .rebuild_symbol_map(gsea_res, data_list$contrast_id)

        display_gene <- .get_display_symbol(current_gene, symbol_map)

        result <- tryCatch(

          extract_boxplot_data(

            GSEAlens_res = gsea_res, contrast_id  = data_list$contrast_id,

            gene_symbol  = display_gene, expr_type = data_list$expression_type),

          error = function(e) {

            shiny::showNotification(

              sprintf("Failed to extract boxplot data: %s", e$message),

              type = "error", duration = 5); NULL })

        if (is.null(result)) return()

        current_boxplot_data(result)

      }

      shiny::showModal(.build_boxplot_export_modal())

    })

    shiny::observeEvent(input$box_exp_dismiss, shiny::removeModal())



    .boxplot_export_code <- function() {

      result <- current_boxplot_data()

      if (is.null(result)) return("")

      data_list <- data_prep_data()

      if (is.null(data_list)) return("")

      lg <- if (is.null(data_list$left_group))  "Left"  else data_list$left_group

      rg <- if (is.null(data_list$right_group)) "Right" else data_list$right_group

      generate_boxplot_image_code(

        box_data_long     = result$long,

        gene_symbol       = result$gene,

        expr_type         = data_list$expression_type,

        left_group        = lg,

        right_group       = rg,

        use_zero_baseline = isTRUE(input$zero_baseline)

      )

    }



    .boxplot_preview_plot <- shiny::reactive({

      code <- tryCatch(.boxplot_export_code(),

        error = function(e) {

          shiny::showNotification(

            sprintf("[boxplot preview] code generation failed: %s", e$message),

            type = "error", duration = 8)

          ""

        })

      if (!nzchar(code)) return(NULL)

      eval_env <- new.env(parent = globalenv())

      eval_env$p <- NULL

      eval_env$print <- base::invisible

      tryCatch(eval(parse(text = code), envir = eval_env),

               error = function(e) {

                 shiny::showNotification(

                   sprintf("[boxplot preview] %s", e$message),

                   type = "error", duration = 8)

                 NULL

               })

      eval_env$p

    })



    output$box_exp_preview <- shiny::renderPlot(

      {

        p <- .boxplot_preview_plot()

        shiny::req(p)

        p

      },

      res  = 100,

      width  = function() {

        d <- .vol_preview_dims(input$box_exp_width, input$box_exp_height); d$width

      },

      height = function() {

        d <- .vol_preview_dims(input$box_exp_width, input$box_exp_height); d$height

      }

    )



    .boxplot_render_to_file <- function(file, fmt) {

      code <- .boxplot_export_code()

      if (!nzchar(code)) return()

      eval_env <- new.env(parent = globalenv())

      eval_env$p <- NULL

      eval_env$print <- base::invisible

      tryCatch(eval(parse(text = code), envir = eval_env),

               error = function(e) shiny::showNotification(

                 sprintf("Render failed: %s", e$message), type = "error", duration = 8))

      if (is.null(eval_env$p)) return()

      ggplot2::ggsave(file, eval_env$p,

                      width  = if (is.null(input$box_exp_width))  5 else input$box_exp_width,

                      height = if (is.null(input$box_exp_height)) 5 else input$box_exp_height,

                      dpi    = if (is.null(input$box_exp_dpi))   300 else input$box_exp_dpi,

                      device = fmt)

    }



    output$box_exp_download_pdf <- shiny::downloadHandler(

      filename = function() sprintf("GSEAlens_boxplot_%s.pdf", Sys.Date()),

      content  = function(file) .boxplot_render_to_file(file, "pdf")

    )

    output$box_exp_download_png <- shiny::downloadHandler(

      filename = function() sprintf("GSEAlens_boxplot_%s.png", Sys.Date()),

      content  = function(file) .boxplot_render_to_file(file, "png")

    )

    output$box_exp_download_other <- shiny::downloadHandler(

      filename = function() sprintf("GSEAlens_boxplot_%s.%s", Sys.Date(),

                                    if (is.null(input$box_exp_format)) "svg" else input$box_exp_format),

      content  = function(file) .boxplot_render_to_file(file, input$box_exp_format)

    )



    shiny::observeEvent(input$box_exp_copy_code, {

      code <- .boxplot_export_code()

      if (!nzchar(code)) {

        shiny::showNotification("Nothing to copy.", type = "warning"); return()

      }

      ok <- tryCatch({ clipr::write_clip(code); TRUE }, error = function(e) FALSE)

      if (isTRUE(ok)) {

        shiny::showNotification("R code copied to clipboard.", type = "message")

      } else {

        shiny::showModal(shiny::modalDialog(

          title = "Reproducible R Code (copy manually)",

          size = "l", easyClose = TRUE,

          shiny::pre(style = "max-height: 60vh; overflow-y: auto; font-size: 11px;", code),

          footer = shiny::modalButton("Close")

        ))

      }

    })



    # ============================================================

    # Pathway Volcano Export Center

    # ------------------------------------------------------------

    # ggplot2 + ggsave (no kaleido/orca dependency), with a reproducible

    # R code generator (generate_volcano_code in 15_code_generator.R).

    # ============================================================



    .build_volcano_export_modal <- function() {

      shiny::modalDialog(

        title = "Export Publication Plot - Pathway Volcano",

        size = "l", footer = NULL, easyClose = TRUE,

        shiny::fluidRow(

          shiny::column(5,

            shiny::h5("Dimensions"),

            shiny::fluidRow(

              shiny::column(6,

                shiny::numericInput(ns("vol_exp_width"),  "Width (inch)",  value = 8, min = 2, max = 24)

              ),

              shiny::column(6,

                shiny::numericInput(ns("vol_exp_height"), "Height (inch)", value = 6, min = 2, max = 24)

              )

            ),

            shiny::numericInput(ns("vol_exp_dpi"), "DPI", value = 300, min = 72, max = 600),

            shiny::hr(),

            shiny::checkboxInput(ns("vol_exp_show_stats"),

              label = "Show stats subtitle (pathways / selected / significant)",

              value = FALSE),

            shiny::hr(),

            shiny::h5("Download"),

            shiny::fluidRow(

              shiny::column(6,

                shiny::downloadButton(ns("vol_exp_download_pdf"),

                  "Download PDF",

                  class = "btn-danger btn-block",

                  icon  = shiny::icon("file-pdf"))

              ),

              shiny::column(6,

                shiny::downloadButton(ns("vol_exp_download_png"),

                  "Download PNG",

                  class = "btn-primary btn-block",

                  icon  = shiny::icon("file-image"))

              )

            ),

            shiny::hr(),

            shiny::fluidRow(

              shiny::column(6,

                shiny::selectInput(ns("vol_exp_format"), "Other format",

                  choices = c("SVG" = "svg", "TIFF" = "tiff"),

                  selected = "svg")

              ),

              shiny::column(6,

                shiny::div(style = "margin-top: 22px;",

                  shiny::downloadButton(ns("vol_exp_download_other"),

                    "Download", class = "btn-default btn-block")

                )

              )

            ),

            shiny::hr(),

            shiny::fluidRow(

              shiny::column(6,

                shiny::actionButton(ns("vol_exp_copy_code"),

                  "Copy R Code",

                  class = "btn-info btn-block",

                  icon = shiny::icon("clipboard"))

              ),

              shiny::column(6,

                shiny::actionButton(ns("vol_exp_dismiss"),

                  "Close",

                  class = "btn-default btn-block")

              )

            )

          ),

          shiny::column(7,

            shiny::h5("Live Preview"),

            shiny::div(style = "background:#f5f5f5; border:1px solid #ddd; border-radius:4px; padding:8px; width:100%; height:520px; display:flex; align-items:center; justify-content:center; overflow:hidden;",

              shiny::plotOutput(ns("vol_exp_preview")) |>

                shinycssloaders::withSpinner(type = 6, color = "#28a745")

            ),

            shiny::tags$small(style = "color: #666; display:block; margin-top:6px;",

              "X axis: NES | Y axis: -log10(pvalue). Preview auto-scaled to fit while keeping the export aspect ratio.")

          )

        )

      )

    }



    volcano_export_df <- shiny::reactive({

      data_list <- data_prep_data()

      shiny::req(data_list)

      data_list$df[, c("ID", "NES", "pvalue", "p.adjust")]

    })



    shiny::observeEvent(input$open_volcano_export, {

      shiny::showModal(.build_volcano_export_modal())

    })

    shiny::observeEvent(input$vol_exp_dismiss, shiny::removeModal())



    # Live preview reactive

    # IRON FIX (2026-06-30): suppress `print(p)` inside the eval'd code so

    # it does not open a separate windows()/pdf() device and conflict with

    # renderPlot's own device. Also surface errors as notifications and add

    # NA guards to width/height/res (numericInput returns NA on cleared).

    #

    # Preview scaling (2026-07-06): the preview pane is a FIXED-SIZE window

    # (~500x480 px). The figure is rendered at the export aspect ratio and

    # scaled to fit inside that window, centered. A tall figure (e.g. 5x20

    # inch) appears as a vertical "pole" thumbnail; a wide figure (16x12)

    # as a horizontal "pole". This mirrors how cairo_pdf/PNG exports look.

    .vol_preview_dims <- function(w_in, h_in, max_w = 500, max_h = 480) {

      if (is.null(w_in) || is.na(w_in) || w_in <= 0) w_in <- 8

      if (is.null(h_in) || is.na(h_in) || h_in <= 0) h_in <- 6

      scale <- min(max_w / w_in, max_h / h_in)

      list(width  = round(w_in * scale),

           height = round(h_in * scale))

    }

    .volcano_preview_plot <- shiny::reactive({

      code <- tryCatch(.volcano_export_code(),

        error = function(e) {

          shiny::showNotification(

            sprintf("[vol preview] code generation failed: %s", e$message),

            type = "error", duration = 8)

          ""

        })

      if (!nzchar(code)) return(NULL)

      eval_env <- new.env(parent = globalenv())

      eval_env$p <- NULL

      eval_env$print <- base::invisible

      tryCatch(eval(parse(text = code), envir = eval_env),

               error = function(e) {

                 shiny::showNotification(

                   sprintf("[vol preview] %s", e$message),

                   type = "error", duration = 8)

                 NULL

               })

      eval_env$p

    })

    output$vol_exp_preview <- shiny::renderPlot(

      {

        p <- .volcano_preview_plot()

        shiny::req(p)

        p

      },

      res  = 100,

      width  = function() {

        d <- .vol_preview_dims(input$vol_exp_width, input$vol_exp_height); d$width

      },

      height = function() {

        d <- .vol_preview_dims(input$vol_exp_width, input$vol_exp_height); d$height

      }

    )



    .volcano_export_code <- function() {

      data_list <- data_prep_data()

      lg <- if (is.null(data_list$left_group))  "Left"  else data_list$left_group

      rg <- if (is.null(data_list$right_group)) "Right" else data_list$right_group

      sel_ids <- selected_pathway_ids()

      n_sel <- length(sel_ids)

      show_stats <- isTRUE(input$vol_exp_show_stats)

      generate_volcano_code(volcano_export_df(),

                            left_group = lg, right_group = rg,

                            n_selected = n_sel,

                            selected_ids = sel_ids,

                            show_annotations = show_stats)

    }



    .volcano_render_to_file <- function(file, fmt) {

      code <- .volcano_export_code()

      if (!nzchar(code)) return()

      # parent=globalenv() so lexical lookup resolves attached packages.

      eval_env <- new.env(parent = globalenv())

      eval_env$p <- NULL

      eval_env$print <- base::invisible   # suppress popup device; ggsave draws itself

      tryCatch(eval(parse(text = code), envir = eval_env),

               error = function(e) shiny::showNotification(

                 sprintf("Render failed: %s", e$message), type = "error", duration = 8))

      if (is.null(eval_env$p)) return()

      ggplot2::ggsave(file, eval_env$p,

                      width  = if (is.null(input$vol_exp_width))  8 else input$vol_exp_width,

                      height = if (is.null(input$vol_exp_height)) 6 else input$vol_exp_height,

                      dpi    = if (is.null(input$vol_exp_dpi))   300 else input$vol_exp_dpi,

                      device = fmt)

    }



    output$vol_exp_download_pdf <- shiny::downloadHandler(

      filename = function() sprintf("GSEAlens_volcano_%s.pdf", Sys.Date()),

      content  = function(file) .volcano_render_to_file(file, "pdf")

    )



    output$vol_exp_download_png <- shiny::downloadHandler(

      filename = function() sprintf("GSEAlens_volcano_%s.png", Sys.Date()),

      content  = function(file) .volcano_render_to_file(file, "png")

    )



    output$vol_exp_download_other <- shiny::downloadHandler(

      filename = function() sprintf("GSEAlens_volcano_%s.%s", Sys.Date(),

                                    if (is.null(input$vol_exp_format)) "svg" else input$vol_exp_format),

      content  = function(file) .volcano_render_to_file(file, input$vol_exp_format)

    )



    shiny::observeEvent(input$vol_exp_copy_code, {

      code <- .volcano_export_code()

      ok <- tryCatch({ clipr::write_clip(code); TRUE }, error = function(e) FALSE)

      if (isTRUE(ok)) {

        shiny::showNotification("R code copied to clipboard.", type = "message")

      } else {

        shiny::showModal(shiny::modalDialog(

          title = "Reproducible R Code (copy manually)",

          size = "l", easyClose = TRUE,

          shiny::pre(style = "max-height: 60vh; overflow-y: auto; font-size: 11px;", code),

          footer = shiny::modalButton("Close")

        ))

      }

    })



    # ============================================================

    # DE Volcano Export Center (gene-level)

    # ------------------------------------------------------------

    # Mirrors the Pathway Volcano export but uses generate_de_volcano_code

    # (six-category classification: both/user/pathway/up/down/ns).

    # ============================================================

    .build_de_volcano_export_modal <- function() {

      shiny::modalDialog(

        title = "Export Publication Plot - DE Volcano",

        size = "l", footer = NULL, easyClose = TRUE,

        shiny::fluidRow(

          shiny::column(5,

            shiny::h5("Dimensions"),

            shiny::fluidRow(

              shiny::column(6,

                shiny::numericInput(ns("de_exp_width"),  "Width (inch)",  value = 8, min = 2, max = 24)

              ),

              shiny::column(6,

                shiny::numericInput(ns("de_exp_height"), "Height (inch)", value = 6, min = 2, max = 24)

              )

            ),

            shiny::numericInput(ns("de_exp_dpi"), "DPI", value = 300, min = 72, max = 600),

            shiny::hr(),

            shiny::h5("Classification thresholds"),

            shiny::numericInput(ns("de_exp_logfc"), "logFC threshold",

                                value = if (is.null(input$volcano_logfc_thresh)) 1

                                        else input$volcano_logfc_thresh,

                                min = 0, max = 22, step = 0.5),

            shiny::numericInput(ns("de_exp_pval"), "P-value threshold",

                                value = if (is.null(input$volcano_pval_thresh)) 0.05

                                        else input$volcano_pval_thresh,

                                min = 0.001, max = 1, step = 0.01),

            shiny::hr(),

            shiny::checkboxInput(ns("de_exp_show_stats"),

              label = "Show stats + direction banners (Up/Down | High-in)",

              value = FALSE),

            shiny::hr(),

            shiny::h5("Download"),

            shiny::fluidRow(

              shiny::column(6,

                shiny::downloadButton(ns("de_exp_download_pdf"),

                  "Download PDF", class = "btn-danger btn-block",

                  icon = shiny::icon("file-pdf"))

              ),

              shiny::column(6,

                shiny::downloadButton(ns("de_exp_download_png"),

                  "Download PNG", class = "btn-primary btn-block",

                  icon = shiny::icon("file-image"))

              )

            ),

            shiny::hr(),

            shiny::fluidRow(

              shiny::column(6,

                shiny::selectInput(ns("de_exp_format"), "Other format",

                  choices = c("SVG" = "svg", "TIFF" = "tiff"), selected = "svg")

              ),

              shiny::column(6,

                shiny::div(style = "margin-top: 22px;",

                  shiny::downloadButton(ns("de_exp_download_other"),

                    "Download", class = "btn-default btn-block")

                )

              )

            ),

            shiny::hr(),

            shiny::fluidRow(

              shiny::column(6,

                shiny::actionButton(ns("de_exp_copy_code"),

                  "Copy R Code", class = "btn-info btn-block",

                  icon = shiny::icon("clipboard"))

              ),

              shiny::column(6,

                shiny::actionButton(ns("de_exp_dismiss"),

                  "Close", class = "btn-default btn-block")

              )

            )

          ),

          shiny::column(7,

            shiny::h5("Live Preview"),

            shiny::div(style = "background:#f5f5f5; border:1px solid #ddd; border-radius:4px; padding:8px; width:100%; height:520px; display:flex; align-items:center; justify-content:center; overflow:hidden;",

              shiny::plotOutput(ns("de_exp_preview")) |>

                shinycssloaders::withSpinner(type = 6, color = "#28a745")

            ),

            shiny::tags$small(style = "color: #666; display:block; margin-top:6px;",

              "Six categories: Both (purple) / Selected (green) / Pathway (orange) / Up (red) / Down (blue) / NS (gray). Preview auto-scaled to fit while keeping the export aspect ratio.")

          )

        )

      )

    }



    shiny::observeEvent(input$open_de_volcano_export, {

      data_list <- data_prep_data()

      de_df <- tryCatch(get_de_table(gsea_res, data_list$contrast_id),

                        error = function(e) NULL)

      if (is.null(de_df) || nrow(de_df) == 0) {

        shiny::showNotification("No DE data to export.", type = "warning"); return()

      }

      shiny::showModal(.build_de_volcano_export_modal())

    })

    shiny::observeEvent(input$de_exp_dismiss, shiny::removeModal())



    .de_volcano_export_df <- shiny::reactive({

      data_list <- data_prep_data()

      shiny::req(data_list)

      de_df <- tryCatch(get_de_table(gsea_res, data_list$contrast_id),

                        error = function(e) NULL)

      shiny::req(de_df)

      if (!"logFC" %in% colnames(de_df)) de_df$logFC <- de_df$log2FoldChange

      if (!"pvalue" %in% colnames(de_df)) de_df$pvalue <- de_df$p.value

      if (!"padj" %in% colnames(de_df)) de_df$padj <- de_df$p.adjust

      de_df <- de_df[!is.na(de_df$logFC) & !is.na(de_df$pvalue), ]

      de_df

    })



    .de_volcano_export_code <- function() {

      de_df <- .de_volcano_export_df()

      if (is.null(de_df) || nrow(de_df) == 0) return("")

      data_list <- data_prep_data()

      if (is.null(data_list)) return("")

      lg <- if (is.null(data_list$left_group))  "Left"  else data_list$left_group

      rg <- if (is.null(data_list$right_group)) "Right" else data_list$right_group

      user_genes <- toupper(highlight_genes_reactive())

      pathway_genes <- toupper(selected_pathway_genes())

      logfc_t <- if (is.null(input$de_exp_logfc)) 1 else input$de_exp_logfc

      pval_t  <- if (is.null(input$de_exp_pval))  0.05 else input$de_exp_pval

      show_stats <- isTRUE(input$de_exp_show_stats)

      generate_de_volcano_code(

        de_df         = de_df,

        user_genes    = user_genes,

        pathway_genes = pathway_genes,

        logfc_thresh  = logfc_t,

        pval_thresh   = pval_t,

        left_group    = lg, right_group = rg,

        show_annotations = show_stats

      )

    }



    .de_volcano_preview_plot <- shiny::reactive({

      code <- tryCatch(.de_volcano_export_code(),

        error = function(e) {

          shiny::showNotification(

            sprintf("[de vol preview] code generation failed: %s", e$message),

            type = "error", duration = 8)

          ""

        })

      if (!nzchar(code)) return(NULL)

      eval_env <- new.env(parent = globalenv())

      eval_env$p <- NULL

      eval_env$print <- base::invisible

      tryCatch(eval(parse(text = code), envir = eval_env),

               error = function(e) {

                 shiny::showNotification(

                   sprintf("[de vol preview] %s", e$message),

                   type = "error", duration = 8)

                 NULL

               })

      eval_env$p

    })



    output$de_exp_preview <- shiny::renderPlot(

      {

        p <- .de_volcano_preview_plot()

        shiny::req(p)

        p

      },

      res  = 100,

      width  = function() {

        d <- .vol_preview_dims(input$de_exp_width, input$de_exp_height); d$width

      },

      height = function() {

        d <- .vol_preview_dims(input$de_exp_width, input$de_exp_height); d$height

      }

    )



    .de_volcano_render_to_file <- function(file, fmt) {

      code <- .de_volcano_export_code()

      if (!nzchar(code)) return()

      eval_env <- new.env(parent = globalenv())

      eval_env$p <- NULL

      eval_env$print <- base::invisible

      tryCatch(eval(parse(text = code), envir = eval_env),

               error = function(e) shiny::showNotification(

                 sprintf("Render failed: %s", e$message), type = "error", duration = 8))

      if (is.null(eval_env$p)) return()

      ggplot2::ggsave(file, eval_env$p,

                      width  = if (is.null(input$de_exp_width))  8 else input$de_exp_width,

                      height = if (is.null(input$de_exp_height)) 6 else input$de_exp_height,

                      dpi    = if (is.null(input$de_exp_dpi))   300 else input$de_exp_dpi,

                      device = fmt)

    }



    output$de_exp_download_pdf <- shiny::downloadHandler(

      filename = function() sprintf("GSEAlens_de_volcano_%s.pdf", Sys.Date()),

      content  = function(file) .de_volcano_render_to_file(file, "pdf")

    )

    output$de_exp_download_png <- shiny::downloadHandler(

      filename = function() sprintf("GSEAlens_de_volcano_%s.png", Sys.Date()),

      content  = function(file) .de_volcano_render_to_file(file, "png")

    )

    output$de_exp_download_other <- shiny::downloadHandler(

      filename = function() sprintf("GSEAlens_de_volcano_%s.%s", Sys.Date(),

                                    if (is.null(input$de_exp_format)) "svg" else input$de_exp_format),

      content  = function(file) .de_volcano_render_to_file(file, input$de_exp_format)

    )



    shiny::observeEvent(input$de_exp_copy_code, {

      code <- .de_volcano_export_code()

      if (!nzchar(code)) {

        shiny::showNotification("Nothing to copy.", type = "warning"); return()

      }

      ok <- tryCatch({ clipr::write_clip(code); TRUE }, error = function(e) FALSE)

      if (isTRUE(ok)) {

        shiny::showNotification("R code copied to clipboard.", type = "message")

      } else {

        shiny::showModal(shiny::modalDialog(

          title = "Reproducible R Code (copy manually)",

          size = "l", easyClose = TRUE,

          shiny::pre(style = "max-height: 60vh; overflow-y: auto; font-size: 11px;", code),

          footer = shiny::modalButton("Close")

        ))

      }

    })

  })

}

