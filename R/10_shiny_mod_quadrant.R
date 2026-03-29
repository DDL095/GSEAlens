#' @title Quadrant Linkage Module

#' @title Quadrant Linkage UI
#' @description UI components for the four-quadrant interactive visualization
#' @param id Module ID
#' @return Shiny UI tagList
#' @keywords internal
mod_quadrant_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      shiny::column(6, shiny::div(class = "white-box",
                                  shiny::h4("1. Pathway Volcano"),
                                  shiny::div(
                                    style = "position: relative;",
                                    plotly::plotlyOutput(ns("volcano_pathway"), height = "400px"),
                                    shiny::div(
                                      style = "position: absolute; top: 10px; left: 10px; background: rgba(255,255,255,0.9); padding: 5px 10px; border-radius: 4px; font-size: 11px; color: #666;",
                                      "Click pathway to highlight | Check 'Joint Plot' in table"
                                    )
                                  ),
                                  shiny::uiOutput(ns("selected_pathways_display"))
      )),
      shiny::column(6, shiny::div(class = "white-box",
                                  shiny::h4("2. Gene Rank Distribution"),
                                  plotly::plotlyOutput(ns("volcano_gene"), height = "400px")))
    ),

    shiny::fluidRow(
      shiny::column(6, shiny::div(class = "white-box",
                                  shiny::h4("3. Differential Expression Volcano"),
                                  shiny::div(
                                    style = "margin-bottom: 10px;",
                                    shiny::checkboxInput(ns("zero_baseline"), label = "Use 0 as baseline", value = FALSE),
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
                                  plotly::plotlyOutput(ns("de_volcano"), height = "400px"))),
      shiny::column(6, shiny::div(class = "white-box",
                                  shiny::h4("4. Full Expression Distribution"),
                                  plotly::plotlyOutput(ns("gene_expr_box"), height = "400px"),
                                  shiny::uiOutput(ns("boxplot_order_status"))))
    ),

    shiny::hr(),

    shiny::fluidRow(
      shiny::column(12, shiny::div(class = "white-box",
                                   shiny::h4("Gene Expression Table (Click View to display boxplot | Click X to remove)"),
                                   shiny::div(
                                     style = "margin-bottom: 10px;",
                                     shiny::actionButton(ns("clear_all_genes_btn_quadrant"), label = "Clear All", class = "btn-warning btn-sm"),
                                     shiny::helpText("Click 'Clear All' to remove all genes from the table")
                                   ),
                                   shiny::div(
                                     style = "margin-bottom: 10px; font-size: 12px; color: #666;",
                                     shiny::HTML("<b>Gene Sources:</b> Interest genes from sidebar | DE volcano clicks | Pathway genes (orange)")
                                   ),
                                   DT::dataTableOutput(ns("gene_expr_table"))
      ))
    ),

    shiny::hr(),
    shiny::fluidRow(
      shiny::column(12, shiny::div(class = "white-box",
                                   shiny::h4("Export Code"),
                                   shiny::actionButton(ns("export_code_btn"), label = "Export Current Plot Code", class = "btn-secondary", style = "width: 100%;"),
                                   shiny::helpText("Generate R code for the currently selected pathways")
      ))
    )
  )
}


#' @title Quadrant Linkage Server
#' @description Server logic with unified gene pool
#' @param id Module ID
#' @param data_prep_list List from data prep module
#' @param gsea_res GseaRes object
#' @keywords internal
mod_quadrant_server <- function(id, data_prep_list, gsea_res) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Color Definitions
    COLOR_LEFT <- "#E41A1C"      # 红色 - 上调
    COLOR_RIGHT <- "#377EB8"     # 蓝色 - 下调
    COLOR_NS <- "#C0C0C0"        # 灰色 - 不显著
    COLOR_USER <- "#4DAF4A"      # 绿色 - 仅用户基因
    COLOR_PATHWAY <- "#FF9800"   # 橙色 - 仅通路基因
    COLOR_BOTH <- "#9C27B0"      # 紫色 - 用户 ∩ 通路

    # Reactive Values
    selected_pathway_ids <- shiny::reactiveVal(character(0))
    selected_pathway_genes <- shiny::reactiveVal(character(0))
    current_boxplot_gene <- shiny::reactiveVal(NULL)

    # Get Reactive References
    data_prep_data <- data_prep_list$data
    highlight_genes_reactive <- data_prep_list$highlight_genes
    boxplot_order_ref <- data_prep_list$boxplot_order

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
        message("[Pathway] Cleared all selected pathways")
      } else {
        current <- selected_pathway_ids()
        new_selection <- setdiff(current, pathway_to_remove)
        selected_pathway_ids(new_selection)

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
                         ifelse(df$NES > 0, COLOR_LEFT, COLOR_RIGHT))
      df$size <- ifelse(df$ID %in% current_selections, 18, 10)
      df$opacity <- ifelse(length(current_selections) == 0, 0.8,
                           ifelse(df$ID %in% current_selections, 1.0, 0.35))
      df$linewidth <- ifelse(df$ID %in% current_selections, 3, 1)

      # Pathway labels
      annotations_list <- list()
      if (length(current_selections) > 0) {
        selected_df <- df[df$ID %in% current_selections, ]
        for (i in 1:nrow(selected_df)) {
          row <- selected_df[i, ]
          ax_offset <- ifelse(i %% 2 == 1, 0, 50)
          annotations_list[[i]] <- list(
            x = row$NES,
            y = -log10(row$p.adjust),
            text = row$ID,
            showarrow = TRUE,
            arrowhead = 2,
            arrowsize = 1,
            arrowwidth = 2,
            arrowcolor = COLOR_PATHWAY,
            ax = ax_offset,
            ay = -30,
            font = list(size = 10, color = COLOR_PATHWAY),
            bgcolor = "rgba(255,255,255,0.95)",
            bordercolor = COLOR_PATHWAY,
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
            text = sprintf("Pathway Volcano: %s vs %s<br><sub>%d pathways | %d selected | %d significant (FDR<0.25)</sub>",
                           data_list$left_group, data_list$right_group,
                           nrow(df), length(current_selections), sum(df$p.adjust < 0.25, na.rm = TRUE)),
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

    # Pathway Volcano Click Event
    shiny::observeEvent(plotly::event_data("plotly_click", source = ns("pathway_volcano")), {
      click <- plotly::event_data("plotly_click", source = ns("pathway_volcano"))
      if (is.null(click) || is.null(click$key)) return()

      clicked_id <- click$key
      current <- selected_pathway_ids()
      data_list <- data_prep_data()

      if (clicked_id %in% current) {
        new_selection <- setdiff(current, clicked_id)
        selected_pathway_ids(new_selection)

        if (length(new_selection) > 0) {
          last_id <- tail(new_selection, 1)
          if (!is.null(data_list) && !is.null(data_list$gsea_res@geneSets[[last_id]])) {
            pathway_genes <- data_list$gsea_res@geneSets[[last_id]]
            selected_pathway_genes(toupper(pathway_genes))
          }
        } else {
          selected_pathway_genes(character(0))
        }
      } else {
        new_selection <- c(current, clicked_id)
        selected_pathway_ids(new_selection)

        if (!is.null(data_list) && !is.null(data_list$gsea_res@geneSets[[clicked_id]])) {
          pathway_genes <- data_list$gsea_res@geneSets[[clicked_id]]
          selected_pathway_genes(toupper(pathway_genes))
        }
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
      ) %>%
        plotly::layout(
          xaxis = list(title = "Gene Rank"),
          yaxis = list(title = "Ranking Metric (Stat)"),
          showlegend = FALSE,
          title = list(
            text = ifelse(length(pathway_genes) > 0,
                          sprintf("Selected: %d pathway genes highlighted", length(pathway_genes)),
                          "Click pathway in volcano above to mark"),
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

      de_df <- tryCatch({
        get_de_table(gsea_res, contrast_id)
      }, error = function(e) NULL)

      shiny::req(de_df)

      if (!"logFC" %in% colnames(de_df)) de_df$logFC <- de_df$log2FoldChange
      if (!"pvalue" %in% colnames(de_df)) de_df$pvalue <- de_df$p.value
      if (!"padj" %in% colnames(de_df)) de_df$padj <- de_df$p.adjust

      de_df <- de_df[!is.na(de_df$logFC) & !is.na(de_df$pvalue), ]
      if (nrow(de_df) == 0) return(NULL)

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
      n_both <- sum(de_df$is_user & de_df$is_pathway, na.rm = TRUE)

      # Six-color logic
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

      title_text <- sprintf(
        "%s vs %s<br><sup>Up %d | Down %d | NS %d | User %d | Both %d | logFC>|%.1f|, p<%.3f</sup>",
        left_group, right_group, n_up, n_down, n_not_sig, n_user, n_both, logfc_thresh, pval_thresh
      )

      annotations_list <- list()
      max_y_val <- max(de_df$y_axis, na.rm = TRUE)

      # High in labels
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

      # User gene labels
      user_genes_df <- de_df[de_df$is_user, ]
      if (nrow(user_genes_df) > 0) {
        for (i in 1:nrow(user_genes_df)) {
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

      # Threshold line labels
      if (logfc_thresh > 0) {
        annotations_list[[length(annotations_list) + 1]] <- list(
          x = logfc_thresh, y = max_y_val * 0.95,
          text = sprintf("logFC=%.1f", logfc_thresh),
          showarrow = FALSE, font = list(size = 10, color = "gray")
        )
        annotations_list[[length(annotations_list) + 1]] <- list(
          x = -logfc_thresh, y = max_y_val * 0.95,
          text = sprintf("-logFC=%.1f", logfc_thresh),
          showarrow = FALSE, font = list(size = 10, color = "gray")
        )
      }

      p <- plotly::plot_ly(
        data = de_df,
        x = de_df$x_axis,
        y = de_df$y_axis,
        type = "scatter",
        mode = "markers",
        marker = list(
          color = de_df$color,
          size = de_df$size,
          opacity = de_df$opacity,
          line = list(color = "white", width = de_df$linewidth)
        ),
        text = ~sprintf("%s<br>logFC: %.2f<br>-log10(p): %.2f<br>FDR: %.2e",
                        gene_symbol, logFC, y_axis, padj),
        hoverinfo = "text",
        key = de_df$gene_upper,
        source = ns("deg_volcano"),
        showlegend = FALSE
      ) %>%
        plotly::layout(
          title = list(text = title_text, font = list(size = 14), x = 0.5, xanchor = "center"),
          xaxis = list(title = "logFC", zeroline = FALSE, showgrid = TRUE, gridcolor = "lightgray"),
          yaxis = list(title = "-log10 (P-value)", zeroline = FALSE, showgrid = TRUE, gridcolor = "lightgray"),
          showlegend = FALSE,
          dragmode = "pan",
          annotations = annotations_list,
          shapes = list(
            list(type = "line", x0 = logfc_thresh, x1 = logfc_thresh,
                 y0 = 0, y1 = max_y_val * 1.05,
                 line = list(color = "gray", dash = "dash", width = 1)),
            list(type = "line", x0 = -logfc_thresh, x1 = -logfc_thresh,
                 y0 = 0, y1 = max_y_val * 1.05,
                 line = list(color = "gray", dash = "dash", width = 1)),
            list(type = "line", x0 = min(de_df$x_axis) * 1.1, x1 = max(de_df$x_axis) * 1.1,
                 y0 = -log10(pval_thresh), y1 = -log10(pval_thresh),
                 line = list(color = "gray", dash = "dash", width = 1))
          )
        )

      p
    })

    # DE Volcano Click Event
    shiny::observeEvent(plotly::event_data("plotly_click", source = ns("deg_volcano")), {
      click <- plotly::event_data("plotly_click", source = ns("deg_volcano"))
      if (is.null(click) || is.null(click$key)) return()

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

      left_grp <- data_list$left_group
      right_grp <- data_list$right_group

      user_genes <- toupper(highlight_genes_reactive())
      pathway_genes <- selected_pathway_genes()

      de_df <- tryCatch({
        get_de_table(gsea_res, data_list$contrast_id)
      }, error = function(e) NULL)

      table_data <- lapply(genes, function(g) {
        gene_upper <- toupper(g)
        is_user <- gene_upper %in% user_genes
        is_pathway <- gene_upper %in% pathway_genes

        color <- if (is_user && is_pathway) COLOR_BOTH else if (is_user) COLOR_USER else COLOR_PATHWAY

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

        pval_display <- ifelse(is.na(pval), "-", if(pval < 0.001) sprintf("%.2e", pval) else sprintf("%.4f", pval))
        padj_display <- ifelse(is.na(padj_val), "-", if(padj_val < 0.001) sprintf("%.2e", padj_val) else sprintf("%.4f", padj_val))

        delete_btn <- sprintf(
          '<button class="btn btn-xs btn-danger" onclick="Shiny.setInputValue(\'%s\', \'%s\', {priority: \'event\'})">X</button>',
          ns("delete_gene_from_table"), g
        )
        view_btn <- sprintf(
          '<button class="btn btn-xs btn-primary" onclick="Shiny.setInputValue(\'%s\', \'%s\', {priority: \'event\'})">View</button>',
          ns("view_boxplot_from_table"), g
        )

        data.frame(
          Gene = g,
          Color = color,
          Log2FC = ifelse(is.na(logfc), "-", sprintf("%.3f", logfc)),
          LinearFC = ifelse(is.na(logfc), "-", sprintf("%.2f", 2^logfc)),
          Pvalue = pval_display,
          Padjust = padj_display,
          UpDown = updown,
          HighIn = high_in,
          IsPathway = ifelse(is_pathway, "Yes", "No"),
          Delete = delete_btn,
          ViewBox = view_btn,
          pval_raw = ifelse(is.na(pval), 1, pval),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      })

      table_df <- do.call(rbind, table_data)

      names(table_df) <- c("Gene", "Color", "Log2FC", "Linear FC", "P-value", "Adj P-value",
                           "Up/Down", "High In", "In Pathway", "Remove", "Boxplot", "pval_sort")

      dt <- DT::datatable(
        table_df,
        escape = FALSE,
        rownames = FALSE,
        selection = "none",
        extensions = c('Scroller'),
        options = list(
          pageLength = 10,
          scrollY = "50vh",
          scroller = TRUE,
          dom = "frtip",
          search = list(caseInsensitive = TRUE),
          columnDefs = list(
            list(visible = FALSE, targets = c(1, 11)),
            list(orderable = FALSE, targets = c(9, 10)),
            list(className = "dt-center", targets = c(2, 3, 4, 5, 6, 7, 8))
          ),
          order = list(list(11, "asc"))
        )
      )

      dt <- dt %>% DT::formatStyle(
        columns = "Up/Down",
        backgroundColor = DT::styleEqual(c("UP", "DOWN", "-"), c("#FFCDD2", "#BBDEFB", "transparent")),
        color = DT::styleEqual(c("UP", "DOWN", "-"), c(COLOR_LEFT, COLOR_RIGHT, "#666")),
        fontWeight = DT::styleEqual(c("UP", "DOWN"), c("bold", "bold"))
      )

      dt <- dt %>% DT::formatStyle(
        columns = "P-value",
        color = DT::styleInterval(0.001, c("red", "black")),
        fontWeight = DT::styleInterval(0.05, c("bold", "normal"))
      )

      dt <- dt %>% DT::formatStyle(
        columns = "Adj P-value",
        color = DT::styleInterval(0.001, c("red", "black")),
        fontWeight = DT::styleInterval(0.05, c("bold", "normal"))
      )

      dt
    })

    # Delete Gene from Table
    shiny::observeEvent(input$delete_gene_from_table, {
      gene_to_remove <- input$delete_gene_from_table
      if (is.null(gene_to_remove)) return()

      current_applied <- highlight_genes_reactive()
      new_applied <- setdiff(toupper(current_applied), toupper(gene_to_remove))
      highlight_genes_reactive(new_applied)

      if (!is.null(current_boxplot_gene()) && toupper(current_boxplot_gene()) == toupper(gene_to_remove)) {
        current_boxplot_gene(NULL)
      }
    })

    # View Boxplot from Table
    shiny::observeEvent(input$view_boxplot_from_table, {
      gene_name <- input$view_boxplot_from_table
      if (is.null(gene_name)) return()
      current_boxplot_gene(gene_name)
    })

    # 5. Expression Boxplot
    output$gene_expr_box <- plotly::renderPlotly({
      current_gene <- current_boxplot_gene()

      if (is.null(current_gene)) {
        return(plotly::plot_ly() %>% plotly::layout(
          title = list(text = "Click gene in DE volcano or table to view expression", font = list(size = 14))
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
          return(plotly::plot_ly() %>% plotly::layout(
            title = sprintf("Gene '%s' not found", current_gene)
          ))
        }

        actual_gene <- rownames(expr_mat)[match_idx[1]]
        display_gene_name <- current_gene

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
          group_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                            "#FF7F00", "#A65628", "#F781BF", "#999999")[1:length(unique_groups)]
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
          ggplot2::geom_jitter(width = 0.2, size = 3, alpha = 0.6,
                               ggplot2::aes(text = sprintf(
                                 "<b>Sample:</b> %s<br><b>Group:</b> %s<br><b>Expression:</b> %.3f",
                                 Sample, Group, Expression
                               ))) +
          ggplot2::scale_fill_manual(values = group_colors) +
          ggplot2::scale_x_discrete(limits = x_categories, drop = FALSE) +
          ggplot2::coord_cartesian(ylim = c(y_min, y_max)) +
          ggplot2::theme_bw(base_size = 12) +
          ggplot2::labs(title = sprintf("Expression: %s", display_gene_name),
                        y = data_list$expression_type, x = NULL) +
          ggplot2::theme(legend.position = "none",
                         axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

        if (use_zero_baseline) {
          p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                                       color = "red", alpha = 0.7, size = 0.8)
        }

        ply <- plotly::ggplotly(p, tooltip = "text")
        ply %>% plotly::layout(
          xaxis = list(categoryorder = "array", categoryarray = x_categories, title = ""),
          dragmode = FALSE
        )

      }, error = function(e) {
        return(plotly::plot_ly() %>% plotly::layout(
          title = sprintf("Error: %s", e$message)
        ))
      })
    })

    # Boxplot Order Status
    output$boxplot_order_status <- shiny::renderUI({
      order_info <- boxplot_order_ref()
      if (is.null(order_info) || order_info == "default") {
        shiny::div(style = "font-size: 12px; color: #666;", "Order: Default")
      } else {
        shiny::div(style = "font-size: 12px; color: #666;",
                   "Order:", gsub(",", " -> ", order_info))
      }
    })

    # 6. Export Code Button
    shiny::observeEvent(input$export_code_btn, {
      data_list <- data_prep_data()
      shiny::req(data_list)

      target_pw <- if (nrow(data_list$df) > 0) {
        data_list$df$ID[1:min(5, nrow(data_list$df))]
      } else {
        character(0)
      }

      shiny::showModal(shiny::modalDialog(
        title = "Generated R Code",
        size = "l",
        easyClose = TRUE,
        shiny::fluidRow(
          shiny::column(12,
                        shiny::div(
                          style = "background: #f5f5f5; padding: 15px; border-radius: 5px; max-height: 500px; overflow: auto;",
                          shiny::tags$pre(
                            shiny::code(
                              generate_pathway_plot_code(
                                GSEAlens_res = gsea_res,
                                contrast_id = data_list$contrast_id,
                                target_pathways = target_pw,
                                user_genes = highlight_genes_reactive(),
                                pathway_genes = selected_pathway_genes(),
                                expr_type = data_list$expression_type
                              )
                            ),
                            style = "font-size: 11px; white-space: pre-wrap; word-break: break-all;"
                          )
                        )
          )
        ),
        footer = shiny::modalButton("Close")
      ))
    })

  })
}
