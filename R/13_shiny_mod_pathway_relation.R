#' @title Pathway Relationship Exploration Module UI

#' @description Network visualization module with dual-mode pathway selection

#' @keywords internal



mod_pathway_relation_ui <- function(id) {

  ns <- shiny::NS(id)

  shiny::tagList(

    shiny::fluidRow(

      # ===== Left Control Panel =====

      shiny::column(

        2,

        shiny::div(

          class = "well",

          style = "padding: 15px;",



          # ---- Mode Selector ----

          shiny::h4("Select Analysis Mode"),

          shiny::radioButtons(

            ns("network_mode"),

            label = NULL,

            choices = c(

              "Top N from Current Set" = "mode_topN",

              "Selected from Main Table" = "mode_select"

            ),

            selected = "mode_topN",

            width = "100%"

          ),

          shiny::hr(),



          # ---- Shared Parameters Panel ----

          shiny::h4("Shared Parameters"),

          shiny::numericInput(

            ns("fdr_threshold"),

            label = "FDR Threshold:",

            min = 0,

            max = 1.0,

            value = 0.25,

            step = 0.01

          ),

          shiny::numericInput(

            ns("min_shared"),

            label = "Min Shared Core Genes:",

            value = 3,

            min = 1,

            max = 9999,

            step = 1

          ),

          shiny::sliderInput(

            ns("max_nodes"),

            label = "Max Nodes:",

            min = 5,

            max = 100,

            value = 30,

            step = 5

          ),

          shiny::selectInput(

            ns("network_layout"),

            label = "Layout Algorithm:",

            choices = c(

              "Fruchterman-Reingold" = "fr",

              "Kamada-Kawai" = "kk",

              "Circle" = "circle"

            ),

            selected = "fr"

          ),

          shiny::numericInput(

            ns("seed"),

            label = "Seed (for reproducibility):",

            value = 42,

            min = 1,

            max = 9999,

            step = 1

          ),

          shiny::helpText(

            style = "color: #666; font-size: 11px;",

            "Seed ensures reproducible layout when parameters change"

          ),

          shiny::hr(),



          # ---- Hover Display Settings ----

          shiny::h4("Hover Display Settings"),

          shiny::sliderInput(

            ns("hover_max_genes"),

            label = "Max Genes in Hover:",

            min = 3,

            max = 30,

            value = 5,

            step = 1

          ),

          shiny::helpText(

            style = "color: #666; font-size: 11px;",

            "Number of shared genes to display in edge hover tooltip"

          ),

          shiny::hr(),

          shiny::conditionalPanel(

            condition = sprintf("input['%s'] == 'mode_topN'", ns("network_mode")),

            shiny::div(

              style = "background: #e3f2fd; padding: 10px; border-radius: 5px;",

              shiny::h5("Top N Configuration"),

              shiny::numericInput(

                ns("topN_count"),

                label = "Top N Count:",

                value = 20,

                min = 1,

                max = 9999,

                step = 1

              )

            )

          ),

          shiny::conditionalPanel(

            condition = sprintf("input['%s'] == 'mode_select'", ns("network_mode")),

            shiny::div(

              style = "background: #f3e5f5; padding: 10px; border-radius: 5px;",

              shiny::h5("Main Table Selection"),

              shiny::helpText(

                style = "color: #6a1b9a; font-size: 11px;",

                "Check pathways in the main table 'Joint Plot' column"

              )

            )

          ),

          shiny::hr(),

          shiny::h4("Pathways to Plot"),

          shiny::uiOutput(ns("pathway_preview_list"))

        )

      ),



      # ===== Right Panel =====

      shiny::column(

        10,

        shiny::tabsetPanel(

          id = ns("active_tab"),

          type = "tabs",

          shiny::tabPanel(

            title = "DotPlot",

            value = "dotplot",

            shiny::div(

              class = "white-box",

              style = "padding: 15px; margin-top: 15px;",

              shiny::uiOutput(ns("dotplot_status")),

              shiny::fluidRow(

                shiny::column(4,

                  shiny::selectInput(

                    ns("dotplot_color_mode"),

                    label = "Color by:",

                    choices = c(

                      "-log10(FDR)" = "padj",

                      "-log10(P-value)" = "pval",

                      "NES" = "nes"

                    ),

                    selected = "padj"

                  )

                ),

                shiny::column(4,

                  shiny::selectInput(

                    ns("dotplot_size_mode"),

                    label = "Size by:",

                    choices = c(

                      "Core Genes Count" = "core_size",

                      "Set Size" = "setsize"

                    ),

                    selected = "core_size"

                  )

                ),

                shiny::column(4,

                  shiny::div(style = "margin-top: 22px;",

                    shiny::actionButton(

                      ns("open_dotplot_export"),

                      label = "Export Publication Plot",

                      class = "btn-success",

                      style = "width: 100%;"

                    )

                  )

                )

              ),

              plotly::plotlyOutput(ns("plot_dotplot"), height = "800px") |>

                shinycssloaders::withSpinner(type = 6, color = "#28a745")

            )

          ),

          shiny::tabPanel(

            title = "Network",

            value = "network",

            shiny::div(

              class = "white-box",

              style = "padding: 15px; margin-top: 15px;",

              shiny::uiOutput(ns("network_status")),

              shiny::fluidRow(

                shiny::column(4,

                  shiny::selectInput(

                    ns("network_edge_width_mode"),

                    label = "Edge Width Mode:",

                    choices = c(

                      "Weight-based (emapplot style)" = "weight",

                      "Rank-based (uniform spacing)"  = "rank"

                    ),

                    selected = "weight"

                  )

                ),

                shiny::column(4,

                  shiny::div(style = "margin-top: 22px;",

                    shiny::actionButton(

                      ns("open_network_export"),

                      label = "Export Publication Plot",

                      class = "btn-success",

                      style = "width: 100%;"

                    )

                  )

                ),

                shiny::column(4,

                  shiny::helpText(

                    style = "margin-top: 22px; color: #666; font-size: 11px;",

                    "Weight-based: faithful to Jaccard values",

                    shiny::br(),

                    "Rank-based: uniform visual spacing"

                  )

                )

              ),

              shiny::div(

                id = ns("selection_panel"),

                style = "background: #e8f4fd; padding: 12px; border-radius: 8px; margin-bottom: 15px;",

                shiny::fluidRow(

                  shiny::column(

                    12,

                    shiny::uiOutput(ns("selection_display"))

                  )

                ),

                shiny::fluidRow(

                  shiny::column(

                    6,

                    shiny::actionButton(

                      ns("show_edge_detail"),

                      label = "Show Edge Detail",

                      class = "btn-primary",

                      style = "width: 100%;",

                      disabled = NA

                    )

                  ),

                  shiny::column(

                    6,

                    shiny::actionButton(

                      ns("clear_selection"),

                      label = "Clear Selection",

                      class = "btn-warning",

                      style = "width: 100%;"

                    )

                  )

                ),

                shiny::div(

                  style = "margin-top: 8px; font-size: 11px; color: #666;",

                  shiny::HTML("<b>Instruction:</b> Click two nodes to select, then click 'Show Edge Detail'")

                )

              ),

              plotly::plotlyOutput(ns("plot_network"), height = "1200px") |>

                shinycssloaders::withSpinner(type = 6, color = "#28a745")

            )

          )

        )

      )

    )

  )

}





#' @title Pathway Relationship Exploration Module Server

#' @description Dual-mode pathway network visualization with defensive rendering

#' @keywords internal



mod_pathway_relation_server <- function(id, data_prep_list, gsea_res, table_result = NULL) {

  shiny::moduleServer(id, function(input, output, session) {

    ns <- session$ns



    # ============================================================

    # 1. Mode State Management

    # ============================================================



    network_mode <- shiny::reactiveVal("mode_topN")



    shiny::observeEvent(input$network_mode, {

      new_mode <- input$network_mode

      if (new_mode != network_mode()) {

        network_mode(new_mode)

      }

    })



    # ============================================================

    # 2. Data Source Reactive

    # ============================================================



    topN_candidates <- shiny::reactive({

      if (network_mode() != "mode_topN") {

        return(character(0))

      }

      data_list <- data_prep_list$data()

      shiny::req(data_list)

      df <- data_list$df

      shiny::req(nrow(df) > 0)

      top_n <- input$topN_count

      if (is.null(top_n)) top_n <- 50

      top_n <- max(1, min(top_n, nrow(df)))

      df[seq_len(top_n), "ID"]

    })



    select_candidates <- shiny::reactive({

      if (network_mode() != "mode_select") {

        return(character(0))

      }

      if (is.null(table_result) || is.null(table_result$selected_pathways)) {

        return(character(0))

      }

      selected <- table_result$selected_pathways()

      if (is.null(selected)) {

        return(character(0))

      }

      return(selected)

    })



    candidate_raw <- shiny::reactive({

      switch(network_mode(),

        "mode_topN" = topN_candidates(),

        "mode_select" = select_candidates(),

        character(0)

      )

    })



    candidate_filtered <- shiny::reactive({

      pathways <- candidate_raw()

      if (length(pathways) == 0) {

        return(character(0))

      }

      fdr_thresh <- input$fdr_threshold

      if (is.null(fdr_thresh)) fdr_thresh <- 0.25

      data_list <- data_prep_list$data()

      shiny::req(data_list)

      df <- data_list$df

      fdr_vec <- df$p.adjust[match(pathways, df$ID)]

      names(fdr_vec) <- pathways

      filtered <- pathways[!is.na(fdr_vec) & fdr_vec < fdr_thresh]

      return(filtered)

    })



    # ============================================================

    # 3. Final Pathway List Management

    # ============================================================



    final_pathways <- shiny::reactiveVal(character(0))



    shiny::observeEvent(candidate_filtered(), {

      new_candidates <- candidate_filtered()

      if (length(new_candidates) > 0) {

        max_n <- input$max_nodes

        if (is.null(max_n)) max_n <- 999

        final <- new_candidates[seq_len(min(length(new_candidates), max_n))]

        final_pathways(final)

      } else {

        final_pathways(character(0))

      }

    })



    # ============================================================

    # 4. Pathway Preview List UI

    # ============================================================



    output$pathway_preview_list <- shiny::renderUI({

      pathways <- final_pathways()

      if (length(pathways) == 0) {

        return(shiny::div(

          style = "background: #fff3cd; padding: 10px; border-radius: 5px; color: #856404;",

          shiny::strong("No pathways available"),

          shiny::br(),

          htmltools::tags$small("Adjust parameters or select pathways in Main Table")

        ))

      }

      data_list <- data_prep_list$data()

      df <- data_list$df

      pathway_info <- lapply(pathways, function(pid) {

        row_idx <- which(df$ID == pid)

        if (length(row_idx) == 0) {

          return(NULL)

        }

        row <- df[row_idx[1], ]

        fdr <- row$p.adjust

        nes <- row$NES

        fdr_str <- if (!is.na(fdr)) sprintf("%.2e", fdr) else "N/A"

        nes_str <- if (!is.na(nes)) sprintf("%.2f", nes) else "N/A"

        shiny::div(

          style = "background: #f8f9fa; padding: 8px; margin-bottom: 5px; border-radius: 4px; border-left: 3px solid #007bff;",

          shiny::div(style = "font-size: 12px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap;", shiny::strong(pid)),

          shiny::div(style = "font-size: 11px; color: #666;", sprintf("FDR: %s | NES: %s", fdr_str, nes_str))

        )

      })

      shiny::tagList(pathway_info)

    })



    # ============================================================

    # 5. DotPlot Status Message

    # ============================================================



    output$dotplot_status <- shiny::renderUI({

      pathways <- final_pathways()

      status_parts <- c(

        sprintf("Mode: %s", ifelse(network_mode() == "mode_topN", "TopN", "Select")),

        sprintf("Pathways: %d", length(pathways)),

        sprintf("FDR < %.2f", input$fdr_threshold)

      )

      if (length(pathways) == 0) {

        div_style <- "background: #fff3cd; padding: 10px; border-radius: 5px; margin-bottom: 10px; color: #856404;"

      } else {

        div_style <- "background: #d4edda; padding: 10px; border-radius: 5px; margin-bottom: 10px; color: #155724;"

      }

      shiny::div(style = div_style, shiny::HTML(paste(status_parts, collapse = " | ")))

    })



    # ============================================================

    # 6. Network Status Message

    # ============================================================



    output$network_status <- shiny::renderUI({

      pathways <- final_pathways()

      status_parts <- c(

        sprintf("Mode: %s", ifelse(network_mode() == "mode_topN", "TopN", "Select")),

        sprintf("Pathways: %d", length(pathways)),

        sprintf("FDR < %.2f", input$fdr_threshold),

        sprintf("min_shared: %d", input$min_shared),

        sprintf("seed: %d", input$seed)

      )

      if (length(pathways) == 0) {

        div_style <- "background: #fff3cd; padding: 10px; border-radius: 5px; margin-bottom: 10px; color: #856404;"

      } else {

        div_style <- "background: #d4edda; padding: 10px; border-radius: 5px; margin-bottom: 10px; color: #155724;"

      }

      shiny::div(style = div_style, shiny::HTML(paste(status_parts, collapse = " | ")))

    })



    # ============================================================

    # 7. DotPlot Rendering Function

    # ============================================================



    output$plot_dotplot <- plotly::renderPlotly({

      pathways <- final_pathways()

      if (length(pathways) == 0) {

        return(plotly::plot_ly() |> plotly::layout(

          title = list(text = "No pathways to display", font = list(size = 14), x = 0.5),

          xaxis = list(showgrid = FALSE), yaxis = list(showgrid = FALSE)

        ))

      }

      data_list <- data_prep_list$data()

      shiny::req(data_list)

      df <- data_list$df

      plot_df <- df[df$ID %in% pathways, ]

      if (nrow(plot_df) == 0) {

        return(plotly::plot_ly() |> plotly::layout(title = list(text = "No matching pathways", font = list(size = 14))))

      }

      task <- list(

        gsea_res = data_list$gsea_res,

        meta = list(

          left_group = data_list$left_group,

          right_group = data_list$right_group,

          contrast_id = data_list$contrast_id

        )

      )

      class(task) <- "GseaTask"

      core_list <- tryCatch(

        {

          get_core_genes_list(task, plot_df$ID)

        },

        error = function(e) NULL

      )

      plot_df$CoreCount <- vapply(plot_df$ID, function(pid) {

        if (is.null(core_list) || is.null(core_list[[pid]])) {

          return(0L)

        }

        length(core_list[[pid]])

      }, integer(1))

      # IRON FIX (2026-06-30): coerce key columns to numeric BEFORE any

      # arithmetic. Upstream GSEA result objects sometimes return p.adjust /

      # pvalue / NES / setSize as character or factor (especially after

      # dplyr subsetting or tibble round-trips). Without this, -log10() or

      # sqrt(pmax()) silently receive character vectors and raise either

      # "non-numeric argument to mathematical function" or

      # "non-numeric argument to binary operator" depending on which operator

      # fires first.

      plot_df$NES       <- suppressWarnings(as.numeric(plot_df$NES))

      plot_df$p.adjust  <- suppressWarnings(as.numeric(plot_df$p.adjust))

      plot_df$pvalue    <- suppressWarnings(as.numeric(plot_df$pvalue))

      plot_df$setSize   <- suppressWarnings(as.numeric(plot_df$setSize))

      plot_df$CoreCount <- suppressWarnings(as.integer(plot_df$CoreCount))

      # Guard against NA values that would propagate to NaN/Inf in downstream

      # log/sqrt operations (these render as blank plotly dots).

      plot_df$NES[is.na(plot_df$NES)]       <- 0

      plot_df$p.adjust[is.na(plot_df$p.adjust)] <- 1

      plot_df$pvalue[is.na(plot_df$pvalue)]     <- 1

      plot_df$setSize[is.na(plot_df$setSize)]   <- 0

      plot_df$CoreCount[is.na(plot_df$CoreCount)] <- 0



      color_mode <- input$dotplot_color_mode

      color_vals <- switch(color_mode,

        "padj" = -log10(plot_df$p.adjust),

        "pval" = -log10(plot_df$pvalue),

        "nes" = abs(plot_df$NES),

        -log10(plot_df$p.adjust)

      )

      color_title <- switch(color_mode,

        "padj" = "-log10(FDR)",

        "pval" = "-log10(P-value)",

        "nes" = "|NES|"

      )

      size_mode <- input$dotplot_size_mode

      size_vals <- switch(size_mode,

        "core_size" = plot_df$CoreCount,

        "setsize" = plot_df$setSize,

        plot_df$CoreCount

      )

      # ----- enrichplot-style size mapping (IRON FIX) -----

      # Previous implementation used per-render min-max normalization, which

      # (a) pinned the smallest/largest pathways to fixed pixel sizes every

      #     render, hiding absolute magnitude differences between CoreCount

      #     (typical 5-50) and setSize (typical 50-500), so switching "Size by"

      #     produced visually identical dot patterns;

      # (b) made dot sizes incomparable across FDR/TopN parameter changes.

      # Fix: use a FIXED domain per variable + sqrt perceptual scale, mirroring

      # ggplot2::scale_size_continuous(range=...) semantics used by

      # clusterProfiler::dotplot / enrichplot::dotplot.

      size_range  <- c(5, 25)

      size_domain <- switch(size_mode,

        "core_size" = c(0, 50),   # CoreCount (leading edge) typical range

        "setsize"   = c(0, 500),  # MSigDB gene-set total size typical range

        c(0, 50)

      )

      sqrt_lo <- sqrt(size_domain[1])

      sqrt_hi <- sqrt(size_domain[2])

      size_scaled <- size_range[1] +

        (sqrt(pmax(size_vals, 0)) - sqrt_lo) / (sqrt_hi - sqrt_lo) *

        (size_range[2] - size_range[1])

      size_scaled <- pmin(pmax(size_scaled, size_range[1]), size_range[2])

      # IRON FIX (2026-06-29): compute sort_idx ONCE on the ORIGINAL plot_df,

      # then reorder plot_df / size_scaled / color_vals by the SAME index.

      #

      # Previous bug:

      #   plot_df <- plot_df[order(plot_df$NES, decreasing = TRUE), ]

      #   size_scaled <- size_scaled[order(plot_df$NES, decreasing = TRUE)]  # WRONG

      #   color_vals  <- color_vals[order(plot_df$NES, decreasing = TRUE)]   # WRONG

      # The second call to order() ran on the ALREADY-SORTED plot_df$NES, so

      # when the source df was not pre-sorted by NES, size_scaled/color_vals

      # got permuted inconsistent with the new plot_df row order. Symptom:

      # "Color by switch has no visible effect" because the misalignment was

      # invariant across color_mode changes.

      sort_idx <- order(plot_df$NES, decreasing = TRUE)

      plot_df    <- plot_df[sort_idx, ]

      size_scaled <- size_scaled[sort_idx]

      color_vals  <- color_vals[sort_idx]

      # Hover shows BOTH size dimensions; bold the one currently driving

      # dot size so the user can correlate visual size with the value.

      core_fmt <- if (size_mode == "core_size") "<b>Core Genes: %d</b>" else "Core Genes: %d"

      set_fmt  <- if (size_mode == "setsize")  "<b>Set Size: %d</b>"   else "Set Size: %d"

      hover_text <- sprintf(

        paste0(

          "<b>%s</b><br>FDR: %.2e<br>P-value: %.2e<br>NES: %.2f<br>",

          core_fmt, " | ", set_fmt

        ),

        plot_df$ID, plot_df$p.adjust, plot_df$pvalue, plot_df$NES,

        plot_df$CoreCount, plot_df$setSize

      )

      plotly::plot_ly(

        data = plot_df,

        x = plot_df$NES,

        y = seq_len(nrow(plot_df)),

        type = "scatter",

        mode = "markers",

        marker = list(

          size = size_scaled,

          color = color_vals,

          # IRON FIX (2026-06-29): unify with the PDF/ggplot2 export, which

          # uses scale_fill_viridis_c(option="D", direction=-1). Previously

          # the interactive DotPlot used "RdYlBu_r", so switching "Color by"

          # changed the value mapping but the visual hue was visibly different

          # from the exported PDF. Viridis-D reversed is perceptually uniform

          # and colorblind-safe (publication-grade).

          colorscale = list(

            list(0.0, "#440154"), list(0.25, "#3B528B"),

            list(0.5,  "#21918C"), list(0.75, "#5EC962"),

            list(1.0,  "#FDE725")

          ),

          cauto = TRUE,

          showscale = TRUE,

          colorbar = list(title = list(text = color_title, font = list(size = 12))),

          line = list(color = "black", width = 1)

        ),

        text = hover_text,

        hoverinfo = "text",

        hovertemplate = "%{text}<extra></extra>"

      ) |>

        plotly::layout(

          title = list(

            text = sprintf(

              "<b>Pathway DotPlot: %s vs %s</b><br><sub style='font-size:10px;color:gray;'>Left (NES&lt;0): Enriched in %s   |   Right (NES&gt;0): Enriched in %s</sub>",

              data_list$left_group, data_list$right_group,

              data_list$right_group, data_list$left_group),

            font = list(size = 12), x = 0.5, xanchor = "center"

          ),

          xaxis = list(title = "NES", zeroline = TRUE),

          yaxis = list(title = "", tickmode = "array", tickvals = seq_len(nrow(plot_df)), ticktext = plot_df$ID, tickfont = list(size = 9)),

          margin = list(l = 250, r = 50, t = 80, b = 50),

          showlegend = FALSE

        ) |>

        plotly::config(displayModeBar = TRUE, displaylogo = FALSE)

    })



    # ============================================================

    # 8. Node Selection State Management

    # ============================================================



    selected_nodes <- shiny::reactiveVal(character(0))



    output$selection_display <- shiny::renderUI({

      sel <- selected_nodes()

      if (length(sel) == 0) {

        shiny::div(style = "color: #666; font-style: italic; text-align: center;", "No nodes selected. Click on network nodes.")

      } else if (length(sel) == 1) {

        shiny::div(style = "color: #007bff; font-weight: bold; text-align: center;", sprintf("Selected: %s", sel[1]))

      } else {

        shiny::div(style = "color: #28a745; font-weight: bold; text-align: center;", sprintf("Selected: %s with %s", sel[1], sel[2]))

      }

    })



    shiny::observe({

      sel <- selected_nodes()

      if (length(sel) == 2) {

        shinyjs::enable("show_edge_detail")

      } else {

        shinyjs::disable("show_edge_detail")

      }

    })



    shiny::observeEvent(input$clear_selection, {

      selected_nodes(character(0))

    })



    shiny::observeEvent(input$show_edge_detail, {

      sel <- selected_nodes()

      if (length(sel) != 2) {

        shiny::showNotification("Please select exactly 2 nodes first", type = "warning")

        return()

      }



      # Read the current edge list from the reactiveVal. This establishes a

      # proper reactive dependency on plot_network's computation. If the user

      # clicks before the network has finished rendering, we degrade gracefully

      # instead of crashing on `NULL$from`.

      edge_list <- edge_list_rv()

      if (is.null(edge_list)) {

        shiny::showNotification(

          "Network is still being computed, please wait a moment and try again.",

          type = "warning",

          duration = 3

        )

        return()

      }



      from_pw <- sel[1]

      to_pw <- sel[2]



      edge_idx <- which(

        (edge_list$from == from_pw & edge_list$to == to_pw) |

          (edge_list$from == to_pw & edge_list$to == from_pw)

      )



      if (length(edge_idx) == 0) {

        shiny::showNotification(

          "No edge found between selected nodes. Try reducing min_shared threshold.",

          type = "warning",

          duration = 5

        )

        return()

      }



      i <- edge_idx[1]

      shared_count <- edge_list$shared[i]

      jaccard <- edge_list$weight[i]

      overlap_coef <- edge_list$overlap_coef[i]

      dice_coef <- edge_list$dice_coef[i]

      shared_genes <- unlist(edge_list$shared_genes[[i]])



      pathway_a_genes <- core_list[[from_pw]]

      pathway_b_genes <- core_list[[to_pw]]



      # ============ Modify start: display all genes in showdetail ============

      # In showdetail, do not limit the number of genes; display all shared genes

      # Use comma and space to separate genes



      # Format all genes for display

      gene_buttons <- vapply(shared_genes, function(g) {

        sprintf(

          '<span style="display:inline-block;background:#e3f2fd;padding:4px 8px;margin:2px;border-radius:4px;font-size:12px;">%s</span>',

          g

        )

      }, character(1))



      # No longer show "+N more" hint

      gene_display <- paste(gene_buttons, collapse = ",")



      # Modal UI

      shiny::showModal(shiny::modalDialog(

        title = HTML(paste0("<b>Edge Detail:</b> ", from_pw, " with ", to_pw)),

        size = "l",

        easyClose = TRUE,

        footer = shiny::modalButton("Close"),



        # Similarity metrics

        HTML("<div style='background:#f8f9fa;padding:15px;border-radius:8px;margin-bottom:15px;'>"),

        HTML("<h4 style='margin-top:0;'>Similarity Metrics</h4>"),

        HTML("<table style='width:100%;text-align:center;'>"),

        HTML("<tr><td><h3>", sprintf("%.3f", jaccard), "</h3><small>Jaccard</small></td>"),

        HTML("<td><h3>", sprintf("%.3f", overlap_coef), "</h3><small>Overlap</small></td>"),

        HTML("<td><h3>", sprintf("%.3f", dice_coef), "</h3><small>Dice</small></td>"),

        HTML("<td><h3>", shared_count, "</h3><small>Shared</small></td></tr>"),

        HTML("</table></div>"),



        # Shared genes

        HTML("<div style='background:#fff;padding:15px;border-radius:8px;border:1px solid #dee2e6;'>"),

        HTML("<h4 style='margin-top:0;'>Shared Core Genes (", shared_count, ")</h4>"),

        HTML(gene_display),

        HTML("</div>")

      ))

    })



    # ============================================================

    # 9. Network Rendering Function

    # ============================================================

    #

    # edge_list and node_df are computed inside output$plot_network and read

    # by observeEvent(input$show_edge_detail). Previously they used `<<-` to

    # mutate module-scope variables, but this had two problems:

    #   1. Bioconductor review: `<<-` is discouraged.

    #   2. Functional bug: `<<-` into a non-reactive variable does NOT trigger

    #      reactive dependency, so observeEvent could read a stale/NULL value

    #      during async re-rendering of plot_network.

    # Using reactiveVal() establishes a proper reactive dependency and lets

    # observeEvent gracefully degrade when the network is still being computed.

    edge_list_rv <- shiny::reactiveVal(NULL)

    node_df_rv <- shiny::reactiveVal(NULL)

    core_list <- NULL



    output$plot_network <- plotly::renderPlotly({

      pathways <- final_pathways()



      if (length(pathways) == 0) {

        return(plotly::plot_ly() |> plotly::layout(

          title = list(text = "No pathways to display", font = list(size = 14), x = 0.5),

          xaxis = list(showgrid = FALSE, showticklabels = FALSE),

          yaxis = list(showgrid = FALSE, showticklabels = FALSE)

        ))

      }



      data_list <- data_prep_list$data()

      shiny::req(data_list)



      task <- list(

        gsea_res = data_list$gsea_res,

        meta = list(

          left_group = data_list$left_group,

          right_group = data_list$right_group,

          contrast_id = data_list$contrast_id

        )

      )

      class(task) <- "GseaTask"



      core_list <- tryCatch(

        {

          get_core_genes_list(task, pathways)

        },

        error = function(e) NULL

      )


      if (is.null(core_list) || length(core_list) == 0) {

        return(plotly::plot_ly() |> plotly::layout(

          title = list(text = "No core genes found", font = list(size = 14)),

          xaxis = list(showgrid = FALSE), yaxis = list(showgrid = FALSE)

        ))

      }



      valid_pathways <- names(core_list)[vapply(core_list, function(x) length(x) > 0, logical(1))]

      if (length(valid_pathways) == 0) {

        return(plotly::plot_ly() |> plotly::layout(

          title = list(text = "All pathways have empty core genes", font = list(size = 14)),

          xaxis = list(showgrid = FALSE), yaxis = list(showgrid = FALSE)

        ))

      }



      min_shared <- input$min_shared

      if (is.null(min_shared)) min_shared <- 3



      # Get hover max genes limit

      hover_max_genes <- input$hover_max_genes

      if (is.null(hover_max_genes)) hover_max_genes <- 10



      # Build edge list (now includes jaccard, overlap_coef, dice_coef).

      # Use a local object for downstream code in this renderPlotly scope,

      # and push it to the module-scope reactiveVal so observeEvent can read

      # the latest value with a proper reactive dependency.

      edge_list <- tryCatch(

        {

          build_edge_list_safely(core_list[valid_pathways], min_shared_genes = min_shared)

        },

        error = function(e) NULL

      )

      edge_list_rv(edge_list)



      if (is.null(edge_list) || nrow(edge_list) == 0) {

        return(plotly::plot_ly() |> plotly::layout(

          title = list(text = paste0("No edges found (min_shared=", min_shared, ")"), font = list(size = 14), x = 0.5),

          xaxis = list(showgrid = FALSE, showticklabels = FALSE),

          yaxis = list(showgrid = FALSE, showticklabels = FALSE)

        ))

      }



      # ============================================================

      # Edge width mapping (user-selectable via input$network_edge_width_mode)

      # ------------------------------------------------------------

      # Two encoding strategies are exposed because they serve different

      # analytical goals:

      #   "weight" (default)  -> edge width proportional to Jaccard value.

      #                          Faithful to the underlying similarity

      #                          magnitude (matches emapplot / ggraph

      #                          conventions). Recommended for publication.

      #   "rank"              -> edge width assigned by Jaccard rank order.

      #                          Guarantees uniform visual spacing across

      #                          edges regardless of absolute weight. Useful

      #                          when the goal is "make every edge visible"

      #                          (e.g. dense networks with small weight

      #                          variance). Original behavior.

      # ============================================================



      # Sort by Jaccard (weight) in descending order

      edge_list <- edge_list[order(edge_list$weight, decreasing = TRUE), ]



      width_mode <- if (is.null(input$network_edge_width_mode)) "weight"

                    else input$network_edge_width_mode



      n_edges <- nrow(edge_list)

      if (width_mode == "rank") {

        # Rank-based: assign 1-5 pixels based on rank (original behavior)

        edge_list$width_rank <- vapply(seq_len(n_edges), function(i) {

          width <- 5 - (i - 1) * (4 / max(1, n_edges - 1))

          max(1, min(5, round(width)))

        }, numeric(1))

      } else {

        # Weight-based: linear mapping of Jaccard in [0,1] to [1,8] pixels

        edge_list$width_rank <- 1 + edge_list$weight * 7

      }



      # Normal edge width range: rank->1-5, weight->1-8 pixels

      edge_list$edge_width_normal <- edge_list$width_rank



      # Selected edge fixed width: 7-8 pixels thicker than normal thickest (5 + 8 = 13)

      edge_list$edge_width_selected <- 13



      # Build igraph object

      g <- tryCatch(

        {

          igraph::graph_from_data_frame(edge_list, directed = FALSE, vertices = valid_pathways)

        },

        error = function(e) NULL

      )



      if (is.null(g)) {

        return(plotly::plot_ly() |> plotly::layout(

          title = list(text = "Graph construction failed", font = list(size = 14)),

          xaxis = list(showgrid = FALSE), yaxis = list(showgrid = FALSE)

        ))

      }



      layout_algo <- input$network_layout

      # Defensive default: if input$seed is NULL (e.g. during fast reactive

      # invalidation before the numericInput has settled), fall back to 42.

      # This guarantees withr::with_seed() always receives an integer.

      seed_val <- if (is.null(input$seed)) 42L else as.integer(input$seed)



      layout_coords <- tryCatch(

        {

          switch(layout_algo,

            "fr" = {

              # Local seed scope: ensures reproducible FR layout across renders

              # without polluting the caller's RNG state.

              # seed_val is controlled by the user via input$seed (Shiny numericInput).

              withr::with_seed(seed_val, igraph::layout_with_fr(g))

            },

            "kk" = igraph::layout_with_kk(g),

            "circle" = igraph::layout_in_circle(g),

            {

              # Default: Fruchterman-Reingold; local seed scope as above.

              withr::with_seed(seed_val, igraph::layout_with_fr(g))

            }

          )

        },

        error = function(e) NULL

      )



      if (is.null(layout_coords)) {

        return(plotly::plot_ly() |> plotly::layout(

          title = list(text = "Layout calculation failed", font = list(size = 14)),

          xaxis = list(showgrid = FALSE), yaxis = list(showgrid = FALSE)

        ))

      }



      node_df <- data.frame(

        name = igraph::V(g)$name,

        x = layout_coords[, 1],

        y = layout_coords[, 2],

        stringsAsFactors = FALSE

      )

      node_df_rv(node_df)



      res_df <- as.data.frame(task$gsea_res@result)

      node_info <- res_df[match(node_df$name, res_df$ID), ]

      node_df$NES <- node_info$NES

      node_df$FDR <- node_info$p.adjust

      node_df$CoreCount <- vapply(node_df$name, function(n) {

        if (is.null(core_list[[n]])) {

          return(0L)

        }

        length(core_list[[n]])

      }, integer(1))

      # IRON FIX (2026-06-30): coerce to numeric to prevent

      # "non-numeric argument to binary operator" when NES/FDR arrive as

      # character from the upstream GSEA result data.frame.

      node_df$NES <- suppressWarnings(as.numeric(node_df$NES))

      node_df$FDR <- suppressWarnings(as.numeric(node_df$FDR))

      node_df$NES[is.na(node_df$NES)] <- 0

      node_df$FDR[is.na(node_df$FDR)] <- 1

      node_df$color_val <- ifelse(is.na(node_df$NES), 0, node_df$NES)



      sel_nodes <- selected_nodes()



      node_hover_text <- sprintf(

        "<b>%s</b><br>NES: %.2f<br>FDR: %.2e<br>Core Genes: %d",

        node_df$name, node_df$NES, node_df$FDR, node_df$CoreCount

      )



      title_text <- sprintf(

        "Pathway Network: %s vs %s<br><sub>%d nodes, %d edges | Width: Jaccard rank (1-5px)</sub>",

        task$meta$left_group, task$meta$right_group,

        igraph::vcount(g), igraph::ecount(g)

      )



      p <- plotly::plot_ly(source = ns("network_plot"))



      # ============================================================

      # Add edges (using Jaccard ranking-mapped thickness)

      # ============================================================



      for (i in seq_len(nrow(edge_list))) {

        from_node <- node_df[node_df$name == edge_list$from[i], ]

        to_node <- node_df[node_df$name == edge_list$to[i], ]

        is_selected <- (edge_list$from[i] %in% sel_nodes) && (edge_list$to[i] %in% sel_nodes)



        if (is_selected) {

          # Selected edge: orange + fixed thickness 13px

          edge_color <- "#FF6600"

          edge_width <- edge_list$edge_width_selected[i]

        } else {

          # Normal edge: gray + Jaccard ranking thickness 1-5px

          edge_color <- "rgba(150, 150, 150, 0.5)"

          edge_width <- edge_list$edge_width_normal[i]

        }



        p <- p |> plotly::add_trace(

          type = "scatter", mode = "lines",

          x = c(from_node$x, to_node$x, NA),

          y = c(from_node$y, to_node$y, NA),

          line = list(color = edge_color, width = edge_width),

          hoverinfo = "skip", showlegend = FALSE, inherit = TRUE

        )

      }



      # ============================================================

      # Add edge hover (enhanced: display shared gene list)

      # ============================================================



      for (i in seq_len(nrow(edge_list))) {

        from_node <- node_df[node_df$name == edge_list$from[i], ]

        to_node <- node_df[node_df$name == edge_list$to[i], ]



        is_selected <- (edge_list$from[i] %in% sel_nodes) && (edge_list$to[i] %in% sel_nodes)



        shared_genes_vec <- edge_list$shared_genes[[i]]

        shared_count <- edge_list$shared[i]



        if (length(shared_genes_vec) > 0) {

          display_genes <- shared_genes_vec[seq_len(min(length(shared_genes_vec), hover_max_genes))]

          genes_display <- paste(display_genes, collapse = ", ")

          if (length(shared_genes_vec) > hover_max_genes) {

            remaining <- length(shared_genes_vec) - hover_max_genes

            genes_display <- paste0(

              genes_display,

              sprintf(" <span style='color:#ccc;'>(+%d more)</span>", remaining)

            )

          }

        } else {

          genes_display <- "(none)"

        }



        # Simplified format to reduce line breaks

        edge_hover <- sprintf(

          "<b>%s with %s</b><br>Shared Genes (%d): %s<br>Jaccard: %.4f | Overlap: %.4f | Dice: %.4f",

          edge_list$from[i], edge_list$to[i], shared_count,

          genes_display,

          edge_list$weight[i], edge_list$overlap_coef[i], edge_list$dice_coef[i]

        )



        hover_bg <- if (is_selected) "#FF8C00" else "#333"



        p <- p |> plotly::add_trace(

          type = "scatter", mode = "markers",

          x = c(mean(c(from_node$x, to_node$x))),

          y = c(mean(c(from_node$y, to_node$y))),

          marker = list(size = 12, opacity = 0, color = "transparent"),

          text = edge_hover,

          hoverinfo = "text",

          hoverlabel = list(

            bgcolor = hover_bg,

            font = list(color = "white", size = 12),

            align = "left", # <- left-aligned

            bordercolor = hover_bg

          ),

          showlegend = FALSE,

          inherit = TRUE

        )

      }



      # ============================================================

      # Add nodes

      # ============================================================



      node_colors <- ifelse(node_df$name %in% sel_nodes, "#FFD700", node_df$color_val)

      node_sizes <- ifelse(node_df$name %in% sel_nodes, 25, 15)

      node_border <- ifelse(node_df$name %in% sel_nodes, 4, 1.5)



      p <- p |> plotly::add_trace(

        type = "scatter", mode = "markers",

        x = node_df$x, y = node_df$y,

        marker = list(

          size = node_sizes,

          color = node_colors,

          colorscale = list(list(0, "#377EB8"), list(0.5, "white"), list(1, "#E41A1C")),

          cauto = FALSE, cmin = -3, cmax = 3, showscale = TRUE,

          colorbar = list(title = list(text = "NES", font = list(size = 12)), len = 0.5, y = 0.5),

          line = list(color = "black", width = node_border)

        ),

        text = node_df$name,

        hovertemplate = paste(node_hover_text, "<extra></extra>"),

        showlegend = FALSE, key = node_df$name, inherit = TRUE

      )



      # ============================================================

      # Add labels (selected nodes are larger and bolder)

      # ------------------------------------------------------------

      # IRON FIX (2026-06-29): Network text rendering improvements.

      #   * Default label size bumped 9 -> 11 (legibility at small zoom).

      #   * White "halo" shadow drawn UNDER the label so dark text remains

      #     readable against clustered nodes / overlapping edges. We achieve

      #     this by drawing TWO text traces: first a white shadow, then the

      #     colored label on top. plotly has no native stroke/fill, so the

      #     two-trace trick is the canonical workaround.

      #   * Vertical offset 0.12 -> 0.18 to avoid covering the node circle.

      #   * Display the pathway short name (strip "HALLMARK_" / dataset

      #     prefix) to reduce horizontal clutter; hover still shows the

      #     full ID.

      # ============================================================



      label_size     <- ifelse(node_df$name %in% sel_nodes, 16, 11)

      label_color    <- ifelse(node_df$name %in% sel_nodes, "#FF6600", "#222")

      label_bold     <- ifelse(node_df$name %in% sel_nodes, "bold",   "normal")

      label_short    <- sub("^[^_]*_", "", node_df$name)

      label_offset   <- 0.18



      # White shadow trace (drawn first -> underneath the colored label)

      p <- p |> plotly::add_trace(

        type = "scatter", mode = "text",

        x = node_df$x, y = node_df$y + label_offset,

        text = label_short,

        textposition = "top center",

        textfont = list(size = label_size + 1, color = "#FFFFFF",

                        family = "Arial Black, sans-serif"),

        hoverinfo = "skip", showlegend = FALSE, inherit = TRUE

      )



      # Colored label trace on top

      p <- p |> plotly::add_trace(

        type = "scatter", mode = "text",

        x = node_df$x, y = node_df$y + label_offset,

        text = label_short,

        textposition = "top center",

        textfont = list(size = label_size, color = label_color,

                        family = "Arial, sans-serif", font = label_bold),

        hoverinfo = "skip", showlegend = FALSE, inherit = TRUE

      )



      x_range <- c(min(layout_coords[, 1]) - 0.5, max(layout_coords[, 1]) + 0.5)

      y_range <- c(min(layout_coords[, 2]) - 0.5, max(layout_coords[, 2]) + 0.5)



      p |>

        plotly::layout(

          title = list(text = title_text, font = list(size = 12), x = 0.5, xanchor = "center"),

          xaxis = list(

            title = "",

            showgrid = FALSE,

            showticklabels = FALSE,

            zeroline = FALSE,

            range = x_range,

            scaleanchor = "x", # Maintain aspect ratio

            scaleratio = 1

          ),

          yaxis = list(

            title = "",

            showgrid = FALSE,

            showticklabels = FALSE,

            zeroline = FALSE,

            range = y_range,

            scaleanchor = "x", # Maintain aspect ratio

            scaleratio = 1

          ),

          hovermode = "closest",

          dragmode = "pan",

          showlegend = FALSE,

          margin = list(l = 50, r = 50, t = 80, b = 50),

          paper_bgcolor = "rgba(0,0,0,0)", # Transparent background

          plot_bgcolor = "rgba(0,0,0,0)", # Transparent canvas

          autosize = TRUE # <- Key: auto-resize

        ) |>

        plotly::config(

          displayModeBar = TRUE,

          displaylogo = FALSE,

          modeBarButtonsToRemove = c("lasso2d", "select2d"),

          responsive = TRUE # <- Key: responsive

        )

    })



    # ============================================================

    # 10. Node Click Event Listener

    # ============================================================



    shiny::observeEvent(plotly::event_data("plotly_click", source = ns("network_plot")), {

      click_data <- plotly::event_data("plotly_click", source = ns("network_plot"))

      if (is.null(click_data) || is.null(click_data$pointNumber)) {

        return()

      }



      n_edges <- nrow(edge_list)

      point_idx <- click_data$pointNumber + 1



      if (point_idx > nrow(node_df) || point_idx < 1) {

        return()

      }



      node_name <- as.character(node_df$name[point_idx])

      current_sel <- selected_nodes()



      if (length(current_sel) < 2) {

        if (!(node_name %in% current_sel)) {

          selected_nodes(c(current_sel, node_name))

        }

      } else {

        selected_nodes(node_name)

      }

    })



    # ============================================================

    # 12. Export Center (DotPlot + Network)

    # ------------------------------------------------------------

    # Provides a unified modal for downloading the current view as a

    # static ggplot2 publication figure (PDF/PNG/SVG/TIFF) and for

    # copying the corresponding reproducible R code.

    #

    # Design decisions:

    #   * Static rendering uses ggplot2 + ggsave (zero external deps,

    #     unlike plotly::save_image which requires kaleido/orca).

    #   * Code is generated via generate_dotplot_code() /

    #     generate_network_code() in R/15_code_generator.R.

    #   * Clipboard copy uses clipr (already in Imports).

    # ============================================================



    # ----- Shared export modal UI builder -----

    # ------------------------------------------------------------

    # Redesign (2026-06-29):

    #   * Live preview pane reflects the current palette/size/width/height.

    #   * Download buttons split into "PDF" and "PNG" (the two formats users

    #     actually use most); SVG/TIFF remain available via the Format

    #     dropdown for advanced users.

    #   * Reactive preview re-runs the generator code in a sandboxed env

    #     (same parent=globalenv() fix used by the download handler).

    # ------------------------------------------------------------

    .build_export_modal <- function(plot_title, show_aesthetics = TRUE) {

      aesthetics_block <- if (isTRUE(show_aesthetics)) {

        shiny::tagList(

          shiny::hr(),

          shiny::h5("Aesthetics"),

          shiny::selectInput(ns("exp_palette"), "Color palette",

            choices = c("Viridis-D (default)" = "D",

                        "Viridis-C"          = "C",

                        "Magma"               = "A",

                        "Inferno"             = "B"),

            selected = "D"),

          shiny::sliderInput(ns("exp_point_range"), "Point size range",

            min = 1, max = 15, value = c(3, 8))

        )

      } else {

        # For Network exports: palette / point range have no effect because

        # the network generator uses a fixed Up/Down color scheme and node

        # sizes are derived from -log10(FDR) internally. Hide these inputs

        # so the UI does not mislead the user.

        shiny::tagList(

          shiny::hr(),

          shiny::div(

            style = "background:#f8f9fa; padding:8px 10px; border-radius:4px; color:#666; font-size:12px;",

            shiny::icon("info-circle"),

            " Network coloring uses a fixed Up (red) / Down (blue) scheme.",

            " Adjust edge width and layout seed in the main panel instead."

          )

        )

      }



      shiny::modalDialog(

        title = sprintf("Export Publication Plot - %s", plot_title),

        size = "l",

        footer = NULL,

        easyClose = TRUE,

        shiny::fluidRow(

          # Left column: parameters

          shiny::column(5,

            shiny::h5("Dimensions"),

            shiny::fluidRow(

              shiny::column(6,

                shiny::numericInput(ns("exp_width"),  "Width (inch)",  value = 9, min = 2, max = 24)

              ),

              shiny::column(6,

                shiny::numericInput(ns("exp_height"), "Height (inch)", value = 7, min = 2, max = 24)

              )

            ),

            shiny::numericInput(ns("exp_dpi"), "DPI", value = 300, min = 72, max = 600),



            shiny::h6("Canvas Margin (pt)"),

            shiny::fluidRow(

              shiny::column(3, shiny::numericInput(ns("exp_margin_top"),    "Top",    value = 18, min = 0, max = 80)),

              shiny::column(3, shiny::numericInput(ns("exp_margin_bottom"), "Bottom", value = 18, min = 0, max = 80)),

              shiny::column(3, shiny::numericInput(ns("exp_margin_left"),   "Left",   value = 18, min = 0, max = 80)),

              shiny::column(3, shiny::numericInput(ns("exp_margin_right"),  "Right",  value = 18, min = 0, max = 80))

            ),

            aesthetics_block,

            shiny::hr(),

            shiny::h5("Download"),

            shiny::fluidRow(

              shiny::column(6,

                shiny::downloadButton(ns("exp_download_pdf"),

                  "Download PDF",

                  class = "btn-danger btn-block",

                  icon  = shiny::icon("file-pdf"))

              ),

              shiny::column(6,

                shiny::downloadButton(ns("exp_download_png"),

                  "Download PNG",

                  class = "btn-primary btn-block",

                  icon  = shiny::icon("file-image"))

              )

            ),

            shiny::hr(),

            shiny::fluidRow(

              shiny::column(6,

                shiny::selectInput(ns("exp_format"), "Other format",

                  choices = c("SVG" = "svg", "TIFF" = "tiff"),

                  selected = "svg")

              ),

              shiny::column(6,

                shiny::div(style = "margin-top: 22px;",

                  shiny::downloadButton(ns("exp_download_other"),

                    "Download",

                    class = "btn-default btn-block")

                )

              )

            ),

            shiny::hr(),

            shiny::fluidRow(

              shiny::column(6,

                shiny::actionButton(ns("exp_copy_code"),

                  "Copy R Code",

                  class = "btn-info btn-block",

                  icon = shiny::icon("clipboard"))

              ),

              shiny::column(6,

                shiny::actionButton(ns("exp_dismiss"),

                  "Close",

                  class = "btn-default btn-block")

              )

            )

          ),

          # Right column: live preview / external canvas

          shiny::column(7,

            shiny::h5("Live Preview"),

            shiny::div(style = "background:#f5f5f5; border:1px solid #ddd; border-radius:4px; padding:8px; width:100%; height:520px; display:flex; align-items:center; justify-content:center; overflow:auto;",

              shiny::plotOutput(ns("exp_preview")) |>

                shinycssloaders::withSpinner(type = 6, color = "#28a745")

            ),

            shiny::tags$small(style = "color: #666; display:block; margin-top:6px;",

              "Preview reflects current Width / Height and is auto-scaled to fit while keeping the export aspect ratio.")

          )

        )

      )

    }



    # ----- Reactive: live preview ggplot object -----

    # IRON FIX (2026-06-30): the generated R code includes a `print(p)` line

    # so the clipboard-copy path renders the plot when the user pastes it

    # into their own R console. But when we eval() the same code inside a

    # Shiny reactive for the live preview, `print(p)` opens a separate

    # graphics device (windows()/pdf()/quartz()) OUTSIDE Shiny's device

    # chain. On ggplot2 >= 4.0 with an active popup device, the subsequent

    # renderPlot print can surface "non-numeric argument to binary operator"

    # from unit/arithmetic ops occurring in the now-stale device context.

    # Fix: override `print` in eval_env to be a no-op so the popup never

    # opens; renderPlot does its own print on its own device.

    .exp_preview_plot <- shiny::reactive({

      tgt <- current_export_target()

      if (is.null(tgt)) return(NULL)

      # Wrap code generation in tryCatch so errors here surface as a Shiny

      # notification (and a clean NULL return) rather than propagating into

      # renderPlot and showing as a grey error pane.

      code <- tryCatch(

        .current_export_code(),

        error = function(e) {

          shiny::showNotification(

            sprintf("[preview] code generation failed: %s", e$message),

            type = "error", duration = 8)

          ""

        })

      if (!nzchar(code)) return(NULL)

      eval_env <- new.env(parent = globalenv())

      eval_env$p <- NULL

      eval_env$print <- base::invisible   # suppress popup device from print(p)
      eval_env$gsea_res <- gsea_res

      tryCatch(

        eval(parse(text = code), envir = eval_env),

        error = function(e) {

          shiny::showNotification(

            sprintf("[preview] %s", e$message),

            type = "error", duration = 8)

          NULL

        }

      )

      eval_env$p

    })



    # Preview plotOutput scales so the on-screen aspect ratio ALWAYS matches

    # the eventual saved figure, with the longer side capped to a max pixel

    # dimension so ultra-wide or ultra-tall specs still fit on screen.

    #

    # Algorithm: take the user's chosen width_in / height_in (inches), find

    # which is the longer side, pin it to max_dim (720 px), and scale the

    # shorter side by the same ratio. This preserves aspect ratio exactly

    # while bounding the preview to a sensible size:

    #   * 16x4 banner  -> 720x180 (wide-and-short)

    #   * 9x7 default  -> 720x560 (matches saved file proportions)

    #   * 4x16 column  -> 180x720 (tall-and-narrow, container scrolls)

    #

    # IRON FIX (2026-06-30): previously width/height were capped

    # INDEPENDENTLY at 720x540, which distorted the aspect ratio when

    # width > 7.2 inch (after that, width was clipped and the preview

    # stopped growing horizontally, showing a squished image).

    # Preview device matching export physical size (2026-07-06 v3):
    # Shiny renderPlot does NOT accept res as a function, so res is fixed
    # at 72 and pixel dims are computed as inch * 72 (with a pixel cap).
    # Physical size = width_px / 72 = w_in inches, matching the export.
    .preview_dims <- function(w_in, h_in, res = 72, max_px = 1600) {

      w_in <- suppressWarnings(as.numeric(w_in[1]))

      h_in <- suppressWarnings(as.numeric(h_in[1]))

      if (length(w_in) == 0 || is.na(w_in) || w_in <= 0) w_in <- 9

      if (length(h_in) == 0 || is.na(h_in) || h_in <= 0) h_in <- 7

      scale <- min(1, max_px / (max(w_in, h_in) * res))

      list(width  = round(w_in * res * scale),

           height = round(h_in * res * scale))

    }



    output$exp_preview <- shiny::renderPlot(

      {

        p <- .exp_preview_plot()

        shiny::req(p)

        p

      },

      res  = 72,

      width  = function() {

        d <- .preview_dims(input$exp_width, input$exp_height); d$width

      },

      height = function() {

        d <- .preview_dims(input$exp_width, input$exp_height); d$height

      }

    )



    # ----- Reactive: current DotPlot data slice -----

    dotplot_export_df <- shiny::reactive({

      pathways <- final_pathways()

      if (length(pathways) == 0) return(NULL)

      data_list <- data_prep_list$data()

      shiny::req(data_list)

      df <- data_list$df

      plot_df <- df[df$ID %in% pathways, ]

      if (nrow(plot_df) == 0) return(NULL)



      task <- list(

        gsea_res = data_list$gsea_res,

        meta = list(left_group = data_list$left_group,

                    right_group = data_list$right_group,

                    contrast_id = data_list$contrast_id)

      )

      class(task) <- "GseaTask"

      core_list <- tryCatch(get_core_genes_list(task, plot_df$ID), error = function(e) NULL)

      plot_df$CoreCount <- vapply(plot_df$ID, function(pid) {

        if (is.null(core_list) || is.null(core_list[[pid]])) 0L else length(core_list[[pid]])

      }, integer(1))

      plot_df

    })



    # ----- Reactive: current Network edges + nodes -----

    network_export_data <- shiny::reactive({

      el <- edge_list_rv(); nd <- node_df_rv()

      if (is.null(el) || is.null(nd) || nrow(el) == 0) return(NULL)

      list(edge_df = el[, c("from", "to", "weight", "shared")],

           node_df = nd)

    })



    # ----- Reactive: ggplot object + code string for the OPEN modal -----

    current_export_target <- shiny::reactiveVal(NULL)



    shiny::observeEvent(input$open_dotplot_export, {

      if (is.null(dotplot_export_df())) {

        shiny::showNotification("No pathways to export.", type = "warning"); return()

      }

      current_export_target("dotplot")

      shiny::showModal(.build_export_modal("Pathway DotPlot"))

    })



    shiny::observeEvent(input$open_network_export, {

      if (is.null(network_export_data())) {

        shiny::showNotification("No network to export.", type = "warning"); return()

      }

      current_export_target("network")

      # Network exports ignore palette/point_range (they use a fixed Up/Down

      # scheme), so hide the Aesthetics block to avoid misleading the user.

      shiny::showModal(.build_export_modal("Pathway Network", show_aesthetics = FALSE))

    })



    shiny::observeEvent(input$exp_dismiss, shiny::removeModal())



    # ----- Shared code builder (depends on current target) -----

    .current_export_code <- function() {

      pal <- if (is.null(input$exp_palette)) "D" else input$exp_palette

      pr  <- if (is.null(input$exp_point_range)) c(3, 8) else input$exp_point_range

      data_list <- data_prep_list$data()

      lg <- if (is.null(data_list$left_group))  "Left"  else data_list$left_group

      rg <- if (is.null(data_list$right_group)) "Right" else data_list$right_group



      tgt <- current_export_target()

      if (is.null(tgt)) return("")



      if (tgt == "dotplot") {

        generate_dotplot_code(

          gsea_res_var = "gsea_res",

          contrast_id  = data_list$contrast_id,

          pathways     = final_pathways(),

          color_mode   = input$dotplot_color_mode,

          size_mode    = input$dotplot_size_mode,

          left_group   = lg, right_group = rg,

          palette      = pal, point_range = pr)

      } else if (tgt == "network") {

        seed_val <- if (is.null(input$seed)) 42L else as.integer(input$seed)

        data_list <- data_prep_list$data()

        generate_network_code(

          gsea_res_var = "gsea_res",

          contrast_id  = data_list$contrast_id,

          pathways     = final_pathways(),

          min_shared_genes = if (is.null(input$min_shared)) 3 else input$min_shared,

          width_mode   = if (is.null(input$network_edge_width_mode)) "weight"

                         else input$network_edge_width_mode,

          seed = seed_val,

          margin_top    = if (is.null(input$exp_margin_top)    || is.na(input$exp_margin_top))    18 else input$exp_margin_top,

          margin_bottom = if (is.null(input$exp_margin_bottom) || is.na(input$exp_margin_bottom)) 18 else input$exp_margin_bottom,

          margin_left   = if (is.null(input$exp_margin_left)   || is.na(input$exp_margin_left))   18 else input$exp_margin_left,

          margin_right  = if (is.null(input$exp_margin_right)  || is.na(input$exp_margin_right))  18 else input$exp_margin_right)

      } else ""

    }



    # ----- Internal helper: render to a specific format -----

    # ------------------------------------------------------------

    # Three download handlers (exp_download_pdf / _png / _other) below all

    # delegate to this helper so behavior is identical modulo format.

    # ------------------------------------------------------------

    .export_render_to_file <- function(file, fmt) {

      code <- .current_export_code()

      if (!nzchar(code)) {

        shiny::showNotification("Nothing to export.", type = "warning"); return()

      }

      # parent=globalenv() so lexical lookup resolves attached packages

      # (see IRON FIX comment in v0.99.11 commit). baseenv() reproduces

      # "could not find function 'ggplot'" and the .htm fallback symptom.

      eval_env <- new.env(parent = globalenv())

      eval_env$p <- NULL

      eval_env$print <- base::invisible   # suppress popup device; ggsave() draws itself
      eval_env$gsea_res <- gsea_res

      tryCatch({

        eval(parse(text = code), envir = eval_env)

      }, error = function(e) {

        shiny::showNotification(sprintf("Render failed: %s", e$message),

                                type = "error", duration = 8)

      })

      if (is.null(eval_env$p)) return()

      ggplot2::ggsave(

        file, eval_env$p,

        width  = if (is.null(input$exp_width))  9 else input$exp_width,

        height = if (is.null(input$exp_height)) 7 else input$exp_height,

        dpi    = if (is.null(input$exp_dpi))   300 else input$exp_dpi,

        device = fmt

      )

    }



    output$exp_download_pdf <- shiny::downloadHandler(

      filename = function() sprintf("GSEAlens_%s_%s.pdf",

                                    current_export_target() %||% "plot", Sys.Date()),

      content  = function(file) .export_render_to_file(file, "pdf")

    )



    output$exp_download_png <- shiny::downloadHandler(

      filename = function() sprintf("GSEAlens_%s_%s.png",

                                    current_export_target() %||% "plot", Sys.Date()),

      content  = function(file) .export_render_to_file(file, "png")

    )



    output$exp_download_other <- shiny::downloadHandler(

      filename = function() sprintf("GSEAlens_%s_%s.%s",

                                    current_export_target() %||% "plot", Sys.Date(),

                                    if (is.null(input$exp_format)) "svg" else input$exp_format),

      content  = function(file) .export_render_to_file(file, input$exp_format)

    )



    # ----- Copy code to clipboard -----

    shiny::observeEvent(input$exp_copy_code, {

      code <- .current_export_code()

      if (!nzchar(code)) {

        shiny::showNotification("Nothing to copy.", type = "warning"); return()

      }

      ok <- tryCatch({ clipr::write_clip(code); TRUE }, error = function(e) FALSE)

      if (isTRUE(ok)) {

        shiny::showNotification("R code copied to clipboard.", type = "message")

      } else {

        # Fallback: show code in a new modal so the user can manually copy

        shiny::showModal(shiny::modalDialog(

          title = "Reproducible R Code (copy manually)",

          size = "l", easyClose = TRUE,

          shiny::pre(style = "max-height: 60vh; overflow-y: auto; font-size: 11px;",

                     code),

          footer = shiny::modalButton("Close")

        ))

      }

    })



    # ============================================================

    # 13. Return Values

    # ============================================================



    return(list(

      final_pathways = final_pathways,

      network_mode = network_mode

    ))

  })

}

