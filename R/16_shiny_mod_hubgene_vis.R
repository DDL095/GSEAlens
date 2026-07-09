#' @title HubGene Network Module UI

#' @description

#' Shiny module UI for visualizing HubGene networks using visNetwork.

#' Provides interactive controls for pathway selection modes, physics

#' simulation parameters, node sizing, and network statistics display.

#'

#' @param id Character string used to namespace the module.

#'

#' @return A Shiny tagList containing the module's UI elements.

#'

#' @keywords internal



mod_hubgene_vis_ui <- function(id) {

    ns <- shiny::NS(id)



    shiny::tagList(

    shiny::fluidRow(

            # ─── Left Control Panel ───

            shiny::column(

        width = 3,

        shiny::div(

                    class = "well",

                    style = "padding: 15px; max-height: 90vh; overflow-y: auto;",

                    shiny::h4("HubGene Network (visNetwork)"),

                    shiny::hr(),



                    # ── Mode Selection ──

                    shiny::h5("Pathway Source Mode"),

                    shiny::radioButtons(

            ns("vis_mode"),

            label = NULL,

            choices = c(

                            "Top N" = "mode_topN",

                            "Selected from Table" = "mode_select"

            ),

            selected = "mode_topN",

            width = "100%"

                    ),

                    shiny::hr(),



                    # ── Top N Configuration ──

                    shiny::conditionalPanel(

            condition = sprintf("input['%s'] == 'mode_topN'", ns("vis_mode")),

            shiny::numericInput(

                            ns("vis_topN"),

                            label = "Top N Count:",

                            value = 5,

                            min = 3,

                            max = 50,

                            step = 1

            )

                    ),

                    shiny::numericInput(

            ns("vis_fdr"),

            label = "FDR Threshold:",

            value = 0.25,

            min = 0,

            max = 1,

            step = 0.01

                    ),

                    shiny::numericInput(

            ns("vis_min_hub"),

            label = "Min Hub Degree:",

            value = 2,

            min = 1,

            max = 20,

            step = 1

                    ),

                    shiny::hr(),



                    # ── Physics Engine Controls ──

                    shiny::h5("Physics & Interaction"),

                    shiny::checkboxInput(

            ns("vis_physics"),

            label = "Enable Physics",

            value = TRUE

                    ),

                    shiny::conditionalPanel(

            condition = sprintf("input['%s'] == true", ns("vis_physics")),

            shiny::div(

                            style = "background: #f8f9fa; padding: 10px; border-radius: 5px; margin-top: 10px;",

                            shiny::sliderInput(

                ns("vis_gravitational"),

                label = "Gravitational Constant:",

                min = -2000,

                max = -50,

                value = -400,

                step = 10

                            ),

                            shiny::sliderInput(

                ns("vis_spring_length"),

                label = "Spring Length:",

                min = 25,

                max = 500,

                value = 150,

                step = 5

                            ),

                            shiny::sliderInput(

                ns("vis_spring_constant"),

                label = "Spring Constant:",

                min = 0,

                max = 1.0,

                value = 0.001,

                step = 0.001

                            ),

                            shiny::sliderInput(

                ns("vis_central_gravity"),

                label = "Central Gravity:",

                min = 0,

                max = 1,

                value = 0.3,

                step = 0.05

                            ),

                            shiny::sliderInput(

                ns("vis_damping"),

                label = "Damping:",

                min = 0.01,

                max = 0.9,

                value = 0.09,

                step = 0.01

                            )

            )

                    ),

                    shiny::checkboxInput(

            ns("vis_nodes_draggable"),

            label = "Allow Node Dragging",

            value = TRUE

                    ),

                    shiny::checkboxInput(

            ns("vis_smooth_edges"),

            label = "Smooth Curved Edges",

            value = TRUE

                    ),

                    shiny::hr(),



                    # ── Node Sizing ──

                    shiny::h5("Node Sizing"),

                    shiny::selectInput(

            ns("vis_pw_size_mode"),

            label = "Pathway Node Size Encoding:",

            choices = c(

                            "Fixed size (slider below)"        = "fixed",

                            "By significance (-log10 FDR)"     = "fdr",

                            "By gene-set size (setSize)"       = "setsize"

            ),

            selected = "setsize"

                    ),

                    shiny::helpText(

            style = "color: #666; font-size: 11px;",

            "setSize mode matches enrichplot::cnetplot convention.",

            "Size range scales around the slider value [0.6x, 1.4x]."

                    ),

                    shiny::sliderInput(

            ns("vis_pw_size"),

            label = "Pathway Node Base Size:",

            min = 15,

            max = 60,

            value = 30,

            step = 1

                    ),

                    shiny::sliderInput(

            ns("vis_gene_size"),

            label = "Gene Node Base Size:",

            min = 5,

            max = 30,

            value = 12,

            step = 1

                    ),

                    shiny::hr(),



                    # ── Random Seed ──

                    shiny::numericInput(

            ns("vis_seed"),

            label = "Random Seed:",

            value = 42,

            min = 1,

            max = 9999,

            step = 1

                    ),

                    shiny::hr(),



                    # ── Statistics ──

                    shiny::h5("Network Statistics"),

                    shiny::verbatimTextOutput(ns("vis_stats")) |>

            shiny::tagAppendAttributes(style = "font-size: 10px; max-height: 100px; overflow-y: auto;"),

                    shiny::hr(),



                    # ── Pathway Preview ──

                    shiny::h5("Pathways Preview"),

                    shiny::uiOutput(ns("vis_pathway_list")) |>

            shiny::tagAppendAttributes(style = "max-height: 150px; overflow-y: auto;")

        )

            ),



            # ─── Main Plotting Area ───

            shiny::column(

        width = 9,

        shiny::div(

                    class = "white-box",

                    style = "min-height: 850px;",



                    # Status Bar

                    shiny::div(

            style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); padding: 12px; border-radius: 5px; margin-bottom: 15px; color: white;",

            shiny::uiOutput(ns("vis_status_text"))

                    ),



                    # visNetwork Output

                    visNetwork::visNetworkOutput(

            ns("vis_network"),

            height = "1200px",

            width = "100%"

                    ) |>

            shinycssloaders::withSpinner(type = 6, color = "#28a745"),

                    shiny::hr(),

                    shiny::div(

            style = "text-align: center;",

            shiny::actionButton(

                            ns("open_hubgene_export"),

                            label = "Export Publication Plot",

                            icon = shiny::icon("file-image"),

                            class = "btn-success",

                            style = "width: 60%;"

            )

                    ),

                    shiny::hr()

        )

            )

    )

    )

}





#' @title HubGene Network (visNetwork) Module Server

#' @description

#' Shiny module server for HubGene network visualization using visNetwork.

#' Handles pathway selection, network building, rendering with physics

#' simulation, and interactive node/edge display.

#'

#' @param id Character string used to namespace the module.

#' @param data_prep_list A reactive expression returning a list containing

#'   GSEA results, data frame, contrast groups, and metadata.

#' @param table_controller A list containing reactive expressions for

#'   table interactions, including selected pathways.

#' @param gsea_res A GSEA result object used for symbol mapping.

#'

#' @return A list containing reactive expressions:

#'   \item{final_pathways}{Reactive character vector of selected pathway IDs}

#'

#' @keywords internal



mod_hubgene_vis_server <- function(id, data_prep_list, table_controller, gsea_res) {

    shiny::moduleServer(id, function(input, output, session) {

    ns <- session$ns



    # ──────────────────────────────────────────────────────────────

    # Debug Log

    # ──────────────────────────────────────────────────────────────

    session$userData$debug_log <- character(0)



    add_debug <- function(msg, level = "INFO") {

            timestamp <- format(Sys.time(), "%H:%M:%S")

            new_entry <- sprintf("[%s] [%s] %s", timestamp, level, msg)

            session$userData$debug_log <- c(new_entry, session$userData$debug_log)[seq_len(min(50, length(session$userData$debug_log) + 1L))]

            message(sprintf("[HubGene-vis] %s", new_entry))

    }



    output$vis_debug <- shiny::renderText({

            paste(session$userData$debug_log, collapse = "\n")

    })



    add_debug("Module started")



    # ──────────────────────────────────────────────────────────────

    # 1. Mode State

    # ──────────────────────────────────────────────────────────────



    vis_mode <- shiny::reactiveVal("mode_topN")



    shiny::observeEvent(input$vis_mode, {

            vis_mode(input$vis_mode)

            add_debug(sprintf("Mode: %s", input$vis_mode))

    })



    # ──────────────────────────────────────────────────────────────

    # 2. Data Source

    # ──────────────────────────────────────────────────────────────



    topN_candidates <- shiny::reactive({

            if (vis_mode() != "mode_topN") {

        return(character(0))

            }

            data_list <- data_prep_list$data()

            shiny::req(data_list)

            df <- data_list$df

            top_n <- input$vis_topN

            if (is.null(top_n)) top_n <- 5

            top_n <- max(3, min(top_n, nrow(df)))

            df[seq_len(top_n), "ID"]

    })



    select_candidates <- shiny::reactive({

            if (vis_mode() != "mode_select") {

        return(character(0))

            }

            sel <- table_controller$selected_pathways()

            if (is.null(sel) || length(sel) == 0) {

        return(character(0))

            }

            return(sel)

    })



    candidate_raw <- shiny::reactive({

            switch(vis_mode(),

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

            fdr_thresh <- input$vis_fdr

            if (is.null(fdr_thresh)) fdr_thresh <- 0.25

            data_list <- data_prep_list$data()

            shiny::req(data_list)

            df <- data_list$df

            fdr_vec <- df$p.adjust[match(pathways, df$ID)]

            names(fdr_vec) <- pathways

            pathways[!is.na(fdr_vec) & fdr_vec < fdr_thresh]

    })



    final_pathways <- shiny::reactiveVal(character(0))



    shiny::observeEvent(candidate_filtered(), {

            final_pathways(candidate_filtered())

            add_debug(sprintf("Pathways: %d", length(candidate_filtered())))

    })



    # ──────────────────────────────────────────────────────────────

    # 3. Network Data

    # ──────────────────────────────────────────────────────────────



    net_data <- shiny::reactive({

            data_list <- data_prep_list$data()

            shiny::req(!is.null(data_list))

            pathway_ids <- final_pathways()

            shiny::req(length(pathway_ids) > 0)



            add_debug(sprintf("Building: %d pathways", length(pathway_ids)))



            net <- build_hubgene_network(

        gsea_task = list(

                    gsea_res = data_list$gsea_res,

                    meta = list(

            left_group = data_list$left_group,

            right_group = data_list$right_group

                    )

        ),

        pathway_ids = pathway_ids,

        min_hub_degree = input$vis_min_hub,

        de_df = NULL,

        res_df = data_list$df,

        seed = input$vis_seed

            )



            if (is.null(net)) {

        add_debug("Build failed", "ERROR")

        return(NULL)

            }



            add_debug(sprintf(

        "Built: %d pw + %d gene + %d edges",

        ifelse(is.null(net$nodes$pathway), 0, nrow(net$nodes$pathway)),

        ifelse(is.null(net$nodes$gene), 0, nrow(net$nodes$gene)),

        ifelse(is.null(net$edges), 0, nrow(net$edges))

            ))

            return(net)

    })



    # ──────────────────────────────────────────────────────────────

    # 4. Render Network

    # ──────────────────────────────────────────────────────────────



    output$vis_network <- visNetwork::renderVisNetwork({

            net <- net_data()

            shiny::req(!is.null(net))



            add_debug("Preparing visNetwork...")



            # Get group information

            data_list <- data_prep_list$data()

            left_group <- data_list$left_group

            right_group <- data_list$right_group



            # ──────────────────────────────────────────────────────────────

            # Fix: Rebuild gene name case mapping

            # ──────────────────────────────────────────────────────────────

            symbol_map <- NULL

            if (!is.null(gsea_res) && !is.null(data_list$contrast_id)) {

        symbol_map <- tryCatch(

                    {

            .rebuild_symbol_map(gsea_res, data_list$contrast_id)

                    },

                    error = function(e) {

            message("[HubGene-vis] symbol_map build failed: ", e$message)

            NULL

                    }

        )

        add_debug(sprintf("symbol_map: %d genes mapped", length(symbol_map)))

            } else {

        add_debug("symbol_map skipped: gsea_res or contrast_id is NULL", "WARN")

            }



            # Node size parameters

            # ------------------------------------------------------------

            # Pathway node size is now driven by a user-selectable encoding

            # (input$vis_pw_size_mode). The slider value (vis_pw_size) acts as

            # the BASE size; the chosen encoding scales around it within

            # [0.6x, 1.4x] to keep visNetwork's force-directed layout stable.

            #   "fixed"   -> constant (legacy behavior)

            #   "fdr"     -> -log10(FDR) normalized to [0.6x, 1.4x] of base

            #   "setsize" -> sqrt(setSize) normalized to [0.6x, 1.4x] of base

            #                (cnetplot convention; sqrt avoids over-weighting

            #                 large gene sets)

            # ------------------------------------------------------------

            pw_base  <- ifelse(is.null(input$vis_pw_size),  30, input$vis_pw_size)

            gene_size <- ifelse(is.null(input$vis_gene_size), 12, input$vis_gene_size)

            pw_size_mode <- if (is.null(input$vis_pw_size_mode)) "setsize"

                                            else input$vis_pw_size_mode



            # Per-pathway scaled size vector (matches row order of net$nodes$pathway)

            #

            # IRON FIX (2026-06-30): build_hubgene_network stores the gene-set

            # size under column "n_core" (NOT "setSize"). The previous code

            # unconditionally accessed net$nodes$pathway$setSize, which returned

            # NULL, breaking sqrt()/range() with

            # "non-numeric argument to mathematical function".

            # Now: detect column name setSize -> n_core -> n_total -> fallback,

            # and coerce all numeric inputs via as.numeric() to defend against

            # character-typed columns from upstream GSEA objects.

            pw_pathway <- net$nodes$pathway

            pw_nrow <- if (is.null(pw_pathway)) 0 else nrow(pw_pathway)



            pw_size <- switch(pw_size_mode,

        "fdr" = {

                    if (pw_nrow > 0) {

            fdr_vals <- suppressWarnings(as.numeric(pw_pathway$FDR))

            fdr_vals[is.na(fdr_vals)] <- 1

            fdr_nl <- -log10(pmax(fdr_vals, 1e-300))

            r <- range(fdr_nl, na.rm = TRUE); if (diff(r) < 1e-6) r <- c(0, 10)

            pw_base * (0.6 + 0.8 * (fdr_nl - r[1]) / diff(r))

                    } else rep(pw_base, 0)

        },

        "setsize" = {

                    if (pw_nrow > 0) {

            # Column-name fallback chain

            sv <- if ("setSize" %in% colnames(pw_pathway)) {

                            suppressWarnings(as.numeric(pw_pathway$setSize))

            } else if ("n_core" %in% colnames(pw_pathway)) {

                            suppressWarnings(as.numeric(pw_pathway$n_core))

            } else if ("n_total" %in% colnames(pw_pathway)) {

                            suppressWarnings(as.numeric(pw_pathway$n_total))

            } else {

                            rep(100, pw_nrow)

            }

            sv[is.na(sv)] <- 100

            r <- range(sv, na.rm = TRUE); if (diff(r) < 1) r <- c(0, 200)

            pw_base * (0.6 + 0.8 * (sqrt(sv) - sqrt(r[1])) /

                                                                    (sqrt(r[2]) - sqrt(r[1])))

                    } else rep(pw_base, 0)

        },

        rep(pw_base, pw_nrow)

            )



            # ──────────────────────────────────────────────────────────────

            # Node Construction (add opacity column)

            # ──────────────────────────────────────────────────────────────



            nodes <- data.frame(

        id = character(),

        label = character(),

        group = character(),

        shape = character(),

        color = character(),

        size = numeric(),

        opacity = numeric(),

        title = character(),

        stringsAsFactors = FALSE

            )



            # Pathway nodes

            if (!is.null(net$nodes$pathway) && nrow(net$nodes$pathway) > 0) {

        # IRON FIX (2026-06-30): build_hubgene_network stores the gene-set

        # size under "n_core" (NOT "setSize"). Accessing row$setSize inside

        # the per-row sprintf returned NULL, which made sprintf produce a

        # zero-length character vector, which in turn made data.frame()

        # fail with "arguments imply differing number of rows: 1, 0".

        # Resolve the size column name once, outside the loop, then index

        # safely with a numeric fallback.

        pw_size_col <- if ("setSize" %in% colnames(net$nodes$pathway)) "setSize"

                                                else if ("n_core" %in% colnames(net$nodes$pathway)) "n_core"

                                                else if ("n_total" %in% colnames(net$nodes$pathway)) "n_total"

                                                else NA_character_

        pw_NES_col  <- if ("NES" %in% colnames(net$nodes$pathway)) "NES" else NA_character_

        pw_FDR_col  <- if ("FDR" %in% colnames(net$nodes$pathway)) "FDR" else NA_character_



        for (i in seq_len(nrow(net$nodes$pathway))) {

                    row <- net$nodes$pathway[i, ]

                    # Safe scalar extraction with NA fallback

                    row_nes <- if (is.na(pw_NES_col)) NA_real_

                                            else suppressWarnings(as.numeric(row[[pw_NES_col]][1]))

                    row_fdr <- if (is.na(pw_FDR_col)) NA_real_

                                            else suppressWarnings(as.numeric(row[[pw_FDR_col]][1]))

                    row_size <- if (is.na(pw_size_col)) NA_integer_

                                            else suppressWarnings(as.integer(row[[pw_size_col]][1]))

                    row_nes[is.na(row_nes)]  <- 0

                    row_fdr[is.na(row_fdr)]  <- 1

                    row_size[is.na(row_size)] <- 0



                    direction <- ifelse(row_nes > 0, left_group, right_group)



                    # Calculate LE connection count

                    le_count <- 0

                    if (!is.null(net$edges) && nrow(net$edges) > 0) {

            le_count <- sum(net$edges$target == row$id & net$edges$is_leading_edge, na.rm = TRUE)

                    }



                    nodes <- rbind(nodes, data.frame(

            id = paste0("pw_", row$id),

            label = gsub("^[^_]*_", "", row$id),

            group = "pathway",

            shape = "diamond",

            color = ifelse(row_nes > 0, "#E41A1C", "#377EB8"),

            size = pw_size[i],

            opacity = 1.0,

            title = sprintf(

                            "<b>%s</b><br>Direction: %s<br>NES: %.3f<br>FDR: %.2e<br>LE Genes: %d<br>Set Size: %d",

                            row$id, direction, row_nes, row_fdr, le_count, row_size

            ),

            stringsAsFactors = FALSE

                    ))

        }

            }



            # Gene nodes (add opacity for distinction)

            if (!is.null(net$nodes$gene) && nrow(net$nodes$gene) > 0) {

        # Resolve gene column names defensively (mirrors the pathway fix

        # above: some upstream paths may rename or drop columns).

        gn_stat_col   <- if ("stat"   %in% colnames(net$nodes$gene)) "stat"   else NA_character_

        gn_degree_col <- if ("degree" %in% colnames(net$nodes$gene)) "degree" else NA_character_



        for (i in seq_len(nrow(net$nodes$gene))) {

                    row <- net$nodes$gene[i, ]

                    # Safe scalar extraction

                    row_stat <- if (is.na(gn_stat_col)) NA_real_

                                            else suppressWarnings(as.numeric(row[[gn_stat_col]][1]))

                    row_degree <- if (is.na(gn_degree_col)) NA_integer_

                        else suppressWarnings(as.integer(row[[gn_degree_col]][1]))

                    row_stat[is.na(row_stat)]    <- 0

                    row_degree[is.na(row_degree)] <- 1



                    # Restore original case

                    gene_label <- .get_display_symbol(row$id, symbol_map)



                    # Check if Leading Edge

                    is_le <- FALSE

                    pw_str <- "None"

                    if (!is.null(net$edges) && nrow(net$edges) > 0) {

            is_le <- any(net$edges$source == row$id & net$edges$is_leading_edge, na.rm = TRUE)

            connected_pws <- unique(net$edges$target[net$edges$source == row$id])

            if (length(connected_pws) > 0) {

                            pw_labels <- gsub("^[^_]*_", "", connected_pws)

                            if (length(pw_labels) <= 5) {

                pw_str <- paste(pw_labels, collapse = ", ")

                            } else {

                pw_str <- paste(pw_labels[seq_len(5)], collapse = ", ")

                pw_str <- paste0(pw_str, sprintf(" (+%d)", length(pw_labels) - 5))

                            }

            }

                    }



                    # Direction

                    if (row_stat > 0) {

            direction <- sprintf("Up in %s", left_group)

                    } else if (row_stat < 0) {

            direction <- sprintf("Up in %s", right_group)

                    } else {

            direction <- "Neutral"

                    }



                    # Opacity: LE=1.0, non-LE=0.5

                    gene_opacity <- ifelse(is_le, 1.0, 0.5)



                    nodes <- rbind(nodes, data.frame(

            id = paste0("gene_", gene_label),

            label = gene_label,

            group = "gene",

            shape = "dot",

            color = ifelse(row_stat > 0, "#E41A1C", ifelse(row_stat < 0, "#377EB8", "#999999")),

            size = gene_size + row_degree * 3,

            opacity = gene_opacity,

            title = sprintf(

                            "<b>%s</b><br>Direction: %s<br>Stat: %.3f<br>Hub Degree: %d<br>LE: %s<br>Connected: %s",

                            gene_label, direction, row_stat, row_degree,

                            ifelse(is_le, "YES", "NO"), pw_str

            ),

            stringsAsFactors = FALSE

                    ))

        }

            }



            # ──────────────────────────────────────────────────────────────

            # Edge Construction (color follows gene stat)

            # ──────────────────────────────────────────────────────────────



            edges <- data.frame(

        from = character(),

        to = character(),

        color = character(),

        width = numeric(),

        dashes = logical(),

        title = character(),

        stringsAsFactors = FALSE

            )



            if (!is.null(net$edges) && nrow(net$edges) > 0) {

        # Resolve edge column names defensively.

        ed_source_col <- if ("source" %in% colnames(net$edges)) "source" else NA_character_

        ed_target_col <- if ("target" %in% colnames(net$edges)) "target" else NA_character_

        ed_le_col     <- if ("is_leading_edge" %in% colnames(net$edges)) "is_leading_edge"

                                                    else NA_character_



        for (i in seq_len(nrow(net$edges))) {

                    row <- net$edges[i, ]



                    ed_source <- if (is.na(ed_source_col)) ""

                                                else as.character(row[[ed_source_col]][1])

                    ed_target <- if (is.na(ed_target_col)) ""

                                                else as.character(row[[ed_target_col]][1])

                    ed_is_le  <- if (is.na(ed_le_col)) FALSE

                                                else isTRUE(as.logical(row[[ed_le_col]][1]))



                    # Get source gene stat value to determine color

                    gene_stat <- 0

                    if (!is.null(net$nodes$gene) && nrow(net$nodes$gene) > 0) {

            if (!is.na(gn_stat_col)) {

                            gene_row <- net$nodes$gene[net$nodes$gene$id == ed_source, ]

                            if (nrow(gene_row) > 0) {

                gene_stat <- suppressWarnings(as.numeric(gene_row[[gn_stat_col]][1]))

                gene_stat <- if (is.na(gene_stat)) 0 else gene_stat

                            }

            }

                    }



                    # Color follows gene stat: light color (with opacity)

                    if (gene_stat > 0) {

            edge_color <- "#E41A1C80"

                    } else if (gene_stat < 0) {

            edge_color <- "#377EB880"

                    } else {

            edge_color <- "#99999980"

                    }



                    # Restore original case

                    gene_display <- .get_display_symbol(ed_source, symbol_map)



                    edges <- rbind(edges, data.frame(

            from = paste0("gene_", gene_display),

            to = paste0("pw_", ed_target),

            color = edge_color,

            width = ifelse(ed_is_le, 3, 1.5),

            dashes = !ed_is_le,

            title = sprintf(

                            "<b>Gene:</b> %s<br><b>Pathway:</b> %s<br><b>Leading Edge:</b> %s",

                            gene_display, ed_target,

                            ifelse(ed_is_le, "YES", "NO")

            ),

            stringsAsFactors = FALSE

                    ))

        }

            }



            add_debug(sprintf("Nodes: %d, Edges: %d", nrow(nodes), nrow(edges)))



            # ──────────────────────────────────────────────────────────────

            # Build visNetwork Object

            # ──────────────────────────────────────────────────────────────



            vis <- visNetwork::visNetwork(nodes, edges, width = "100%", height = "650px")



            # Physics Engine

            physics_enabled <- isTRUE(input$vis_physics)

            grav <- ifelse(is.null(input$vis_gravitational), -400, input$vis_gravitational)



            spring_length <- ifelse(is.null(input$vis_spring_length), 150, input$vis_spring_length)

            spring_constant <- ifelse(is.null(input$vis_spring_constant), 0.01, input$vis_spring_constant)

            damping <- ifelse(is.null(input$vis_damping), 0.09, input$vis_damping)

            central_gravity <- ifelse(is.null(input$vis_central_gravity), 0.3, input$vis_central_gravity)



            if (physics_enabled) {

        vis <- vis |> visNetwork::visPhysics(

                    enabled = TRUE,

                    solver = "barnesHut",

                    barnesHut = list(

            gravitationalConstant = grav,

            centralGravity = central_gravity,

            springLength = spring_length,

            springConstant = spring_constant,

            damping = damping

                    )

        )

        add_debug(sprintf("Physics: ON, solver=barnesHut, grav=%.0f", grav))

            } else {

        vis <- vis |> visNetwork::visPhysics(enabled = FALSE)

        add_debug("Physics: OFF")

            }



            # Layout

            seed <- ifelse(is.null(input$vis_seed), 42, input$vis_seed)

            vis <- vis |> visNetwork::visLayout(randomSeed = seed)



            # Interaction

            vis <- vis |> visNetwork::visInteraction(

        dragNodes = TRUE,

        dragView = TRUE,

        zoomView = TRUE,

        hover = TRUE,

        tooltipDelay = 300,

        navigationButtons = TRUE

            )



            # Click to highlight related nodes and edges

            vis <- vis |> visNetwork::visOptions(

        highlightNearest = TRUE,

        nodesIdSelection = TRUE

            )



            # Click Event

            vis <- vis |> visNetwork::visEvents(

        click = sprintf("function(props) {

            var n = props.nodes[0];

            if(n) Shiny.setInputValue('%s', {id: n, ts: Date.now()}, {priority:'event'});

    }", ns("vis_click"))

            )



            add_debug("Render complete")

            vis

    })



    # ──────────────────────────────────────────────────────────────

    # 5. Events

    # ──────────────────────────────────────────────────────────────



    shiny::observeEvent(input$vis_click, {

            add_debug(sprintf("Click: %s", input$vis_click$id), "CLICK")

    })



    # ──────────────────────────────────────────────────────────────

    # 6. Statistics

    # ──────────────────────────────────────────────────────────────



    output$vis_stats <- shiny::renderPrint({

            net <- tryCatch(net_data(), error = function(e) NULL)

            if (is.null(net)) {

        cat("No data\n")

        return()

            }

            pw <- ifelse(is.null(net$nodes$pathway), 0, nrow(net$nodes$pathway))

            gene <- ifelse(is.null(net$nodes$gene), 0, nrow(net$nodes$gene))

            edge <- ifelse(is.null(net$edges), 0, nrow(net$edges))

            le_edges <- ifelse(is.null(net$edges), 0, sum(net$edges$is_leading_edge, na.rm = TRUE))

            cat(sprintf("Pathways: %d\n", pw))

            cat(sprintf("Hub Genes: %d\n", gene))

            cat(sprintf("Total Edges: %d\n", edge))

            cat(sprintf("LE Edges: %d\n", le_edges))

    })



    # ──────────────────────────────────────────────────────────────

    # 7. Status Bar

    # ──────────────────────────────────────────────────────────────



    output$vis_status_text <- shiny::renderUI({

            net <- tryCatch(net_data(), error = function(e) NULL)

            if (is.null(net)) {

        return(shiny::span("No data", style = "color: #ffcccc;"))

            }

            pw <- ifelse(is.null(net$nodes$pathway), 0, nrow(net$nodes$pathway))

            gene <- ifelse(is.null(net$nodes$gene), 0, nrow(net$nodes$gene))

            edge <- ifelse(is.null(net$edges), 0, nrow(net$edges))

            le_edges <- ifelse(is.null(net$edges), 0, sum(net$edges$is_leading_edge, na.rm = TRUE))

            physics <- if (isTRUE(input$vis_physics)) "ON" else "OFF"



            shiny::tagList(

        shiny::strong(sprintf(

                    "HubGene Network | %d Pathways + %d Hub Genes + %d Edges (LE: %d)",

                    pw, gene, edge, le_edges

        )),

        htmltools::tags$br(),

        htmltools::tags$small(sprintf(

                    "Physics: %s (Barnes-Hut) | FDR < %.2f | tooltipDelay: 300ms ",

                    physics, input$vis_fdr %||% 0.25

        ))

            )

    })



    # ──────────────────────────────────────────────────────────────

    # 8. Pathway Preview

    # ──────────────────────────────────────────────────────────────



    output$vis_pathway_list <- shiny::renderUI({

            pathways <- final_pathways()

            net <- tryCatch(net_data(), error = function(e) NULL)



            if (length(pathways) == 0) {

        return(shiny::div("No pathways", style = "color: #856404;"))

            }

            data_list <- data_prep_list$data()

            df <- data_list$df



            lapply(pathways, function(pid) {

        row_idx <- which(df$ID == pid)

        if (length(row_idx) == 0) {

                    return(NULL)

        }

        row <- df[row_idx[1], ]

        le_count <- 0

        if (!is.null(net) && !is.null(net$edges)) {

                    le_count <- sum(net$edges$target == pid & net$edges$is_leading_edge, na.rm = TRUE)

        }

        shiny::div(

                    style = "background: #f8f9fa; padding: 5px; margin-bottom: 3px; border-left: 3px solid #007bff; font-size: 10px;",

                    shiny::strong(gsub("^[^_]*_", "", pid)),

                    htmltools::tags$br(),

                    sprintf("NES: %.2f | LE: %d", as.numeric(row$NES), le_count)

        )

            })

    })



    # ──────────────────────────────────────────────────────────────

    # 9. Detail Tables

    # ──────────────────────────────────────────────────────────────



    output$vis_node_detail <- DT::renderDataTable({

            net <- tryCatch(net_data(), error = function(e) NULL)

            if (is.null(net) || (is.null(net$nodes$pathway) && is.null(net$nodes$gene))) {

        return(DT::datatable(data.frame(Message = "No nodes"), options = list(dom = "t")))

            }



            # Get symbol_map for displaying original case

            symbol_map <- NULL

            data_list <- tryCatch(data_prep_list$data(), error = function(e) NULL)

            if (!is.null(data_list) && !is.null(data_prep_list$gsea_res) && !is.null(data_list$contrast_id)) {

        symbol_map <- tryCatch(

                    {

            .rebuild_symbol_map(data_prep_list$gsea_res, data_list$contrast_id)

                    },

                    error = function(e) NULL

        )

            }



            rows <- list()



            # Pathway nodes

            if (!is.null(net$nodes$pathway) && nrow(net$nodes$pathway) > 0) {

        for (i in seq_len(nrow(net$nodes$pathway))) {

                    r <- net$nodes$pathway[i, ]

                    le_count <- sum(net$edges$target == r$id & net$edges$is_leading_edge, na.rm = TRUE)

                    total_count <- sum(net$edges$target == r$id, na.rm = TRUE)

                    rows[[length(rows) + 1]] <- data.frame(

            ID = r$id,

            Type = "Pathway",

            NES = round(r$NES, 3),

            Direction = ifelse(r$NES > 0, "Up", "Down"),

            Degree = NA_integer_,

            LE = le_count,

            Total = total_count,

            stringsAsFactors = FALSE

                    )

        }

            }



            # Gene nodes

            if (!is.null(net$nodes$gene) && nrow(net$nodes$gene) > 0) {

        for (i in seq_len(nrow(net$nodes$gene))) {

                    r <- net$nodes$gene[i, ]

                    is_le <- any(net$edges$source == r$id & net$edges$is_leading_edge, na.rm = TRUE)

                    gene_display <- .get_display_symbol(r$id, symbol_map)

                    rows[[length(rows) + 1]] <- data.frame(

            ID = gene_display,

            Type = "Gene",

            NES = round(r$stat, 3),

            Direction = ifelse(r$stat > 0, "Up", ifelse(r$stat < 0, "Down", "Neutral")),

            Degree = r$degree,

            LE = ifelse(is_le, "YES", "NO"),

            Total = NA_integer_,

            stringsAsFactors = FALSE

                    )

        }

            }



            df <- do.call(rbind, rows)

            DT::datatable(df, rownames = FALSE, options = list(pageLength = 10, dom = "t"))

    })



    output$vis_edge_detail <- DT::renderDataTable({

            net <- tryCatch(net_data(), error = function(e) NULL)

            if (is.null(net) || is.null(net$edges) || nrow(net$edges) == 0) {

        return(DT::datatable(data.frame(Message = "No edges"), options = list(dom = "t")))

            }



            # Get symbol_map for displaying original case

            symbol_map <- NULL

            data_list <- tryCatch(data_prep_list$data(), error = function(e) NULL)

            if (!is.null(data_list) && !is.null(data_prep_list$gsea_res) && !is.null(data_list$contrast_id)) {

        symbol_map <- tryCatch(

                    {

            .rebuild_symbol_map(data_prep_list$gsea_res, data_list$contrast_id)

                    },

                    error = function(e) NULL

        )

            }



            rows <- lapply(seq_len(nrow(net$edges)), function(i) {

        r <- net$edges[i, ]

        gene_display <- .get_display_symbol(r$source, symbol_map)

        data.frame(

                    Gene = gene_display,

                    Pathway = r$target,

                    LeadingEdge = ifelse(r$is_leading_edge, "YES", "NO"),

                    Width = ifelse(r$is_leading_edge, 3, 1.5),

                    Style = ifelse(r$is_leading_edge, "Solid", "Dashed"),

                    stringsAsFactors = FALSE

        )

            })



            df <- do.call(rbind, rows)

            DT::datatable(df, rownames = FALSE, options = list(pageLength = 10, dom = "t"))

    })



    add_debug("Ready")



    # ============================================================

    # HubGene Export Center

    # ------------------------------------------------------------

    # Static ggplot2 reproduction of the visNetwork bipartite graph.

    # Reuses generate_hubgene_code() from 15_code_generator.R.

    # ============================================================



    .build_hubgene_export_modal <- function() {

            shiny::modalDialog(

        title = "Export Publication Plot - HubGene Network",

        size = "l", footer = NULL, easyClose = TRUE,

        shiny::fluidRow(

                    shiny::column(5,

            shiny::h5("Dimensions"),

            shiny::fluidRow(

                            shiny::column(6,

                shiny::numericInput(ns("hub_exp_width"),  "Width (inch)",  value = 10, min = 2, max = 24)

                            ),

                            shiny::column(6,

                shiny::numericInput(ns("hub_exp_height"), "Height (inch)", value = 8, min = 2, max = 24)

                            )

            ),

            shiny::numericInput(ns("hub_exp_dpi"), "DPI", value = 300, min = 72, max = 600),



            shiny::h6("Canvas Margin (pt)"),

            shiny::fluidRow(

                            shiny::column(3, shiny::numericInput(ns("hub_exp_margin_top"),    "Top",    value = 18, min = 0, max = 80)),

                            shiny::column(3, shiny::numericInput(ns("hub_exp_margin_bottom"), "Bottom", value = 18, min = 0, max = 80)),

                            shiny::column(3, shiny::numericInput(ns("hub_exp_margin_left"),   "Left",   value = 18, min = 0, max = 80)),

                            shiny::column(3, shiny::numericInput(ns("hub_exp_margin_right"),  "Right",  value = 18, min = 0, max = 80))

            ),

            shiny::hr(),

            shiny::h5("Download"),

            shiny::fluidRow(

                            shiny::column(6,

                shiny::downloadButton(ns("hub_exp_download_pdf"),

                                    "Download PDF",

                                    class = "btn-danger btn-block",

                                    icon  = shiny::icon("file-pdf"))

                            ),

                            shiny::column(6,

                shiny::downloadButton(ns("hub_exp_download_png"),

                                    "Download PNG",

                                    class = "btn-primary btn-block",

                                    icon  = shiny::icon("file-image"))

                            )

            ),

            shiny::hr(),

            shiny::fluidRow(

                            shiny::column(6,

                shiny::selectInput(ns("hub_exp_format"), "Other format",

                                    choices = c("SVG" = "svg", "TIFF" = "tiff"),

                                    selected = "svg")

                            ),

                            shiny::column(6,

                shiny::div(style = "margin-top: 22px;",

                                    shiny::downloadButton(ns("hub_exp_download_other"),

                    "Download", class = "btn-default btn-block")

                )

                            )

            ),

            shiny::hr(),

            shiny::fluidRow(

                            shiny::column(6,

                shiny::actionButton(ns("hub_exp_copy_code"),

                                    "Copy R Code",

                                    class = "btn-info btn-block",

                                    icon = shiny::icon("clipboard"))

                            ),

                            shiny::column(6,

                shiny::actionButton(ns("hub_exp_dismiss"),

                                    "Close",

                                    class = "btn-default btn-block")

                            )

            )

                    ),

                    shiny::column(7,

            shiny::h5("Live Preview"),

            shiny::div(style = "background:#f5f5f5; border:1px solid #ddd; border-radius:4px; padding:8px; width:100%; height:520px; display:flex; align-items:center; justify-content:center; overflow:auto;",

                            shiny::plotOutput(ns("hub_exp_preview")) |>

                shinycssloaders::withSpinner(type = 6, color = "#28a745")

            ),

            shiny::tags$small(style = "color: #666; display:block; margin-top:6px;",

                            "Pathway nodes: diamonds; Gene nodes: circles. Preview auto-scaled to fit while keeping the export aspect ratio.")

                    )

        )

            )

    }



    shiny::observeEvent(input$open_hubgene_export, {

            net <- tryCatch(net_data(), error = function(e) NULL)

            if (is.null(net) || (is.null(net$nodes$pathway) && is.null(net$nodes$gene))) {

        shiny::showNotification("No network to export.", type = "warning"); return()

            }

            shiny::showModal(.build_hubgene_export_modal())

    })

    shiny::observeEvent(input$hub_exp_dismiss, shiny::removeModal())



    # Live preview reactive

    # IRON FIX (2026-06-30): same pattern as pathway/quadrant modules —

    # suppress `print(p)` from opening a separate graphics device, surface

    # errors as notifications, and guard against NA from cleared numericInput.

    #

    # Preview scaling (2026-07-06): the preview pane is a FIXED-SIZE window

    # (~500x480 px). The figure is rendered at the export aspect ratio and

    # scaled to fit inside that window, centered. A tall figure appears as a

    # vertical "pole" thumbnail; a wide figure as a horizontal "pole".

    # Preview device matching export physical size (2026-07-06 v3):
    # Shiny renderPlot does NOT accept res as a function, so res is fixed
    # at 72 and pixel dims are computed as inch * 72 (with a pixel cap).
    # Physical size = width_px / 72 = w_in inches, matching the export.
    .hub_preview_dims <- function(w_in, h_in, res = 72, max_px = 1600) {

            w_in <- suppressWarnings(as.numeric(w_in[1]))

            h_in <- suppressWarnings(as.numeric(h_in[1]))

            if (length(w_in) == 0 || is.na(w_in) || w_in <= 0) w_in <- 10

            if (length(h_in) == 0 || is.na(h_in) || h_in <= 0) h_in <- 8

            scale <- min(1, max_px / (max(w_in, h_in) * res))

            list(width  = round(w_in * res * scale),

                        height = round(h_in * res * scale))

    }

    .hubgene_preview_plot <- shiny::reactive({

            code <- tryCatch(.hubgene_export_code(),

        error = function(e) {

                    shiny::showNotification(

            sprintf("[hub preview] code generation failed: %s", e$message),

            type = "error", duration = 8)

                    ""

        })

            if (!nzchar(code)) return(NULL)

            eval_env <- new.env(parent = globalenv())

            eval_env$p <- NULL

            eval_env$print <- base::invisible
            eval_env$gsea_res <- gsea_res

            tryCatch(eval(parse(text = code), envir = eval_env),

                                error = function(e) {

                                    shiny::showNotification(

                                        sprintf("[hub preview] %s", e$message),

                                        type = "error", duration = 8)

                                    NULL

                                })

            eval_env$p

    })

    output$hub_exp_preview <- shiny::renderPlot(

            {

        p <- .hubgene_preview_plot()

        shiny::req(p)

        p

            },

            res  = 72,

            width  = function() {

        d <- .hub_preview_dims(input$hub_exp_width, input$hub_exp_height); d$width

            },

            height = function() {

        d <- .hub_preview_dims(input$hub_exp_width, input$hub_exp_height); d$height

            }

    )



    .hubgene_export_code <- function() {

            data_list <- data_prep_list$data()

            if (is.null(data_list)) return("")

            seed_val <- if (is.null(input$vis_seed)) 42L else as.integer(input$vis_seed)

            pw_mode  <- if (is.null(input$vis_pw_size_mode)) "setsize" else input$vis_pw_size_mode

            generate_hubgene_code(gsea_res_var = "gsea_res",

                            contrast_id  = data_list$contrast_id,

                            pathways     = final_pathways(),

                            min_hub_degree = if (is.null(input$vis_min_hub)) 2 else input$vis_min_hub,

                            pw_size_mode = pw_mode,

                            seed = seed_val,

                            target_collection = if (is.null(data_list$collections)) "ALL" else data_list$collections,

                            margin_top    = if (is.null(input$hub_exp_margin_top)    || is.na(input$hub_exp_margin_top))    18 else input$hub_exp_margin_top,

                            margin_bottom = if (is.null(input$hub_exp_margin_bottom) || is.na(input$hub_exp_margin_bottom)) 18 else input$hub_exp_margin_bottom,

                            margin_left   = if (is.null(input$hub_exp_margin_left)   || is.na(input$hub_exp_margin_left))   18 else input$hub_exp_margin_left,

                            margin_right  = if (is.null(input$hub_exp_margin_right)  || is.na(input$hub_exp_margin_right))  18 else input$hub_exp_margin_right)

    }



    .hubgene_render_to_file <- function(file, fmt) {

            code <- .hubgene_export_code()

            if (!nzchar(code)) return()

            # parent=globalenv() so lexical lookup resolves attached packages.

            # (baseenv() reproduces the v0.99.11 "could not find function" bug +

            # the .htm fallback symptom; see IRON FIX comment in commit log.)

            eval_env <- new.env(parent = globalenv())

            eval_env$p <- NULL

            eval_env$print <- base::invisible   # suppress popup device; ggsave draws itself
            eval_env$gsea_res <- gsea_res

            tryCatch(eval(parse(text = code), envir = eval_env),

                                error = function(e) shiny::showNotification(

                                    sprintf("Render failed: %s", e$message), type = "error", duration = 8))

            if (is.null(eval_env$p)) return()

            ggplot2::ggsave(file, eval_env$p,

                                            width  = if (is.null(input$hub_exp_width))  10 else input$hub_exp_width,

                                            height = if (is.null(input$hub_exp_height))  8 else input$hub_exp_height,

                                            dpi    = if (is.null(input$hub_exp_dpi))   300 else input$hub_exp_dpi,

                                            device = fmt)

    }



    output$hub_exp_download_pdf <- shiny::downloadHandler(

            filename = function() sprintf("GSEAlens_hubgene_%s.pdf", Sys.Date()),

            content  = function(file) .hubgene_render_to_file(file, "pdf")

    )



    output$hub_exp_download_png <- shiny::downloadHandler(

            filename = function() sprintf("GSEAlens_hubgene_%s.png", Sys.Date()),

            content  = function(file) .hubgene_render_to_file(file, "png")

    )



    output$hub_exp_download_other <- shiny::downloadHandler(

            filename = function() sprintf("GSEAlens_hubgene_%s.%s", Sys.Date(),

                                    if (is.null(input$hub_exp_format)) "svg" else input$hub_exp_format),

            content  = function(file) .hubgene_render_to_file(file, input$hub_exp_format)

    )



    shiny::observeEvent(input$hub_exp_copy_code, {

            code <- .hubgene_export_code()

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



    return(list(final_pathways = final_pathways))

    })

}

