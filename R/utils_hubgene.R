# =============================================================================
# HubGene Network Utility Functions
# =============================================================================
#' @title Extract Hub Genes from Selected Pathways
#' @description Extract hub genes that appear in multiple selected pathways,
#'   with their connection information and direction data.
#' @param gsea_task GseaTask object
#' @param pathway_ids Character vector of pathway IDs to analyze
#' @param min_degree Minimum number of pathways a gene must appear in (default: 2)
#' @param de_df Optional differential expression data frame with logFC
#' @return data.frame with columns: gene, degree, pathways, is_hub, log2FC
#' @export
#' @examples
#' \dontrun{
#' hub_genes <- extract_hub_genes(task, c("HALLMARK_APOPTOSIS", "KEGG_CELL_CYCLE"))
#' }
extract_hub_genes <- function(gsea_task, pathway_ids, min_degree = 2, de_df = NULL) {

  # Validate inputs
  if (is.null(pathway_ids) || length(pathway_ids) == 0) {
    return(NULL)
  }

  # Extract gene sets from GSEA result
  gene_sets <- gsea_task$gsea_res@geneSets[pathway_ids]

  if (length(gene_sets) == 0) {
    return(NULL)
  }

  # Count gene occurrences across pathways
  all_genes <- unlist(gene_sets)
  gene_counts <- table(toupper(all_genes))

  # Create result data frame
  gene_info <- data.frame(
    gene = names(gene_counts),
    degree = as.integer(gene_counts),
    stringsAsFactors = FALSE
  )

  # Determine which pathways each gene appears in
  gene_info$pathways <- sapply(gene_info$gene, function(g) {
    g_upper <- toupper(g)
    pws <- names(gene_sets)[sapply(gene_sets, function(x) g_upper %in% toupper(x))]
    paste(pws, collapse = ", ")
  })

  # Mark hub genes
  gene_info$is_hub <- gene_info$degree >= min_degree

  # Add log2FC if provided
  if (!is.null(de_df) && "logFC" %in% colnames(de_df)) {
    fc_map <- setNames(de_df$logFC, toupper(de_df$gene_symbol))
    gene_info$log2FC <- fc_map[gene_info$gene]
    gene_info$log2FC[is.na(gene_info$log2FC)] <- 0
  } else {
    gene_info$log2FC <- NA_real_
  }

  # Sort by degree descending
  gene_info <- gene_info[order(gene_info$degree, decreasing = TRUE), ]

  return(gene_info)
}


#' @title Build HubGene Network Data Structure
#' @description Build complete network data (nodes and edges) for HubGene Network
#' @param gsea_task GseaTask object
#' @param pathway_ids Selected pathway IDs
#' @param min_hub_degree Minimum degree for hub genes
#' @param de_df Optional DE data frame (will extract stat if no log2FC)
#' @param res_df Optional pathway result data frame
#' @return list with nodes and edges data frames
#' @export
build_hubgene_network <- function(gsea_task, pathway_ids,
                                  min_hub_degree = 2,
                                  de_df = NULL,
                                  res_df = NULL) {

  # Extract hub genes
  hub_df <- extract_hub_genes(gsea_task, pathway_ids, min_hub_degree, de_df)

  if (is.null(hub_df) || nrow(hub_df) == 0) {
    return(list(nodes = NULL, edges = NULL))
  }

  # Extract gene sets
  gene_sets <- gsea_task$gsea_res@geneSets[pathway_ids]

  # Get pathway result info if not provided
  if (is.null(res_df)) {
    res_df <- as.data.frame(gsea_task$gsea_res@result)
  }

  # Build pathway nodes
  pathway_nodes <- data.frame(
    id = pathway_ids,
    type = "pathway",
    stringsAsFactors = FALSE
  )

  # Add pathway statistics
  pw_stats <- res_df[match(pathway_ids, res_df$ID), ]
  pathway_nodes$NES <- pw_stats$NES
  pathway_nodes$FDR <- pw_stats$p.adjust
  pathway_nodes$pvalue <- pw_stats$pvalue
  pathway_nodes$n_core <- pw_stats$setSize
  pathway_nodes$direction <- ifelse(pathway_nodes$NES > 0, "up", "down")

  # Count total genes per pathway
  pathway_nodes$n_total <- sapply(pathway_ids, function(pw) {
    length(gene_sets[[pw]])
  })

  # Calculate core ratio
  pathway_nodes$core_ratio <- pathway_nodes$n_core / pathway_nodes$n_total

  # ─── 修改：使用 stat 值代替 log2FC ───

  # 如果没有提供 de_df，尝试从 gsea_task 获取
  if (is.null(de_df)) {
    # 尝试从 gsea_res@geneList 获取排序信息
    gene_list <- gsea_task$gsea_res@geneList
    if (!is.null(gene_list)) {
      stat_vec <- as.numeric(gene_list)
      names(stat_vec) <- names(gene_list)

      # 为 hub_df 添加 stat 值
      hub_df$stat <- sapply(hub_df$gene, function(g) {
        idx <- which(toupper(names(stat_vec)) == toupper(g))
        if (length(idx) > 0) stat_vec[idx[1]] else NA_real_
      })

      # 如果 stat 也没有，使用 gene degree 作为方向（不太理想但至少能用）
      if (all(is.na(hub_df$stat))) {
        hub_df$stat <- hub_df$degree
      }
    } else {
      hub_df$stat <- hub_df$degree
    }
  } else {
    # 从 de_df 获取 stat 值
    if ("stat" %in% colnames(de_df)) {
      stat_map <- setNames(de_df$stat, toupper(de_df$gene_symbol))
      hub_df$stat <- stat_map[hub_df$gene]
    } else if ("logFC" %in% colnames(de_df)) {
      stat_map <- setNames(de_df$logFC, toupper(de_df$gene_symbol))
      hub_df$stat <- stat_map[hub_df$gene]
    } else if ("log2FoldChange" %in% colnames(de_df)) {
      stat_map <- setNames(de_df$log2FoldChange, toupper(de_df$gene_symbol))
      hub_df$stat <- stat_map[hub_df$gene]
    } else {
      # 最后尝试 t 统计量
      stat_col <- grep("^t$|^t\\.", colnames(de_df), value = TRUE, ignore.case = TRUE)[1]
      if (!is.na(stat_col)) {
        stat_map <- setNames(de_df[[stat_col]], toupper(de_df$gene_symbol))
        hub_df$stat <- stat_map[hub_df$gene]
      } else {
        hub_df$stat <- hub_df$degree
      }
    }
  }

  # 填充 NA
  hub_df$stat[is.na(hub_df$stat)] <- 0

  # Build gene nodes (only hub genes)
  hub_genes <- hub_df[hub_df$is_hub, ]
  gene_nodes <- data.frame(
    id = hub_genes$gene,
    type = "gene",
    degree = hub_genes$degree,
    stat = hub_genes$stat,
    pathways = hub_genes$pathways,
    stringsAsFactors = FALSE
  )

  # Gene direction based on stat
  gene_nodes$direction <- ifelse(gene_nodes$stat > 0, "up",
                                 ifelse(gene_nodes$stat < 0, "down", "neutral"))

  # Build edges (hub genes -> pathways)
  edges_list <- lapply(hub_genes$gene, function(g) {
    g_upper <- toupper(g)
    connected_pws <- names(gene_sets)[sapply(gene_sets, function(x) g_upper %in% toupper(x))]

    if (length(connected_pws) > 0) {
      data.frame(
        source = g,
        target = connected_pws,
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
  })

  edges <- dplyr::bind_rows(edges_list)

  return(list(
    nodes = list(pathway = pathway_nodes, gene = gene_nodes),
    edges = edges,
    hub_df = hub_df
  ))
}


#' @title Prepare HubGene Nodes for Plotting
#' @description Prepare node data with positions for plotly
#' @param network_data Output from build_hubgene_network()
#' @param layout Layout algorithm ("fr", "kk", "circle")
#' @param seed Random seed for reproducibility
#' @return Named list with node data frames
#' @export
prepare_hubgene_nodes <- function(network_data, layout = "fr", seed = 42) {

  if (is.null(network_data) || is.null(network_data$nodes)) {
    return(NULL)
  }

  nodes <- network_data$nodes
  edges <- network_data$edges

  if (nrow(edges) == 0) {
    return(NULL)
  }

  # Combine all nodes for igraph
  all_nodes <- dplyr::bind_rows(nodes$pathway, nodes$gene)
  all_nodes$id <- as.character(all_nodes$id)

  # Build igraph graph
  g <- igraph::graph_from_data_frame(
    edges,
    vertices = all_nodes,
    directed = FALSE
  )

  # Calculate layout
  set.seed(seed)
  if (layout == "fr") {
    coords <- igraph::layout_with_fr(g)
  } else if (layout == "kk") {
    coords <- igraph::layout_with_kk(g)
  } else if (layout == "circle") {
    coords <- igraph::layout_in_circle(g)
  } else {
    coords <- igraph::layout_with_fr(g)
  }

  # Add coordinates to nodes
  all_nodes$x <- coords[, 1]
  all_nodes$y <- coords[, 2]

  # Calculate node sizes
  # Pathway nodes: based on -log10(FDR) * scale
  if ("FDR" %in% colnames(all_nodes)) {
    all_nodes$size_pathway <- -log10(all_nodes$FDR + 1e-10) * 8 + 15
    all_nodes$size_pathway <- pmin(all_nodes$size_pathway, 40)
  } else {
    all_nodes$size_pathway <- 25
  }

  # Gene nodes: based on degree * scale
  if ("degree" %in% colnames(all_nodes)) {
    all_nodes$size_gene <- all_nodes$degree * 5 + 10
    all_nodes$size_gene <- pmin(all_nodes$size_gene, 30)
  } else {
    all_nodes$size_gene <- 15
  }

  # Assign final size based on type
  all_nodes$size <- ifelse(all_nodes$type == "pathway",
                           all_nodes$size_pathway,
                           all_nodes$size_gene)

  # Split back into pathway and gene nodes
  pathway_nodes_out <- all_nodes[all_nodes$type == "pathway", ]
  gene_nodes_out <- all_nodes[all_nodes$type == "gene", ]

  return(list(
    pathway = pathway_nodes_out,
    gene = gene_nodes_out,
    all = all_nodes
  ))
}


#' @title Color Nodes by Direction
#' @description Assign colors based on NES (pathways) or log2FC (genes)
#' @param node_data Node data frame with type and direction columns
#' @param color_mode Color mode: "logFC", "pathway", or "uniform"
#' @param left_group Left group name for label
#' @param right_group Right group name for label
#' @return Node data with added color column
#' @export
color_by_direction <- function(node_data, color_mode = "logFC",
                               left_group = "A", right_group = "B") {

  # Define color palette (consistent with volcano/NES plots)
  color_up <- "#E41A1C"      # Red for up-regulated
  color_down <- "#377EB8"    # Blue for down-regulated
  color_neutral <- "#999999" # Gray for neutral
  color_hub <- "#FFD700"     # Gold for hub genes

  if (color_mode == "logFC") {
    # For gene nodes: color by log2FC
    # For pathway nodes: color by NES

    node_data$color <- sapply(seq_len(nrow(node_data)), function(i) {
      row <- node_data[i, ]

      if (row$type == "pathway") {
        # Pathway: color by NES direction
        if (!is.na(row$NES)) {
          if (row$NES > 0) return(color_up)
          else if (row$NES < 0) return(color_down)
        }
        return(color_neutral)

      } else if (row$type == "gene") {
        # Gene: color by log2FC
        if (!is.na(row$log2FC)) {
          if (row$log2FC > 0) return(color_up)
          else if (row$log2FC < 0) return(color_down)
        }
        return(color_neutral)
      }

      return(color_neutral)
    })

  } else if (color_mode == "pathway") {
    # Color genes by which pathway they connect to most
    node_data$color <- color_hub  # Default hub color

  } else {
    # Uniform color
    node_data$color <- color_neutral
  }

  # Add color legend labels
  node_data$color_label <- sapply(seq_len(nrow(node_data)), function(i) {
    row <- node_data[i, ]

    if (row$type == "pathway") {
      if (!is.na(row$NES)) {
        if (row$NES > 0) return(paste0("Up in ", left_group))
        else return(paste0("Up in ", right_group))
      }
      return("Neutral")
    } else {
      if (!is.na(row$log2FC)) {
        if (row$log2FC > 0) return(paste0("Up in ", left_group))
        else if (row$log2FC < 0) return(paste0("Up in ", right_group))
      }
      return("Neutral")
    }
  })

  return(node_data)
}


#' @title Generate HubGene Hover Text
#' @description Generate HTML text for plotly hover
#' @param node Node data row (as list with named elements)
#' @param node_type "pathway" or "gene"
#' @param left_group Left group name
#' @param right_group Right group name
#' @return Character string with HTML for hover
#' @export
generate_hubgene_hover_text <- function(node, node_type,
                                        left_group = "A",
                                        right_group = "B") {

  # Safe type conversion - ensure all numeric fields are properly converted
  safe_num <- function(x, default = 0) {
    if (is.null(x) || is.na(x)) return(default)
    x <- as.numeric(x)
    if (is.na(x)) return(default)
    return(x)
  }

  safe_char <- function(x, default = "") {
    if (is.null(x) || is.na(x) || is.factor(x)) return(default)
    return(as.character(x))
  }

  if (node_type == "pathway") {
    # Pathway hover text
    nes <- safe_num(node$NES)
    fdr <- safe_num(node$FDR)
    pval <- safe_num(node$pvalue)
    n_core <- safe_num(node$n_core, 0)
    n_total <- safe_num(node$n_total, 1)
    core_ratio <- if (n_total > 0) n_core / n_total else 0

    direction_text <- if (nes > 0) {
      paste0("▲ Up in ", left_group)
    } else {
      paste0("▼ Up in ", right_group)
    }

    id <- safe_char(node$id, "Unknown")

    hover <- sprintf(
      "<b style='font-size:14px;'>%s</b><br><br>
      <table style='font-size:12px;'>
        <tr><td><b>NES:</b></td><td>%.3f</td></tr>
        <tr><td><b>FDR:</b></td><td>%.2e</td></tr>
        <tr><td><b>p-value:</b></td><td>%.2e</td></tr>
        <tr><td><b>Direction:</b></td><td>%s</td></tr>
        <tr><td><b>Core Genes:</b></td><td>%d / %d (%.1f%%)</td></tr>
      </table>",
      id,
      nes,
      fdr,
      pval,
      direction_text,
      as.integer(n_core),
      as.integer(n_total),
      core_ratio * 100
    )

  } else if (node_type == "gene") {
    # Gene hover text
    id <- safe_char(node$id, "Unknown")
    log2fc <- safe_num(node$log2FC)
    degree <- safe_num(node$degree, 0)
    pathways <- safe_char(node$pathways, "None")
    color <- safe_char(node$color, "#999999")

    if (log2fc > 0) {
      direction_text <- paste0("▲ Up in ", left_group)
    } else if (log2fc < 0) {
      direction_text <- paste0("▼ Up in ", right_group)
    } else {
      direction_text <- "— Neutral"
    }

    n_pathways <- length(unlist(strsplit(pathways, ", ")))

    hover <- sprintf(
      "<b style='font-size:14px; color:%s;'>%s</b><br><br>
      <table style='font-size:12px;'>
        <tr><td><b>log2FC:</b></td><td>%.3f</td></tr>
        <tr><td><b>Direction:</b></td><td>%s</td></tr>
        <tr><td><b>Hub Degree:</b></td><td>%d</td></tr>
        <tr><td><b>Pathways:</b></td><td>%d</td></tr>
      </table><br>
      <b>Appears in:</b><br>
      <small>%s</small>",
      color,
      id,
      log2fc,
      direction_text,
      as.integer(degree),
      n_pathways,
      pathways
    )
  }

  return(hover)
}


#' @title Get Color Legend Configuration
#' @description Generate color legend data for plot
#' @param left_group Left group name
#' @param right_group Right group name
#' @return Named list with legend colors and labels
#' @export
get_hubgene_legend <- function(left_group = "A", right_group = "B") {
  list(
    pathway_up = "#E41A1C",
    pathway_down = "#377EB8",
    gene_up = "#E41A1C",
    gene_down = "#377EB8",
    labels = list(
      pathway_up = paste0("Pathway: Up in ", left_group),
      pathway_down = paste0("Pathway: Up in ", right_group),
      gene_up = paste0("Gene: Up in ", left_group),
      gene_down = paste0("Gene: Up in ", right_group),
      hub_gene = "Hub Gene (degree ≥ 2)"
    )
  )
}
