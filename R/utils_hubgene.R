# ==============================================================================
# File: R/utils_hubgene.R
# Purpose: HubGene network construction and visualization utilities
# ==============================================================================

# ------------------------------------------------------------------------------
# Modified: build_hubgene_network() function - Added de_df parameter passing
# ------------------------------------------------------------------------------

#' Build HubGene Network Data Structure (Enhanced with Leading Edge)
#'
#' Build complete network data with leading edge information for edges.
#' Now properly passes de_df to extract_hub_genes for stat extraction.
#'
#' @param gsea_task GseaTask object
#' @param pathway_ids Character vector of selected pathway IDs
#' @param min_hub_degree Minimum degree threshold for hub genes (default: 2)
#' @param de_df Optional DE data frame containing stat column
#' @param res_df Optional pathway result data frame
#' @param seed Random seed for reproducibility (default: 123)
#'
#' @return A list with three components:
#'   \item{nodes}{A list containing pathway and gene data frames}
#'   \item{edges}{Edge data frame with source, target, and leading edge info}
#'   \item{hub_df}{Hub gene data frame from extraction}
#'
#' @export

#' @examples
#' if(interactive()){
#' # Placeholder for function example
#' }
build_hubgene_network <- function(gsea_task, pathway_ids,
                                  min_hub_degree = 2,
                                  de_df = NULL,
                                  res_df = NULL,
                                  seed = 123) {
  # Extract hub genes (now includes leading edge info)
  # Fixed: Properly pass de_df parameter
  hub_df <- extract_hub_genes(gsea_task, pathway_ids, min_hub_degree, de_df = de_df)

  if (is.null(hub_df) || nrow(hub_df) == 0) {
    return(list(nodes = NULL, edges = NULL, hub_df = NULL))
  }

  # Extract gene sets from GSEA result
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
    ifelse(gene_nodes$stat < 0, "down", "neutral")
  )

  # Build edges with leading edge info
  edges_list <- lapply(seq_len(nrow(hub_genes)), function(i) {
    g <- hub_genes$gene[i]
    g_upper <- toupper(g)
    leading_edges <- hub_genes$pathway_leading_edges[[i]]

    connected_pws <- names(gene_sets)[sapply(gene_sets, function(x) g_upper %in% toupper(x))]

    if (length(connected_pws) > 0) {
      edge_info <- lapply(connected_pws, function(pw) {
        is_leading <- if (pw %in% names(leading_edges)) {
          leading_edges[[pw]]
        } else {
          FALSE
        }

        data.frame(
          source = g,
          target = pw,
          is_leading_edge = is_leading,
          edge_width = ifelse(is_leading, 3, 1),
          stringsAsFactors = FALSE
        )
      })

      return(do.call(rbind, edge_info))
    } else {
      return(NULL)
    }
  })

  if (length(edges_list) > 0 && !all(sapply(edges_list, is.null))) {
    edges <- do.call(plyr::rbind.fill, edges_list)
  } else {
    edges <- NULL
  }

  return(list(
    nodes = list(pathway = pathway_nodes, gene = gene_nodes),
    edges = edges,
    hub_df = hub_df
  ))
}


# ------------------------------------------------------------------------------
# Modified: extract_hub_genes() function - Added debug output for diagnosis
# ------------------------------------------------------------------------------

#' Extract Hub Genes from Selected Pathways (Enhanced with Leading Edge)
#'
#' Extract hub genes with their leading edge status for each pathway.
#' A gene is considered a hub if it appears in at least min_degree pathways.
#'
#' @param gsea_task GseaTask object
#' @param pathway_ids Character vector of pathway IDs to analyze
#' @param min_degree Minimum number of pathways a gene must appear in (default: 2)
#' @param de_df Optional differential expression data frame with stat column
#'
#' @return A data frame with columns:
#'   \item{gene}{Gene symbol}
#'   \item{degree}{Number of pathways the gene appears in}
#'   \item{pathways}{Comma-separated list of pathway IDs}
#'   \item{is_hub}{Logical indicating if gene is a hub (degree >= min_degree)}
#'   \item{stat}{Statistic value from de_df or geneList}
#'   \item{pathway_leading_edges}{List column with per-pathway leading edge status}
#'
#' @export

#' @examples
#' if(interactive()){
#' # Placeholder for function example
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

  # ===========================================
  # Extract Leading Edge gene lists
  # ===========================================
  leading_edge_genes <- list()
  res_df <- as.data.frame(gsea_task$gsea_res@result)

  for (pw_id in pathway_ids) {
    row_idx <- which(res_df$ID == pw_id)
    if (length(row_idx) > 0) {
      core_str <- as.character(res_df$core_enrichment[row_idx[1]])
      if (!is.na(core_str) && core_str != "") {
        core_genes <- unlist(strsplit(core_str, "/"))
        core_genes <- toupper(trimws(core_genes))
        core_genes <- core_genes[core_genes != ""]
        leading_edge_genes[[pw_id]] <- unique(core_genes)
      } else {
        leading_edge_genes[[pw_id]] <- character(0)
      }
    } else {
      leading_edge_genes[[pw_id]] <- character(0)
    }
  }

  # ===========================================
  # Source filtering: keep only expressed genes
  # ===========================================

  valid_genes <- NULL

  if (!is.null(de_df)) {
    # Preferentially use genes with stat values from DE table
    if ("gene_symbol" %in% colnames(de_df) && "stat" %in% colnames(de_df)) {
      valid_genes <- de_df$gene_symbol[!is.na(de_df$stat)]
      valid_genes <- toupper(trimws(as.character(valid_genes)))
      valid_genes <- valid_genes[!is.na(valid_genes) & valid_genes != ""]
      valid_genes <- unique(valid_genes)
      message("[extract_hub_genes] Using stat from de_df, ", length(valid_genes), " genes with valid stat")
    } else if ("gene_symbol" %in% colnames(de_df) && "logFC" %in% colnames(de_df)) {
      valid_genes <- de_df$gene_symbol[!is.na(de_df$logFC)]
      valid_genes <- toupper(trimws(as.character(valid_genes)))
      valid_genes <- valid_genes[!is.na(valid_genes) & valid_genes != ""]
      valid_genes <- unique(valid_genes)
      message("[extract_hub_genes] Using logFC from de_df, ", length(valid_genes), " genes with valid logFC")
    }
  }

  # Fallback to geneList if DE table is unavailable
  if (is.null(valid_genes) || length(valid_genes) == 0) {
    if (!is.null(gsea_task$gsea_res@geneList)) {
      valid_genes <- toupper(names(gsea_task$gsea_res@geneList))
      valid_genes <- unique(valid_genes)
      message("[extract_hub_genes] Falling back to geneList, ", length(valid_genes), " genes")
    }
  }

  if (is.null(valid_genes) || length(valid_genes) == 0) {
    warning("[extract_hub_genes] Cannot determine valid genes, returning NULL")
    return(NULL)
  }

  # Filter gene_sets
  gene_sets <- lapply(gene_sets, function(genes) {
    genes_upper <- toupper(trimws(as.character(genes)))
    genes_upper <- genes_upper[!is.na(genes_upper) & genes_upper != ""]
    intersect(genes_upper, valid_genes)
  })

  gene_sets <- gene_sets[sapply(gene_sets, length) > 0]

  if (length(gene_sets) == 0) {
    return(NULL)
  }

  # ===========================================
  # Count gene occurrences across pathways
  # ===========================================

  all_genes <- unlist(gene_sets)
  gene_counts <- table(toupper(all_genes))

  gene_info <- data.frame(
    gene = names(gene_counts),
    degree = as.integer(gene_counts),
    stringsAsFactors = FALSE
  )

  gene_info$pathways <- sapply(gene_info$gene, function(g) {
    g_upper <- toupper(g)
    pws <- names(gene_sets)[sapply(gene_sets, function(x) g_upper %in% toupper(x))]
    paste(pws, collapse = ", ")
  })

  # Determine leading edge status for each pathway
  gene_info$pathway_leading_edges <- lapply(gene_info$gene, function(g) {
    g_upper <- toupper(g)
    leading_status <- sapply(pathway_ids, function(pw_id) {
      if (pw_id %in% names(leading_edge_genes)) {
        g_upper %in% toupper(leading_edge_genes[[pw_id]])
      } else {
        FALSE
      }
    })
    names(leading_status) <- pathway_ids
    return(leading_status)
  })

  gene_info$is_hub <- gene_info$degree >= min_degree

  # ===========================================
  # Critical fix: Properly read stat information
  # ===========================================

  if (!is.null(de_df) && is.data.frame(de_df)) {
    # Ensure gene_symbol column exists
    if ("gene_symbol" %in% colnames(de_df)) {
      # Create gene name to stat mapping
      de_df$gene_upper <- toupper(as.character(de_df$gene_symbol))

      # Prefer stat column
      if ("stat" %in% colnames(de_df)) {
        stat_map <- setNames(de_df$stat, de_df$gene_upper)
        gene_info$stat <- as.numeric(stat_map[gene_info$gene])
        message("[extract_hub_genes] Mapped stat from de_df$stat")
      }
      # Use logFC if no stat column
      else if ("logFC" %in% colnames(de_df)) {
        stat_map <- setNames(de_df$logFC, de_df$gene_upper)
        gene_info$stat <- as.numeric(stat_map[gene_info$gene])
        message("[extract_hub_genes] Mapped stat from de_df$logFC")
      }
      # Use log2FoldChange if no logFC column (DESeq2)
      else if ("log2FoldChange" %in% colnames(de_df)) {
        stat_map <- setNames(de_df$log2FoldChange, de_df$gene_upper)
        gene_info$stat <- as.numeric(stat_map[gene_info$gene])
        message("[extract_hub_genes] Mapped stat from de_df$log2FoldChange")
      } else {
        gene_info$stat <- 0
        warning("[extract_hub_genes] de_df has no stat/logFC/log2FoldChange column")
      }

      # Handle NA values
      gene_info$stat[is.na(gene_info$stat)] <- 0
    } else {
      gene_info$stat <- 0
      warning("[extract_hub_genes] de_df has no gene_symbol column")
    }
  } else {
    # No de_df available, use geneList values as fallback
    gene_list <- gsea_task$gsea_res@geneList
    if (!is.null(gene_list)) {
      stat_map <- setNames(as.numeric(gene_list), toupper(names(gene_list)))
      gene_info$stat <- as.numeric(stat_map[gene_info$gene])
      gene_info$stat[is.na(gene_info$stat)] <- 0
      message("[extract_hub_genes] Mapped stat from geneList (fallback)")
    } else {
      gene_info$stat <- 0
    }
  }

  gene_info$log2FC <- gene_info$stat

  # Sort by degree descending
  gene_info <- gene_info[order(gene_info$degree, decreasing = TRUE), ]

  return(gene_info)
}


#' Prepare HubGene Nodes for Plotting
#'
#' Prepare node data with layout positions for plotly visualization.
#'
#' @param network_data Output from \code{\link{build_hubgene_network}}
#' @param layout Layout algorithm: "fr" (Fruchterman-Reingold),
#'   "kk" (Kamada-Kawai), or "circle"
#' @param seed Random seed for reproducibility (default: 42)
#'
#' @return A named list with three components:
#'   \item{pathway}{Pathway node data frame with coordinates}
#'   \item{gene}{Gene node data frame with coordinates}
#'   \item{all}{Combined node data frame}
#'
#' @export

#' @examples
#' if(interactive()){
#' # Placeholder for function example
#' }
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
  # TODO: BiocCheck: set.seed accepted for reproducibility
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
    all_nodes$size_gene
  )

  # Split back into pathway and gene nodes
  pathway_nodes_out <- all_nodes[all_nodes$type == "pathway", ]
  gene_nodes_out <- all_nodes[all_nodes$type == "gene", ]

  return(list(
    pathway = pathway_nodes_out,
    gene = gene_nodes_out,
    all = all_nodes
  ))
}


#' Color Nodes by Direction
#'
#' Assign colors based on NES (pathways) or log2FC (genes).
#'
#' @param node_data Node data frame with type and direction columns
#' @param color_mode Color mode: "logFC", "pathway", or "uniform"
#' @param left_group Left group name for legend label
#' @param right_group Right group name for legend label
#'
#' @return Node data with added color and color_label columns
#'
#' @export

#' @examples
#' if(interactive()){
#' # Placeholder for function example
#' }
color_by_direction <- function(node_data, color_mode = "logFC",
                               left_group = "A", right_group = "B") {
  # Define color palette (consistent with volcano/NES plots)
  color_up <- "#E41A1C" # Red for up-regulated
  color_down <- "#377EB8" # Blue for down-regulated
  color_neutral <- "#999999" # Gray for neutral
  color_hub <- "#FFD700" # Gold for hub genes

  if (color_mode == "logFC") {
    # For gene nodes: color by log2FC
    # For pathway nodes: color by NES

    node_data$color <- sapply(seq_len(nrow(node_data)), function(i) {
      row <- node_data[i, ]

      if (row$type == "pathway") {
        # Pathway: color by NES direction
        if (!is.na(row$NES)) {
          if (row$NES > 0) {
            return(color_up)
          } else if (row$NES < 0) {
            return(color_down)
          }
        }
        return(color_neutral)
      } else if (row$type == "gene") {
        # Gene: color by log2FC
        if (!is.na(row$log2FC)) {
          if (row$log2FC > 0) {
            return(color_up)
          } else if (row$log2FC < 0) {
            return(color_down)
          }
        }
        return(color_neutral)
      }

      return(color_neutral)
    })
  } else if (color_mode == "pathway") {
    # Color genes by which pathway they connect to most
    node_data$color <- color_hub # Default hub color
  } else {
    # Uniform color
    node_data$color <- color_neutral
  }

  # Add color legend labels
  node_data$color_label <- sapply(seq_len(nrow(node_data)), function(i) {
    row <- node_data[i, ]

    if (row$type == "pathway") {
      if (!is.na(row$NES)) {
        if (row$NES > 0) {
          return(paste0("Up in ", left_group))
        } else {
          return(paste0("Up in ", right_group))
        }
      }
      return("Neutral")
    } else {
      if (!is.na(row$log2FC)) {
        if (row$log2FC > 0) {
          return(paste0("Up in ", left_group))
        } else if (row$log2FC < 0) {
          return(paste0("Up in ", right_group))
        }
      }
      return("Neutral")
    }
  })

  return(node_data)
}


#' Generate HubGene Hover Text
#'
#' Generate HTML text for plotly hover tooltips.
#'
#' @param node Node data row (as list with named elements)
#' @param node_type Node type: "pathway" or "gene"
#' @param left_group Left group name
#' @param right_group Right group name
#'
#' @return Character string with HTML content for hover tooltip
#'
#' @export

#' @examples
#' if(interactive()){
#' # Placeholder for function example
#' }
generate_hubgene_hover_text <- function(node, node_type,
                                        left_group = "A",
                                        right_group = "B") {
  # Safe type conversion - ensure all numeric fields are properly converted
  safe_num <- function(x, default = 0) {
    if (is.null(x) || is.na(x)) {
      return(default)
    }
    x <- as.numeric(x)
    if (is.na(x)) {
      return(default)
    }
    return(x)
  }

  safe_char <- function(x, default = "") {
    if (is.null(x) || is.na(x) || is.factor(x)) {
      return(default)
    }
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
      paste0("Up in ", left_group)
    } else {
      paste0("Up in ", right_group)
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
      direction_text <- paste0("Up in ", left_group)
    } else if (log2fc < 0) {
      direction_text <- paste0("Up in ", right_group)
    } else {
      direction_text <- "Neutral"
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


#' Get Color Legend Configuration
#'
#' Generate color legend data for HubGene network plot.
#'
#' @param left_group Left group name (e.g., treatment or case)
#' @param right_group Right group name (e.g., control or reference)
#'
#' @return A named list containing:
#'   \item{pathway_up}{Color for up-regulated pathways}
#'   \item{pathway_down}{Color for down-regulated pathways}
#'   \item{gene_up}{Color for up-regulated genes}
#'   \item{gene_down}{Color for down-regulated genes}
#'   \item{labels}{Named list of legend labels}
#'
#' @export

#' @examples
#' if(interactive()){
#' # get_hubgene_legend example
#' }
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
      hub_gene = "Hub Gene (degree >= 2)"
    )
  )
}
