# ==============================================================================
# File: R/vis-gsea-core.R
# Purpose: Internal GSEA visualization core (GseaVis-independent)
# ==============================================================================

#' @title Internal GSEA Data Extractor
#' @description Extract running score data from gseaResult object.
#'   Adapted from GseaVis::gsInfo, using enrichit::gseaScores internally.
#' @param object A gseaResult object.
#' @param geneSetID Character, gene set ID.
#' @return data.frame with runningScore, position, geneList, etc.
#' @keywords internal
.gs_info <- function(object, geneSetID) {
  gseaScores <- enrichit::gseaScores

  geneList <- object@geneList

  if (is.numeric(geneSetID)) {
    geneSetID <- object@result[geneSetID, "ID"]
  }

  geneSet <- object@geneSets[[geneSetID]]
  if (is.null(geneSet)) {
    stop(sprintf("Gene set '%s' not found in the GSEA result object. It may belong to a different gene set collection.", geneSetID))
  }
  exponent <- object@params[["exponent"]]

  df <- gseaScores(geneList, geneSet, exponent, fortify = TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore)) / 20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  df$Description <- object@result[geneSetID, "Description"]

  return(df)
}


#' @title Internal GSEA Native Plotting Engine
#' @description Core plotting logic migrated from GseaVis::gseaNb,
#'   stripped of unused features (newGsea, filePath, KEGG, geneExpHt, etc.).
#'   Returns a list of ggplot objects for downstream assembly via patchwork.
#' @param object gseaResult object.
#' @param geneSetID Character vector of pathway IDs.
#' @param subPlot Integer, 1/2/3.
#' @param curveCol Curve color vector.
#' @param htCol Heatmap gradient colors.
#' @param lineSize Line width for enrichment curve.
#' @param addPval Logical, add NES/Pvalue text annotation (single pathway only).
#' @param pvalX X position for pval text (relative fraction).
#' @param pvalY Y position for pval text (relative fraction).
#' @param pvalSize Text size for pval annotation.
#' @param pCol Color for pval text.
#' @param nesDigit NES rounding digits.
#' @param pDigit P-value rounding digits.
#' @param show_contrast_in_axis Logical, if TRUE shows "Left vs Right" on X-axis.
#' @param left_group Left group label.
#' @param right_group Right group label.
#' @param base_size Theme base size.
#' @param subRatio Relative heights for 3 panels.
#' @param rankSeq X-axis break interval for rank plot.
#' @param htHeight Relative height of heatmap strip inside Panel 2.
#' @param htAlpha Alpha for heatmap rectangles.
#' @return Named list: p1 (curve), p2 (heatmap strip), p3 (rank distribution),
#'   data_ga, gsdata, glist.
#' @keywords internal
.gsea_nb_core <- function(object,
                          geneSetID,
                          subPlot = 3,
                          curveCol = c("#76BA99", "#EB4747", "#996699"),
                          htCol = c("#08519C", "#A50F15"),
                          lineSize = 0.8,
                          addPval = FALSE,
                          pvalX = 0.9,
                          pvalY = 0.9,
                          pvalSize = 4,
                          pCol = "grey0",
                          nesDigit = 2,
                          pDigit = 2,
                          show_contrast_in_axis = FALSE,
                          left_group = "Left",
                          right_group = "Right",
                          base_size = 14,
                          subRatio = c(0.5, 0.2, 0.3),
                          rankSeq = 5000,
                          htHeight = 0.3,
                          htAlpha = 0.8) {

  # ============================================================================
  # 1. Data preparation
  # ============================================================================
  glist <- object@geneList

  # Running score data for all pathways
  gsdata_list <- lapply(geneSetID, function(setid) {
    tmp <- .gs_info(object, geneSetID = setid)
    tmp$id <- setid
    tmp
  })
  gsdata <- dplyr::bind_rows(gsdata_list)

  # In-pathway gene positions (only used for tick marks in Panel 2)
  gsdata1_list <- lapply(unique(gsdata$Description), function(setid) {
    tmp <- gsdata[gsdata$Description == setid, ]
    tmp[tmp$position == 1, ]
  })
  gsdata1 <- dplyr::bind_rows(gsdata1_list)

  # Pathway statistics
  data_ga <- as.data.frame(object)
  data_ga <- data_ga[data_ga$ID %in% geneSetID, ]
  data_ga <- data_ga[match(unique(gsdata$id), data_ga$ID), ]

  n_pathways <- length(geneSetID)

  # ============================================================================
  # 2. Panel 1: Enrichment Curve
  # ============================================================================
  if (n_pathways == 1) {
    if (length(curveCol) == 1) {
      p1 <- ggplot2::ggplot(gsdata, ggplot2::aes(x = .data$x, y = .data$runningScore)) +
        ggplot2::geom_line(color = curveCol[1], linewidth = lineSize)
    } else {
      p1 <- ggplot2::ggplot(gsdata, ggplot2::aes(x = .data$x, y = .data$runningScore)) +
        ggplot2::geom_line(ggplot2::aes(color = .data$runningScore), linewidth = lineSize) +
        ggplot2::scale_color_gradient(low = curveCol[1], high = curveCol[2])
    }
  } else {
    mulcol <- curveCol[seq_len(min(n_pathways, length(curveCol)))]
    names(mulcol) <- geneSetID
    gsdata$id <- factor(gsdata$id, levels = geneSetID)

    p1 <- ggplot2::ggplot(gsdata, ggplot2::aes(x = .data$x, y = .data$runningScore)) +
      ggplot2::geom_line(ggplot2::aes(color = .data$id), linewidth = lineSize) +
      ggplot2::scale_color_manual(values = mulcol, name = "Term Name")
  }

  p1 <- p1 +
    ggplot2::geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::theme(
      legend.position = if (n_pathways == 1) "none" else "right",
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 0.2, r = 0.2, b = 0, l = 0.2, unit = "cm")
    ) +
    ggplot2::ylab("Running Enrichment Score")

  # P-value annotation (single pathway only)
  if (isTRUE(addPval) && n_pathways == 1 && nrow(data_ga) > 0) {
    p_label <- paste0(
      "NES: ", round(data_ga$NES[1], digits = nesDigit), "\n",
      "Pvalue: ", ifelse(data_ga$pvalue[1] < 0.001, "< 0.001", round(data_ga$pvalue[1], digits = pDigit)), "\n",
      "Adjusted Pvalue: ", ifelse(data_ga$p.adjust[1] < 0.001, "< 0.001", round(data_ga$p.adjust[1], digits = pDigit))
    )

    px <- pvalX * nrow(gsdata[gsdata$id == geneSetID[1], ])
    py <- pvalY * sum(abs(range(gsdata$runningScore))) + min(gsdata$runningScore)

    p1 <- p1 +
      ggplot2::annotate(
        "text",
        x = px,
        y = py,
        label = p_label,
        size = pvalSize,
        color = pCol,
        fontface = "italic",
        hjust = 1
      )
  }

  # ============================================================================
  # 3. Panel 2: Heatmap Strip + Gene Ticks
  # ============================================================================
  p2 <- NULL

  if (subPlot >= 2) {
    # Bin geneList into 10 intervals for heatmap
    df_rank <- data.frame(pos = seq_along(glist), fc = glist)
    qt <- stats::quantile(df_rank$fc, probs = c(0.1, 0.9))

    df_rank <- df_rank |>
      dplyr::mutate(fc = dplyr::case_when(
        .data$fc >= qt[2] ~ qt[2],
        .data$fc <= qt[1] ~ qt[1],
        .default = .data$fc
      ))

    break_intervals <- cut(df_rank$fc, breaks = 10, include.lowest = TRUE)
    intervals <- levels(break_intervals)

    interval_bounds <- data.frame(do.call(rbind, strsplit(gsub("[()\\[\\]]", "", intervals), ",")))
    interval_bounds$X1 <- as.numeric(vapply(strsplit(interval_bounds$X1, split = "\\[|\\("), "[", character(1), 2))
    interval_bounds$X2 <- as.numeric(vapply(strsplit(interval_bounds$X2, split = "\\]|\\)"), "[", character(1), 1))

    start_positions <- vapply(interval_bounds[, 1], function(bound) which.min(abs(df_rank$fc - bound)), integer(1))
    start_positions[1] <- nrow(df_rank)
    end_positions <- vapply(interval_bounds[, 2], function(bound) which.min(abs(df_rank$fc - bound)), integer(1))
    interval_means <- tapply(df_rank$fc, break_intervals, mean)

    result_df <- data.frame(
      Interval = seq_len(10),
      Start = interval_bounds[, 1],
      End = interval_bounds[, 2],
      xmax = start_positions,
      xmin = end_positions,
      ymin = 0,
      ymax = htHeight,
      Mean_LogFC = as.numeric(interval_means)
    )

    # Replicate for each pathway
    result_df_all <- dplyr::bind_rows(lapply(unique(gsdata$id), function(setid) {
      result_df$id <- setid
      result_df
    }))

    # Factor ordering
    if (gsdata$id[1] == gsdata$Description[1]) {
      result_df_all$id <- factor(result_df_all$id, levels = geneSetID)
      gsdata1$id <- factor(gsdata1$id, levels = geneSetID)
    } else {
      result_df_all$id <- factor(result_df_all$id, levels = data_ga$ID)
      gsdata1$id <- factor(gsdata1$id, levels = data_ga$ID)
    }

    # Build plot
    p2 <- ggplot2::ggplot() +
      ggplot2::geom_rect(
        data = result_df_all,
        ggplot2::aes(
          xmin = .data$xmin, xmax = .data$xmax,
          ymin = .data$ymin, ymax = .data$ymax,
          fill = .data$Mean_LogFC
        ),
        color = NA, inherit.aes = FALSE, alpha = htAlpha, show.legend = FALSE
      ) +
      ggplot2::scale_fill_gradient2(low = htCol[1], mid = "white", high = htCol[2], midpoint = 0)

    if (n_pathways > 1) {
      p2 <- p2 +
        ggplot2::geom_segment(
          data = gsdata1,
          ggplot2::aes(x = .data$x, xend = .data$x, y = 0, yend = 1, color = .data$id),
          show.legend = FALSE
        ) +
        ggplot2::scale_color_manual(values = mulcol)
    } else {
      p2 <- p2 +
        ggplot2::geom_segment(
          data = gsdata1,
          ggplot2::aes(x = .data$x, xend = .data$x, y = 0, yend = 1),
          color = "black",
          show.legend = FALSE
        )
    }

    p2 <- p2 +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(
        axis.ticks = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_blank(),
        panel.spacing = ggplot2::unit(0.1, "cm"),
        plot.margin = ggplot2::margin(t = 0, r = 0.2, b = 0.2, l = 0.2, unit = "cm")
      ) +
      ggplot2::xlab("Rank in Ordered Dataset")

    if (n_pathways > 1) {
      p2 <- p2 + ggplot2::facet_wrap(~id, ncol = 1)
    }

    if (subPlot > 2) {
      p2 <- p2 + ggplot2::theme(axis.title.x = ggplot2::element_blank())
    }
  }

  # ============================================================================
  # 4. Panel 3: Ranked List Distribution (Native Red-Blue)
  # ============================================================================
  p3 <- NULL

  if (subPlot == 3) {
    df_rank_plot <- data.frame(x = seq_along(glist), y = as.numeric(glist))

    x_axis_label <- if (isTRUE(show_contrast_in_axis)) {
      sprintf("%s vs %s", left_group, right_group)
    } else {
      "Rank in Ordered Dataset"
    }

    p3 <- ggplot2::ggplot(df_rank_plot, ggplot2::aes(x = .data$x, y = .data$y)) +
      ggplot2::geom_col(ggplot2::aes(fill = .data$y), width = 1, color = NA, show.legend = FALSE) +
      ggplot2::scale_fill_gradient2(low = "#08519C", mid = "white", high = "#A50F15", midpoint = 0) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.5, color = "black", linetype = "dashed") +
      ggplot2::scale_x_continuous(breaks = seq(0, length(glist), rankSeq)) +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text = ggplot2::element_text(colour = "black"),
        plot.margin = ggplot2::margin(t = -0.1, r = 0.2, b = 0.2, l = 0.2, unit = "cm")
      ) +
      ggplot2::coord_cartesian(expand = 0) +
      ggplot2::labs(x = x_axis_label, y = "Ranked List")

    # Joint canvas mode: bold X-axis label
    if (isTRUE(show_contrast_in_axis)) {
      p3 <- p3 +
        ggplot2::theme(axis.title.x = ggplot2::element_text(face = "bold", size = 18))
    }

    # Zero cross and group annotations
    z_cross <- sum(glist > 0)
    m_rank <- length(glist)

    p3 <- p3 +
      ggplot2::geom_vline(xintercept = z_cross, linetype = "dashed", color = "grey50") +
      ggplot2::annotate(
        "text",
        x = z_cross, y = 0,
        label = paste0("Zero cross at ", z_cross),
        vjust = 1.5, hjust = -0.05, size = 3, color = "grey30"
      ) +
      ggplot2::annotate(
        "text",
        x = m_rank * 0.01, y = max(glist) * 0.85,
        label = sprintf("'%s' (pos)", left_group),
        color = "#A50F15", hjust = 0, size = 4, fontface = "italic"
      ) +
      ggplot2::annotate(
        "text",
        x = m_rank * 0.99, y = min(glist) * 0.85,
        label = sprintf("'%s' (neg)", right_group),
        color = "#08519C", hjust = 1, size = 4, fontface = "italic"
      )
  }

  # ============================================================================
  # 5. Return
  # ============================================================================
  return(list(
    p1 = p1,
    p2 = p2,
    p3 = p3,
    data_ga = data_ga,
    gsdata = gsdata,
    glist = glist
  ))
}
