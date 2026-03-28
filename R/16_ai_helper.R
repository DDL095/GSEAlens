# =============================================================================
# AI Helper Functions (Pre-reservation for Future Integration)
# =============================================================================

#' @title Prepare Pathway Summary for AI Analysis
#' @description Generate structured summary of GSEA results suitable for AI
#'   interpretation. Includes pathway statistics, quality assessment, and
#'   interpretation guidelines.
#' @param pathways Data frame of GSEA results (must contain columns:
#'   ID, Description, NES, pvalue, p.adjust, setSize, core_enrichment)
#' @param include_leading_edge Logical. Whether to include leading edge genes.
#'   Default TRUE.
#' @param format Character. Output format: "markdown", "json", or "prompt".
#'   Default "markdown".
#' @param quality_threshold List. Thresholds for quality assessment with named
#'   elements: nes (default 1.0), pval (default 0.01), padj (default 0.25).
#' @return A character string in the specified format
#' @export
#' @examples
#' \dontrun{
#' # Extract significant pathways
#' sig_pathways <- df[df$p.adjust < 0.25, ]
#'
#' # Generate markdown report
#' prompt <- prepare_pathway_prompt(sig_pathways)
#' cat(prompt)
#'
#' # Generate JSON for API
#' json_data <- prepare_pathway_prompt(sig_pathways, format = "json")
#' }
prepare_pathway_prompt <- function(pathways,
                                   include_leading_edge = TRUE,
                                   format = "markdown",
                                   quality_threshold = list(nes = 1.0, pval = 0.01, padj = 0.25)) {

  # Validate input
  required_cols <- c("ID", "NES", "pvalue", "p.adjust")
  missing_cols <- setdiff(required_cols, colnames(pathways))
  if (length(missing_cols) > 0) {
    stop("pathways must contain columns: ", paste(required_cols, collapse = ", "))
  }

  # Filter pathways with valid statistics
  pathways_filtered <- pathways %>%
    dplyr::filter(!is.na(NES), !is.na(pvalue), !is.na(p.adjust))

  # Sort by absolute NES
  pathways_filtered <- pathways_filtered %>%
    dplyr::arrange(dplyr::desc(abs(NES)))

  # Build summary list
  summary_list <- lapply(seq_len(nrow(pathways_filtered)), function(i) {
    row <- pathways_filtered[i, ]

    # Quality assessment
    quality <- assess_pathway_quality(
      nes = row$NES,
      pvalue = row$pvalue,
      padj = row$p.adjust,
      thresholds = quality_threshold
    )

    # Leading edge genes
    leading_edge <- NULL
    if (include_leading_edge && "core_enrichment" %in% colnames(row)) {
      core_str <- row$core_enrichment
      if (!is.na(core_str) && core_str != "") {
        leading_edge <- unlist(strsplit(as.character(core_str), "/"))
        leading_edge <- leading_edge[leading_edge != ""]
      }
    }

    # Description fallback
    description <- if ("Description" %in% colnames(row) && !is.na(row$Description)) {
      as.character(row$Description)
    } else {
      as.character(row$ID)
    }

    list(
      rank = i,
      id = as.character(row$ID),
      description = description,
      nes = round(row$NES, 3),
      abs_nes = round(abs(row$NES), 3),
      pvalue = row$pvalue,
      padj = row$p.adjust,
      set_size = if ("setSize" %in% colnames(row)) row$setSize else NA,
      enriched_in = if (row$NES > 0) "up" else "down",
      quality = quality,
      leading_edge = leading_edge
    )
  })

  # Format output
  switch(format,
         markdown = .format_as_markdown(summary_list),
         json = .format_as_json(summary_list),
         prompt = .build_ai_prompt(summary_list),
         stop("Unsupported format: ", format, ". Use 'markdown', 'json', or 'prompt'")
  )
}

#' @title Assess Pathway Quality
#' @description Evaluate GSEA pathway quality based on statistical criteria.
#'   Returns a character string describing the pathway quality with specific issues.
#' @param nes Numeric. Normalized Enrichment Score
#' @param pvalue Numeric. P-value
#' @param padj Numeric. Adjusted P-value (FDR)
#' @param thresholds List. Named thresholds for assessment with elements:
#'   \itemize{
#'     \item nes: Minimum |NES| threshold (default 1.0)
#'     \item pval: Maximum p-value threshold (default 0.01)
#'     \item padj: Maximum FDR threshold (default 0.25)
#'   }
#' @return Character string describing pathway quality
#' @export
#' @examples
#' \dontrun{
#' assess_pathway_quality(2.5, 0.001, 0.01)
#' # Returns: "High Quality: NES=2.50, FDR=0.010, suitable for detailed analysis"
#'
#' assess_pathway_quality(0.8, 0.05, 0.35)
#' # Returns: "Low Quality: NES<1 (0.80), FDR>0.25 (0.35)"
#' }
assess_pathway_quality <- function(nes, pvalue, padj,
                                   thresholds = list(nes = 1.0, pval = 0.01, padj = 0.25)) {

  issues <- character(0)

  # Check NES
  if (abs(nes) < thresholds$nes) {
    issues <- c(issues, sprintf("NES=%.2f below threshold (|NES|<%.1f)", abs(nes), thresholds$nes))
  }

  # Check FDR
  if (padj > thresholds$padj) {
    issues <- c(issues, sprintf("FDR=%.3f exceeds threshold (FDR>%.2f)", padj, thresholds$padj))
  }

  # Check p-value
  if (pvalue > thresholds$pval) {
    issues <- c(issues, sprintf("P=%.3f exceeds threshold (P>%.2f)", pvalue, thresholds$pval))
  }

  # Determine quality level
  if (length(issues) == 0) {
    quality_level <- "High Quality"
    recommendation <- "Recommended for detailed analysis"
    emoji <- ""
  } else if (length(issues) == 1) {
    quality_level <- "Moderate Quality"
    recommendation <- "Use with caution"
    emoji <- ""
  } else {
    quality_level <- "Low Quality"
    recommendation <- "Interpret carefully or exclude"
    emoji <- ""
  }

  # Format output
  if (length(issues) == 0) {
    sprintf("%s [%s]: NES=%.2f, FDR=%.3f, P=%.2e | %s",
            emoji, quality_level, nes, padj, pvalue, recommendation)
  } else {
    sprintf("%s [%s]: %s | %s",
            emoji, quality_level,
            paste(issues, collapse = "; "),
            recommendation)
  }
}

# =============================================================================
# Internal Helper Functions
# =============================================================================

.format_as_markdown <- function(summary_list) {
  lines <- character(0)

  # Header
  lines <- c(lines, "# GSEA Pathway Summary Report")
  lines <- c(lines, "")
  lines <- c(lines, sprintf("Generated: %s", Sys.time()))
  lines <- c(lines, sprintf("Total Pathways: %d", length(summary_list)))
  lines <- c(lines, "")
  lines <- c(lines, "---")
  lines <- c(lines, "")
  lines <- c(lines, "## Pathway Table")
  lines <- c(lines, "")
  lines <- c(lines, "| Rank | Pathway ID | Description | NES | FDR | Quality |")
  lines <- c(lines, "|------|-------------|-------------|-----|-----|---------|")

  for (item in summary_list) {
    desc <- substr(item$description, 1, 40)
    if (nchar(item$description) > 40) desc <- paste0(desc, "...")

    lines <- c(lines, sprintf("| %d | %s | %s | %.2f | %.3f | %s |",
                              item$rank,
                              item$id,
                              desc,
                              item$nes,
                              item$padj,
                              item$quality))
  }

  lines <- c(lines, "")
  lines <- c(lines, "---")
  lines <- c(lines, "")
  lines <- c(lines, "## Quality Legend")
  lines <- c(lines, "")
  lines <- c(lines, "- **High Quality**: |NES|>=1, FDR<=0.25, P<=0.01 | Recommended")
  lines <- c(lines, "- **Moderate Quality**: One criterion not met | Use with caution")
  lines <- c(lines, "- **Low Quality**: Multiple criteria not met | Interpret carefully")
  lines <- c(lines, "")

  # Leading edge genes for top pathways
  top_pathways <- head(summary_list, 5)
  top_with_le <- Filter(function(x) !is.null(x$leading_edge), top_pathways)

  if (length(top_with_le) > 0) {
    lines <- c(lines, "---")
    lines <- c(lines, "")
    lines <- c(lines, "## Top Pathway Leading Edge Genes")
    lines <- c(lines, "")

    for (item in top_with_le) {
      genes <- paste(head(item$leading_edge, 20), collapse = ", ")
      if (length(item$leading_edge) > 20) {
        genes <- paste0(genes, sprintf(" ... (+%d more)", length(item$leading_edge) - 20))
      }
      lines <- c(lines, sprintf("**%s** (%s): %s",
                                item$id, item$quality, genes))
      lines <- c(lines, "")
    }
  }

  paste(lines, collapse = "\n")
}

.format_as_json <- function(summary_list) {
  # Convert to JSON-friendly format
  json_list <- lapply(summary_list, function(item) {
    if (is.null(item$leading_edge) || length(item$leading_edge) == 0) {
      item$leading_edge <- NA
    }
    item
  })

  jsonlite::toJSON(json_list, auto_unbox = TRUE, pretty = TRUE)
}

.build_ai_prompt <- function(summary_list) {
  lines <- character(0)

  lines <- c(lines, "# Role Definition")
  lines <- c(lines, "")
  lines <- c(lines, "You are a bioinformatics expert specializing in transcriptomics analysis.")
  lines <- c(lines, "Your task is to help users understand the biological significance of GSEA pathway enrichment results.")
  lines <- c(lines, "")
  lines <- c(lines, "# Input Data Description")
  lines <- c(lines, "")
  lines <- c(lines, "The following is a Gene Set Enrichment Analysis (GSEA) result from RNA-seq data analysis.")
  lines <- c(lines, "Each pathway contains: pathway ID, description, NES value, p-value, FDR, gene set size, and enrichment direction.")
  lines <- c(lines, "")
  lines <- c(lines, "# Interpretation Guidelines")
  lines <- c(lines, "")
  lines <- c(lines, "## NES (Normalized Enrichment Score)")
  lines <- c(lines, "- |NES| > 1: Strong enrichment signal")
  lines <- c(lines, "- |NES| < 1: Weak enrichment signal")
  lines <- c(lines, "- NES > 0: Pathway enriched at the top of the ranked list (typically up-regulated)")
  lines <- c(lines, "- NES < 0: Pathway enriched at the bottom of the ranked list (typically down-regulated)")
  lines <- c(lines, "")
  lines <- c(lines, "## Statistical Significance")
  lines <- c(lines, "- FDR < 0.25: MSigDB standard threshold")
  lines <- c(lines, "- FDR < 0.05: Stringent threshold")
  lines <- c(lines, "- Consider both raw p-value and FDR together")
  lines <- c(lines, "")
  lines <- c(lines, "## Leading Edge Genes")
  lines <- c(lines, "- Core genes contributing most to the enrichment")
  lines <- c(lines, "- Critical for understanding pathway regulatory mechanisms")
  lines <- c(lines, "")
  lines <- c(lines, "## Quality Assessment")
  lines <- c(lines, "- High Quality: All criteria met | Recommended for detailed analysis")
  lines <- c(lines, "- Moderate Quality: One criterion not met | Use with caution")
  lines <- c(lines, "- Low Quality: Multiple criteria not met | Interpret carefully")
  lines <- c(lines, "")
  lines <- c(lines, "# Output Requirements")
  lines <- c(lines, "Please provide a systematic analysis of the following pathways:")
  lines <- c(lines, "1. Identify the most significant pathways")
  lines <- c(lines, "2. Explain their biological significance")
  lines <- c(lines, "3. Point out potential connections between pathways")
  lines <- c(lines, "4. Pay special attention to high-quality pathways (green/yellow markers)")
  lines <- c(lines, "5. Provide cautious evaluation for low-quality pathways (red markers)")
  lines <- c(lines, "")
  lines <- c(lines, "# Pathway Data")
  lines <- c(lines, "")

  # Add markdown table
  lines <- c(lines, .format_as_markdown(summary_list))

  paste(lines, collapse = "\n")
}
