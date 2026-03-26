# Phase 0: Core Genes Extraction and Data Access Layer (Unified Interface for 03/04 Modules)


#' @title Extract Pathway Core Genes (Leading Edge)
#' @description Parse the core_enrichment field from GSEA results to extract gene list.
#'   This is the fundamental data interface for all visualization modules in 03/04.
#' @param gsea_res_obj GseaRes object or single GSEA result (list(status=, data=, ...))
#' @param pathway_id Character, pathway ID
#' @return Character vector, core genes (uppercase normalized)
#' @export

#' @examples
#' \dontrun{
#' core_genes <- get_core_genes_for_pathway(gsea_res, "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
#' }
get_core_genes_for_pathway <- function(gsea_res_obj, pathway_id) {
  # 统一处理输入：可能是GseaRes对象或已提取的task结果
  if (inherits(gsea_res_obj, "GseaRes")) {
    # 需要指定contrast_id，这里简化处理，假设传入的是单个对比结果
    stop("For GseaRes object, please use extract_gsea_task to extract specific contrast first")
  }

  # 处理提取后的GseaTask对象或直接的结果列表
  if (inherits(gsea_res_obj, "GseaTask")) {
    gsea_result <- gsea_res_obj$gsea_res
  } else if (is.list(gsea_res_obj) && !is.null(gsea_res_obj$data)) {
    # 来自results列表的直接元素
    if (gsea_res_obj$status != "Success") return(character(0))
    gsea_result <- gsea_res_obj$data
  } else {
    gsea_result <- gsea_res_obj
  }

  # 从S4对象中提取结果表
  if (methods::is(gsea_result, "gseaResult")) {
    res_df <- as.data.frame(gsea_result@result)
  } else {
    return(character(0))
  }

  # 匹配通路ID
  row_idx <- which(res_df$ID == pathway_id)
  if (length(row_idx) == 0) return(character(0))

  # 解析core_enrichment字段（以/分隔的基因字符串）
  core_str <- as.character(res_df$core_enrichment[row_idx[1]])
  if (is.na(core_str) || core_str == "") return(character(0))

  # 分割并清洗
  genes <- unlist(strsplit(core_str, "/"))
  genes <- toupper(trimws(genes))
  genes <- genes[genes != ""]

  return(unique(genes))
}

#' @title Batch Extract Core Genes for Multiple Pathways
#' @description Batch data preparation for 03 module Network/UpSet/Chord visualizations
#' @param gsea_task_obj GseaTask object (single contrast)
#' @param pathway_ids Character vector, pathway ID list
#' @return Named list, names=pathway_id, values=core_genes vector
#' @export

get_core_genes_list <- function(gsea_task_obj, pathway_ids) {
  if (!inherits(gsea_task_obj, "GseaTask")) {
    stop("Must pass a GseaTask object")
  }

  result <- lapply(pathway_ids, function(pid) {
    get_core_genes_for_pathway(gsea_task_obj, pid)
  })
  names(result) <- pathway_ids
  return(result)
}

#' @title Calculate ORA Ratio (for DotPlot)
#' @description Calculate overlap ratio between pathway genes and differentially expressed genes
#' @param pathway_genes Gene set of pathway (defined by TERM2GENE)
#' @param de_genes Differentially expressed genes (filtered by pvalue threshold)
#' @param ratio_mode Character, "ora" (intersection/pathway size) or "leading" (intersection/DE size)
#' @return Numeric, ratio value
#' @export

calculate_overlap_ratio <- function(pathway_genes, de_genes, ratio_mode = c("ora", "leading")) {
  ratio_mode <- match.arg(ratio_mode)

  pathway_genes <- toupper(pathway_genes)
  de_genes <- toupper(de_genes)

  overlap <- intersect(pathway_genes, de_genes)

  if (ratio_mode == "ora") {
    # 经典ORA：交集 / 通路总基因数
    return(length(overlap) / length(pathway_genes))
  } else {
    # Leading Edge模式：交集 / DE基因总数
    return(length(overlap) / length(de_genes))
  }
}

#' @title Get Full Pathway Gene Set (TERM2GENE)
#' @description Extract complete gene list for a pathway from geneset_info
#' @param gsea_res GseaRes object
#' @param pathway_id Pathway ID
#' @return Character vector
#' @export

get_term_genes <- function(gsea_res, pathway_id) {
  term2gene <- gsea_res$geneset_info$term2gene
  genes <- term2gene$gene_symbol[term2gene$gs_name == pathway_id]
  return(toupper(unique(genes)))
}

#' @title Safe Parameter Validation with Fallback (Global Strategy)
#' @description Unified handling of invalid parameters (e.g., ncol<=0) to ensure 03/04 module stability
#' @param value Input value
#' @param default Default value
#' @param min_val Minimum value (inclusive)
#' @param max_val Maximum value (inclusive, optional)
#' @param param_name Parameter name (for warning messages)
#' @return Validated value
#' @export

validate_param <- function(value, default, min_val = 1, max_val = NULL, param_name = "parameter") {
  # 处理NULL或NA
  if (is.null(value) || is.na(value)) {
    message(sprintf("[Param Check] %s is NULL/NA, fallback to default: %s", param_name, default))
    return(default)
  }

  # 转换为数值
  val <- suppressWarnings(as.numeric(value))
  if (is.na(val)) {
    message(sprintf("[Param Check] %s conversion failed, fallback to: %s", param_name, default))
    return(default)
  }

  # 检查范围
  if (val < min_val) {
    message(sprintf("[Param Check] %s=%s < min(%s), fallback to: %s", param_name, val, min_val, default))
    return(default)
  }

  if (!is.null(max_val) && val > max_val) {
    message(sprintf("[Param Check] %s=%s > max(%s), capped to: %s", param_name, val, max_val, max_val))
    return(max_val)
  }

  return(as.integer(val))
}
