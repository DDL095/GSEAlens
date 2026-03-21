
# Phase 0: Core Genes 提取与数据访问层（为03/04模块提供统一接口）


#' @title 提取通路的Core Genes（Leading Edge）
#' @description 从GSEA结果中解析core_enrichment字段，提取基因列表。
#'   这是03/04所有可视化模块的基础数据接口。
#' @param gsea_res_obj GseaRes对象或单个GSEA结果（list(status=, data=, ...)）
#' @param pathway_id 字符串，通路ID
#' @return 字符向量，core genes（大写标准化）
#' @export
#' @examples
#' \dontrun{
#' core_genes <- get_core_genes_for_pathway(gsea_res, "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
#' }
get_core_genes_for_pathway <- function(gsea_res_obj, pathway_id) {
  # 统一处理输入：可能是GseaRes对象或已提取的task结果
  if (inherits(gsea_res_obj, "GseaRes")) {
    # 需要指定contrast_id，这里简化处理，假设传入的是单个对比结果
    stop("对于GseaRes对象，请先使用extract_gsea_task提取特定对比组")
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

#' @title 批量提取多个通路的Core Genes
#' @description 为03模块的Network/UpSet/Chord提供批量数据准备
#' @param gsea_task_obj GseaTask对象（单个对比组）
#' @param pathway_ids 字符向量，通路ID列表
#' @return 命名列表，names=pathway_id, values=core_genes向量
#' @export
get_core_genes_list <- function(gsea_task_obj, pathway_ids) {
  if (!inherits(gsea_task_obj, "GseaTask")) {
    stop("必须传入GseaTask对象")
  }

  result <- lapply(pathway_ids, function(pid) {
    get_core_genes_for_pathway(gsea_task_obj, pid)
  })
  names(result) <- pathway_ids
  return(result)
}

#' @title 计算ORA Ratio（用于DotPlot）
#' @description 计算通路基因与差异表达基因的 overlap ratio
#' @param pathway_genes 通路的基因集（TERM2GENE定义）
#' @param de_genes 差异表达基因（根据pvalue阈值筛选）
#' @param ratio_mode 字符串，"ora"（交集/通路大小）或"leading"（交集/DE大小）
#' @return 数值，ratio值
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

#' @title 获取通路的完整基因集（TERM2GENE）
#' @description 从geneset_info中提取通路定义的完整基因列表
#' @param gsea_res GseaRes对象
#' @param pathway_id 通路ID
#' @return 字符向量
#' @export
get_term_genes <- function(gsea_res, pathway_id) {
  term2gene <- gsea_res$geneset_info$term2gene
  genes <- term2gene$gene_symbol[term2gene$gs_name == pathway_id]
  return(toupper(unique(genes)))
}

#' @title 安全参数校验与回退（全局策略）
#' @description 统一处理非法参数（如ncol<=0等），确保03/04模块稳定性
#' @param value 输入值
#' @param default 默认值
#' @param min_val 最小值（含）
#' @param max_val 最大值（含，可选）
#' @param param_name 参数名（用于警告信息）
#' @return 校验后的值
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
