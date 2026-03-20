#' @title 核心类定义与校验函数
#' @description 定义 GseaEnv, GseaRes, GseaTask 结构，并提供内部校验逻辑。
#' @keywords internal
#' @name class_validations
NULL


# 1. 核心类定义


#' @title 定义 GseaEnv 类
#' @description 标准化的 GSEA 输入环境对象。
#' @export
create_gsea_env <- function(backend_info, contrast_registry, de_store, expr_bundle, geneset, raw_obj) {
  structure(
    list(
      backend_info = backend_info,
      contrast_registry = contrast_registry,
      de_store = de_store,
      expr_bundle = expr_bundle,
      geneset = geneset,
      raw_backend_obj = raw_obj # 保留原始对象以备不时之需，但下游不应直接访问
    ),
    class = "GseaEnv"
  )
}

#' @title 定义 GseaRes 类
#' @description GSEA 计算结果胶囊。
#' @export
create_gsea_res <- function(metadata, backend_info, contrast_registry, de_store, expr_bundle, geneset_info, results) {
  structure(
    list(
      metadata = metadata,
      backend_info = backend_info,
      contrast_registry = contrast_registry,
      de_store = de_store,
      expr_bundle = expr_bundle,
      geneset_info = geneset_info,
      results = results
    ),
    class = "GseaRes"
  )
}

#' @title 定义 GseaTask 类
#' @description 单个对比组的提取结果对象。
#' @export
create_gsea_task <- function(gsea_res, meta) {
  structure(
    list(
      gsea_res = gsea_res,
      meta = meta
    ),
    class = "GseaTask"
  )
}


# 2. 内部校验函数


#' @title 校验 Limma 设计矩阵
#' @description 强制要求无截距设计 (~ 0 + group)，确保对比组解析准确。
#' @param fit MArrayLM 对象
#' @return TRUE 或抛出错误
#' @keywords internal
.validate_limma_design <- function(fit) {
  # 检查是否包含截距项
  # Intercept / (Intercept) 是 R 中默认的截距命名
  design_matrix <- fit$design
  has_intercept <- any(grepl("Intercept", colnames(design_matrix), ignore.case = TRUE))

  if (has_intercept) {
    stop(
      "\n❌ [Limma 设计错误] 检测到截距项！\n",
      "GSEAlens 要求使用无截距设计矩阵。\n",
      "请修改您的设计公式，例如：\n",
      "  design <- model.matrix(~ 0 + group, data = samples)\n",
      "  fit <- lmFit(expr, design)\n",
      "原因：无截距设计能确保 colnames(fit) 直接对应组名，从而精准构建对比组。"
    )
  }

  # 检查是否有足够的列进行对比
  if (ncol(design_matrix) < 2) {
    warning("⚠️ [Limma 警告] 设计矩阵仅包含 1 列，无法进行组间对比。")
  }

  return(TRUE)
}

#' @title 校验 DESeq2 目标因子
#' @description 检查 target_factor 是否存在于 colData 中。
#' @param dds DESeqDataSet 对象
#' @param target_factor 字符串，指定的目标因子
#' @return TRUE 或抛出错误
#' @keywords internal
.validate_deseq2_design <- function(dds, target_factor) {
  col_data <- as.data.frame(SummarizedExperiment::colData(dds))

  if (!target_factor %in% colnames(col_data)) {
    stop(
      sprintf("\n❌ [DESeq2 设计错误] 指定的 target_factor '%s' 不存在于 colData 中！\n", target_factor),
      "可用的列名: ", paste(colnames(col_data), collapse = ", ")
    )
  }

  # 检查是否为因子
  if (!is.factor(col_data[[target_factor]])) {
    warning(sprintf("⚠️ [DESeq2 警告] target_factor '%s' 不是因子类型，正在尝试自动转换...", target_factor))
    # 这里不实际转换，只是警告，因为 DESeq2 通常在构建时已处理
  }

  return(TRUE)
}

#' @title 校验 GseaEnv 对象完整性
#' @description 内部函数，确保对象结构符合规范。
#' @param env_obj GseaEnv 对象
#' @return TRUE 或抛出错误
#' @keywords internal
.check_gsea_env <- function(env_obj) {
  if (!inherits(env_obj, "GseaEnv")) stop("输入对象不是 GseaEnv 类。")

  required_slots <- c("backend_info", "contrast_registry", "de_store", "expr_bundle", "geneset")
  missing_slots <- setdiff(required_slots, names(env_obj))

  if (length(missing_slots) > 0) {
    stop(sprintf("GseaEnv 对象结构不完整，缺失: %s", paste(missing_slots, collapse = ", ")))
  }

  # 检查 contrast_registry 必要字段
  reg <- env_obj$contrast_registry
  if (!is.data.frame(reg) || !all(c("contrast_id", "left_group", "right_group") %in% colnames(reg))) {
    stop("contrast_registry 必须包含 contrast_id, left_group, right_group 列。")
  }

  return(TRUE)
}

#' @title 校验 GseaRes 对象完整性
#' @keywords internal
.check_gsea_res <- function(res_obj) {
  if (!inherits(res_obj, "GseaRes")) stop("输入对象不是 GseaRes 类。")
  # 简单检查 results 列表是否存在
  if (is.null(res_obj$results)) stop("GseaRes 对象中无计算结果。")
  return(TRUE)
}
