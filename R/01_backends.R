#' @title 后端数据提取器
#' @description 从 limma 或 DESeq2 对象中提取标准化数据结构。
#' @keywords internal
#' @name backends
#' @noRd
NULL


# 1. Limma 后端提取器


#' @title 提取 Limma-Voom 数据
#' @description 从 MArrayLM 对象提取对比组和差异分析结果。
#' @param fit MArrayLM 对象 (必须经过 eBayes)
#' @param expr_data 可选。DGEList 或表达矩阵。
#' @return 包含 contrast_registry, de_store, expr_bundle 的列表
#' @keywords internal
.extract_limma_data <- function(fit, expr_data = NULL) {

  # 1. 强制校验：无截距设计
  .validate_limma_design(fit)

  # 2. 解析对比组
  coef_names <- colnames(fit)

  # 检查是否包含对比符号 " - " (标准 limma contrast 命名)
  is_contrast_obj <- !is.null(fit$contrasts) || any(grepl(" - ", coef_names))

  if (!is_contrast_obj) {
    stop(
      "\n❌ [Limma 输入错误] 检测到设计矩阵列名为组名 (如 'GroupA', 'GroupB')，而非对比 (如 'GroupA - GroupB')。\n",
      "GSEAlens 要求传入已经定义好对比的 fit 对象。\n",
      "请使用 makeContrasts 和 contrasts.fit 定义您的比较组。"
    )
  }

  # 3. 构建 contrast_registry
  parsed_contrasts <- lapply(coef_names, function(name) {
    parts <- strsplit(name, " - ")[[1]]
    if (length(parts) == 2) {
      return(list(left = trimws(parts[1]), right = trimws(parts[2])))
    } else {
      return(list(left = name, right = "Background"))
    }
  })

  contrast_registry <- tibble::tibble(
    contrast_id = sapply(parsed_contrasts, function(x) paste(x$left, x$right, sep = "_vs_")),
    left_group = sapply(parsed_contrasts, `[[`, "left"),
    right_group = sapply(parsed_contrasts, `[[`, "right"),
    source_name = coef_names,
    backend = "limma_voom"
  )

  # 4. 提取 DE 结果
  de_store <- list()

  for (i in seq_len(nrow(contrast_registry))) {
    reg_row <- contrast_registry[i, ]
    tt <- limma::topTable(fit, coef = reg_row$source_name, number = Inf, sort.by = "none")
    de_store[[reg_row$contrast_id]] <- .standardize_de_columns(df = tt, backend = "limma_voom")
  }

  # 5. 构建 expr_bundle
  expr_bundle <- .build_expr_bundle(expr_data, backend = "limma_voom")

  return(list(
    contrast_registry = contrast_registry,
    de_store = de_store,
    expr_bundle = expr_bundle
  ))
}


# 2. DESeq2 后端提取器


#' @title 提取 DESeq2 数据
#' @description 从 DESeqDataSet 对象提取对比组和差异分析结果。
#' @param dds DESeqDataSet 对象 (必须已运行 DESeq())
#' @param target_factor 字符串。目标因子。若为 NULL，自动推断为设计公式最后一项。
#' @return 包含 contrast_registry, de_store, expr_bundle 的列表
#' @keywords internal
#' @noRd
.extract_deseq2_data <- function(dds, target_factor = NULL) {

  # 1. 确定 target_factor
  design_formula <- DESeq2::design(dds)
  design_terms <- attr(terms(design_formula), "term.labels")

  if (is.null(target_factor)) {
    target_factor <- utils::tail(design_terms, 1)
    message(sprintf("🔍 [DESeq2] 未指定 target_factor，自动推断为: '%s'", target_factor))
  }

  # 校验 target_factor
  .validate_deseq2_design(dds, target_factor)

  # 2. 获取所有 levels
  col_data <- as.data.frame(SummarizedExperiment::colData(dds))
  factor_levels <- levels(col_data[[target_factor]])

  if (length(factor_levels) < 2) {
    stop(sprintf("因子 '%s' 的水平数少于 2，无法进行比较。", target_factor))
  }

  # 3. 生成所有成对比较
  combos <- combn(factor_levels, 2, simplify = FALSE)

  contrast_registry_list <- list()
  de_store_list <- list() # 定义列表变量

  for (combo in combos) {
    left <- combo[1]
    right <- combo[2]
    contrast_id <- paste(left, right, sep = "_vs_")

    # DESeq2 contrast 格式: c(target_factor, numerator, denominator)
    contrast_vec <- c(target_factor, left, right)

    # 提取结果
    res <- tryCatch({
      DESeq2::results(dds, contrast = contrast_vec)
    }, error = function(e) {
      warning(sprintf("提取对比 %s 时出错: %s", contrast_id, e$message))
      return(NULL)
    })

    if (!is.null(res)) {
      # 添加到 registry
      contrast_registry_list[[contrast_id]] <- tibble::tibble(
        contrast_id = contrast_id,
        left_group = left,
        right_group = right,
        source_name = paste(left, "vs", right, sep = " "),
        backend = "deseq2"
      )

      # 标准化列名
      de_store_list[[contrast_id]] <- .standardize_de_columns(
        df = as.data.frame(res),
        backend = "deseq2"
      )
    }
  }

  contrast_registry <- dplyr::bind_rows(contrast_registry_list)

  # 4. 构建 expr_bundle
  expr_bundle <- .build_expr_bundle(dds, backend = "deseq2")

  # 🌟 修复点：返回 de_store_list 而非 de_store
  return(list(
    contrast_registry = contrast_registry,
    de_store = de_store_list,
    expr_bundle = expr_bundle
  ))
}


# 3. 辅助函数


#' @title 标准化差异分析表列名
#' @description 将不同后端的列名统一为: gene_symbol, logFC, stat, pvalue, padj
#' @keywords internal
#' @noRd
.standardize_de_columns <- function(df, backend) {

  # 添加基因名列
  if ("gene_symbol" %in% colnames(df)) {
    # pass
  } else if ("SYMBOL" %in% colnames(df)) {
    df$gene_symbol <- df$SYMBOL
  } else {
    df$gene_symbol <- rownames(df)
  }

  if (backend == "limma_voom") {
    df$stat <- df$t
    df$logFC <- df$logFC
    df$pvalue <- df$P.Value
    df$padj <- df$adj.P.Val

  } else if (backend == "deseq2") {
    df$stat <- df$stat
    df$logFC <- df$log2FoldChange
    df$pvalue <- df$pvalue
    df$padj <- df$padj
  }

  # 保留核心列
  core_cols <- c("gene_symbol", "logFC", "stat", "pvalue", "padj")

  # 确保列存在
  missing <- setdiff(core_cols, colnames(df))
  if (length(missing) > 0) stop(sprintf("标准化失败，缺失列: %s", paste(missing, collapse=", ")))

  return(df[, c(core_cols, setdiff(colnames(df), core_cols))])
}

#' @title 构建表达数据包
#' @description 统一封装表达矩阵和元数据
#' @keywords internal
#' @noRd
.build_expr_bundle <- function(obj, backend) {

  if (is.null(obj)) {
    return(list(
      raw_counts = NULL,
      display_expr = NULL,
      sample_meta = NULL,
      gene_meta = NULL
    ))
  }

  if (backend == "limma_voom") {
    if (inherits(obj, "DGEList")) {
      raw_counts <- obj$counts
      sample_meta <- obj$samples
      gene_meta <- obj$genes
      display_expr <- edgeR::cpm(obj, log = TRUE)
    } else {
      raw_counts <- as.matrix(obj)
      sample_meta <- data.frame(row.names = colnames(obj))
      gene_meta <- NULL
      display_expr <- log2(raw_counts + 1)
    }

  } else if (backend == "deseq2") {
    raw_counts <- counts(obj, normalized = FALSE)
    sample_meta <- as.data.frame(colData(obj))
    gene_meta <- as.data.frame(rowData(obj))

    # 默认展示 log2(normalized counts + 1)
    norm_counts <- counts(obj, normalized = TRUE)
    display_expr <- log2(norm_counts + 1)
  }

  return(list(
    raw_counts = raw_counts,
    display_expr = display_expr,
    sample_meta = sample_meta,
    gene_meta = gene_meta,
    dge_list = if(backend == "limma_voom" && inherits(obj, "DGEList")) obj else NULL,
    dds_obj = if(backend == "deseq2") obj else NULL
  ))
}
