#' @title Data Accessor Utilities
#' @description Provides unified data extraction interface, isolating underlying object structure differences. Fixes DESeq2 sample name matching issues.
#' @keywords internal
#' @name utils-accessors
NULL


# 1. 表达矩阵访问器


#' @title Get Expression Matrix
#' @description Extract expression matrix from GseaEnv or GseaRes, supporting multiple normalization methods.
#' @param obj GseaEnv or GseaRes object
#' @param type Character. Options: "default" (default display), "raw", "cpm", "logcpm", "vst", "fpkm", "tpm".
#' @param ... Additional arguments
#' @return Expression matrix (genes x samples)
#' @export
#'
get_expr_matrix <- function(obj, type = "default", ...) {
  UseMethod("get_expr_matrix")
}

#' @export

get_expr_matrix.GseaEnv <- function(obj, type = "default", ...) {
  .get_expr_internal(obj$expr_bundle, obj$backend_info, type, ...)
}

#' @export

get_expr_matrix.GseaRes <- function(obj, type = "default", ...) {
  .get_expr_internal(obj$expr_bundle, obj$backend_info, type, ...)
}

#' @title Internal Expression Matrix Extraction Logic
#' @description Fixes DESeq2 backend expression matrix and sample metadata matching issues.
#' @keywords internal


.get_expr_internal <- function(expr_bundle, backend_info, type = "default", ...) {

  # 1. 标准化类型名称
  type_normalized <- tolower(type)
  type_aliases <- c(
    "log2cpm" = "logcpm",
    "cpm" = "cpm",
    "log2tpm" = "logtpm",
    "tpm" = "tpm",
    "vst" = "vst",
    "log2fpkm" = "logfpkm",
    "fpkm" = "fpkm",
    "lognorm" = "lognorm",
    "default" = "default",
    "raw" = "raw"
  )

  if (type_normalized %in% names(type_aliases)) {
    type_normalized <- type_aliases[[type_normalized]]
  }

  backend <- backend_info$backend

  # 2. 如果请求默认展示矩阵
  if (type_normalized == "default" && !is.null(expr_bundle$display_expr)) {
    expr_mat <- expr_bundle$display_expr
    if (!is.null(expr_bundle$sample_meta)) {
      if (!is.null(rownames(expr_bundle$sample_meta))) {
        common_samples <- intersect(colnames(expr_mat), rownames(expr_bundle$sample_meta))
        if (length(common_samples) > 0) {
          expr_mat <- expr_mat[, common_samples, drop = FALSE]
        }
      }
    }
    return(expr_mat)
  }

  # 3. 获取原始数据
  raw_counts <- expr_bundle$raw_counts
  sample_meta <- expr_bundle$sample_meta
  gene_lengths <- expr_bundle$gene_meta$length  # 假设有基因长度信息

  if (is.null(raw_counts)) {
    warning("Raw count matrix not available in object; cannot compute expression values dynamically.")
    return(NULL)
  }

  # 4. 严格根据后端类型分发
  res <- switch(type_normalized,
                "raw" = raw_counts,

                # CPM 相关（通用）
                "cpm" = {
                  if (backend == "limma_voom" && !is.null(expr_bundle$dge_list)) {
                    edgeR::cpm(expr_bundle$dge_list, log = FALSE)
                  } else if (backend == "deseq2" && !is.null(expr_bundle$dds_obj)) {
                    counts <- DESeq2::counts(expr_bundle$dds_obj, normalized = FALSE)
                    t(t(counts) / colSums(counts)) * 1e6
                  } else {
                    t(t(raw_counts) / colSums(raw_counts)) * 1e6
                  }
                },

                "logcpm" = {
                  if (backend == "limma_voom" && !is.null(expr_bundle$dge_list)) {
                    edgeR::cpm(expr_bundle$dge_list, log = TRUE)
                  } else if (backend == "deseq2" && !is.null(expr_bundle$dds_obj)) {
                    counts <- DESeq2::counts(expr_bundle$dds_obj, normalized = FALSE)
                    log2(t(t(counts) / colSums(counts)) * 1e6 + 1)
                  } else {
                    log2(t(t(raw_counts) / colSums(raw_counts)) * 1e6 + 1)
                  }
                },

                # 🔧 FPKM计算（需要基因长度）
                "fpkm" = {
                  if (is.null(gene_lengths)) {
                    warning("FPKM calculation requires gene length information (gene_meta$length); falling back to CPM.")
                    t(t(raw_counts) / colSums(raw_counts)) * 1e6
                  } else {
                    # FPKM = (counts / gene_length_kb) / (total_counts / 1e6)
                    gene_lengths_kb <- gene_lengths / 1000
                    rpm <- t(t(raw_counts) / colSums(raw_counts)) * 1e6
                    rpm / gene_lengths_kb
                  }
                },

                "logfpkm" = {
                  if (is.null(gene_lengths)) {
                    warning("logFPKM calculation requires gene length information; falling back to logCPM.")
                    log2(t(t(raw_counts) / colSums(raw_counts)) * 1e6 + 1)
                  } else {
                    gene_lengths_kb <- gene_lengths / 1000
                    rpm <- t(t(raw_counts) / colSums(raw_counts)) * 1e6
                    fpkm <- rpm / gene_lengths_kb
                    log2(fpkm + 1)
                  }
                },

                # VST (DESeq2 专属，严格检查)
                "vst" = {
                  if (backend != "deseq2") {
                    # 🔧 关键修复：Limma流程中如果请求VST，给出明确警告并回退
                    warning(sprintf("VST (variance stabilizing transformation) is a DESeq2-exclusive method; current backend is '%s', falling back to logCPM.", backend))
                    if (backend == "limma_voom" && !is.null(expr_bundle$dge_list)) {
                      edgeR::cpm(expr_bundle$dge_list, log = TRUE)
                    } else {
                      log2(t(t(raw_counts) / colSums(raw_counts)) * 1e6 + 1)
                    }
                  } else if (!is.null(expr_bundle$vst_matrix)) {
                    expr_bundle$vst_matrix
                  } else if (!is.null(expr_bundle$dds_obj)) {
                    tryCatch({
                      SummarizedExperiment::assay(DESeq2::vst(expr_bundle$dds_obj, blind = FALSE))
                    }, error = function(e) {
                      warning("VST calculation failed, falling back to logCPM: ", e$message)
                      counts <- DESeq2::counts(expr_bundle$dds_obj, normalized = FALSE)
                      log2(t(t(counts) / colSums(counts)) * 1e6 + 1)
                    })
                  } else {
                    stop("VST matrix was not pre-computed and cannot be computed from dds_obj.")
                  }
                },

                # log2 Normalized counts (DESeq2)
                "lognorm" = {
                  if (backend == "deseq2" && !is.null(expr_bundle$dds_obj)) {
                    norm_counts <- DESeq2::counts(expr_bundle$dds_obj, normalized = TRUE)
                    log2(norm_counts + 1)
                  } else {
                    # Limma中也支持：使用logCPM作为替代
                    if (backend == "limma_voom" && !is.null(expr_bundle$dge_list)) {
                      edgeR::cpm(expr_bundle$dge_list, log = TRUE)
                    } else {
                      log2(t(t(raw_counts) / colSums(raw_counts)) * 1e6 + 1)
                    }
                  }
                },

                # 未知类型
                stop(sprintf("Unsupported expression value type: %s", type))
  )

  # 确保表达矩阵列名与sample_meta行名匹配
  if (!is.null(sample_meta) && !is.null(rownames(sample_meta))) {
    common_samples <- intersect(colnames(res), rownames(sample_meta))
    if (length(common_samples) == 0) {
      warning("Expression matrix column names do not match sample metadata row names! Please check sample identifiers.")
    } else {
      res <- res[, common_samples, drop = FALSE]
    }
  }

  return(res)
}


# 2. 差异分析表访问器


#' @title Get Differential Expression Table
#' @description Extract results for a specific contrast from de_store.
#' @param obj GseaEnv or GseaRes object
#' @param contrast_id Character. Contrast ID (e.g., "A_vs_B").
#' @return data.frame (containing gene_symbol, logFC, stat, pvalue, padj)
#' @export
get_de_table <- function(obj, contrast_id) {
  UseMethod("get_de_table")
}
#' @export
get_de_table.GseaRes <- function(obj, contrast_id) {
  # 🔧 修复3：检查是否为反向对比（如 PREA_vs_IBAA）
  # 如果是反向，找到正向对比并翻转logFC符号

  if (contrast_id %in% names(obj$de_store)) {
    # 正向对比直接返回
    return(obj$de_store[[contrast_id]])
  }

  # 🔧 新增：尝试解析反向对比
  parts <- strsplit(contrast_id, "_vs_")[[1]]
  if (length(parts) == 2) {
    # 构建反向ID（交换左右）
    reverse_id <- paste(parts[2], parts[1], sep = "_vs_")

    if (reverse_id %in% names(obj$de_store)) {
      # 找到正向数据，翻转符号
      de_df <- obj$de_store[[reverse_id]]

      # 🆕 新增：翻转logFC和stat符号（关键！）
      if ("logFC" %in% colnames(de_df)) {
        de_df$logFC <- -de_df$logFC
      }
      if ("stat" %in% colnames(de_df)) {
        de_df$stat <- -de_df$stat
      }
      if ("t" %in% colnames(de_df)) {  # limma的t统计量
        de_df$t <- -de_df$t
      }

      message(sprintf("Auto-mapped reverse contrast: %s -> %s (logFC sign flipped)",
                      contrast_id, reverse_id))
      return(de_df)
    }
  }

  # 如果都找不到，报错
  stop(sprintf("Contrast '%s' and its reverse '%s' not found in de_store. Available: %s",
               contrast_id,
               ifelse(length(parts)==2, paste(parts[2], parts[1], sep="_vs_"), "N/A"),
               paste(names(obj$de_store), collapse = ", ")))
}


# 3. 元数据访问器（关键修复区域）


#' @title Get Sample Metadata
#' @description Fixes DESeq2 backend: unifies group column name to 'group', ensuring row names match expression matrix.
#' @export
get_sample_meta <- function(obj) {
  UseMethod("get_sample_meta")
}

#' @export

get_sample_meta.GseaEnv <- function(obj) {
  .process_sample_meta(obj$expr_bundle$sample_meta, obj$expr_bundle)
}

#' @export

get_sample_meta.GseaRes <- function(obj) {
  .process_sample_meta(obj$expr_bundle$sample_meta, obj$expr_bundle)
}

#' @title Internal Sample Metadata Processing
#' @description Unified sample metadata processing: fixes row names, unifies group column names.
#' @keywords internal

.process_sample_meta <- function(sample_meta, expr_bundle) {
  if (is.null(sample_meta)) return(NULL)

  # 转换为data.frame（如果是tibble或其他）
  if (!is.data.frame(sample_meta)) {
    sample_meta <- as.data.frame(sample_meta)
  }

  # 关键修复1：确保有rownames（来自原始colData的行名）
  if (is.null(rownames(sample_meta))) {
    # 尝试从raw_counts的colnames获取
    if (!is.null(expr_bundle$raw_counts)) {
      if (ncol(expr_bundle$raw_counts) == nrow(sample_meta)) {
        rownames(sample_meta) <- colnames(expr_bundle$raw_counts)
        message("[Accessor] Sample metadata row names restored from expression matrix column names")
      }
    }
  }

  # 关键修复2：统一分组列名（DESeq2使用"分组"，limma使用"group"）
  # 优先使用"group"，如果不存在则查找其他可能的列名
  if (!"group" %in% colnames(sample_meta)) {
    # 常见的中文分组列名
    alt_names <- c("分组", "Group", "condition", "Condition", "treatment", "Treatment")
    found_name <- NULL
    for (name in alt_names) {
      if (name %in% colnames(sample_meta)) {
        found_name <- name
        break
      }
    }

    if (!is.null(found_name)) {
      sample_meta$group <- sample_meta[[found_name]]
      message(sprintf("[Accessor] Mapped group column '%s' to 'group'", found_name))
    } else {
      # 如果没有找到，尝试推断第一个factor列作为分组
      factor_cols <- names(sample_meta)[sapply(sample_meta, is.factor)]
      if (length(factor_cols) > 0) {
        sample_meta$group <- sample_meta[[factor_cols[1]]]
        warning(sprintf("[Accessor] Standard group column not found; using '%s' as group", factor_cols[1]))
      }
    }
  }

  # 关键修复3：确保group列是factor类型
  if ("group" %in% colnames(sample_meta) && !is.factor(sample_meta$group)) {
    sample_meta$group <- as.factor(sample_meta$group)
  }

  return(sample_meta)
}

#' @title Get Contrast Registry
#' @export
get_contrast_registry <- function(obj) {
  UseMethod("get_contrast_registry")
}

#' @export

get_contrast_registry.GseaEnv <- function(obj) {
  obj$contrast_registry
}

#' @export

get_contrast_registry.GseaRes <- function(obj) {
  obj$contrast_registry
}

#' @title Get Gene Set Information
#' @export
get_geneset_info <- function(obj) {
  UseMethod("get_geneset_info")
}

#' @export

get_geneset_info.GseaRes <- function(obj) {
  obj$geneset_info
}

#' @export

get_geneset_info.GseaEnv <- function(obj) {
  obj$geneset
}
