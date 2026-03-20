#' @title 数据访问器工具函数
#' @description 提供统一的数据提取接口，隔离底层对象结构差异。修复DESeq2样本名匹配问题。
#' @keywords internal
#' @name utils-accessors
NULL

# ==============================================================================
# 1. 表达矩阵访问器
# ==============================================================================

#' @title 获取表达矩阵
#' @description 从 GseaEnv 或 GseaRes 中提取表达矩阵，支持多种标准化方式。
#' @param obj GseaEnv 或 GseaRes 对象
#' @param type 字符串。可选值："default" (默认展示), "raw", "cpm", "logcpm", "vst", "fpkm", "tpm"。
#' @param ... 额外参数
#' @return 表达矩阵 (基因 x 样本)
#' @export
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

#' @title 内部表达矩阵提取逻辑
#' @description 修复DESeq2后端表达矩阵与样本元数据匹配问题
#' @keywords internal
.get_expr_internal <- function(expr_bundle, backend_info, type = "default", ...) {

  # 1. 标准化类型名称（处理别名）
  type_normalized <- tolower(type)
  type_aliases <- c(
    "log2cpm" = "logcpm",
    "cpm" = "cpm",
    "log2tpm" = "logtpm",
    "tpm" = "tpm",
    "vst" = "vst",
    "log2fpkm" = "logfpkm",
    "fpkm" = "fpkm",
    "default" = "default",
    "raw" = "raw"
  )

  if (type_normalized %in% names(type_aliases)) {
    type_normalized <- type_aliases[[type_normalized]]
  }

  # 2. 如果请求默认展示矩阵，且已存在，直接返回
  if (type_normalized == "default" && !is.null(expr_bundle$display_expr)) {
    # 确保列名与sample_meta行名匹配（关键修复）
    expr_mat <- expr_bundle$display_expr
    if (!is.null(expr_bundle$sample_meta)) {
      # 如果sample_meta有行名，确保表达矩阵列名与之匹配
      if (!is.null(rownames(expr_bundle$sample_meta))) {
        common_samples <- intersect(colnames(expr_mat), rownames(expr_bundle$sample_meta))
        if (length(common_samples) > 0) {
          expr_mat <- expr_mat[, common_samples, drop = FALSE]
        }
      }
    }
    return(expr_mat)
  }

  # 3. 否则根据类型动态计算
  raw_counts <- expr_bundle$raw_counts
  sample_meta <- expr_bundle$sample_meta

  # 如果 raw_counts 为空
  if (is.null(raw_counts)) {
    warning("对象中未包含原始计数矩阵，无法动态计算表达量。")
    return(NULL)
  }

  # 4. 根据后端类型分发计算逻辑
  backend <- backend_info$backend

  res <- switch(type_normalized,
                "raw" = raw_counts,

                # CPM 相关
                "cpm" = {
                  if (backend == "limma_voom" && !is.null(expr_bundle$dge_list)) {
                    edgeR::cpm(expr_bundle$dge_list, log = FALSE)
                  } else if (backend == "deseq2" && !is.null(expr_bundle$dds_obj)) {
                    # DESeq2: 手动计算CPM
                    counts <- DESeq2::counts(expr_bundle$dds_obj, normalized = FALSE)
                    t(t(counts) / colSums(counts)) * 1e6
                  } else {
                    # 通用计算
                    t(t(raw_counts) / colSums(raw_counts)) * 1e6
                  }
                },

                "logcpm" = {
                  if (backend == "limma_voom" && !is.null(expr_bundle$dge_list)) {
                    edgeR::cpm(expr_bundle$dge_list, log = TRUE)
                  } else if (backend == "deseq2" && !is.null(expr_bundle$dds_obj)) {
                    # DESeq2: log2(CPM+1)
                    counts <- DESeq2::counts(expr_bundle$dds_obj, normalized = FALSE)
                    log2(t(t(counts) / colSums(counts)) * 1e6 + 1)
                  } else {
                    log2(t(t(raw_counts) / colSums(raw_counts)) * 1e6 + 1)
                  }
                },

                # VST (DESeq2 专属)
                "vst" = {
                  if (backend != "deseq2") {
                    warning("VST 仅支持 DESeq2 后端，回退到 logCPM。")
                    log2(t(t(raw_counts) / colSums(raw_counts)) * 1e6 + 1)
                  } else if (!is.null(expr_bundle$vst_matrix)) {
                    expr_bundle$vst_matrix
                  } else if (!is.null(expr_bundle$dds_obj)) {
                    # 尝试从 dds_obj 计算 VST
                    tryCatch({
                      SummarizedExperiment::assay(DESeq2::vst(expr_bundle$dds_obj, blind = FALSE))
                    }, error = function(e) {
                      warning("VST 计算失败，回退到 logCPM: ", e$message)
                      log2(t(t(raw_counts) / colSums(raw_counts)) * 1e6 + 1)
                    })
                  } else {
                    stop("VST 矩阵未预计算且无法从 dds_obj 计算。")
                  }
                },

                # TPM/FPKM (需要基因长度)
                "tpm" = , "logtpm" = , "fpkm" = , "logfpkm" = {
                  stop(sprintf("%s 计算需要基因长度信息，当前未实现。", type))
                },

                # 未知类型
                stop(sprintf("不支持的表达量类型: %s", type))
  )

  # 关键修复：确保表达矩阵列名与sample_meta行名匹配（取交集）
  if (!is.null(sample_meta) && !is.null(rownames(sample_meta))) {
    common_samples <- intersect(colnames(res), rownames(sample_meta))
    if (length(common_samples) == 0) {
      warning("表达矩阵列名与样本元数据行名无匹配！请检查样本标识。")
    } else {
      res <- res[, common_samples, drop = FALSE]
    }
  }

  return(res)
}

# ==============================================================================
# 2. 差异分析表访问器
# ==============================================================================

#' @title 获取差异分析表
#' @description 从 de_store 中提取指定对比组的结果。
#' @param obj GseaEnv 或 GseaRes 对象
#' @param contrast_id 字符串。对比组 ID (如 "A_vs_B")。
#' @return data.frame (包含 gene_symbol, logFC, stat, pvalue, padj)
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

      message(sprintf("🔄 自动映射反向对比: %s -> %s (已翻转logFC符号)",
                      contrast_id, reverse_id))
      return(de_df)
    }
  }

  # 如果都找不到，报错
  stop(sprintf("对比组 '%s' 及其反向 '%s' 都不存在于 de_store 中。可用: %s",
               contrast_id,
               ifelse(length(parts)==2, paste(parts[2], parts[1], sep="_vs_"), "N/A"),
               paste(names(obj$de_store), collapse = ", ")))
}

# ==============================================================================
# 3. 元数据访问器（关键修复区域）
# ==============================================================================

#' @title 获取样本元数据
#' @description 修复DESeq2后端：统一分组列名为'group'，确保行名与表达矩阵匹配
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

#' @title 内部样本元数据处理
#' @description 统一处理样本元数据：修复行名、统一分组列名
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
        message("[Accessor] 已从表达矩阵列名恢复样本元数据行名")
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
      message(sprintf("[Accessor] 已将分组列 '%s' 映射为 'group'", found_name))
    } else {
      # 如果没有找到，尝试推断第一个factor列作为分组
      factor_cols <- names(sample_meta)[sapply(sample_meta, is.factor)]
      if (length(factor_cols) > 0) {
        sample_meta$group <- sample_meta[[factor_cols[1]]]
        warning(sprintf("[Accessor] 未找到标准分组列，使用 '%s' 作为分组", factor_cols[1]))
      }
    }
  }

  # 关键修复3：确保group列是factor类型
  if ("group" %in% colnames(sample_meta) && !is.factor(sample_meta$group)) {
    sample_meta$group <- as.factor(sample_meta$group)
  }

  return(sample_meta)
}

#' @title 获取对比组注册表
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

#' @title 获取基因集信息
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
