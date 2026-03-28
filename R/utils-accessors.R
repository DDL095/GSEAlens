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

                "fpkm" = {
                  if (is.null(gene_lengths)) {
                    warning("FPKM calculation requires gene length information (gene_meta$length); falling back to CPM.")
                    t(t(raw_counts) / colSums(raw_counts)) * 1e6
                  } else {
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
  # 如果是反向，找到正向对比并翻转logFC符号

  if (contrast_id %in% names(obj$de_store)) {
    # 正向对比直接返回
    return(obj$de_store[[contrast_id]])
  }

  parts <- strsplit(contrast_id, "_vs_")[[1]]
  if (length(parts) == 2) {
    # 构建反向ID（交换左右）
    reverse_id <- paste(parts[2], parts[1], sep = "_vs_")

    if (reverse_id %in% names(obj$de_store)) {
      # 找到正向数据，翻转符号
      de_df <- obj$de_store[[reverse_id]]

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


# 3. 元数据访问器


#' @title Get Sample Metadata
#' @description Fixes DESeq2 backend: handles DFrame properly, unifies group column names
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
#' @keywords internal
.process_sample_meta <- function(sample_meta, expr_bundle) {
  if (is.null(sample_meta)) {
    message("[SampleMeta] Input is NULL")
    return(NULL)
  }

  if (inherits(sample_meta, "DFrame") || inherits(sample_meta, "DataFrame")) {
    message("[SampleMeta] Converting DFrame to data.frame")
    # 保存原始行名和维度信息
    original_rownames <- rownames(sample_meta)
    original_ncol <- ncol(sample_meta)

    # 转换为 data.frame
    sample_meta <- as.data.frame(sample_meta, stringsAsFactors = FALSE)

    # 如果行名丢失但原始有行名，恢复它们
    if ((is.null(rownames(sample_meta)) || all(rownames(sample_meta) == ""))
        && !is.null(original_rownames)) {
      rownames(sample_meta) <- original_rownames
      message(sprintf("[SampleMeta] Restored %d rownames after DFrame conversion",
                      length(original_rownames)))
    }

    # 确保列数一致（DFrame 转换有时会出问题）
    if (ncol(sample_meta) != original_ncol) {
      warning(sprintf("[SampleMeta] Column count mismatch after conversion: %d vs %d",
                      ncol(sample_meta), original_ncol))
    }
  } else if (!is.data.frame(sample_meta)) {
    # 其他非 data.frame 类型也转换
    sample_meta <- as.data.frame(sample_meta, stringsAsFactors = FALSE)
  }

  # 确保有 rownames（样本ID）
  if (is.null(rownames(sample_meta)) || all(rownames(sample_meta) == "")) {
    message("[SampleMeta] Row names missing, attempting recovery...")

    # 尝试从 raw_counts 恢复
    if (!is.null(expr_bundle$raw_counts)) {
      if (ncol(expr_bundle$raw_counts) == nrow(sample_meta)) {
        rownames(sample_meta) <- colnames(expr_bundle$raw_counts)
        message("[SampleMeta] Recovered rownames from raw_counts column names")
      } else {
        message(sprintf("[SampleMeta] Dimension mismatch: sample_meta (%d) vs raw_counts (%d)",
                        nrow(sample_meta), ncol(expr_bundle$raw_counts)))
      }
    }

    # 尝试从 dds_obj 恢复（DESeq2）
    if ((is.null(rownames(sample_meta)) || all(rownames(sample_meta) == ""))
        && !is.null(expr_bundle$dds_obj)) {
      dds_colnames <- colnames(expr_bundle$dds_obj)
      if (length(dds_colnames) == nrow(sample_meta)) {
        rownames(sample_meta) <- dds_colnames
        message("[SampleMeta] Recovered rownames from dds_obj")
      }
    }
  }

  # 统一分组列名
  if (!"group" %in% colnames(sample_meta)) {
    alt_names <- c("group", "Group", "condition", "Condition",
                   "treatment", "Treatment", "grp")
    found_name <- NULL

    for (name in alt_names) {
      if (name %in% colnames(sample_meta)) {
        found_name <- name
        break
      }
    }

    if (!is.null(found_name)) {
      sample_meta$group <- sample_meta[[found_name]]
      message(sprintf("[SampleMeta] Mapped column '%s' to 'group'", found_name))
    } else {
      # 尝试推断第一个因子列
      factor_cols <- names(sample_meta)[sapply(sample_meta, is.factor)]
      if (length(factor_cols) > 0) {
        # 选择水平数大于1的因子
        valid_factors <- factor_cols[sapply(sample_meta[factor_cols], function(x) length(levels(x)) > 1)]
        if (length(valid_factors) > 0) {
          sample_meta$group <- as.character(sample_meta[[valid_factors[1]]])
          warning(sprintf("[SampleMeta] Using first valid factor '%s' as group", valid_factors[1]))
        }
      }
    }
  }

  # 处理数值型分组编码
  if ("group" %in% colnames(sample_meta)) {
    if (is.numeric(sample_meta$group) || is.integer(sample_meta$group)) {
      message("[SampleMeta] Converting numeric group codes to labels...")

      # 尝试从 dds_obj 获取原始因子信息
      converted <- FALSE
      if (!is.null(expr_bundle$dds_obj)) {
        tryCatch({
          original_colData <- as.data.frame(SummarizedExperiment::colData(expr_bundle$dds_obj))
          factor_cols <- names(original_colData)[sapply(original_colData, is.factor)]

          for (fc in factor_cols) {
            if (length(levels(original_colData[[fc]])) > 1) {
              group_levels <- levels(original_colData[[fc]])
              # 匹配样本
              sample_match <- match(rownames(sample_meta), rownames(original_colData))
              if (any(!is.na(sample_match))) {
                group_codes <- as.numeric(original_colData[[fc]])[sample_match]
                mapped_groups <- group_levels[group_codes]
                if (any(!is.na(mapped_groups))) {
                  sample_meta$group <- mapped_groups
                  message(sprintf("[SampleMeta] Converted using dds_obj$colData[['%s']]", fc))
                  converted <- TRUE
                  break
                }
              }
            }
          }
        }, error = function(e) {
          message(sprintf("[SampleMeta] Error accessing dds_obj: %s", e$message))
        })
      }

      if (!converted) {
        warning("[SampleMeta] Could not convert numeric codes, converting to character as-is")
        sample_meta$group <- as.character(sample_meta$group)
      }
    }

    # 最终确保 group 是字符型
    if (is.factor(sample_meta$group)) {
      sample_meta$group <- as.character(sample_meta$group)
    }
  }

  # 最终验证
  if (is.null(rownames(sample_meta)) || all(rownames(sample_meta) == "")) {
    warning("[SampleMeta] CRITICAL: Failed to recover rownames, sample matching will fail!")
  } else {
    message(sprintf("[SampleMeta] Final: %d samples, groups: %s",
                    nrow(sample_meta),
                    paste(unique(sample_meta$group), collapse = ", ")))
  }

  return(sample_meta)
}

#' @title .process_sample_meta_simple
#' @description Unified sample metadata processing: fixes row names, unifies group column names,
#'   and converts numeric group codes to character labels.
#' @keywords internal

.process_sample_meta_simple <- function(sample_meta, expr_bundle) {
  if (is.null(sample_meta)) return(NULL)

  # 转换为data.frame（如果是tibble或其他）
  if (!is.data.frame(sample_meta)) {
    sample_meta <- as.data.frame(sample_meta)
  }

  # 关键修复1：确保有rownames（来自原始colData的行名）
  if (is.null(rownames(sample_meta))) {
    if (!is.null(expr_bundle$raw_counts)) {
      if (ncol(expr_bundle$raw_counts) == nrow(sample_meta)) {
        rownames(sample_meta) <- colnames(expr_bundle$raw_counts)
        message("[Accessor] Sample metadata row names restored from expression matrix column names")
      }
    }
  }

  # 关键修复2：统一分组列名（DESeq2使用"分组"，limma使用"group"）
  if (!"group" %in% colnames(sample_meta)) {
    alt_names <- c("group", "Group", "condition", "Condition", "treatment", "Treatment")
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

  if ("group" %in% colnames(sample_meta)) {
    # 检查 group 列是否为数值型
    if (is.numeric(sample_meta$group) || is.integer(sample_meta$group)) {
      message("[Accessor] Detected numeric group codes, attempting to convert to character labels...")

      # 尝试从原始 colData 中获取因子的 levels
      # DESeq2 的 colData 中，分组通常是 factor 类型
      if (!is.null(expr_bundle$dds_obj)) {
        # 从 DDS 对象获取原始分组信息
        original_colData <- as.data.frame(SummarizedExperiment::colData(expr_bundle$dds_obj))
        # 查找可能的分组列（因子类型）
        factor_cols_in_coldata <- names(original_colData)[sapply(original_colData, is.factor)]

        if (length(factor_cols_in_coldata) > 0) {
          # 使用第一个因子列作为分组
          for (fc in factor_cols_in_coldata) {
            if (length(levels(original_colData[[fc]])) > 1) {
              # 获取因子 levels 并映射
              group_levels <- levels(original_colData[[fc]])
              group_codes <- as.numeric(original_colData[[fc]])
              # 创建映射：数值编码 → 字符标签
              mapped_groups <- group_levels[group_codes]
              names(mapped_groups) <- rownames(original_colData)
              # 应用到 sample_meta（按行名匹配）
              sample_meta$group <- mapped_groups[rownames(sample_meta)]
              message(sprintf("[Accessor] Successfully converted numeric codes to character labels using '%s'", fc))
              break
            }
          }
        }
      }

      # 如果上述方法失败，尝试从 dge_list 获取
      if (is.numeric(sample_meta$group) && !is.null(expr_bundle$dge_list)) {
        original_samples <- expr_bundle$dge_list$samples
        if ("group" %in% colnames(original_samples) && is.factor(original_samples$group)) {
          group_levels <- levels(original_samples$group)
          group_codes <- as.numeric(original_samples$group)
          mapped_groups <- group_levels[group_codes]
          names(mapped_groups) <- rownames(original_samples)
          sample_meta$group <- mapped_groups[rownames(sample_meta)]
          message("[Accessor] Converted numeric codes using dge_list$group")
        }
      }

      # 最终检查：如果仍然是数值，给出警告
      if (is.numeric(sample_meta$group)) {
        warning("[Accessor] Could not convert numeric group codes to character labels. Please check your colData.")
      }
    }

    # 确保group列是factor类型
    if (!is.factor(sample_meta$group)) {
      sample_meta$group <- as.factor(sample_meta$group)
    }
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
