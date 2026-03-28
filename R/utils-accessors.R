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

#' @title Internal Sample Metadata Processing (DESeq2 Target Factor Aware)
#' @description Now automatically infers target_factor from dds_obj@design for DESeq2 workflows.
#' @keywords internal
.process_sample_meta <- function(sample_meta, expr_bundle) {
  if (is.null(sample_meta)) {
    message("[SampleMeta] Input is NULL")
    return(NULL)
  }

  if (inherits(sample_meta, "DFrame") || inherits(sample_meta, "DataFrame")) {
    message("[SampleMeta] Converting DFrame to data.frame")
    original_rownames <- rownames(sample_meta)
    original_ncol <- ncol(sample_meta)
    sample_meta <- as.data.frame(sample_meta, stringsAsFactors = FALSE)

    if ((is.null(rownames(sample_meta)) || all(rownames(sample_meta) == ""))
        && !is.null(original_rownames)) {
      rownames(sample_meta) <- original_rownames
      message(sprintf("[SampleMeta] Restored %d rownames after DFrame conversion",
                      length(original_rownames)))
    }

    if (ncol(sample_meta) != original_ncol) {
      warning(sprintf("[SampleMeta] Column count mismatch after conversion: %d vs %d",
                      ncol(sample_meta), original_ncol))
    }
  } else if (!is.data.frame(sample_meta)) {
    sample_meta <- as.data.frame(sample_meta, stringsAsFactors = FALSE)
  }
  if (is.null(rownames(sample_meta)) || all(rownames(sample_meta) == "")) {
    message("[SampleMeta] Row names missing, attempting recovery...")

    if (!is.null(expr_bundle$raw_counts)) {
      if (ncol(expr_bundle$raw_counts) == nrow(sample_meta)) {
        rownames(sample_meta) <- colnames(expr_bundle$raw_counts)
        message("[SampleMeta] Recovered rownames from raw_counts column names")
      }
    }

    if ((is.null(rownames(sample_meta)) || all(rownames(sample_meta) == ""))
        && !is.null(expr_bundle$dds_obj)) {
      dds_colnames <- colnames(expr_bundle$dds_obj)
      if (length(dds_colnames) == nrow(sample_meta)) {
        rownames(sample_meta) <- dds_colnames
        message("[SampleMeta] Recovered rownames from dds_obj")
      }
    }
  }

  target_factor <- NULL

  # 策略A: 从 dds_obj@design 推断（DESeq2 流程，最高优先级）
  if (!is.null(expr_bundle$dds_obj)) {
    tryCatch({
      design_formula <- DESeq2::design(expr_bundle$dds_obj)
      if (inherits(design_formula, "formula")) {
        design_terms <- attr(terms(design_formula), "term.labels")
        if (length(design_terms) > 0) {
          # 取设计公式的最后一个变量（通常是目标分组变量）
          inferred_factor <- tail(design_terms, 1)
          # 验证该变量确实存在于 sample_meta 中
          if (inferred_factor %in% colnames(sample_meta)) {
            target_factor <- inferred_factor
            message(sprintf("[SampleMeta] Inferred target_factor from dds_obj@design: '%s'",
                            target_factor))
          }
        }
      }
    }, error = function(e) {
      message(sprintf("[SampleMeta] Failed to infer from dds_obj@design: %s", e$message))
    })
  }

  # 策略B: 从 expr_bundle$target_factor 读取（如果之前已存储）
  if (is.null(target_factor) && !is.null(expr_bundle$target_factor)) {
    if (expr_bundle$target_factor %in% colnames(sample_meta)) {
      target_factor <- expr_bundle$target_factor
      message(sprintf("[SampleMeta] Using stored target_factor: '%s'", target_factor))
    }
  }

  # 策略C: 回退到原有逻辑（选择第一个有效因子）
  if (is.null(target_factor)) {
    warning("[SampleMeta] Could not determine target_factor from dds_obj, using first valid factor")
    factor_cols <- names(sample_meta)[sapply(sample_meta, is.factor)]
    valid_factors <- factor_cols[sapply(sample_meta[factor_cols], function(x) length(levels(x)) > 1)]
    if (length(valid_factors) > 0) {
      target_factor <- valid_factors[1]
    }
  }

  if (!"group" %in% colnames(sample_meta)) {
    alt_names <- c("group", "Group", "condition", "Condition",
                   "treatment", "Treatment", "grp")
    found_name <- NULL

    # 优先使用推断的 target_factor
    if (!is.null(target_factor) && target_factor %in% colnames(sample_meta)) {
      found_name <- target_factor
    } else {
      for (name in alt_names) {
        if (name %in% colnames(sample_meta)) {
          found_name <- name
          break
        }
      }
    }

    if (!is.null(found_name)) {
      sample_meta$group <- sample_meta[[found_name]]
      message(sprintf("[SampleMeta] Mapped column '%s' to 'group'", found_name))
    } else {
      # 最终回退：使用推断的 target_factor（即使不在 alt_names 中）
      if (!is.null(target_factor) && target_factor %in% colnames(sample_meta)) {
        sample_meta$group <- sample_meta[[target_factor]]
        message(sprintf("[SampleMeta] Using inferred factor '%s' as 'group'", target_factor))
      } else {
        factor_cols <- names(sample_meta)[sapply(sample_meta, is.factor)]
        if (length(factor_cols) > 0) {
          valid_factors <- factor_cols[sapply(sample_meta[factor_cols], function(x) length(levels(x)) > 1)]
          if (length(valid_factors) > 0) {
            sample_meta$group <- as.character(sample_meta[[valid_factors[1]]])
            warning(sprintf("[SampleMeta] Using first valid factor '%s' as group", valid_factors[1]))
          }
        }
      }
    }
  }


  if ("group" %in% colnames(sample_meta)) {
    if (is.numeric(sample_meta$group) || is.integer(sample_meta$group)) {
      message("[SampleMeta] Converting numeric group codes to character labels...")

      converted <- FALSE
      if (!is.null(expr_bundle$dds_obj)) {
        tryCatch({
          original_colData <- as.data.frame(SummarizedExperiment::colData(expr_bundle$dds_obj))
          factor_cols <- names(original_colData)[sapply(original_colData, is.factor)]

          for (fc in factor_cols) {
            if (length(levels(original_colData[[fc]])) > 1) {
              group_levels <- levels(original_colData[[fc]])
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

      if (!converted && !is.null(expr_bundle$dge_list)) {
        original_samples <- expr_bundle$dge_list$samples
        if ("group" %in% colnames(original_samples) && is.factor(original_samples$group)) {
          group_levels <- levels(original_samples$group)
          group_codes <- as.numeric(original_samples$group)
          mapped_groups <- group_levels[group_codes]
          names(mapped_groups) <- rownames(original_samples)
          sample_meta$group <- mapped_groups[rownames(sample_meta)]
          message("[SampleMeta] Converted using dge_list$group")
        }
      }

      if (is.numeric(sample_meta$group)) {
        warning("[SampleMeta] Could not convert numeric group codes to character labels.")
      }
    }

    if (!is.factor(sample_meta$group)) {
      sample_meta$group <- as.factor(sample_meta$group)
    }
  }


  if (is.null(rownames(sample_meta)) || all(rownames(sample_meta) == "")) {
    warning("[SampleMeta] CRITICAL: Failed to recover rownames, sample matching will fail!")
  } else {
    unique_groups <- if (!is.null(sample_meta$group)) {
      paste(unique(sample_meta$group), collapse = ", ")
    } else {
      "N/A"
    }
    message(sprintf("[SampleMeta] Final: %d samples, groups: %s",
                    nrow(sample_meta), unique_groups))
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



# Addition Data Helper Functions


#' @title Read Addition Data
#' @description Read and validate addition_data from CSV or RDS file.
#'   Supports compressed formats (.gz). The data frame must contain an 'ID'
#'   column as the primary key for joining with pathway results.
#' @param file_path Character. File path. Supports formats:
#'   \itemize{
#'     \item .csv / .csv.gz - Comma-separated values (optionally gzip compressed)
#'     \item .rds / .rds.gz - R serialized data (optionally gzip compressed)
#'   }
#' @return A data frame with validated ID column. Rows are deduplicated by ID.
#' @export
#' @examples
#' \dontrun{
#' # Read from CSV
#' add_data <- read_addition_data("pathway_annotations.csv")
#'
#' # Read from compressed CSV
#' add_data <- read_addition_data("pathway_annotations.csv.gz")
#'
#' # Read from RDS
#' add_data <- read_addition_data("pathway_annotations.rds")
#' }
read_addition_data <- function(file_path) {

  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }

  # Detect file extension
  ext <- tools::file_ext(file_path)

  # Read based on extension
  if (grepl("rds", ext, ignore.case = TRUE)) {
    # RDS format
    data <- readRDS(file_path)
  } else if (grepl("csv", ext, ignore.case = TRUE) || ext == "gz") {
    # CSV or compressed CSV
    if (grepl("\\.gz$", file_path)) {
      # Gzip compressed CSV
      con <- gzfile(file_path, "rt", encoding = "UTF-8")
      data <- tryCatch({
        read.csv(con, stringsAsFactors = FALSE, check.names = FALSE)
      }, finally = {
        close(con)
      })
    } else {
      # Plain CSV
      data <- read.csv(file_path, stringsAsFactors = FALSE,
                       check.names = FALSE, encoding = "UTF-8")
    }
  } else {
    stop("Unsupported file format: ", ext,
         "\nSupported formats: .csv, .csv.gz, .rds, .rds.gz")
  }

  # Validate: must be data.frame
  if (!is.data.frame(data)) {
    stop("Input file must contain a data.frame, got: ", class(data)[1])
  }

  # Validate: must have ID column
  if (!"ID" %in% colnames(data)) {
    stop("Data frame must contain 'ID' column as the primary key. ",
         "Available columns: ", paste(colnames(data), collapse = ", "))
  }

  # Ensure ID column is character type
  data$ID <- as.character(data$ID)

  # Remove duplicates, keeping first occurrence
  if (any(duplicated(data$ID))) {
    n_dup <- sum(duplicated(data$ID))
    message("[read_addition_data] Found ", n_dup,
            " duplicate ID(s), keeping first occurrence")
    data <- data[!duplicated(data$ID), ]
  }

  message("[read_addition_data] Loaded: ", nrow(data), " rows x ",
          ncol(data), " columns from ", basename(file_path))

  return(data)
}

#' @title Create Addition Data RDS File
#' @description Convert a CSV file to standardized gzip-compressed RDS format
#'   with validation. The output file will be saved to the specified directory
#'   with a standardized naming convention.
#' @param csv_path Character. Input CSV file path.
#' @param output_dir Character. Output directory path.
#'   Default: current working directory.
#' @param output_name Character. Output filename.
#'   Default: "addition_data_gsealens.rds"
#' @param overwrite Logical. Whether to overwrite existing file.
#'   Default: FALSE.
#' @return Invisible TRUE on success. The RDS file path is printed to console.
#' @export
#' @examples
#' \dontrun{
#' # Convert CSV to RDS in current directory
#' creat_addition_data_rdsfile("my_annotations.csv")
#'
#' # Specify output directory and filename
#' creat_addition_data_rdsfile("annotations.csv",
#'                             output_dir = "/path/to/project",
#'                             output_name = "custom_name.rds")
#' }
creat_addition_data_rdsfile <- function(csv_path,
                                        output_dir = getwd(),
                                        output_name = "addition_data_gsealens.rds") {

  # Validate input CSV exists
  if (!file.exists(csv_path)) {
    stop("CSV file not found: ", csv_path)
  }

  # Read CSV
  message("[creat_addition_data_rdsfile] Reading: ", csv_path)
  data <- read.csv(csv_path, stringsAsFactors = FALSE,
                   check.names = FALSE, encoding = "UTF-8")

  # Validate: must be data.frame
  if (!is.data.frame(data)) {
    stop("CSV content must be a data.frame, got: ", class(data)[1])
  }

  # Validate: must have ID column
  if (!"ID" %in% colnames(data)) {
    stop("CSV must contain 'ID' column as the primary key. ",
         "Available columns: ", paste(colnames(data), collapse = ", "))
  }

  # Ensure ID column is character type
  data$ID <- as.character(data$ID)

  # Remove duplicates
  if (any(duplicated(data$ID))) {
    n_dup <- sum(duplicated(data$ID))
    message("[creat_addition_data_rdsfile] Removed ", n_dup, " duplicate ID(s)")
    data <- data[!duplicated(data$ID), ]
  }

  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Construct output path
  output_path <- file.path(normalizePath(output_dir), output_name)

  # Check for existing file
  if (file.exists(output_path)) {
    if (!overwrite) {
      stop("Output file already exists: ", output_path,
           "\nUse overwrite = TRUE to replace")
    }
    message("[creat_addition_data_rdsfile] Overwriting existing file...")
  }

  # Save as gzip-compressed RDS
  saveRDS(data, output_path, compress = "gzip")

  message("[creat_addition_data_rdsfile] Saved: ", output_path)
  message("   Dimensions: ", nrow(data), " rows x ", ncol(data), " columns")
  message("   Compression: gzip")
  message("   File size: ", format(file.info(output_path)$size / 1024, digits = 2), " KB")

  invisible(TRUE)
}

#' @title Create Addition Data Template
#' @description Generate a template CSV file containing all pathway IDs from
#'   the GseaRes object. Users can edit this template to add custom annotations
#'   such as Brief_Description, Custom_Category, Custom_Tag, and User_Notes.
#' @param gsea_res GseaRes object (generated by batch_calc_gsea)
#' @param output_path Character. Output template file path.
#'   Default: "addition_template_gsealens.csv"
#' @param include_cols Character vector. Additional columns to include besides ID.
#'   Default includes standard annotation columns.
#' @return Invisible TRUE on success.
#' @export
#' @examples
#' \dontrun{
#' # Create template with all pathway IDs
#' create_addition_template(gsea_res, "my_template.csv")
#'
#' # After editing, convert to RDS
#' creat_addition_data_rdsfile("my_template.csv")
#' }
create_addition_template <- function(gsea_res,
                                     output_path = "addition_template_gsealens.csv",
                                     include_cols = c("Brief_Description",
                                                      "Custom_Category",
                                                      "Custom_Tag",
                                                      "User_Notes")) {

  # Validate GseaRes object
  if (!inherits(gsea_res, "GseaRes")) {
    stop("Input must be a GseaRes object")
  }

  # Extract all pathway IDs from geneset_info
  if (is.null(gsea_res$geneset_info$meta_dict)) {
    stop("GseaRes object missing geneset_info$meta_dict")
  }

  meta_dict <- gsea_res$geneset_info$meta_dict

  # Get unique pathway IDs
  all_ids <- unique(meta_dict$ID)

  # Create template data frame
  template <- data.frame(
    ID = sort(all_ids),
    stringsAsFactors = FALSE
  )

  # Add optional annotation columns
  for (col in include_cols) {
    template[[col]] <- NA_character_
  }

  # Ensure output directory exists
  output_dir <- dirname(output_path)
  if (output_dir != "" && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Write CSV (not compressed for easy editing)
  write.csv(template, output_path, row.names = FALSE, quote = TRUE,
            fileEncoding = "UTF-8")

  message("[create_addition_template] Template created: ", output_path)
  message("   Contains ", nrow(template), " pathway IDs")
  message("   Columns: ", paste(colnames(template), collapse = ", "))
  message("")
  message("Instructions:")
  message("  1. Edit the CSV file to fill in annotations")
  message("  2. Convert to RDS: creat_addition_data_rdsfile('", basename(output_path), "')")
  message("  3. Load in GSEAlens: launch_gsea_app(gsea_res, '",
          gsub("\\.csv$", ".rds", basename(output_path)), "')")

  invisible(TRUE)
}
