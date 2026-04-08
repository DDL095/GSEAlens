# Section: Data Accessor Utilities ----

#' @title Data Accessor Utilities
#' @description Provides unified data extraction interface, isolating underlying
#'   object structure differences. Fixes DESeq2 sample name matching issues.
#' @keywords internal
#' @name utils-accessors
NULL



# Section: Expression Matrix Accessors ----


#' @title Get Expression Matrix
#' @description Extract expression matrix from GseaEnv or GseaRes object, supporting
#'   multiple normalization methods including raw counts, CPM, logCPM, VST, FPKM,
#'   TPM, and log-normalized counts.
#' @param obj A GseaEnv or GseaRes object.
#' @param type Character string specifying the type of expression values to return.
#'   Options include: \code{"default"} (default display matrix), \code{"raw"}
#'   (raw counts), \code{"cpm"} (counts per million), \code{"logcpm"}
#'   (log2-transformed CPM), \code{"vst"} (variance stabilizing transformation,
#'   DESeq2 only), \code{"fpkm"} (fragments per kilobase per million),
#'   \code{"tpm"} (transcripts per million), and \code{"lognorm"} (log2 of
#'   normalized counts).
#' @param ... Additional arguments passed to backend-specific methods.
#' @return A numeric matrix with genes as rows and samples as columns.
#' @export
#'
#' @examples
#' \dontrun{
#' # Get default expression matrix
#' expr <- get_expr_matrix(gsea_res)
#'
#' # Get log2 CPM values
#' expr <- get_expr_matrix(gsea_res, type = "logcpm")
#'
#' # Get variance stabilizing transformed values (DESeq2 only)
#' expr <- get_expr_matrix(gsea_res, type = "vst")
#' }
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
#' @description Core logic for extracting and computing expression matrices.
#'   Handles backend-specific computations (DESeq2, limma-voom) and normalizes
#'   sample ordering between expression matrices and metadata.
#' @keywords internal
#'
.get_expr_internal <- function(expr_bundle, backend_info, type = "default", ...) {

  # Normalize type names
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

  # Return default display matrix if requested
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

  # Retrieve raw data
  raw_counts <- expr_bundle$raw_counts
  sample_meta <- expr_bundle$sample_meta
  gene_lengths <- expr_bundle$gene_meta$length

  if (is.null(raw_counts)) {
    warning("Raw count matrix not available in object; cannot compute expression values dynamically.")
    return(NULL)
  }

  # Compute expression values based on type and backend
  res <- switch(type_normalized,
                "raw" = raw_counts,

                # CPM (generic)
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

                # VST (DESeq2-exclusive)
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

                # Log2 normalized counts (DESeq2)
                "lognorm" = {
                  if (backend == "deseq2" && !is.null(expr_bundle$dds_obj)) {
                    norm_counts <- DESeq2::counts(expr_bundle$dds_obj, normalized = TRUE)
                    log2(norm_counts + 1)
                  } else {
                    if (backend == "limma_voom" && !is.null(expr_bundle$dge_list)) {
                      edgeR::cpm(expr_bundle$dge_list, log = TRUE)
                    } else {
                      log2(t(t(raw_counts) / colSums(raw_counts)) * 1e6 + 1)
                    }
                  }
                },

                # Unknown type
                stop(sprintf("Unsupported expression value type: %s", type))
  )

  # Align expression matrix column order with sample metadata
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



# Section: Differential Expression Table Accessors ----


#' @title Get Differential Expression Table
#' @description Extract differential expression results for a specific contrast
#'   from the stored results. Automatically handles reverse contrasts (e.g.,
#'   B_vs_A) by retrieving the forward contrast and flipping the sign of logFC.
#' @param obj A GseaRes object containing differential expression results.
#' @param contrast_id Character string specifying the contrast ID (e.g., "A_vs_B").
#'   If the reverse contrast exists in storage, it will be retrieved and
#'   logFC/stat signs will be inverted.
#' @return A data frame containing at minimum: gene_symbol, logFC, and p-value
#'   columns. Additional columns (stat, pvalue, padj, t) may be present depending
#'   on the analysis backend.
#' @export
#'
#' @examples
#' \dontrun{
#' # Get results for a specific contrast
#' de_result <- get_de_table(gsea_res, contrast_id = "Treatment_vs_Control")
#'
#' # Reverse contrasts are automatically handled
#' de_result <- get_de_table(gsea_res, contrast_id = "Control_vs_Treatment")
#' # logFC and stat values will have inverted signs
#' }
get_de_table <- function(obj, contrast_id) {
  UseMethod("get_de_table")
}

#' @export
get_de_table.GseaRes <- function(obj, contrast_id) {

  if (contrast_id %in% names(obj$de_store)) {
    return(obj$de_store[[contrast_id]])
  }

  parts <- strsplit(contrast_id, "_vs_")[[1]]
  if (length(parts) == 2) {
    reverse_id <- paste(parts[2], parts[1], sep = "_vs_")

    if (reverse_id %in% names(obj$de_store)) {
      de_df <- obj$de_store[[reverse_id]]

      if ("logFC" %in% colnames(de_df)) {
        de_df$logFC <- -de_df$logFC
      }
      if ("stat" %in% colnames(de_df)) {
        de_df$stat <- -de_df$stat
      }
      if ("t" %in% colnames(de_df)) {
        de_df$t <- -de_df$t
      }

      message(sprintf("Auto-mapped reverse contrast: %s -> %s (logFC sign flipped)",
                      contrast_id, reverse_id))
      return(de_df)
    }
  }

  stop(sprintf("Contrast '%s' and its reverse '%s' not found in de_store. Available: %s",
               contrast_id,
               ifelse(length(parts)==2, paste(parts[2], parts[1], sep="_vs_"), "N/A"),
               paste(names(obj$de_store), collapse = ", ")))
}



# Section: Metadata Accessors ----


#' @title Get Sample Metadata
#' @description Extract and process sample metadata, handling DFrame conversion
#'   from DESeq2 and unifying group column naming conventions across backends.
#' @return A data frame with sample metadata. Includes a standardized 'group'
#'   column derived from the target factor (inferred from dds_obj@design for
#'   DESeq2 workflows). Sample row names are preserved for downstream matching.
#' @export
#'
#' @examples
#' \dontrun{
#' # Get processed sample metadata
#' meta <- get_sample_meta(gsea_res)
#' head(meta)
#' }
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
#' @description Automatically infers target_factor from dds_obj@design for
#'   DESeq2 workflows, handles DFrame conversion, and creates a unified
#'   'group' column for consistent downstream processing.
#' @param sample_meta Raw sample metadata (possibly DFrame).
#' @param expr_bundle Expression bundle containing raw_counts, dds_obj, and
#'   other components for row name recovery.
#' @keywords internal
#'
.process_sample_meta <- function(sample_meta, expr_bundle) {
  if (is.null(sample_meta)) {
    message("[SampleMeta] Input is NULL")
    return(NULL)
  }

  # Convert DFrame/DataFrame to data.frame
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

  # Recover row names if missing
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

  # Strategy A: Infer from dds_obj@design (DESeq2 workflow, highest priority)
  if (!is.null(expr_bundle$dds_obj)) {
    tryCatch({
      design_formula <- DESeq2::design(expr_bundle$dds_obj)
      if (inherits(design_formula, "formula")) {
        design_terms <- attr(terms(design_formula), "term.labels")
        if (length(design_terms) > 0) {
          inferred_factor <- tail(design_terms, 1)
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

  # Strategy B: Use stored target_factor
  if (is.null(target_factor) && !is.null(expr_bundle$target_factor)) {
    if (expr_bundle$target_factor %in% colnames(sample_meta)) {
      target_factor <- expr_bundle$target_factor
      message(sprintf("[SampleMeta] Using stored target_factor: '%s'", target_factor))
    }
  }

  # Strategy C: Fallback to first valid factor
  if (is.null(target_factor)) {
    is_limma_backend <- !is.null(expr_bundle$dge_list) && is.null(expr_bundle$dds_obj)

    if (!is_limma_backend) {
      warning("[SampleMeta] Could not determine target_factor from dds_obj, using first valid factor")
    }

    factor_cols <- names(sample_meta)[sapply(sample_meta, is.factor)]
    valid_factors <- factor_cols[sapply(sample_meta[factor_cols], function(x) length(levels(x)) > 1)]
    if (length(valid_factors) > 0) {
      target_factor <- valid_factors[1]
    }
  }

  # Create unified 'group' column
  if (!"group" %in% colnames(sample_meta)) {
    alt_names <- c("group", "Group", "condition", "Condition",
                   "treatment", "Treatment", "grp")
    found_name <- NULL

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

  # Convert numeric group codes to character labels
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

  # Final validation and reporting
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
#' @description Retrieve the registered contrasts used for differential expression
#'   analysis, including contrast definitions and metadata.
#' @param obj A GseaEnv or GseaRes object.
#' @return The contrast registry object stored in the input object.
#' @export
#'
#' @examples
#' \dontrun{
#' registry <- get_contrast_registry(gsea_res)
#' names(registry)
#' }
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
#' @description Retrieve gene set collection metadata including pathway definitions
#'   and annotation dictionary.
#' @param obj A GseaEnv or GseaRes object.
#' @return An object containing geneset_info with pathway definitions and
#'   metadata dictionary.
#' @export
#'
#' @examples
#' \dontrun{
#' gs_info <- get_geneset_info(gsea_res)
#' names(gs_info)
#' }
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



# Section: Additional Data Helper Functions ----


#' @title Read Additional Data
#' @description Read and validate additional data from CSV or RDS file. Supports
#'   gzip-compressed formats (.gz). The data frame must contain an 'ID' column
#'   as the primary key for joining with pathway results.
#' @param file_path Character string specifying the file path. Supported formats:
#'   \itemize{
#'     \item \code{.csv} / \code{.csv.gz} - Comma-separated values (optionally
#'       gzip compressed)
#'     \item \code{.rds} / \code{.rds.gz} - R serialized data (optionally gzip
#'       compressed)
#'   }
#' @return A data frame with validated ID column. Rows are deduplicated by ID,
#'   keeping the first occurrence of each unique ID.
#' @export
#'
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
#'
#' # Use with GseaRes
#' add_data <- read_addition_data("my_annotations.rds")
#' gsea_app <- launch_gsea_app(gsea_res, add_data = add_data)
#' }
read_addition_data <- function(file_path) {

  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }

  ext <- tools::file_ext(file_path)

  if (grepl("rds", ext, ignore.case = TRUE)) {
    data <- readRDS(file_path)
  } else if (grepl("csv", ext, ignore.case = TRUE) || ext == "gz") {
    if (grepl("\\.gz$", file_path)) {
      con <- gzfile(file_path, "rt", encoding = "UTF-8")
      data <- tryCatch({
        read.csv(con, stringsAsFactors = FALSE, check.names = FALSE)
      }, finally = {
        close(con)
      })
    } else {
      data <- read.csv(file_path, stringsAsFactors = FALSE,
                       check.names = FALSE, encoding = "UTF-8")
    }
  } else {
    stop("Unsupported file format: ", ext,
         "\nSupported formats: .csv, .csv.gz, .rds, .rds.gz")
  }

  if (!is.data.frame(data)) {
    stop("Input file must contain a data.frame, got: ", class(data)[1])
  }

  if (!"ID" %in% colnames(data)) {
    stop("Data frame must contain 'ID' column as the primary key. ",
         "Available columns: ", paste(colnames(data), collapse = ", "))
  }

  data$ID <- as.character(data$ID)

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

#' @title Create Additional Data RDS File
#' @description Convert a CSV file to standardized gzip-compressed RDS format
#'   with validation. The output file will be saved to the specified directory
#'   with a standardized naming convention.
#' @param csv_path Character string specifying the input CSV file path.
#' @param output_dir Character string specifying the output directory path.
#'   Default: current working directory.
#' @param output_name Character string specifying the output filename.
#'   Default: \code{"addition_data_gsealens.rds"}
#' @param overwrite Logical scalar indicating whether to overwrite an existing
#'   file. Default: \code{FALSE}.
#' @return Invisible \code{TRUE} on success. The RDS file path is printed to
#'   the console.
#' @export
#'
#' @examples
#' \dontrun{
#' # Convert CSV to RDS in current directory
#' creat_addition_data_rdsfile("my_annotations.csv")
#'
#' # Specify output directory and filename
#' creat_addition_data_rdsfile("annotations.csv",
#'                             output_dir = "/path/to/project",
#'                             output_name = "custom_name.rds")
#'
#' # Overwrite existing file
#' creat_addition_data_rdsfile("updated_annotations.csv",
#'                             overwrite = TRUE)
#' }
creat_addition_data_rdsfile <- function(csv_path,
                                        output_dir = getwd(),
                                        output_name = "addition_data_gsealens.rds",
                                        overwrite = FALSE) {

  if (!file.exists(csv_path)) {
    stop("CSV file not found: ", csv_path)
  }

  message("[creat_addition_data_rdsfile] Reading: ", csv_path)
  data <- read.csv(csv_path, stringsAsFactors = FALSE,
                   check.names = FALSE, encoding = "UTF-8")

  if (!is.data.frame(data)) {
    stop("CSV content must be a data.frame, got: ", class(data)[1])
  }

  if (!"ID" %in% colnames(data)) {
    stop("CSV must contain 'ID' column as the primary key. ",
         "Available columns: ", paste(colnames(data), collapse = ", "))
  }

  data$ID <- as.character(data$ID)

  if (any(duplicated(data$ID))) {
    n_dup <- sum(duplicated(data$ID))
    message("[creat_addition_data_rdsfile] Removed ", n_dup, " duplicate ID(s)")
    data <- data[!duplicated(data$ID), ]
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  output_path <- file.path(normalizePath(output_dir), output_name)

  if (file.exists(output_path)) {
    if (!overwrite) {
      stop("Output file already exists: ", output_path,
           "\nUse overwrite = TRUE to replace")
    }
    message("[creat_addition_data_rdsfile] Overwriting existing file...")
  }

  saveRDS(data, output_path, compress = "gzip")

  message("[creat_addition_data_rdsfile] Saved: ", output_path)
  message("   Dimensions: ", nrow(data), " rows x ", ncol(data), " columns")
  message("   Compression: gzip")
  message("   File size: ", format(file.info(output_path)$size / 1024, digits = 2), " KB")

  invisible(TRUE)
}

#' @title Create Additional Data Template
#' @description Generate a template CSV file containing all pathway IDs from
#'   the GseaRes object. Users can edit this template to add custom annotations
#'   such as Brief_Description, Custom_Category, Custom_Tag, and User_Notes.
#' @param gsea_res A GseaRes object (generated by \code{batch_calc_gsea}).
#' @param output_path Character string specifying the output template file path.
#'   Default: \code{"addition_template_gsealens.csv"}.
#' @param include_cols Character vector specifying additional columns to include
#'   besides ID. Default includes standard annotation columns.
#' @return Invisible \code{TRUE} on success.
#' @export
#'
#' @examples
#' \dontrun{
#' # Create template with all pathway IDs
#' create_addition_template(gsea_res, "my_template.csv")
#'
#' # After editing the CSV with custom annotations, convert to RDS
#' creat_addition_data_rdsfile("my_template.csv")
#'
#' # Then load in GSEAlens
#' gsea_app <- launch_gsea_app(gsea_res, "my_template.rds")
#' }
create_addition_template <- function(gsea_res,
                                     output_path = "addition_template_gsealens.csv",
                                     include_cols = c("Brief_Description",
                                                      "Custom_Category",
                                                      "Custom_Tag",
                                                      "User_Notes")) {

  if (!inherits(gsea_res, "GseaRes")) {
    stop("Input must be a GseaRes object")
  }

  if (is.null(gsea_res$geneset_info$meta_dict)) {
    stop("GseaRes object missing geneset_info$meta_dict")
  }

  meta_dict <- gsea_res$geneset_info$meta_dict
  all_ids <- unique(meta_dict$ID)

  template <- data.frame(
    ID = sort(all_ids),
    stringsAsFactors = FALSE
  )

  for (col in include_cols) {
    template[[col]] <- NA_character_
  }

  output_dir <- dirname(output_path)
  if (output_dir != "" && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

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



# Section: Gene Symbol Processing Utilities ----
# Cross-species, cross-backend gene name handling utilities
# Compatible with both limma-voom and DESeq2 pipelines

# Housekeeping gene probes for detection (uppercase for matching)
DETECTION_PROBES <- c(
  # Animal housekeeping genes
  "GAPDH", "ACTB", "ACTIN", "TUBA", "TUBB", "UBC", "UBB", "RPLP0", "RPS18",
  "RPL13A", "PPIA", "YWHAZ", "HPRT1", "B2M",
  # Plant housekeeping genes
  "ACTIN2", "ACTIN7", "TUB8", "EF1A", "UBQ10", "UBQ11",
  # Eukaryotic common
  "18S", "28S", "RPS5", "RPSA"
)


#' @title Detect Gene Symbol Column
#' @description Automatically detect the gene symbol column in gene metadata.
#'   Uses housekeeping gene probe matching for robust cross-species detection.
#' @param gene_meta A data frame containing gene metadata.
#' @return The detected column name if successful, or \code{NULL} if detection
#'   fails.
#' @keywords internal
#'
.detect_gene_symbol_column <- function(gene_meta) {
  if (is.null(gene_meta) || !is.data.frame(gene_meta) || nrow(gene_meta) == 0) {
    return(NULL)
  }

  # Strategy 1: Probe matching (most robust)
  for (col_name in colnames(gene_meta)) {
    col_data <- gene_meta[[col_name]]
    if (!is.character(col_data) && !is.factor(col_data)) next
    if (length(col_data) == 0) next

    col_upper <- toupper(as.character(col_data))
    matched <- sum(col_upper %in% DETECTION_PROBES, na.rm = TRUE)

    if (matched >= 2) {
      message(sprintf("[GeneDetector] Detected gene column '%s' by probe matching (matched %d housekeeping genes)",
                      col_name, matched))
      return(col_name)
    }
  }

  # Strategy 2: Fallback to common column names
  fallback_names <- c(
    "SYMBOL", "symbol", "Gene", "gene_name", "gene",
    "GeneSymbol", "gene_symbol", "gene_id", "GeneID",
    "Gene_Name", "geneID"
  )
  for (name in fallback_names) {
    if (name %in% colnames(gene_meta)) {
      message(sprintf("[GeneDetector] Using fallback column '%s'", name))
      return(name)
    }
  }

  warning("[GeneDetector] Failed to detect gene symbol column in gene_meta")
  return(NULL)
}


#' @title Rebuild Gene Symbol Case Mapping
#' @description Reconstruct gene symbol case mapping from differential expression
#'   results for display purposes.
#' @param gsea_res A GseaRes object.
#' @param contrast_id Character string specifying the contrast ID.
#' @return A named character vector where names are uppercase gene symbols and
#'   values are the original case. Returns \code{NULL} if mapping cannot be
#'   constructed.
#' @keywords internal
#'
.rebuild_symbol_map <- function(gsea_res, contrast_id) {
  de_df <- gsea_res$de_store[[contrast_id]]
  if (is.null(de_df) || !"gene_symbol" %in% colnames(de_df)) {
    warning(sprintf("[SymbolMap] de_store not available for contrast '%s'", contrast_id))
    return(NULL)
  }

  gene_syms <- de_df$gene_symbol
  gene_syms <- gene_syms[!is.na(gene_syms) & gene_syms != ""]

  if (length(gene_syms) == 0) {
    warning("[SymbolMap] No valid gene symbols found in de_store")
    return(NULL)
  }

  symbol_map <- setNames(gene_syms, toupper(gene_syms))
  message(sprintf("[SymbolMap] Built mapping for %d genes (contrast: %s)",
                  length(symbol_map), contrast_id))
  return(symbol_map)
}


#' @title Get Display Gene Symbol
#' @description Retrieve the original case gene symbol based on the case mapping.
#' @param gene A gene identifier (case-insensitive).
#' @param symbol_map A case mapping vector generated by \code{.rebuild_symbol_map}.
#' @return The original case gene symbol, or the input if no mapping exists.
#' @keywords internal
#'
.get_display_symbol <- function(gene, symbol_map) {
  if (is.null(symbol_map) || is.null(gene) || is.na(gene) || gene == "") {
    return(gene)
  }
  upper_gene <- toupper(as.character(gene))
  if (upper_gene %in% names(symbol_map)) {
    return(symbol_map[[upper_gene]])
  }
  return(gene)
}


#' @title Batch Get Display Gene Symbols
#' @param genes A character vector of gene identifiers.
#' @param symbol_map A case mapping vector generated by \code{.rebuild_symbol_map}.
#' @return A character vector of gene symbols in their original case.
#' @keywords internal
#'
.get_display_symbols <- function(genes, symbol_map) {
  sapply(genes, function(g) .get_display_symbol(g, symbol_map), USE.NAMES = FALSE)
}


#' @title Convert Expression Matrix Row Names to Gene Symbols
#' @description Automatically detect and convert expression matrix row names,
#'   supporting Ensembl ID to gene symbol conversion.
#' @param expr_mat An expression matrix.
#' @param gene_meta A gene metadata data frame containing both ID and symbol columns.
#' @param gsea_res A GseaRes object (used to access backend information).
#' @return The expression matrix with updated row names.
#' @keywords internal
#'
.convert_expr_rownames <- function(expr_mat, gene_meta, gsea_res) {
  if (is.null(expr_mat) || is.null(rownames(expr_mat))) {
    return(expr_mat)
  }

  current_rownames <- rownames(expr_mat)

  # Check if row names are Ensembl IDs
  is_ensembl <- any(grepl("^ENS[A-Z]{0,3}G[0-9]+", current_rownames, ignore.case = TRUE))

  if (!is_ensembl) {
    detected_col <- .detect_gene_symbol_column(gene_meta)
    if (!is.null(detected_col)) {
      message(sprintf("[ExprConvert] Row names already appear to be gene symbols"))
    }
    return(expr_mat)
  }

  message("[ExprConvert] Detected Ensembl IDs as row names, attempting conversion...")

  gene_symbol_col <- .detect_gene_symbol_column(gene_meta)

  if (is.null(gene_symbol_col)) {
    stop("Ensembl IDs detected but gene symbol column not found in gene_meta. ",
         "Cannot convert row names for heatmap plotting.")
  }

  gm <- gene_meta

  geneid_col <- NULL
  for (name in c("GeneID", "gene_id", "geneID", "ENSEMBL", "ensembl_id")) {
    if (name %in% colnames(gm)) {
      geneid_col <- name
      break
    }
  }

  if (is.null(geneid_col)) {
    stop("Ensembl IDs detected but GeneID column not found in gene_meta. ",
         "Available columns: ", paste(colnames(gm), collapse = ", "))
  }

  ensembl_ids <- gm[[geneid_col]]
  gene_symbols <- gm[[gene_symbol_col]]
  names(gene_symbols) <- ensembl_ids

  gene_symbols <- gene_symbols[!duplicated(names(gene_symbols))]

  converted <- gene_symbols[current_rownames]
  valid_mask <- !is.na(converted) & converted != ""

  n_converted <- sum(valid_mask)
  n_total <- length(current_rownames)

  if (n_converted == 0) {
    stop("Conversion failed: no genes matched between expr_mat rownames and gene_meta$",
         geneid_col)
  }

  rownames(expr_mat)[valid_mask] <- converted[valid_mask]

  message(sprintf("[ExprConvert] Successfully converted %d/%d genes (%.1f%%)",
                  n_converted, n_total, 100*n_converted/n_total))

  return(expr_mat)
}
