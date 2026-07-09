# Use rlang::`%||%` (imported via utils-imports.R @importFrom rlang %||%).
# This avoids the package-local semantics divergence (length-0 / all-NA fallback)
# and ensures consistent NA handling across the whole codebase. (B-4 IRON FIX, 2026-07-01)

#' @title Backend Data Extractor

#' @description Extract standardized data structures from limma or DESeq2 objects.

#' @keywords internal

#' @name backends



NULL





# 1. Limma 后端提取器





#' @title Extract Limma-Voom Data

#' @description Extract contrast groups and differential analysis results from MArrayLM object.

#' @param fit MArrayLM object (must have been processed with eBayes)

#' @param expr_data Optional. DGEList or expression matrix.

#' @return A list containing contrast_registry, de_store, and expr_bundle

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

            "\n[Limma Input Error] Design matrix column names detected as group names (e.g., 'GroupA', 'GroupB') instead of contrasts (e.g., 'GroupA - GroupB').\nGSEAlens requires a fit object with contrasts already defined.\nPlease use makeContrasts and contrasts.fit to define your comparison groups."

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

    contrast_id = vapply(parsed_contrasts, function(x) paste(x$left, x$right, sep = "_vs_"), character(1)),

    left_group = vapply(parsed_contrasts, `[[`, character(1), "left"),

    right_group = vapply(parsed_contrasts, `[[`, character(1), "right"),

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





#' @title Extract DESeq2 Data

#' @description Extract contrast groups and differential analysis results from DESeqDataSet object.

#' @param dds DESeqDataSet object (must have been processed with DESeq())

#' @param target_factor Character. Target factor. If NULL, automatically inferred as the last term in the design formula.

#' @return A list containing contrast_registry, de_store, and expr_bundle

#' @keywords internal



.extract_deseq2_data <- function(dds, target_factor = NULL) {

    # 1. 确定 target_factor

    design_formula <- DESeq2::design(dds)

    design_terms <- attr(terms(design_formula), "term.labels")



    if (is.null(target_factor)) {

    target_factor <- utils::tail(design_terms, 1)

    message(sprintf("[DESeq2] target_factor not specified, automatically inferred as: '%s'", target_factor))

    }



    # 校验 target_factor

    .validate_deseq2_design(dds, target_factor)



    # 2. 获取所有 levels

    col_data <- as.data.frame(SummarizedExperiment::colData(dds))

    factor_levels <- levels(col_data[[target_factor]])



    if (length(factor_levels) < 2) {

    stop(sprintf("Factor '%s' has fewer than 2 levels, cannot perform comparison.", target_factor))

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

    res <- tryCatch(

            {

        DESeq2::results(dds, contrast = contrast_vec)

            },

            error = function(e) {

        warning(sprintf("Error extracting contrast %s: %s", contrast_id, e$message))

        return(NULL)

            }

    )



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





#' @title Standardize Differential Expression Table Column Names

#' @description Unify column names from different backends to: gene_symbol, logFC, stat, pvalue, padj

#' @keywords internal

.standardize_de_columns <- function(df, backend) {

    # 添加基因名列（保持大小写敏感，但确保存在）

    if ("gene_symbol" %in% colnames(df)) {

    # 保持原样

    } else if ("SYMBOL" %in% colnames(df)) {

    df$gene_symbol <- df$SYMBOL

    } else if ("row.names" %in% colnames(df)) {

    df$gene_symbol <- df$row.names

    } else {

    df$gene_symbol <- rownames(df)

    }



    # 标准化统计列（根据后端类型）

    if (backend == "deseq2") {

    # DESeq2: stat, log2FoldChange, pvalue, padj

    df$stat <- df$stat %||% df$WaldStatistic %||% NA_real_

    df$logFC <- df$log2FoldChange %||% df$logFC %||% NA_real_

    df$pvalue <- df$pvalue %||% df$p.value %||% NA_real_

    df$padj <- df$padj %||% df$adj.P.Val %||% NA_real_

    } else if (backend == "limma_voom") {

    # Limma-Voom: t统计量作为stat, P.Value作为pvalue, adj.P.Val作为padj

    df$stat <- df$stat %||% df$t %||% NA_real_

    df$logFC <- df$logFC %||% NA_real_

    df$pvalue <- df$pvalue %||% df$P.Value %||% NA_real_

    df$padj <- df$padj %||% df$adj.P.Val %||% NA_real_

    }



    # 核心列检查

    core_cols <- c("gene_symbol", "logFC", "stat", "pvalue", "padj")

    missing <- setdiff(core_cols, colnames(df))

    if (length(missing) > 0) {

    stop(sprintf("Standardization failed, missing columns: %s", paste(missing, collapse = ", ")))

    }



    # 保留所有原始列，但确保核心列在前

    other_cols <- setdiff(colnames(df), core_cols)

    return(df[, c(core_cols, other_cols), drop = FALSE])

}



#' @title Build Expression Data Bundle

#' @description Uniformly encapsulate expression matrix and metadata

#' @keywords internal



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

    # IRON FIX (2026-07-01, D-4): namespace counts() to DESeq2::counts() so
    # we do not collide with Matrix::counts() or any other package export.
    raw_counts <- DESeq2::counts(obj, normalized = FALSE)

    gene_meta <- as.data.frame(rowData(obj))

    sample_meta <- as.data.frame(colData(obj))

    # 默认展示 log2(normalized counts + 1)

    norm_counts <- DESeq2::counts(obj, normalized = TRUE)

    display_expr <- log2(norm_counts + 1)

    }



    return(list(

    raw_counts = raw_counts,

    display_expr = display_expr,

    sample_meta = sample_meta,

    gene_meta = gene_meta,

    dge_list = if (backend == "limma_voom" && inherits(obj, "DGEList")) obj else NULL,

    dds_obj = if (backend == "deseq2") obj else NULL

    ))

}

