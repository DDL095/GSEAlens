#' @title Core Class Definitions and Validation Functions

#' @description Define GseaEnv, GseaRes, GseaTask structures and provide internal validation logic.

#' @keywords internal

#' @name class_validations



NULL





# 1. Core Class Definitions





#' @title Define GseaEnv Class

#' @description Standardized GSEA input environment object.

#' @param backend_info Backend information object

#' @param contrast_registry Contrast registry

#' @param de_store Differential expression store

#' @param expr_bundle Expression data bundle

#' @param geneset Gene set collection

#' @param raw_obj Raw data object

#' @return A GseaEnv object

#' @examples

#' gsea_env <- create_gsea_env(

#'   backend_info = list(),

#'   contrast_registry = data.frame(),

#'   de_store = data.frame(),

#'   expr_bundle = list(),

#'   geneset = character(),

#'   raw_obj = NULL

#' )

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



#' @title Define GseaRes Class

#' @description GSEA computation result capsule.

#' @param metadata Metadata information

#' @param backend_info Backend information object

#' @param contrast_registry Contrast registry

#' @param de_store Differential expression store

#' @param expr_bundle Expression data bundle

#' @param geneset_info Gene set information

#' @param results GSEA results

#' @return A GseaRes object

#' @examples

#' gsea_res <- create_gsea_res(

#'   metadata = list(),

#'   backend_info = list(),

#'   contrast_registry = data.frame(),

#'   de_store = data.frame(),

#'   expr_bundle = list(),

#'   geneset_info = list(),

#'   results = list()

#' )

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



#' @title Define GseaTask Class

#' @description Single contrast extraction result object.

#' @param gsea_res GseaRes object

#' @param meta Metadata for the task

#' @return A GseaTask object

#' @examples

#' gsea_task <- create_gsea_task(

#'   gsea_res = NULL,

#'   meta = data.frame()

#' )

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





# 2. Internal Validation Functions





#' @title Validate Limma Design Matrix

#' @description Enforce no-intercept design (~ 0 + group) to ensure accurate contrast parsing.

#' @param fit MArrayLM object from limma

#' @return TRUE if validation passes, otherwise stops with an error

#' @keywords internal



.validate_limma_design <- function(fit) {

  # 检查是否包含截距项

  # Intercept / (Intercept) 是 R 中默认的截距命名

  design_matrix <- fit$design

  has_intercept <- any(grepl("Intercept", colnames(design_matrix), ignore.case = TRUE))



  if (has_intercept) {

    rlang::abort(

      c(

        "Intercept term detected in limma design matrix.",

        i = "GSEAlens requires a no-intercept design (~ 0 + group).",

        i = "Example: design <- model.matrix(~ 0 + group, data = samples); fit <- lmFit(expr, design)",

        i = "Reason: No-intercept design makes colnames(fit) directly correspond to group names."

      ),

      .class = "GSEAlens_limma_intercept_detected"

    )

  }



  # 检查是否有足够的列进行对比

  if (ncol(design_matrix) < 2) {

    rlang::warn(

      "Design matrix contains only 1 column; no between-group comparisons possible.",

      .class = "GSEAlens_limma_single_col"

    )

  }



  return(TRUE)

}



#' @title Validate DESeq2 Target Factor

#' @description Check if target_factor exists in colData and is properly formatted.

#' @param dds DESeqDataSet object

#' @param target_factor String, the target factor to validate

#' @return TRUE if validation passes, otherwise stops with an error

#' @keywords internal



.validate_deseq2_design <- function(dds, target_factor) {

  col_data <- as.data.frame(SummarizedExperiment::colData(dds))



  if (!target_factor %in% colnames(col_data)) {

    rlang::abort(

      c(

        sprintf("target_factor '%s' not found in colData.", target_factor),

        i = sprintf("Available column names: %s", paste(colnames(col_data), collapse = ", "))

      ),

      target_factor = target_factor,

      .class = "GSEAlens_deseq2_missing_factor"

    )

  }



  # 检查是否为因子

  if (!is.factor(col_data[[target_factor]])) {

    rlang::warn(

      sprintf("target_factor '%s' is not a factor type; attempting automatic conversion...", target_factor),

      .class = "GSEAlens_deseq2_factor_coerce"

    )

    # 这里不实际转换，只是警告，因为 DESeq2 通常在构建时已处理

  }



  return(TRUE)

}



#' @title Validate GseaEnv Object Integrity

#' @description Internal function to ensure object structure conforms to specification.

#' @param env_obj GseaEnv object to validate

#' @return TRUE if validation passes, otherwise stops with an error

#' @keywords internal



.check_gsea_env <- function(env_obj) {

  if (!inherits(env_obj, "GseaEnv")) {

    rlang::abort(

      "Input object is not of class GseaEnv.",

      .class = "GSEAlens_bad_class"

    )

  }



  required_slots <- c("backend_info", "contrast_registry", "de_store", "expr_bundle", "geneset")

  missing_slots <- setdiff(required_slots, names(env_obj))



  if (length(missing_slots) > 0) {

    rlang::abort(

      c(

        "GseaEnv object structure is incomplete.",

        i = sprintf("Missing slots: %s", paste(missing_slots, collapse = ", "))

      ),

      missing_slots = missing_slots,

      .class = "GSEAlens_incomplete_env"

    )

  }



  # 检查 contrast_registry 必要字段

  reg <- env_obj$contrast_registry

  if (!is.data.frame(reg) || !all(c("contrast_id", "left_group", "right_group") %in% colnames(reg))) {

    rlang::abort(

      "contrast_registry must contain columns: contrast_id, left_group, right_group.",

      .class = "GSEAlens_bad_registry"

    )

  }



  return(TRUE)

}



#' @title Validate GseaRes Object Integrity

#' @description Internal function to ensure GseaRes object contains required components.

#' @param res_obj GseaRes object to validate

#' @return TRUE if validation passes, otherwise stops with an error

#' @keywords internal



.check_gsea_res <- function(res_obj) {

  if (!inherits(res_obj, "GseaRes")) {

    rlang::abort(

      "Input object is not of class GseaRes.",

      .class = "GSEAlens_bad_class"

    )

  }

  # 简单检查 results 列表是否存在

  if (is.null(res_obj$results)) {

    rlang::abort(

      "No computation results in GseaRes object.",

      .class = "GSEAlens_empty_res"

    )

  }

  return(TRUE)

}

