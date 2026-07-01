#' @title Extract Pathway Core Genes (Leading Edge)

#' @description Parse the core_enrichment field from GSEA results to extract the gene list.

#'   This is the fundamental data interface for all visualization modules in 03/04.

#' @param gsea_res_obj GseaRes object or single GSEA result (list(status=, data=, ...))

#' @param pathway_id Character, pathway ID

#' @return Character vector, core genes (uppercase normalized)

#' @export

#'

#' @examples

#' # Load a pre-computed GseaRes and extract a task

#' gsea_res <- readRDS(system.file(

#'   "extdata", "precomputed_gseares.rds", package = "GSEAlens"

#' ))

#' task <- extract_gsea_task(gsea_res, contrast_id = "untrt_vs_trt")

#' res_df <- as.data.frame(task$gsea_res)

#' pw_id <- res_df$ID[1]

#' core_genes <- get_core_genes_for_pathway(task, pw_id)

#' head(core_genes)

get_core_genes_for_pathway <- function(gsea_res_obj, pathway_id) {

  # 缁熶竴澶勭悊杈撳叆锛氬彲鑳芥槸GseaRes瀵硅薄鎴栧凡鎻愬彇鐨則ask缁撴灉

  if (inherits(gsea_res_obj, "GseaRes")) {

    # 闇€瑕佹寚瀹歝ontrast_id锛岃繖閲岀畝鍖栧鐞嗭紝鍋囪浼犲叆鐨勬槸鍗曚釜瀵规瘮缁撴灉

    stop("For GseaRes object, please use extract_gsea_task to extract specific contrast first")

  }



  # 澶勭悊鎻愬彇鍚庣殑GseaTask瀵硅薄鎴栫洿鎺ョ殑缁撴灉鍒楄〃

  if (inherits(gsea_res_obj, "GseaTask")) {

    gsea_result <- gsea_res_obj$gsea_res

  } else if (is.list(gsea_res_obj) && !is.null(gsea_res_obj$data)) {

    # 鏉ヨ嚜results鍒楄〃鐨勭洿鎺ュ厓绱?

    if (gsea_res_obj$status != "Success") {

      return(character(0))

    }

    gsea_result <- gsea_res_obj$data

  } else {

    gsea_result <- gsea_res_obj

  }



  # 浠嶴4瀵硅薄涓彁鍙栫粨鏋滆〃

  if (methods::is(gsea_result, "gseaResult")) {

    res_df <- as.data.frame(gsea_result@result)

  } else {

    return(character(0))

  }



  # 鍖归厤閫氳矾ID

  row_idx <- which(res_df$ID == pathway_id)

  if (length(row_idx) == 0) {

    return(character(0))

  }



  # 瑙ｆ瀽core_enrichment瀛楁锛堜互/鍒嗛殧鐨勫熀鍥犲瓧绗︿覆锛?

  core_str <- as.character(res_df$core_enrichment[row_idx[1]])

  if (is.na(core_str) || core_str == "") {

    return(character(0))

  }



  # 鍒嗗壊骞舵竻娲?

  genes <- unlist(strsplit(core_str, "/"))

  genes <- toupper(trimws(genes))

  genes <- genes[genes != ""]



  return(unique(genes))

}



# ==============================================================================

# Section ----

# ==============================================================================



## Subsection ----



### Sub-subsection ----



#' @title Batch Extract Core Genes for Multiple Pathways

#' @description Batch data preparation for 03 module Network/UpSet/Chord visualizations.

#' @param gsea_task_obj GseaTask object (single contrast)

#' @param pathway_ids Character vector, pathway ID list

#' @return Named list, names = pathway_id, values = core_genes vector

#' @export

#'



#' @examples

#' # Load a pre-computed GseaRes and extract a task

#' gsea_res <- readRDS(system.file(

#'   "extdata", "precomputed_gseares.rds", package = "GSEAlens"

#' ))

#' task <- extract_gsea_task(gsea_res, contrast_id = "untrt_vs_trt")

#' res_df <- as.data.frame(task$gsea_res)

#' pw_ids <- head(res_df$ID, 3)

#' genes_list <- get_core_genes_list(task, pw_ids)

#' length(genes_list)

get_core_genes_list <- function(gsea_task_obj, pathway_ids) {

  if (!inherits(gsea_task_obj, "GseaTask")) {

    stop("Must pass a GseaTask object")

  }



  result <- lapply(pathway_ids, function(pid) {

    get_core_genes_for_pathway(gsea_task_obj, pid)

  })

  names(result) <- pathway_ids

  return(result)

}



# ==============================================================================

# Section ----

# ==============================================================================



## Subsection ----



### Sub-subsection ----



#' @title Calculate ORA Ratio (for DotPlot)

#' @description Calculate overlap ratio between pathway genes and differentially

#'   expressed genes.

#' @param pathway_genes Gene set of pathway (defined by TERM2GENE)

#' @param de_genes Differentially expressed genes (filtered by p-value threshold)

#' @param ratio_mode Character, "ora" (intersection/pathway size) or "leading"

#'   (intersection/DE size)

#' @return Numeric, ratio value

#' @export

#'



#' @examples

#' ratio <- calculate_overlap_ratio(c("GAPDH", "TP53", "MYC"), c("TP53", "MYC", "EGFR"))

calculate_overlap_ratio <- function(pathway_genes, de_genes, ratio_mode = c("ora", "leading")) {

  ratio_mode <- match.arg(ratio_mode)



  pathway_genes <- toupper(pathway_genes)

  de_genes <- toupper(de_genes)



  overlap <- intersect(pathway_genes, de_genes)



  if (ratio_mode == "ora") {

    # 缁忓吀ORA锛氫氦闆?/ 閫氳矾鎬诲熀鍥犳暟

    return(length(overlap) / length(pathway_genes))

  } else {

    # Leading Edge妯″紡锛氫氦闆?/ DE鍩哄洜鎬绘暟

    return(length(overlap) / length(de_genes))

  }

}



# ==============================================================================

# Section ----

# ==============================================================================



## Subsection ----



### Sub-subsection ----



#' @title Get Full Pathway Gene Set (TERM2GENE)

#' @description Extract complete gene list for a pathway from geneset_info.

#' @param gsea_res GseaRes object

#' @param pathway_id Pathway ID

#' @return Character vector

#' @export

#'



#' @examples

#' # Load a pre-computed GseaRes shipped with the package

#' gsea_res <- readRDS(system.file(

#'   "extdata", "precomputed_gseares.rds", package = "GSEAlens"

#' ))

#' pw_id <- gsea_res$geneset_info$term2gene$gs_name[1]

#' genes <- get_term_genes(gsea_res, pw_id)

#' head(genes)

get_term_genes <- function(gsea_res, pathway_id) {

  term2gene <- gsea_res$geneset_info$term2gene

  genes <- term2gene$gene_symbol[term2gene$gs_name == pathway_id]

  return(toupper(unique(genes)))

}



# ==============================================================================

# Section ----

# ==============================================================================



## Subsection ----



### Sub-subsection ----



#' @title Safe Parameter Validation with Fallback (Global Strategy)

#' @description Unified handling of invalid parameters (e.g., ncol <= 0) to ensure

#'   03/04 module stability.

#' @param value Input value

#' @param default Default value

#' @param min_val Minimum value (inclusive)

#' @param max_val Maximum value (inclusive, optional)

#' @param param_name Parameter name (for warning messages)

#' @return Validated value

#' @export

#'



#' @examples

#' workers <- validate_param(4, default = 2, min_val = 1, max_val = 8, param_name = "workers")

validate_param <- function(value, default, min_val = 1, max_val = NULL, param_name = "parameter") {

  # 澶勭悊NULL鎴朜A

  if (is.null(value) || (length(value) == 1 && is.na(value))) {

    message(sprintf("[Param Check] %s is NULL/NA, fallback to default: %s", param_name, default))

    return(default)

  }



  # 优先 is.numeric 判断，避免触发 as.numeric 的强制转换警告

  # 仅对非 numeric 输入（如字符串）才走 suppressWarnings 兜底分支

  val <- if (is.numeric(value)) {

    as.numeric(value[1])

  } else {

    parsed <- suppressWarnings(as.numeric(value[1]))

    if (is.na(parsed)) {

      message(sprintf("[Param Check] %s conversion failed, fallback to: %s", param_name, default))

      return(default)

    }

    parsed

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





# ==============================================================================

# Section ----

# ==============================================================================



## Subsection ----



### Sub-subsection ----



#' @title Build Pathway Similarity Network Edges (Safe Version)

#' @description Safely construct edge list for pathway network based on Jaccard

#'   similarity of core genes. Includes defensive checks for edge cases.

#'   Now includes Overlap Coefficient and Dice Coefficient for comprehensive

#'   similarity assessment.

#' @param core_genes_list Named list: pathway_id -> core_gene_vector

#' @param min_shared_genes Minimum number of shared genes between pathways to form an edge

#' @return Data frame with columns: from, to, shared, weight (Jaccard index),

#'   overlap_coef, dice_coef, shared_genes

#' @export

#'



#' @examples

#' # Build edges from two pathways sharing 2 core genes

#' core_genes <- list(

#'   PATHWAY_A = c("GENE_A", "GENE_B", "GENE_C"),

#'   PATHWAY_B = c("GENE_B", "GENE_C", "GENE_D")

#' )

#' edges <- build_edge_list_safely(core_genes, min_shared_genes = 2)

#' print(edges)

build_edge_list_safely <- function(core_genes_list, min_shared_genes = 2) {

  # ==== 闃插尽鎬ф鏌?====



  # 闃插尽1锛氭鏌ヨ緭鍏ユ槸鍚︿负 NULL 鎴栫┖鍒楄〃

  if (is.null(core_genes_list) || length(core_genes_list) == 0) {

    message("[build_edge_list_safely] Input core_genes_list is NULL or empty")

    return(NULL)

  }



  # 闃插尽2锛氱Щ闄ょ┖鐨勬牳蹇冨熀鍥犲悜閲?

  valid_idx <- which(vapply(core_genes_list, function(x) {

    length(x) > 0 && !all(is.na(x))

  }, logical(1)))



  if (length(valid_idx) == 0) {

    message("[build_edge_list_safely] All pathways have empty core genes")

    return(NULL)

  }



  core_genes_list <- core_genes_list[valid_idx]



  pathway_ids <- names(core_genes_list)

  n <- length(pathway_ids)



  # 闃插尽3锛氬鏋滃彧鏈?0 鎴?1 涓?pathway锛屾棤娉曞舰鎴愯竟

  if (n <= 1) {

    message(sprintf("[build_edge_list_safely] Only %d pathway(s), cannot form edges", n))

    return(NULL)

  }



  # 闃插尽4锛氬鏋?min_shared_genes <= 0锛岃缃粯璁ゅ€?

  if (is.null(min_shared_genes) || min_shared_genes < 1) {

    message("[build_edge_list_safely] min_shared_genes must be >= 1, using default 2")

    min_shared_genes <- 2

  }



  # ==== 杈规瀯寤?====



  # 棰勫垎閰嶈竟瀛樺偍锛堜娇鐢ㄥ垪琛ㄩ伩鍏嶅姩鎬佹墿灞曪級

  edges_list <- vector("list", length = n * (n - 1) / 2)

  edge_count <- 0



  # 鍙屽眰寰幆璁＄畻 Jaccard 鐩镐技搴?

  for (i in seq_len(n - 1)) {

    for (j in (i + 1):n) {

      p1 <- pathway_ids[i]

      p2 <- pathway_ids[j]



      genes1 <- core_genes_list[[p1]]

      genes2 <- core_genes_list[[p2]]



      # 璺宠繃浠讳綍鍖呭惈 NA 鐨勫熀鍥犲悜閲?

      if (any(is.na(genes1)) || any(is.na(genes2))) next



      # 澶у皬鍐欎笉鏁忔劅鐨勪氦闆?骞堕泦璁＄畻

      genes1_upper <- toupper(as.character(genes1))

      genes2_upper <- toupper(as.character(genes2))



      # 璁＄畻鍏变韩鍩哄洜鏁?

      shared_count <- length(intersect(genes1_upper, genes2_upper))



      # 妫€鏌ユ槸鍚︽弧瓒虫渶灏忓叡浜槇鍊?

      if (shared_count < min_shared_genes) next



      # 璁＄畻 Jaccard 鐩镐技搴?

      union_count <- length(union(genes1_upper, genes2_upper))

      jaccard <- if (union_count > 0) shared_count / union_count else 0



      # 璁＄畻 Overlap Coefficient (Simpson Index)

      # 鍏紡锛殀A鈭〣| / min(|A|, |B|)

      # 閫傚悎妫€娴嬪瓙闆嗗叧绯?

      min_size <- min(length(genes1_upper), length(genes2_upper))

      overlap_coef <- if (min_size > 0) shared_count / min_size else 0



      # 璁＄畻 Dice Coefficient

      # 鍏紡锛?|A鈭〣| / (|A| + |B|)

      # 瀵逛氦闆嗘洿鏁忔劅锛岃寖鍥?[0, 1]

      sum_size <- length(genes1_upper) + length(genes2_upper)

      dice_coef <- if (sum_size > 0) (2 * shared_count) / sum_size else 0



      # 鑾峰彇鍏变韩鍩哄洜鍒楄〃

      shared_genes_list <- intersect(genes1_upper, genes2_upper)



      edge_count <- edge_count + 1

      edges_list[[edge_count]] <- data.frame(

        from = p1,

        to = p2,

        shared = shared_count,

        weight = jaccard, # Jaccard index (primary weight)

        overlap_coef = overlap_coef, # Overlap coefficient

        dice_coef = dice_coef, # Dice coefficient

        shared_genes = I(list(shared_genes_list)),

        stringsAsFactors = FALSE

      )

    }

  }



  # ==== 杩斿洖缁撴灉 ====



  if (edge_count == 0) {

    message(

      "[build_edge_list_safely] No edges formed with current threshold (min_shared=",

      min_shared_genes, ")"

    )

    return(NULL)

  }



  # 鍚堝苟杈瑰垪琛?

  edge_df <- do.call(rbind, edges_list[seq_len(edge_count)])

  rownames(edge_df) <- NULL



  message(sprintf(

    "[build_edge_list_safely] Built %d edges from %d pathways",

    nrow(edge_df), n

  ))



  return(edge_df)

}

