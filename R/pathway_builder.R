#' @title Build Custom Mixed MSigDB Dataframe (V10 修复版)
#' @export
build_custom_msigdb <- function(collection_list = list(c("C2", "CP:KEGG_LEGACY"),
                                                       c("C2", "CP:KEGG_MEDICUS"),
                                                       c("C5", "GO:BP"),
                                                       c("C5", "GO:CC"),
                                                       c("C5", "GO:MF")),
                                species = "Homo sapiens") {
  message("\U0001F527 Building custom mixed MSigDB database...")

  df_list <- lapply(collection_list, function(x) {
    cat <- x[1]
    subcat <- ifelse(length(x) > 1, x[2], "")
    message(sprintf("  -> Fetching %s : %s", cat, subcat))

    # 🌟 修复：参数更新为 collection 和 subcollection
    tmp <- msigdbr::msigdbr(species = species, collection = cat, subcollection = subcat)
    tmp %>% dplyr::select(gs_name, gene_symbol, gs_description, gs_url, gs_collection = gs_collection, gs_subcollection = gs_subcat)
  })

  final_df <- do.call(rbind, df_list)
  message(sprintf("✅ Built successfully! Total pathways: %d", length(unique(final_df$gs_name))))
  return(final_df)
}

#' @title Setup GSEA Environment (Polymorphic Engine - V10 修复版)
#' @export
setup_gsea_env <- function(fit, custom_df = NULL, gmt_file = NULL, custom_set_name = "Default_Set",
                           category = "H", subcategory = NULL, species = "Homo sapiens", expr_data = NULL) {

  # 复用自动对比组提取 (直接沿用您之前的提取方式或提取引擎)
  if (!inherits(fit, "MArrayLM")) stop("[Error] Input must be a limma 'fit' object.")
  c_names <- colnames(fit)
  parsed <- lapply(c_names, function(x) {
    p <- strsplit(x, "\\s*-\\s*")[[1]]
    if(length(p) == 2) c(p[1], p[2]) else c(x, "Unknown")
  })
  parsed_df <- do.call(rbind, parsed)
  contrasts_df <- tibble::tibble(ID = seq_along(c_names), Contrast_Name = c_names, Num = parsed_df[,1], Den = parsed_df[,2])

  term2gene <- NULL; meta_dict <- NULL; set_identifier <- custom_set_name

  if (!is.null(gmt_file)) {
    # (外部 GMT 处理，保持原样)
    message("\U0001F4E5 Mode: External GMT File detected.")
    lines <- readLines(gmt_file)
    gmt_list <- lapply(lines, function(l) strsplit(l, "\t")[[1]])
    meta_dict <- data.frame(ID = sapply(gmt_list, `[`, 1), Description = sapply(gmt_list, `[`, 2), URL = NA, Collection = "Custom_GMT", stringsAsFactors = FALSE)
    t2g_list <- lapply(gmt_list, function(x) data.frame(gs_name = x[1], gene_symbol = x[-(1:2)], stringsAsFactors = FALSE))
    term2gene <- do.call(rbind, t2g_list)

  } else if (!is.null(custom_df)) {
    # (自定义 df 处理，保持原样)
    message("\U0001F9EC Mode: Mixed Custom DataFrame detected.")
    term2gene <- custom_df %>% dplyr::select(gs_name, gene_symbol)
    if (!"gs_description" %in% colnames(custom_df)) custom_df$gs_description <- custom_df$gs_name
    if (!"gs_url" %in% colnames(custom_df)) custom_df$gs_url <- NA
    if (!"gs_collection" %in% colnames(custom_df)) custom_df$gs_collection <- "Mixed_Custom"
    meta_dict <- custom_df %>% dplyr::select(ID = gs_name, Description = gs_description, URL = gs_url, Collection = gs_collection) %>% dplyr::distinct(ID, .keep_all = TRUE)
    if (custom_set_name == "Default_Set") set_identifier <- "Mixed_Custom_Set"

  } else {
    # 🌟 原生模式：消除 category 警告
    if (!is.null(subcategory) && (is.na(subcategory) || nchar(trimws(subcategory)) == 0)) subcategory <- NULL

    if (is.null(subcategory)) {
      message(sprintf("\U0001F310 Mode: Native MSigDB (Collection: %s | No subcollection)", category))
      set_identifier <- category
    } else {
      message(sprintf("\U0001F310 Mode: Native MSigDB (Collection: %s | Subcollection: %s)", category, subcategory))
      set_identifier <- paste(category, subcategory, sep = "_")
    }

    # 🌟 修复：底层使用 collection 替代 category
    df <- msigdbr::msigdbr(species = species, collection = category, subcollection = subcategory)
    term2gene <- df %>% dplyr::select(gs_name, gene_symbol)
    meta_dict <- df %>% dplyr::select(ID = gs_name, Description = gs_description, URL = gs_url, Collection = gs_collection) %>% dplyr::distinct(ID, .keep_all = TRUE)
  }

  env_obj <- list(fit = fit, contrasts = contrasts_df, geneset = list(name = set_identifier, term2gene = term2gene, meta_dict = meta_dict), expr_data = expr_data)
  class(env_obj) <- "GseaEnv"
  return(env_obj)
}









#' @title 检查并打印 GSEA 环境状态
#' @description 提供一个精美格式化的控制台面板，展示胶囊内部状态。
#' @param gsea_env 由 setup_gsea_env() 创建的对象
#' @export
inspect_gsea_env <- function(gsea_env) {
  if (!inherits(gsea_env, "GseaEnv")) {
    stop("❌ [错误] 传入的对象不是标准的 GseaEnv 胶囊。")
  }

  gs_name <- gsea_env$geneset$name
  total_pw <- nrow(gsea_env$geneset$meta_dict)
  total_genes <- length(unique(gsea_env$geneset$term2gene$gene_symbol))
  c_df <- gsea_env$contrasts

  cat(rep("=", 65), "\n", sep = "")
  cat(sprintf("%-20s %s %-20s\n", "", "\U0001F9EC GSEA Environment Summary", ""))
  cat(rep("=", 65), "\n\n", sep = "")

  cat("\U0001F4E6 [1] GENE SET DATABASE\n")
  cat(sprintf("  • Set Name       : %s\n", gs_name))
  cat(sprintf("  • Total Pathways : %s\n", format(total_pw, big.mark = ",")))
  cat(sprintf("  • Unique Genes   : %s\n\n", format(total_genes, big.mark = ",")))

  cat("⚖️  [2] AVAILABLE CONTRASTS (Auto-extracted)\n")
  header <- sprintf("  %-4s | %-20s | %-15s | %-15s", "ID", "Contrast_Name", "POSITIVE (Left)", "NEGATIVE (Right)")
  cat(header, "\n")
  cat("  ", rep("-", nchar(header)-2), "\n", sep="")

  for (i in 1:nrow(c_df)) {
    cat(sprintf("  %-4s | %-20s | %-15s | %-15s\n",
                c_df$ID[i], c_df$Contrast_Name[i], c_df$Num[i], c_df$Den[i]))
  }
  cat("\n")

  cat("\U0001F321️  [3] EXPRESSION MATRIX STATUS\n")
  if (is.null(gsea_env$expr_data)) {
    cat("  ⚠️ 未检测到表达矩阵 (HTML 报告将只出折线图，不出热图)\n\n")
  } else {
    cat("  ✅ 表达矩阵已成功冻结入囊，随时可调取绘制热图！\n\n")
  }

  cat("\U0001F680 [4] NEXT STEP GUIDE\n")
  cat("  -> 启动全局计算与渲染引擎:\n")
  cat("     batch_parallel_gsea(gsea_env, ...)\n")
  cat(rep("=", 65), "\n", sep = "")

  return(invisible(gsea_env))
}

#' @export
print.GseaEnv <- function(x, ...) {
  inspect_gsea_env(x)
}
