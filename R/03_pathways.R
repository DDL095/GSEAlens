#' @title Build GSEA Pathway Database with Interactive Selection
#' @description Provides an interactive or automated gene set selection interface
#'   based on msigdbr, generating standardized pathway objects with intelligent
#'   semantic tagging. Supports dynamic handling of unknown future collections
#'   and includes C9 perturbation signature compatibility.
#' @param species Character string specifying the species, default is "Homo sapiens".
#'   Supports all species available in msigdbr.
#' @param auto_select Optional parameter for automated selection. Can be:
#'   \itemize{
#'     \item NULL: Launch interactive menu (default)
#'     \item "ALL": Load entire MSigDB database
#'     \item Integer vector: Index numbers (e.g., c(17, 26))
#'     \item Character vector: Collection names (e.g., c("H", "C2:CP:KEGG_LEGACY"))
#'   }
#' @return A list object containing:
#'   \itemize{
#'     \item TERM2GENE: Data frame mapping pathway IDs to gene symbols
#'     \item meta_dict: Metadata dictionary with descriptions, URLs, and collection info
#'     \item SuperTag: Intelligent batch identifier for the selected collections
#'     \item collections_used: Data frame of selected collection metadata
#'   }
#' @export
#' @examples
#' \dontrun{
#' # Interactive mode
#' pathways <- build_gsea_pathways()
#'
#' # Automated selection by name
#' pathways <- build_gsea_pathways(auto_select = c("H", "C5:GO:BP"))
#'
#' # Load all collections
#' pathways <- build_gsea_pathways(auto_select = "ALL")
#' }
#' @importFrom dplyr arrange left_join select distinct mutate all_of bind_rows filter
#' @importFrom msigdbr msigdbr_collections msigdbr
build_gsea_pathways <- function(species = "Homo sapiens", auto_select = NULL) {

  # 动态拉取可用的基因集集合列表
  avail_colls <- msigdbr::msigdbr_collections() %>%
    dplyr::arrange(gs_collection, gs_subcollection)

  # 白名单字典：包含标准集合与C9修复
  dict <- data.frame(
    gs_collection = c(
      "C1", "C2", "C2", "C2", "C2", "C2", "C2", "C2", "C2",
      "C3", "C3", "C3", "C3", "C4", "C4", "C4",
      "C5", "C5", "C5", "C5", "C6", "C7", "C7", "C8", "C9", "H"
    ),
    gs_subcollection = c(
      "", "CGP", "CP", "CP:BIOCARTA", "CP:KEGG_LEGACY", "CP:KEGG_MEDICUS",
      "CP:PID", "CP:REACTOME", "CP:WIKIPATHWAYS", "MIR:MIRDB", "MIR:MIR_LEGACY",
      "TFT:GTRD", "TFT:TFT_LEGACY", "3CA", "CGN", "CM",
      "GO:BP", "GO:CC", "GO:MF", "HPO", "", "IMMUNESIGDB", "VAX", "", "", ""
    ),
    short_tag = c(
      "C1Pos", "CGP", "CP", "BioC", "KeggL", "KeggM",
      "PID", "Reac", "Wiki", "MirDB", "MirL",
      "TFTG", "TFTL", "3CA", "CGN", "CM",
      "GoBP", "GoCC", "GoMF", "HPO", "C6Onc", "ImmS", "VAX", "C8Cell", "C9Pert", "Hal"
    ),
    description = c(
      "C1 Positional Gene Sets",
      "C2 Chemical and Genetic Perturbations",
      "C2 Canonical Pathways (All)",
      "C2 BioCarta Pathways",
      "C2 KEGG Legacy",
      "C2 KEGG Medicus",
      "C2 PID Pathways",
      "C2 Reactome",
      "C2 WikiPathways",
      "C3 miRNA Targets (MIRDB)",
      "C3 miRNA Legacy",
      "C3 TF Targets (GTRD)",
      "C3 TF Legacy",
      "C4 3CA Cancer Features",
      "C4 Cancer Gene Neighborhoods",
      "C4 Cancer Modules",
      "C5 GO Biological Process",
      "C5 GO Cellular Component",
      "C5 GO Molecular Function",
      "C5 Human Phenotype Ontology",
      "C6 Oncogenic Signatures",
      "C7 Immunologic Signatures",
      "C7 Vaccine Responses",
      "C8 Cell Type Signatures",
      "C9 Computational Perturbation Signatures",
      "H Hallmark Gene Sets"
    ),
    stringsAsFactors = FALSE
  )

  # 数据清洗：处理 NO_SUB 标记并合并字典
  avail_colls$gs_subcollection_clean <- gsub(".*:NO_SUB", "", avail_colls$gs_subcollection)
  menu_df <- dplyr::left_join(avail_colls, dict,
                              by = c("gs_collection", "gs_subcollection_clean" = "gs_subcollection"))

  # 动态截断机制：为未知的新集合生成临时标签
  missing_idx <- is.na(menu_df$short_tag)
  if (any(missing_idx)) {
    menu_df$short_tag[missing_idx] <- paste0(
      menu_df$gs_collection[missing_idx], "_",
      substr(gsub("[^A-Za-z]", "", menu_df$gs_subcollection[missing_idx]), 1, 4)
    )
    menu_df$description[missing_idx] <- paste0("New Collection: ", menu_df$gs_subcollection[missing_idx])
  }

  # 生成组合名称用于匹配
  menu_df$combo_name <- paste0(
    menu_df$gs_collection,
    ifelse(menu_df$gs_subcollection_clean == "", "",
           paste0(":", menu_df$gs_subcollection_clean))
  )

  # 交互界面逻辑
  if (is.null(auto_select)) {
    message("\n", rep("=", 70))
    message(sprintf("Welcome to GSEAlens Pathway Wizard (%s)", species))
    message(rep("=", 70))
    message("Available Gene Set Collections:")
    message(sprintf("%-4s %-8s %-20s %s", "ID", "Tag", "Collection", "Description"))
    message(rep("-", 70))

    for (i in 1:nrow(menu_df)) {
      cat(sprintf("[%2d] %-6s | %-18s | %s\n",
                  i, menu_df$short_tag[i], menu_df$combo_name[i], menu_df$description[i]))
    }

    user_input <- readline(prompt = "\nEnter selection (comma-separated, e.g., 17,26): ")
    selected_idx <- as.integer(unlist(strsplit(user_input, "[, ]+")))
    selected_idx <- selected_idx[!is.na(selected_idx) &
                                   selected_idx >= 1 &
                                   selected_idx <= nrow(menu_df)]

    if (length(selected_idx) == 0) {
      stop("Invalid input: No valid selection provided.")
    }
  } else {
    # 自动选择模式
    if (length(auto_select) == 1 && toupper(auto_select) == "ALL") {
      message("Loading complete MSigDB database (this may include thousands of pathways)...")
      selected_idx <- 1:nrow(menu_df)
    } else if (is.character(auto_select)) {
      selected_idx <- match(auto_select, menu_df$combo_name)
      if (any(is.na(selected_idx))) {
        missing <- auto_select[is.na(selected_idx)]
        stop(sprintf("Collection not found: %s. Please check spelling.",
                     paste(missing, collapse = ", ")))
      }
    } else if (is.numeric(auto_select)) {
      selected_idx <- auto_select
      if (any(selected_idx < 1 | selected_idx > nrow(menu_df))) {
        stop(sprintf("Index out of range. Valid range: 1-%d", nrow(menu_df)))
      }
    } else {
      stop("Invalid auto_select parameter type.")
    }
  }

  selected_rows <- menu_df[selected_idx, ]

  # 智能批次标签生成
  if (length(auto_select) == 1 && toupper(auto_select) == "ALL") {
    super_tag <- "ALL_MSigDB_Global"
  } else {
    super_tag <- ifelse(
      nrow(selected_rows) <= 4,
      paste(selected_rows$short_tag, collapse = "_"),
      sprintf("Mix%d_%s", nrow(selected_rows), selected_rows$short_tag[1])
    )
  }

  message(sprintf("\nSelected %d collection(s). Batch Tag: [%s]",
                  nrow(selected_rows), super_tag))

  # 核心提取：使用 collection/subcollection 参数（废弃旧的 category 系统）
  pathway_list <- lapply(1:nrow(selected_rows), function(i) {
    c_coll <- selected_rows$gs_collection[i]
    c_sub <- selected_rows$gs_subcollection[i]

    # 判断是否为无子集合的情况
    if (grepl("NO_SUB", c_sub) || c_sub == "") {
      msigdbr::msigdbr(species = species, collection = c_coll)
    } else {
      msigdbr::msigdbr(species = species, collection = c_coll, subcollection = c_sub)
    }
  })

  all_pathways <- dplyr::bind_rows(pathway_list)

  # 动态列名检测：兼容新旧版 msigdbr（gs_cat vs gs_collection）
  cat_col_name <- if ("gs_cat" %in% colnames(all_pathways)) "gs_cat" else "gs_collection"
  subcat_col_name <- if ("gs_subcat" %in% colnames(all_pathways)) "gs_subcat" else "gs_subcollection"

  # 构建 TERM2GENE 映射表
  TERM2GENE <- all_pathways %>%
    dplyr::select(gs_name, gene_symbol)

  # 构建元数据字典
  TERM2NAME <- all_pathways %>%
    dplyr::select(
      ID = gs_name,
      Description = gs_description,
      URL = gs_url,
      Collection = dplyr::all_of(cat_col_name),
      Subcollection = dplyr::all_of(subcat_col_name)
    ) %>%
    dplyr::distinct(ID, .keep_all = TRUE) %>%
    dplyr::mutate(
      # 生成标准格式的组合名称（如 C2:CP:KEGG_LEGACY）
      Combo_Name = ifelse(
        Subcollection == "" | is.na(Subcollection),
        Collection,
        paste0(Collection, ":", Subcollection)
      )
    )

  # 返回结果列表
  return(list(
    TERM2GENE = TERM2GENE,
    meta_dict = TERM2NAME,
    SuperTag = super_tag,
    collections_used = selected_rows
  ))
}
