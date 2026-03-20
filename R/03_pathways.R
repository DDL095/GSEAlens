#' @title 构建基因集数据库
#' @description 基于 msigdbr 提供交互式或自动化的基因集选择，生成标准化的基因集对象。
#' 支持动态截断未知的新集合，并生成智能语义化标签。
#' @param species 物种名称，默认 "Homo sapiens"。
#' @param auto_select 可选。可以是序号向量 (如 c(1, 5)) 或名称向量 (如 c("H", "C2:CP:KEGG_LEGACY"))。
#' 若为 NULL，则启动交互式菜单。若为 "ALL"，则载入全库。
#' @return 包含 TERM2GENE, meta_dict, SuperTag 的列表对象。
#' @export
build_gsea_pathways <- function(species = "Homo sapiens", auto_select = NULL) {

  # 动态拉取库
  avail_colls <- msigdbr::msigdbr_collections() %>%
    dplyr::arrange(gs_collection, gs_subcollection)

  # 白名单字典 (包含 C9 修复)
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
      "C1 位置基因集", "C2 化学和遗传微扰", "C2 经典通路(总)", "C2 BioCarta", "C2 KEGG(旧)", "C2 KEGG(医学)",
      "C2 PID 通路", "C2 Reactome", "C2 WikiPathways", "C3 miRNA 靶标", "C3 miRNA(旧)",
      "C3 TF 靶标(GTRD)", "C3 TF(旧)", "C4 3CA 癌症", "C4 癌症邻居", "C4 癌症模块",
      "C5 GO 生物过程(BP)", "C5 GO 细胞组分(CC)", "C5 GO 分子功能(MF)", "C5 人类表型(HPO)", "C6 肿瘤特征", "C7 免疫特征", "C7 疫苗特征", "C8 细胞类型", "C9 计算扰动特征", "H 癌症标志"
    ),
    stringsAsFactors = FALSE
  )

  # 数据清洗与合并
  avail_colls$gs_subcollection_clean <- gsub(".*:NO_SUB", "", avail_colls$gs_subcollection)
  menu_df <- dplyr::left_join(avail_colls, dict,
                              by = c("gs_collection", "gs_subcollection_clean" = "gs_subcollection"))

  # 动态截断机制 (应对未来未知集合)
  missing_idx <- is.na(menu_df$short_tag)
  if (any(missing_idx)) {
    menu_df$short_tag[missing_idx] <- paste0(
      menu_df$gs_collection[missing_idx], "_",
      substr(gsub("[^A-Za-z]", "", menu_df$gs_subcollection[missing_idx]), 1, 4)
    )
    menu_df$description[missing_idx] <- paste0("新集合: ", menu_df$gs_subcollection[missing_idx])
  }

  menu_df$combo_name <- paste0(menu_df$gs_collection,
                               ifelse(menu_df$gs_subcollection_clean == "", "", paste0(":", menu_df$gs_subcollection_clean)))

  # 交互界面逻辑
  if (is.null(auto_select)) {
    message("\n", rep("=", 60))
    message(sprintf("🌟 欢迎使用 GSEAlens 基因集向导 (%s)", species))
    message(rep("=", 60))
    for (i in 1:nrow(menu_df)) {
      cat(sprintf("[%2d] %-6s | %-16s | %s\n",
                  i, menu_df$short_tag[i], menu_df$combo_name[i], menu_df$description[i]))
    }
    user_input <- readline(prompt = "\n👉 请输入编号 (用逗号分隔, 如 17,26): ")
    selected_idx <- as.integer(unlist(strsplit(user_input, "[, ]+")))
    selected_idx <- selected_idx[!is.na(selected_idx) & selected_idx >= 1 & selected_idx <= nrow(menu_df)]
    if (length(selected_idx) == 0) stop("❌ 无效的输入！")
  } else {
    if (length(auto_select) == 1 && toupper(auto_select) == "ALL") {
      message("🚨 检测到 [ALL] 指令：正在载入 MSigDB 全库...")
      selected_idx <- 1:nrow(menu_df)
    } else if (is.character(auto_select)) {
      selected_idx <- match(auto_select, menu_df$combo_name)
      if (any(is.na(selected_idx))) stop("❌ 找不到对应的集合名称，请检查拼写。")
    } else {
      selected_idx <- auto_select
    }
  }

  selected_rows <- menu_df[selected_idx, ]

  # 智能命名
  if (length(auto_select) == 1 && toupper(auto_select) == "ALL") {
    super_tag <- "ALL_MSigDB_Global"
  } else {
    super_tag <- ifelse(nrow(selected_rows) <= 4,
                        paste(selected_rows$short_tag, collapse = "_"),
                        sprintf("Mix%d_%s", nrow(selected_rows), selected_rows$short_tag[1]))
  }

  message(sprintf("\n✅ 已选定 %d 个集合。智能批次 Tag: [%s]", nrow(selected_rows), super_tag))

  # 核心提取：使用 collection 参数
  pathway_list <- lapply(1:nrow(selected_rows), function(i) {
    c_coll <- selected_rows$gs_collection[i]
    c_sub <- selected_rows$gs_subcollection[i]
    if (grepl("NO_SUB", c_sub) || c_sub == "") {
      msigdbr::msigdbr(species = species, collection = c_coll)
    } else {
      msigdbr::msigdbr(species = species, collection = c_coll, subcollection = c_sub)
    }
  })

  all_pathways <- dplyr::bind_rows(pathway_list)

  # 动态检测列名
  cat_col_name <- if("gs_cat" %in% colnames(all_pathways)) "gs_cat" else "gs_collection"
  subcat_col_name <- if("gs_subcat" %in% colnames(all_pathways)) "gs_subcat" else "gs_subcollection"

  TERM2GENE <- all_pathways %>% dplyr::select(gs_name, gene_symbol)

  TERM2NAME <- all_pathways %>%
    dplyr::select(ID = gs_name,
                  Description = gs_description,
                  URL = gs_url,
                  Collection = dplyr::all_of(cat_col_name),
                  Subcollection = dplyr::all_of(subcat_col_name)) %>%
    dplyr::distinct(ID, .keep_all = TRUE) %>%
    dplyr::mutate(
      Combo_Name = ifelse(Subcollection == "",
                          Collection,
                          paste0(Collection, ":", Subcollection))
    )

  return(list(TERM2GENE = TERM2GENE, meta_dict = TERM2NAME, SuperTag = super_tag, collections_used = selected_rows))
}
