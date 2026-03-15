#' @title 交互式构建 GSEA 基因集与超级标签 (Pro 引擎 - 兼容 V10 与 C9 修复版)
#' @description 基于 msigdbr 提供交互式菜单，选择基因集并生成智能语义化标签。支持未知未来子集的动态截断。
#' @param species 物种，默认 "Homo sapiens"
#' @param auto_select 可选。可传入序号向量 c(17, 26)；或者更安全的名称向量 c("C5:GO:BP", "H")
#' @return 返回包含 TERM2GENE、TERM2NAME、SuperTag 的列表对象
#' @export
build_gsea_pathways_pro <- function(species = "Homo sapiens", auto_select = NULL) {

  suppressPackageStartupMessages(require(msigdbr))
  suppressPackageStartupMessages(require(dplyr))

  # 动态拉取库
  avail_colls <- msigdbr::msigdbr_collections() %>%
    dplyr::arrange(gs_collection, gs_subcollection)

  # 🌟 白名单字典 (已正式收编 C9 集合！)
  dict <- data.frame(
    gs_collection = c(
      "C1", "C2", "C2", "C2", "C2", "C2", "C2", "C2", "C2",
      "C3", "C3", "C3", "C3", "C4", "C4", "C4",
      "C5", "C5", "C5", "C5", "C6", "C7", "C7", "C8", "C9", "H" # 加入了 C9
    ),
    gs_subcollection = c(
      "", "CGP", "CP", "CP:BIOCARTA", "CP:KEGG_LEGACY", "CP:KEGG_MEDICUS",
      "CP:PID", "CP:REACTOME", "CP:WIKIPATHWAYS", "MIR:MIRDB", "MIR:MIR_LEGACY",
      "TFT:GTRD", "TFT:TFT_LEGACY", "3CA", "CGN", "CM",
      "GO:BP", "GO:CC", "GO:MF", "HPO", "", "IMMUNESIGDB", "VAX", "", "", "" # C9 无子库，留空
    ),
    short_tag = c(
      "C1Pos", "CGP", "CP", "BioC", "KeggL", "KeggM",
      "PID", "Reac", "Wiki", "MirDB", "MirL",
      "TFTG", "TFTL", "3CA", "CGN", "CM",
      "GoBP", "GoCC", "GoMF", "HPO", "C6Onc", "ImmS", "VAX", "C8Cell", "C9Pert", "Hal" # C9 的短标签
    ),
    description = c(
      "C1 位置基因集", "C2 化学和遗传微扰", "C2 经典通路(总)", "C2 BioCarta", "C2 KEGG(旧)", "C2 KEGG(医学)",
      "C2 PID 通路", "C2 Reactome", "C2 WikiPathways", "C3 miRNA 靶标", "C3 miRNA(旧)",
      "C3 TF 靶标(GTRD)", "C3 TF(旧)", "C4 3CA 癌症", "C4 癌症邻居", "C4 癌症模块",
      "C5 GO 生物过程(BP)", "C5 GO 细胞组分(CC)", "C5 GO 分子功能(MF)", "C5 人类表型(HPO)", "C6 肿瘤特征", "C7 免疫特征", "C7 疫苗特征", "C8 细胞类型", "C9 计算扰动特征(Perturbation)", "H 癌症标志(Hallmark)"
    ),
    stringsAsFactors = FALSE
  )

  avail_colls$gs_subcollection_clean <- gsub(".*:NO_SUB", "", avail_colls$gs_subcollection)
  menu_df <- dplyr::left_join(avail_colls, dict,
                              by = c("gs_collection", "gs_subcollection_clean" = "gs_subcollection"))

  # 🔮 动态截断机制 (应对未来更多未知集合)
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

  # 交互界面
  if (is.null(auto_select)) {
    message("\n", rep("=", 60))
    message(sprintf("🌟 欢迎使用 GSEAlens PRO 基因集向导 (%s)", species))
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
    if (is.character(auto_select)) {
      selected_idx <- match(auto_select, menu_df$combo_name)
      if (any(is.na(selected_idx))) stop(sprintf("❌ 找不到集合名: %s", paste(auto_select[is.na(selected_idx)], collapse = ", ")))
    } else {
      selected_idx <- auto_select
    }
  }

  selected_rows <- menu_df[selected_idx, ]
  super_tag <- ifelse(nrow(selected_rows) <= 4,
                      paste(selected_rows$short_tag, collapse = "_"),
                      sprintf("Mix%d_%s", nrow(selected_rows), selected_rows$short_tag[1]))

  message(sprintf("\n✅ 已选定 %d 个集合。智能批次 Tag: [%s]", nrow(selected_rows), super_tag))

  # 🌟 核心修复：彻底废弃 category，使用 collection
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
  TERM2GENE <- all_pathways %>% dplyr::select(gs_name, gene_symbol)
  TERM2NAME <- all_pathways %>% dplyr::select(ID = gs_name, Description = gs_description, URL = gs_url, Collection = gs_cat) %>% dplyr::distinct(ID, .keep_all = TRUE)

  return(list(TERM2GENE = TERM2GENE, meta_dict = TERM2NAME, SuperTag = super_tag, collections_used = selected_rows))
}





#' @title 组装计算胶囊 (Pro 引擎)
#' @description 将差异结果(fit)、基因集与表达矩阵完美焊死在一起，实现一次打包，终身复现。
#' @param fit limma 分析得到的 MArrayLM 对象
#' @param pathway_obj build_gsea_pathways_pro() 返回的基因集对象
#' @param expr_data 你的 DGEList 或者标准化后的表达矩阵 (用于画热图，可为 NULL)
#' @export
setup_gsea_env_pro <- function(fit, pathway_obj, expr_data = NULL) {
  if (!inherits(fit, "MArrayLM")) stop("❌ fit 必须是 limma 的对象！")

  c_names <- colnames(fit)
  parsed <- lapply(c_names, function(x) {
    p <- strsplit(x, "\\s*-\\s*")[[1]]
    if(length(p) == 2) c(p[1], p[2]) else c(x, "Unknown")
  })
  parsed_df <- do.call(rbind, parsed)
  contrasts_df <- tibble::tibble(ID = seq_along(c_names), Contrast_Name = c_names, Num = parsed_df[,1], Den = parsed_df[,2])

  env_obj <- list(
    fit = fit,
    contrasts = contrasts_df,
    geneset = list(
      name = pathway_obj$SuperTag,
      term2gene = pathway_obj$TERM2GENE,
      meta_dict = pathway_obj$meta_dict,
      used_collections = pathway_obj$collections_used
    ),
    expr_data = expr_data
  )

  class(env_obj) <- "GseaEnvPro"
  message(sprintf("✅ 胶囊封装完毕！对比组数: %d | Tag: [%s]", nrow(contrasts_df), pathway_obj$SuperTag))
  return(env_obj)
}

#' @title 检查计算胶囊 (Pro 引擎)
#' @param env_pro setup_gsea_env_pro 创建的对象
#' @export
inspect_gsea_env_pro <- function(env_pro) {
  # (直接复制 txt 中的 inspect_gsea_env_pro)
}
