#' @title 构建 GSEA 基因集与超级标签
#' @description 基于 msigdbr 提供交互式菜单，选择基因集并生成智能语义化标签。
#'   支持未知未来子集的动态截断，完全兼容 MSigDB V10+ 及 C9 集合。
#' @param species 字符型，物种名称，默认 "Homo sapiens"
#' @param auto_select 可选参数。可传入：
#'   - 序号向量如 c(17, 26)
#'   - 名称向量如 c("C5:GO:BP", "H")
#'   - 字符串 "ALL" 获取全库
#' @return 返回包含以下字段的列表：
#'   \item{TERM2GENE}{两列数据框：gs_name, gene_symbol}
#'   \item{meta_dict}{通路元数据：ID, Description, URL, Collection, Subcollection, Combo_Name}
#'   \item{SuperTag}{智能生成的批次标签}
#'   \item{collections_used}{选用的原始集合信息}
#' @export
#' @importFrom msigdbr msigdbr msigdbr_collections
#' @importFrom dplyr arrange left_join bind_rows select distinct mutate filter
build_gsea_pathways <- function(species = "Homo sapiens", auto_select = NULL) {

  suppressPackageStartupMessages(requireNamespace("msigdbr"))
  suppressPackageStartupMessages(requireNamespace("dplyr"))

  # 动态拉取可用集合
  avail_colls <- msigdbr::msigdbr_collections() %>%
    dplyr::arrange(gs_collection, gs_subcollection)

  # 白名单字典（含 C9 集合）
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
      "C1 位置基因集", "C2 化学和遗传微扰", "C2 经典通路(总)", "C2 BioCarta",
      "C2 KEGG(旧)", "C2 KEGG(医学)", "C2 PID 通路", "C2 Reactome", "C2 WikiPathways",
      "C3 miRNA 靶标", "C3 miRNA(旧)", "C3 TF 靶标(GTRD)", "C3 TF(旧)",
      "C4 3CA 癌症", "C4 癌症邻居", "C4 癌症模块",
      "C5 GO 生物过程(BP)", "C5 GO 细胞组分(CC)", "C5 GO 分子功能(MF)", "C5 人类表型(HPO)",
      "C6 肿瘤特征", "C7 免疫特征", "C7 疫苗特征", "C8 细胞类型",
      "C9 计算扰动特征(Perturbation)", "H 癌症标志(Hallmark)"
    ),
    stringsAsFactors = FALSE
  )

  # 处理子集名称匹配（去除 NO_SUB 后缀）
  avail_colls$gs_subcollection_clean <- gsub(".*:NO_SUB", "", avail_colls$gs_subcollection)
  menu_df <- dplyr::left_join(
    avail_colls, dict,
    by = c("gs_collection", "gs_subcollection_clean" = "gs_subcollection")
  )

  # 动态截断机制：为未知未来集合生成临时标签
  missing_idx <- is.na(menu_df$short_tag)
  if (any(missing_idx)) {
    menu_df$short_tag[missing_idx] <- paste0(
      menu_df$gs_collection[missing_idx], "_",
      substr(gsub("[^A-Za-z]", "", menu_df$gs_subcollection[missing_idx]), 1, 4)
    )
    menu_df$description[missing_idx] <- paste0("新集合: ", menu_df$gs_subcollection[missing_idx])
  }

  menu_df$combo_name <- paste0(
    menu_df$gs_collection,
    ifelse(menu_df$gs_subcollection_clean == "", "",
           paste0(":", menu_df$gs_subcollection_clean))
  )

  # 交互式选择逻辑
  if (is.null(auto_select)) {
    message("\n", rep("=", 60))
    message(sprintf("欢迎使用 GSEAlens 基因集向导 (%s)", species))
    message(rep("=", 60))
    for (i in seq_len(nrow(menu_df))) {
      cat(sprintf("[%2d] %-6s | %-16s | %s\n",
                  i, menu_df$short_tag[i], menu_df$combo_name[i], menu_df$description[i]))
    }
    user_input <- readline(prompt = "\n请输入编号 (用逗号分隔, 如 17,26): ")
    selected_idx <- as.integer(unlist(strsplit(user_input, "[, ]+")))
    selected_idx <- selected_idx[!is.na(selected_idx) & selected_idx >= 1 &
                                   selected_idx <= nrow(menu_df)]
    if (length(selected_idx) == 0) stop("无效的输入！")
  } else {
    if (length(auto_select) == 1 && toupper(auto_select) == "ALL") {
      message("检测到 [ALL] 指令：正在载入 MSigDB 全库...")
      selected_idx <- seq_len(nrow(menu_df))
    } else if (is.character(auto_select)) {
      selected_idx <- match(auto_select, menu_df$combo_name)
      if (any(is.na(selected_idx))) stop("找不到对应的集合名称，请检查拼写。")
    } else {
      selected_idx <- auto_select
    }
  }

  selected_rows <- menu_df[selected_idx, , drop = FALSE]

  # 智能命名策略
  if (length(auto_select) == 1 && toupper(auto_select) == "ALL") {
    super_tag <- "ALL_MSigDB_Global"
  } else {
    super_tag <- ifelse(
      nrow(selected_rows) <= 4,
      paste(selected_rows$short_tag, collapse = "_"),
      sprintf("Mix%d_%s", nrow(selected_rows), selected_rows$short_tag[1])
    )
  }

  message(sprintf("\n已选定 %d 个集合。智能批次 Tag: [%s]",
                  nrow(selected_rows), super_tag))

  # 核心：使用 collection 参数（兼容 V10）
  pathway_list <- lapply(seq_len(nrow(selected_rows)), function(i) {
    c_coll <- selected_rows$gs_collection[i]
    c_sub <- selected_rows$gs_subcollection[i]
    if (grepl("NO_SUB", c_sub) || c_sub == "") {
      msigdbr::msigdbr(species = species, collection = c_coll)
    } else {
      msigdbr::msigdbr(species = species, collection = c_coll, subcollection = c_sub)
    }
  })

  all_pathways <- dplyr::bind_rows(pathway_list)

  # 动态列名检测（兼容不同 msigdbr 版本）
  cat_col_name <- if ("gs_cat" %in% colnames(all_pathways)) "gs_cat" else "gs_collection"
  subcat_col_name <- if ("gs_subcat" %in% colnames(all_pathways)) "gs_subcat" else "gs_subcollection"

  TERM2GENE <- all_pathways %>%
    dplyr::select(gs_name, gene_symbol)

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
      Combo_Name = ifelse(Subcollection == "", Collection,
                          paste0(Collection, ":", Subcollection))
    )

  list(
    TERM2GENE = TERM2GENE,
    meta_dict = TERM2NAME,
    SuperTag = super_tag,
    collections_used = selected_rows
  )
}


#' @title 组装 GSEA 计算环境（多态入口）
#' @description 统一入口函数，自动识别输入对象类型（MArrayLM 或 DESeqDataSet），
#'   抽提标准化信息后封装为统一的 GseaEnv 对象。这是重构后架构的唯一入口。
#' @param x 输入对象，必须是 MArrayLM（limma-voom）或 DESeqDataSet（DESeq2）
#' @param pathway_obj 基因集对象，由 build_gsea_pathways() 生成
#' @param expr_data 可选的表达数据对象。若为 NULL，尝试从 x 中提取（DGEList 或 DESeqDataSet）
#' @param target_factor 字符型，目标比较因子名称。DESeq2 多因素设计时必须提供；
#'   limma 中可选，用于元数据记录。默认为 NULL
#' @param ... 额外参数传递给具体 backend 的提取函数
#' @return 返回标准化的 GseaEnv 对象（S3 类），包含：
#'   \item{backend_info}{backend 类型、设计公式等元数据}
#'   \item{contrast_registry}{统一对比注册表}
#'   \item{de_store}{标准化差异表达结果列表}
#'   \item{expr_bundle}{标准化表达数据包}
#'   \item{geneset}{基因集信息}
#'   \item{raw_backend_obj}{原始对象存档（用于 debug）}
#' @export
#' @seealso build_gsea_pathways, validate_gsea_env
setup_gsea_env <- function(x, pathway_obj, expr_data = NULL, target_factor = NULL, ...) {

  # 严格类型检查
  if (!inherits(x, c("MArrayLM", "DESeqDataSet"))) {
    stop("x 必须是 MArrayLM（limma）或 DESeqDataSet（DESeq2）对象")
  }

  # 基因集对象验证
  if (is.null(pathway_obj$TERM2GENE) || is.null(pathway_obj$meta_dict)) {
    stop("pathway_obj 结构缺失！请确保使用 build_gsea_pathways() 生成")
  }

  # S3 方法分派
  UseMethod("setup_gsea_env", x)
}


#' @title 组装 GSEA 计算环境（limma-voom 专用方法）
#' @description MArrayLM 对象的专用提取逻辑，生成标准化 GseaEnv。
#'   自动从 fit$coefficients 提取对比矩阵，构建 contrast_registry。
#' @param x MArrayLM 对象，必须包含 contrast（colnames 存在）
#' @param pathway_obj 基因集对象
#' @param expr_data 可选的 DGEList 对象。若为 NULL，尝试从 x$genes 等位置推断
#' @param target_factor 可选的目标因子名称（仅作记录）
#' @param ... 忽略
#' @return GseaEnv 对象
#' @keywords internal
#' @export
setup_gsea_env.MArrayLM <- function(x, pathway_obj, expr_data = NULL,
                                    target_factor = NULL, ...) {

  message("检测到 limma-voom backend，开始标准化提取...")

  fit_object <- x

  # 验证 contrast 存在
  coef_names <- colnames(fit_object)
  if (is.null(coef_names)) {
    stop("MArrayLM 对象缺少 colnames（对比组名），请检查 limma 设计矩阵")
  }

  # 推断 target_factor（若未提供）
  inferred_target_factor <- target_factor
  if (is.null(inferred_target_factor)) {
    # 尝试从 design 矩阵列名推断（简单启发式：找包含对比组名的因子）
    # limma 中通常无法自动推断，保持 NULL
    inferred_target_factor <- NULL
  }

  # 构建 contrast_registry
  contrast_reg <- create_contrast_registry_limma(
    coef_names = coef_names,
    target_factor = inferred_target_factor
  )

  # 构建 de_store（预先提取所有 contrast 的 topTable）
  de_store <- create_de_store_limma(
    fit = fit_object,
    contrast_registry = contrast_reg
  )

  # 构建 expr_bundle
  expr_bundle <- create_expr_bundle(
    expr_data = expr_data,
    backend = "limma_voom",
    target_factor = inferred_target_factor
  )

  # backend 元信息
  backend_info <- list(
    backend = "limma_voom",
    input_class = "MArrayLM",
    design_formula = NA_character_,  # limma 对象中难以反向解析公式
    target_factor = inferred_target_factor,
    supported_mode = "pairwise_factor",
    n_contrasts = nrow(contrast_reg),
    n_genes = nrow(fit_object$coefficients)
  )

  # 组装终极胶囊
  env_obj <- list(
    backend_info = backend_info,
    contrast_registry = contrast_reg,
    de_store = de_store,
    expr_bundle = expr_bundle,
    geneset = list(
      name = pathway_obj$SuperTag,
      term2gene = pathway_obj$TERM2GENE,
      meta_dict = pathway_obj$meta_dict,
      used_collections = pathway_obj$collections_used
    ),
    raw_backend_obj = list(fit = fit_object)  # 仅存档，下游不应直接访问
  )

  class(env_obj) <- c("GseaEnv", "list")

  # 验证与反馈
  validate_gsea_env(env_obj)

  message("\n", rep("-", 60))
  message("胶囊 [GseaEnv] 封装完毕（limma-voom backend）")
  message(sprintf("   基因集 Tag : [%s]", pathway_obj$SuperTag))
  message(sprintf("   通路总数量 : %d 条", nrow(pathway_obj$meta_dict)))
  message(sprintf("   对比组数量 : %d 个", nrow(contrast_reg)))
  if (pathway_obj$SuperTag == "ALL_MSigDB_Global") {
    message("\n[注意] 已装载 MSigDB 全库，后续计算建议设置 pvalueCutoff = 1")
  }
  message(rep("-", 60), "\n")

  return(env_obj)
}
