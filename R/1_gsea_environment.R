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
    if (length(auto_select) == 1 && toupper(auto_select) == "ALL") {
      # 💥 触发全库获取模式！
      message("🚨 检测到 [ALL] 指令：正在载入 MSigDB 全库 (可能包含上万条通路)...")
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

  # 🌟 核心防呆修复：动态检测列名，彻底解决 gs_cat/gs_subcat 消失的报错
  cat_col_name <- if("gs_cat" %in% colnames(all_pathways)) "gs_cat" else "gs_collection"
  subcat_col_name <- if("gs_subcat" %in% colnames(all_pathways)) "gs_subcat" else "gs_subcollection"

  TERM2GENE <- all_pathways %>% dplyr::select(gs_name, gene_symbol)

  TERM2NAME <- all_pathways %>%
    dplyr::select(ID = gs_name,
                  Description = gs_description,
                  URL = gs_url,
                  Collection = dplyr::all_of(cat_col_name),
                  Subcollection = dplyr::all_of(subcat_col_name)) %>%  # <--- 加入了子集！
    dplyr::distinct(ID, .keep_all = TRUE) %>%
    dplyr::mutate(
      # 🌟 新增：拼接成 C2:CP:KEGG_LEGACY 的标准格式，方便后续精准切片！
      Combo_Name = ifelse(Subcollection == "",
                          Collection,
                          paste0(Collection, ":", Subcollection))
    )

  return(list(TERM2GENE = TERM2GENE, meta_dict = TERM2NAME, SuperTag = super_tag, collections_used = selected_rows))
}





#' @title 组装计算胶囊 (Pro 引擎 - 终极版)
#' @description 将差异结果(fit)、基因集字典(带亚组血统)与表达矩阵完美焊死在一起，实现一次打包，终身复现。
#' @param fit limma 分析得到的 MArrayLM 对象 (必须包含 contrast)
#' @param pathway_obj 必须是 build_gsea_pathways_pro() 返回的完整列表对象
#' @param expr_data 你的 DGEList 或者标准化后的表达矩阵 (用于后续动态热图和火山图，可为 NULL)
#' @return 返回一个类为 "GseaEnvPro" 的复合计算胶囊
#' @export
setup_gsea_env_pro <- function(fit, pathway_obj, expr_data = NULL) {

  # 1. 严格拦截与防呆检查
  if (!inherits(fit, "MArrayLM")) {
    stop("❌ 严重错误: fit 必须是 limma 流程中生成的 MArrayLM 对象！")
  }
  if (is.null(pathway_obj$TERM2GENE) || is.null(pathway_obj$meta_dict)) {
    stop("❌ 严重错误: pathway_obj 结构缺失！请确保使用的是 build_gsea_pathways_pro() 生成的对象。")
  }

  # 2. 智能解析对比组 (Contrasts)
  c_names <- colnames(fit)
  if (is.null(c_names)) stop("❌ fit 对象中找不到 colnames(对比组名)，请检查您的 limma 设计矩阵。")

  parsed <- lapply(c_names, function(x) {
    # 兼容 "Treat - Control" 或 "Treat-Control"
    p <- strsplit(x, "\\s*-\\s*")[[1]]
    if (length(p) == 2) {
      return(c(p[1], p[2]))
    } else {
      return(c(x, "Background")) # 如果解析失败，提供安全兜底
    }
  })
  parsed_df <- do.call(rbind, parsed)

  contrasts_df <- tibble::tibble(
    ID = seq_along(c_names),
    Contrast_Name = c_names,
    Num = parsed_df[, 1],
    Den = parsed_df[, 2]
  )

  # 3. 组装终极胶囊 (完美保留亚组血统 meta_dict)
  env_obj <- list(
    fit = fit,
    contrasts = contrasts_df,
    geneset = list(
      name = pathway_obj$SuperTag,
      term2gene = pathway_obj$TERM2GENE,
      meta_dict = pathway_obj$meta_dict,           # 核心：保留了 Collection 供后续按需切片！
      used_collections = pathway_obj$collections_used
    ),
    expr_data = expr_data
  )

  class(env_obj) <- "GseaEnvPro"

  # 4. 动态反馈与预警机制
  message("\n", rep("-", 60))
  message("✅ 胶囊 [GseaEnvPro] 封装完毕！")
  message(sprintf("   🏷️  基因集 Tag : [%s]", pathway_obj$SuperTag))
  message(sprintf("   🧬  通路总数量 : %d 条", nrow(pathway_obj$meta_dict)))
  message(sprintf("   ⚖️  发现对比组 : %d 个", nrow(contrasts_df)))

  # ALL 模式专属提醒
  if (pathway_obj$SuperTag == "ALL_MSigDB_Global") {
    message("\n🚨 [高能预警] 您已装载 MSigDB 全库！")
    message("   下一步扔进 batch_calc_gsea_pro() 时，请务必设置 pvalueCutoff = 1")
    message("   这需要 1~5 分钟的计算时间，请保持耐心，让子弹飞一会儿~")
  }
  message(rep("-", 60), "\n")

  return(env_obj)
}




#' @title 载入并智能归位 GSEA 计算胶囊 (档案管理员)
#' @description 安全载入计算胶囊。若文件脱离原始项目路径（比如被拷贝到了桌面），引擎将自动在当前工作目录复原标准文件夹，并将胶囊护送归位。
#' @param file_path 字符型。胶囊 rds 文件的绝对或相对路径。
#' @param auto_relocate 逻辑值。是否自动根据内置血统在当前目录下修复文件夹并归位？默认 TRUE。
#' @return 返回解析后的 GseaResPro 或 GseaEnv 胶囊对象。
#' @export
import_gsea_capsule <- function(file_path, auto_relocate = TRUE) {
  if (!file.exists(file_path)) stop(sprintf("❌ 文件不存在: %s", file_path))

  message("📦 正在唤醒计算胶囊...")
  capsule <- readRDS(file_path)

  if (inherits(capsule, c("GseaEnvPro", "GseaEnv"))) {
    message("✅ 成功载入 [GseaEnv] 环境封装胶囊 (尚未进行并行计算)。")
    return(invisible(capsule))
  }

  if (!inherits(capsule, "GseaResPro")) {
    warning("⚠️ 该文件似乎不是标准的 GSEAlens 胶囊！")
    return(invisible(capsule))
  }

  info <- capsule$metadata$project_info
  if (is.null(info)) {
    message("⚠️ 该胶囊为旧版本生成，缺乏项目血统记忆。直接载入。")
  } else {
    current_abs <- normalizePath(file_path, winslash = "/", mustWork = FALSE)
    expected_abs <- normalizePath(info$rds_path, winslash = "/", mustWork = FALSE)

    # 路径不一致，且开启了自动归位
    if (current_abs != expected_abs && auto_relocate) {
      message(sprintf("🚨 [血统警报] 胶囊当前处于非标准路径: %s", current_abs))
      message(sprintf("   原籍隶属于项目 : [%s]", info$custom_series_name))

      # 根据当前工作目录(WD)重建标准档案库
      local_series_dir <- file.path(getwd(), "GSEA_Output", info$custom_series_name)
      if (!dir.exists(local_series_dir)) {
        dir.create(local_series_dir, recursive = TRUE)
      }

      target_file <- file.path(local_series_dir, basename(file_path))
      if (!file.exists(target_file) || target_file == current_abs) {
        file.copy(from = file_path, to = target_file, overwrite = TRUE)
        message(sprintf("   🚑 已自动将胶囊遣返归位至标准档案库: %s", local_series_dir))

        # 覆写胶囊里的地址，防止下次原位读取时再次报警告
        capsule$metadata$project_info$output_dir <- file.path(getwd(), "GSEA_Output")
        capsule$metadata$project_info$series_dir <- local_series_dir
        capsule$metadata$project_info$rds_path <- target_file
      } else {
        message("   ✅ 标准档案库中已有该文件备份。")
      }
    }
  }

  # 调用探针函数打印概况 (如果环境中加载了的话)
  if (exists("inspect_gsea_res_pro", mode = "function")) {
    inspect_gsea_res_pro(capsule)
  } else {
    message("✅ 成功载入 [GseaResPro] 计算完成结果胶囊！")
  }

  return(invisible(capsule))
}
