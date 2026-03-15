#' @title 并行计算 GSEA 核心引擎 (Pro 智能缓存版 - 带血统记忆)
#' @description 自动化执行 limma 结果清洗、去重，并下发多核 GSEA 计算，最终生成带有项目路径记忆的胶囊。
#' @param gsea_env 标准 GseaEnv 或 GseaEnvPro 胶囊对象
#' @param custom_series_name 字符串。分析的系列名称（将作为文件夹名称和血统名称）
#' @param output_dir 字符串。输出的基础路径，默认 "./GSEA_Output"
#' @param force 逻辑值。是否强制重新计算而不使用缓存，默认 FALSE
#' @param bidirectional 逻辑值。是否进行双向对比（自动生成正反双向），默认 TRUE
#' @param workers 并行核心数。若为 NULL，自动保留 4 核，其余全用。
#' @param minGSSize 基因集最小包含基因数，默认 10
#' @param maxGSSize 基因集最大包含基因数，默认 500
#' @param pvalueCutoff GSEA 的 P 值阈值（强制默认 1，保留全貌用于后续精细过滤）
#' @return 返回 GseaResPro 对象
#' @export
batch_calc_gsea_pro <- function(gsea_env, custom_series_name = "Auto_Analysis", output_dir = "./GSEA_Output",
                                force = FALSE, bidirectional = TRUE, workers = 20,
                                minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1) {

  # 兼容您新旧版本的对象命名
  if (!inherits(gsea_env, c("GseaEnvPro", "GseaEnv"))) {
    stop("❌ 请传入标准 GseaEnv 或 GseaEnvPro 胶囊对象！")
  }

  super_tag <- gsea_env$geneset$name

  # 🔴 1. 智能缓存拦截系统
  series_dir <- file.path(output_dir, custom_series_name)
  if (!dir.exists(series_dir)) dir.create(series_dir, recursive = TRUE, showWarnings = FALSE)

  rds_name <- sprintf("GSEA_Capsule_[%s]_[%s].rds", custom_series_name, super_tag)
  rds_path <- file.path(series_dir, rds_name)

  if (file.exists(rds_path) && force == FALSE) {
    message(sprintf("✅ 命中缓存！检测到已存在的 GSEA 胶囊: %s", rds_name))
    message("   触发极速载入模式，直接跳过计算环节...")
    return(readRDS(rds_path))
  }

  # 🟡 2. 硬件调度与数据清洗系统
  total_cores <- parallel::detectCores(logical = TRUE)
  use_cores <- if (is.null(workers)) max(1, total_cores - 4) else min(total_cores, max(1, workers))
  message(sprintf("🖥️ 硬件侦测: 发现 %d 个逻辑核心。调度 %d 核执行并行计算！", total_cores, use_cores))

  # 修复了原始代码中遗漏的 future 显存控制，防止大数据崩溃
  options(future.globals.maxSize = 32000 * 1024^2)
  future::plan(future::multisession, workers = use_cores)

  fit <- gsea_env$fit
  contrasts <- gsea_env$contrasts
  tasks <- list()

  message("🔍 正在提取 topTable 并进行去重与大小写清洗...")
  for (i in 1:nrow(contrasts)) {
    c_name <- contrasts$Contrast_Name[i]
    num <- contrasts$Num[i]; den <- contrasts$Den[i]

    tt <- limma::topTable(fit, coef = c_name, number = Inf) %>% as.data.frame()
    if (!"SYMBOL" %in% colnames(tt)) tt$SYMBOL <- rownames(tt)

    # 完美的内置清洗逻辑保留
    genelist_base <- tt %>%
      dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
      dplyr::mutate(SYMBOL = toupper(SYMBOL)) %>%
      dplyr::arrange(dplyr::desc(abs(t))) %>%
      dplyr::distinct(SYMBOL, .keep_all = TRUE) %>%
      dplyr::select(SYMBOL, t) %>%
      tibble::deframe()

    tasks[[paste0(num, "_vs_", den)]] <- sort(genelist_base, decreasing = TRUE)
    if (bidirectional) {
      tasks[[paste0(den, "_vs_", num)]] <- sort(-genelist_base, decreasing = TRUE)
    }
  }

  message(sprintf("🚀 洗练完毕，生成 %d 个独立对比 genelist。发射多核置换计算...", length(tasks)))

  term2gene <- gsea_env$geneset$term2gene
  res_list <- future.apply::future_lapply(names(tasks), function(task_name) {
    genelist <- tasks[[task_name]]
    set.seed(123)
    gsea_res <- tryCatch({
      clusterProfiler::GSEA(
        geneList = genelist, TERM2GENE = term2gene, minGSSize = minGSSize,
        maxGSSize = maxGSSize, pvalueCutoff = pvalueCutoff, pAdjustMethod = "BH", verbose = FALSE, seed = 123,
        eps = 0 #为了屏蔽报错
      )
    }, error = function(e) NULL)

    status <- if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) "Success" else "Failed/NoEnrich"
    return(list(name = task_name, status = status, data = gsea_res, genelist = genelist))
  }, future.seed = TRUE)

  names(res_list) <- names(tasks)
  future::plan(future::sequential)

  # 🌟 终极打包：注入 Project Info 血统记忆
  final_obj <- list(
    metadata = list(
      run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      workers_used = use_cores,
      parameters = list(bidirectional = bidirectional, minGSSize = minGSSize),
      project_info = list(                      # <--- 新增记忆中枢
        custom_series_name = custom_series_name,
        output_dir = normalizePath(output_dir, mustWork = FALSE),
        series_dir = normalizePath(series_dir, mustWork = FALSE),
        rds_path = normalizePath(rds_path, mustWork = FALSE)
      )
    ),
    geneset_info = gsea_env$geneset,
    results = res_list,
    expr_data = gsea_env$expr_data
  )
  class(final_obj) <- "GseaResPro"

  # 🟢 3. 资产落袋保护
  saveRDS(final_obj, rds_path)

  example_task <- names(tasks)[1]
  cat(sprintf("\n=============================================================\n"))
  cat(sprintf("✅ [数据资产已存档] 计算胶囊已安全保存，防崩溃保护已启动！\n"))
  cat(sprintf("📂 文件路径: %s\n\n", rds_path))
  cat(sprintf("💡 未来如需写文章画单图，请复制以下 R 代码直接调取：\n"))
  cat(sprintf("-------------------------------------------------------------\n"))
  cat(sprintf("my_capsule <- GSEAlens::import_gsea_capsule(\"%s\")\n", rds_path)) # 推荐使用新版智能载入
  cat(sprintf("my_task <- GSEAlens::extract_gsea_task_pro(my_capsule, \"%s\")\n", example_task))
  cat(sprintf("GSEAlens::plot_directional_gsea(my_task, target_pathways = c(\"ID_1\"))\n"))
  cat(sprintf("=============================================================\n\n"))

  return(final_obj)
}





#' @title 展示 GSEA 胶囊中包含的所有子集 (科学排序)
#' @description 提取并按 H -> C1 -> C2...-> C8 的生物学逻辑严格排序，全量展示子集及其通路数量。
#' @param meta_dict 胶囊中的 meta_dict 数据框
#' @return 隐式返回排序后的统计数据框
#' @export
show_gsea_subsets <- function(meta_dict) {
  if(is.null(meta_dict) || !("Combo_Name" %in% colnames(meta_dict))) {
    stop("❌ 字典格式错误或缺失 Combo_Name 列！")
  }

  # 1. 统计数量
  combo_counts <- table(meta_dict$Combo_Name)
  combo_df <- as.data.frame(combo_counts, stringsAsFactors = FALSE)
  colnames(combo_df) <- c("Combo_Name", "Count")

  # 2. 提取主集前缀 (如 "C2:CP:KEGG" 提取出 "C2")
  combo_df$Base <- sapply(strsplit(combo_df$Combo_Name, ":"), `[`, 1)

  # 3. 设定科学排序因子 (H 第一，C1-C8 随后，其他垫底)
  base_levels <- c("H", paste0("C", 1:8))
  other_bases <- setdiff(unique(combo_df$Base), base_levels)
  combo_df$Base_Factor <- factor(combo_df$Base, levels = c(base_levels, sort(other_bases)))

  # 4. 终极排序：先按主集，再按子集字母表排序
  combo_df <- combo_df[order(combo_df$Base_Factor, combo_df$Combo_Name), ]

  # 5. 全量优雅打印
  message(sprintf("📦 本胶囊包含 %d 个可用子集 (共 %d 条通路，已按主集 H->C8 科学排序):",
                  nrow(combo_df), sum(combo_df$Count)))

  for (i in seq_len(nrow(combo_df))) {
    # 使用 %-25s 保证冒号对齐，极其舒适
    message(sprintf("   - %-25s : %5d 条", combo_df$Combo_Name[i], combo_df$Count[i]))
  }

  invisible(combo_df)
}




#' @title 检查并概览 GSEA 计算胶囊 (Pro 引擎 - 终极版)
#' @description 快速查看胶囊计算状态，调用探针展示全量排序子集，并提供主/子集混合提取的保姆级代码。
#' @param res_capsule batch_calc_gsea_pro 返回的计算胶囊
#' @export
inspect_gsea_res_pro <- function(res_capsule) {
  if (!is.list(res_capsule) || is.null(res_capsule$results)) {
    stop("❌ 无效的计算胶囊！请确保传入的是 batch_calc_gsea_pro() 的结果。")
  }

  message("\n", rep("=", 70))
  message("🌟 GSEA 计算胶囊概览 (GSEAlens PRO)")
  message(rep("=", 70))

  meta_dict <- res_capsule$geneset_info$meta_dict
  tag <- res_capsule$geneset_info$name
  message(sprintf("🏷️  全局 Tag: [%s]", tag))
  message(sprintf("🧬  总通路数: %d 条", nrow(meta_dict)))

  message("\n📊 对比组计算状态:")
  for (c_name in names(res_capsule$results)) {
    status <- res_capsule$results[[c_name]]$status
    if (status == "Success") {
      n_sig <- sum(res_capsule$results[[c_name]]$data@result$p.adjust < 0.05, na.rm = TRUE)
      message(sprintf("   ✅ [%s] -> 成功 (未切片全局 FDR < 0.05: %d 条)", c_name, n_sig))
    } else {
      message(sprintf("   ❌ [%s] -> 失败或无富集", c_name))
    }
  }

  message("\n" , rep("-", 70))
  message("💡 下游切片提取指南 (extract_gsea_task_pro):")
  message("您可以基于以下 [子集标签] 进行精准提取，引擎将自动重算并拯救真实的 FDR (p.adjust)。\n")

  # 🌟 直接调用刚才写的探针函数，全量、科学排序展示！
  show_gsea_subsets(meta_dict)

  message("\n👨‍💻 多重组合提取代码示例 (支持主集、子集、或任意自由组合):")
  example_contrast <- names(res_capsule$results)[1]

  message("\n  # 🎯 场景1：单独提取某个特定的子集 (例如只看 KEGG)")
  message(sprintf('  res_kegg <- extract_gsea_task_pro('))
  message(sprintf('    gsea_capsule = my_res_capsule, task_name = "%s",', example_contrast))
  message(sprintf('    target_collection = "C2:CP:KEGG_LEGACY"'))
  message(sprintf('  )'))

  message("\n  # 🎯 场景2：提取多个主集的全部内容 (例如 C2 和 C5 全要)")
  message(sprintf('  res_c2_c5 <- extract_gsea_task_pro('))
  message(sprintf('    gsea_capsule = my_res_capsule, task_name = "%s",', example_contrast))
  message(sprintf('    target_collection = c("C2", "C5") # <--- 直接传主集名即可！'))
  message(sprintf('  )'))

  message("\n  # 🎯 场景3：主集 + 子集 混合点杀 (例如 Hallmark 全要 + GOBP)")
  message(sprintf('  res_mix <- extract_gsea_task_pro('))
  message(sprintf('    gsea_capsule = my_res_capsule, task_name = "%s",', example_contrast))
  message(sprintf('    target_collection = c("H", "C5:GO:BP") # <--- 底层将合并并重新校正 FDR！'))
  message(sprintf('  )'))

  message(rep("=", 70), "\n")
  invisible(res_capsule)
}









#' @title Pro 数据提取器 (文章精细绘图与按需切片专用)
#' @description 提取对比组计算结果。如果使用 ALL 库进行计算，支持传入 target_collection 进行动态切片并自动重算 FDR。
#' @param gsea_capsule GseaResPro 计算胶囊
#' @param task_name 想提取的组别名 (如 "Treat_vs_Control")
#' @param target_collection 想要切片的亚组名(如 "H" 或 "C2:CP:KEGG_LEGACY")。默认为 "ALL" (不切片)。
#' @return DirectionalGSEA 对象，完美兼容底层画图，并将矩阵等所有附加信息收拢于 meta。
#' @export
extract_gsea_task_pro <- function(gsea_capsule, task_name, target_collection = "ALL") {

  if (!inherits(gsea_capsule, "GseaResPro")) stop("❌ 请传入标准 GseaResPro 胶囊对象！")
  if (!(task_name %in% names(gsea_capsule$results))) stop(sprintf("❌ 任务 '%s' 不存在于胶囊中！", task_name))

  task <- gsea_capsule$results[[task_name]]
  if (task$status != "Success") stop("❌ 该组无成功计算结果！")

  # 1. 解析对比组名称 (兼容处理兜底)
  parts <- strsplit(task_name, "_vs_")[[1]]
  left_grp <- if(length(parts) >= 1) parts[1] else task_name
  right_grp <- if(length(parts) >= 2) parts[2] else "Background"

  # 2. 拿到原始计算对象和全局字典
  gsea_res <- task$data
  master_meta <- gsea_capsule$geneset_info$meta_dict
  res_df <- as.data.frame(gsea_res@result)

  # 3. 🛡️ 血统贴回：将字典中的 Collection / Combo_Name 贴到结果表中
  matched_idx <- match(res_df$ID, master_meta$ID)

  # 防呆提取列
  desc_col <- if("Description" %in% colnames(master_meta)) master_meta$Description else master_meta$ID
  url_col <- if("URL" %in% colnames(master_meta)) master_meta$URL else NA
  coll_col <- if("Collection" %in% colnames(master_meta)) master_meta$Collection else "Unknown"
  combo_col <- if("Combo_Name" %in% colnames(master_meta)) master_meta$Combo_Name else coll_col

  res_df$Collection <- coll_col[matched_idx]
  res_df$Combo_Name <- combo_col[matched_idx]


  # 🌟 4. 核心：动态切片与 FDR (p.adjust) 重算！

  is_slice_mode <- length(target_collection) != 1 || toupper(target_collection[1]) != "ALL"

  if (is_slice_mode) {
    # 允许按 Collection (如 "H") 或 Combo_Name (如 "C2:CP:KEGG_LEGACY") 切片
    slice_idx <- res_df$Collection %in% target_collection | res_df$Combo_Name %in% target_collection
    res_df <- res_df[slice_idx, , drop = FALSE]

    if (nrow(res_df) == 0) {
      stop(sprintf("❌ 切片失败！组别 [%s] 中没有发现属于 '%s' 的通路。",
                   task_name, paste(target_collection, collapse = ", ")))
    }

    # 💥 消除全库背景的惩罚，针对小子集重新校正多重假设检验 💥
    res_df$p.adjust <- p.adjust(res_df$pvalue, method = "BH")
    # 如果底层引擎生成了 qvalue，为了严谨性同步校正为 p.adjust (或设为 NA 避免误导)
    if ("qvalue" %in% colnames(res_df)) res_df$qvalue <- res_df$p.adjust

    message(sprintf("✂️ [智能切片] 成功抽出 %d 条 '%s' 通路 | 🌟 已重新校正 FDR (p.adjust)！",
                    nrow(res_df), paste(target_collection, collapse = ", ")))
  } else {
    message(sprintf("📦 [全库提取] 提取出 %d 条通路 (未进行切片与 FDR 重算)", nrow(res_df)))
  }

  # 把修改完的 dataframe 塞回 S4 对象中
  gsea_res@result <- res_df

  # 5. 重构轻量级字典用于前端交互
  meta_dict_for_res_obj <- data.frame(
    ID = res_df$ID,
    long_description_for_html = desc_col[match(res_df$ID, master_meta$ID)],
    URL = url_col[match(res_df$ID, master_meta$ID)],
    Collection = res_df$Collection,
    Combo_Name = res_df$Combo_Name,
    stringsAsFactors = FALSE
  )

  # 🌟 6. 极简打包原则：万物皆入 meta！
  res_obj <- list(
    gsea_res = gsea_res,
    meta = list(
      left_group = left_grp,
      right_group = right_grp,
      geneset_name = ifelse(is_slice_mode,
                            paste(target_collection, collapse = "_"),
                            gsea_capsule$geneset_info$name),
      meta_dict = meta_dict_for_res_obj,
      expr_data = gsea_capsule$expr_data,               # <--- 将表达矩阵纳入 meta
      project_info = gsea_capsule$metadata$project_info # <--- 将项目血统纳入 meta
    )
  )
  class(res_obj) <- "DirectionalGSEA"

  return(res_obj)
}
