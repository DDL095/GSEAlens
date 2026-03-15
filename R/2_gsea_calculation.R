
# 模块：GSEA Calculation Engine (核心计算与并行引擎)
# 职责：承接环境胶囊，执行高速并行计算、智能缓存拦截、以及独立结果提取


#' @title 并行计算 GSEA 核心引擎 (Pro 智能缓存版)
#' @description 自动感知硬件多核算力，对封装在胶囊中的所有对比组提取 topTable 进行去重并真实双向计算。
#' @param gsea_env 由 setup_gsea_env_pro() 或 setup_gsea_env() 创建的计算胶囊
#' @param custom_series_name 自定义分析系列名称，用于创建专属文件夹和缓存名
#' @param output_dir 输出的总文件夹路径
#' @param force 逻辑值。是否强制无视已存在的缓存 RDS 文件，重新计算？默认 FALSE。
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

  final_obj <- list(
    metadata = list(
      run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      workers_used = use_cores,
      parameters = list(bidirectional = bidirectional, minGSSize = minGSSize)
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
  cat(sprintf("my_capsule <- readRDS(\"%s\")\n", rds_path))
  # 修复了这里的包名前缀
  cat(sprintf("my_task <- GSEAlens::extract_gsea_task_pro(my_capsule, \"%s\")\n", example_task))
  cat(sprintf("GSEAlens::plot_directional_gsea(my_task, target_pathways = c(\"ID_1\"))\n"))
  cat(sprintf("=============================================================\n\n"))

  return(final_obj)
}


#' @title 查看 GSEA 计算胶囊结果
#' @description 格式化打印 GseaResPro 对象的元数据和各个对比组的富集状态，避免直接 View 导致卡顿。
#' @param res_pro batch_calc_gsea_pro() 生成的对象
#' @export
inspect_gsea_res_pro <- function(res_pro) {
  if (!inherits(res_pro, "GseaResPro")) stop("❌ 必须传入 GseaResPro 对象！")

  cat("\n", rep("=", 65), "\n", sep = "")
  cat(sprintf("%-18s %s %-18s\n", "", "🌟 GSEA Pro 运行报告 (Capsule)", ""))
  cat(rep("=", 65), "\n\n", sep = "")

  cat(sprintf("📅 运行时间 : %s\n", res_pro$metadata$run_time))
  cat(sprintf("🖥️ 调度核数 : %s Threads 并行\n", res_pro$metadata$workers_used))
  cat(sprintf("🧬 基因集合 : [%s]\n", res_pro$geneset_info$name)) # 修复了字段名调用

  expr_status <- if (!is.null(res_pro$expr_data)) "✅ 已封装 (支持自动渲染热图)" else "⚠️ 未封装 (仅出折线图)"
  cat(sprintf("🌡️ 表达矩阵 : %s\n", expr_status))
  cat("\n", rep("-", 65), "\n", sep = "")

  cat(sprintf("📊 组别运行状态 (共 %d 组):\n\n", length(res_pro$results)))
  for (i in seq_along(res_pro$results)) {
    task_name <- names(res_pro$results)[i]
    task_info <- res_pro$results[[i]]

    if (task_info$status == "Success") {
      gsea_obj <- task_info$data
      total_count <- nrow(gsea_obj@result)
      sig_count <- sum(gsea_obj@result$p.adjust < 0.05, na.rm = TRUE)
      cat(sprintf("  [%02d] ✅ %-18s | 成功 (总富集: %4d 条, p.adj < 0.05: %3d 条)\n", i, task_name, total_count, sig_count))
    } else {
      cat(sprintf("  [%02d] ❌ %-18s | 失败 (无显著结果或抛出异常)\n", i, task_name))
    }
  }

  cat("\n", rep("=", 65), "\n", sep = "")
}


#' @title Pro 数据提取器 (文章精细绘图专用)
#' @description 直接输入对比组名称，提取并包装成原生 DirectionalGSEA 对象，完美兼容底层画图。
#' @param gsea_capsule GseaResPro 计算胶囊
#' @param task_name 想提取的组别名 (如 "LTA_vs_POSTA")
#' @return DirectionalGSEA 对象。
#' @export
extract_gsea_task_pro <- function(gsea_capsule, task_name) {
  if (!inherits(gsea_capsule, "GseaResPro")) stop("❌ 请传入标准 GseaResPro 胶囊对象！")
  if (!(task_name %in% names(gsea_capsule$results))) stop(sprintf("❌ 任务 '%s' 不存在于胶囊中！", task_name))

  task <- gsea_capsule$results[[task_name]]
  if(task$status != "Success") stop("❌ 该组无成功计算结果！")

  parts <- strsplit(task_name, "_vs_")[[1]]
  left_grp <- parts[1]
  right_grp <- parts[2]

  master_geneset_meta <- gsea_capsule$geneset_info$meta_dict
  current_gsea_ids <- task$data@result$ID
  matched_indices <- match(current_gsea_ids, master_geneset_meta$ID)

  # 防呆：如果有的通路因为意外没有元数据（如自定义 gmt 导致不匹配）
  desc_col <- if("Description" %in% colnames(master_geneset_meta)) master_geneset_meta$Description else master_geneset_meta$ID
  url_col <- if("URL" %in% colnames(master_geneset_meta)) master_geneset_meta$URL else NA
  coll_col <- if("Collection" %in% colnames(master_geneset_meta)) master_geneset_meta$Collection else "Unknown"

  meta_dict_for_res_obj <- data.frame(
    ID = current_gsea_ids,
    long_description_for_html = desc_col[matched_indices],
    URL = url_col[matched_indices],
    Collection = coll_col[matched_indices],
    stringsAsFactors = FALSE
  )

  res_obj <- list(
    gsea_res = task$data,
    meta = list(
      left_group = left_grp,
      right_group = right_grp,
      geneset_name = gsea_capsule$geneset_info$name,
      meta_dict = meta_dict_for_res_obj
    )
  )
  class(res_obj) <- "DirectionalGSEA"

  return(res_obj)
}
