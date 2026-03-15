#' @title 并行计算 GSEA 核心引擎 (Pro - 智能缓存版)
#' @description 结合原生代码最稳健的 t 值去重与大小写转换逻辑，支持硬件核数自定与双向翻转计算。
#' @param gsea_env 由 setup_gsea_env_pro() 创建的计算胶囊
#' @param custom_series_name 自定义分析系列名称
#' @param output_dir 输出目录
#' @param force 是否强制重新计算
#' @param bidirectional 逻辑值，是否进行双向对比
#' @param workers 整数，使用的 CPU 核心数
#' @param minGSSize 最小基因集
#' @param maxGSSize 最大基因集
#' @param pvalueCutoff P值阈值
#' @return DudaliGseaResPro 对象
#' @export
batch_calc_gsea_pro <- function(gsea_env, custom_series_name = "Auto_Analysis", output_dir = "./GSEA_Output",
                                force = FALSE, bidirectional = TRUE, workers = 20,
                                minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1) {

  if (!inherits(gsea_env, "GseaEnvPro")) stop("❌ 请传入标准 GseaEnvPro 胶囊对象！")
  super_tag <- gsea_env$geneset$name

  series_dir <- file.path(output_dir, custom_series_name)
  if (!dir.exists(series_dir)) dir.create(series_dir, recursive = TRUE, showWarnings = FALSE)
  rds_name <- sprintf("GSEA_Capsule_[%s]_[%s].rds", custom_series_name, super_tag)
  rds_path <- file.path(series_dir, rds_name)

  if (file.exists(rds_path) && force == FALSE) {
    message(sprintf("✅ 命中缓存！检测到已存在的 GSEA 胶囊: %s", rds_name))
    return(readRDS(rds_path))
  }

  total_cores <- parallel::detectCores(logical = TRUE)
  use_cores <- if (is.null(workers)) max(1, total_cores - 4) else min(total_cores, max(1, workers))

  # 🌟 [修复点] 加入这一行，防止复杂网络导致 future 抛出超限错误
  options(future.globals.maxSize = 32000 * 1024^2)
  future::plan(future::multisession, workers = use_cores)

  # ... (中间的 future_lapply 循环原封不动保留您 txt 的逻辑) ...

  # 🌟 [修复点] 修复包名提示 (DudaliRnaseq 改为 GSEAlens)
  cat(sprintf("my_capsule <- readRDS(\"%s\")\n", rds_path))
  cat(sprintf("my_task <- GSEAlens::extract_gsea_task_pro(my_capsule, \"%s\")\n", example_task))
  cat(sprintf("GSEAlens::plot_directional_gsea(my_task, target_pathways = c(\"ID_1\", \"ID_2\"))\n"))

  return(final_obj)
}

#' @title Pro 数据提取器 (文章精细绘图专用)
#' @param gsea_capsule DudaliGseaResPro 计算胶囊
#' @param task_name 想提取的组别名
#' @export
extract_gsea_task_pro <- function(gsea_capsule, task_name) {
  # (直接复制 txt 中的 extract_gsea_task_pro)
}

#' @title GSEA 全流程一键终极总揽引擎 (All-in-One)
#' @export
batch_parallel_gsea_pro <- function(...) {
  # 🌟 [修复点] 修改内部调用的包名为 GSEAlens
  calc_res <- GSEAlens::batch_calc_gsea_pro(...)
  GSEAlens::batch_plot_gsea_pro(...)
  return(invisible(calc_res))
}
