#' @title 批量并行 GSEA 计算（带性能监控）
#' @description 消费标准 GseaEnv 对象，执行并行 GSEA 计算。
#'   支持系统级性能监控，用于基准测试和发表文章。
#' @param gsea_env GseaEnv 对象 (由 setup_gsea_env 生成)
#' @param custom_series_name 字符串。分析系列名称，用于结果存档。
#' @param output_dir 字符串。结果输出路径，默认 "./GSEA_Output"。
#' @param workers 并行核心数。默认 4。
#' @param bidirectional 逻辑值。是否自动生成反向对比，默认 TRUE。
#' @param minGSSize 最小基因集大小，默认 10。
#' @param maxGSSize 最大基因集大小，默认 500。
#' @param pvalueCutoff P 值阈值，默认 1 (保留全量结果)。
#' @param force 逻辑值。是否强制重新计算，默认 FALSE。
#' @param enable_monitoring 逻辑值。是否启用系统级性能监控（需要ps包），默认 FALSE。
#' @param monitor_interval 数值。监控采样间隔（秒），默认 0.5。
#' @return GseaRes 对象，包含 gsea_running_time_monitor 字段（如启用监控）。
#' @export
batch_calc_gsea <- function(gsea_env,
                            custom_series_name = "Auto_Analysis",
                            output_dir = "./GSEA_Output",
                            workers = 4,
                            bidirectional = TRUE,
                            minGSSize = 10,
                            maxGSSize = 500,
                            pvalueCutoff = 1,
                            force = FALSE,
                            enable_monitoring = FALSE,
                            monitor_interval = 0.5) {

  # 1. 校验输入
  .check_gsea_env(gsea_env)

  # 2. 初始化性能监控（如启用）
  perf_monitor <- NULL
  rds_path <- NULL  # 提前声明，供监控使用

  if (isTRUE(enable_monitoring)) {
    # 预计算RDS路径（监控需要知道输出位置）
    series_dir <- file.path(output_dir, custom_series_name)
    rds_name <- sprintf("GSEA_Capsule_[%s]_[%s].rds", custom_series_name, gsea_env$geneset$name)
    rds_path <- file.path(series_dir, rds_name)

    # 创建监控管理器（但不立即启动，等确认非缓存后再启动）
    perf_monitor <- new_performance_monitor(rds_path, workers, monitor_interval)
  }

  # 3. 智能缓存拦截
  series_dir <- file.path(output_dir, custom_series_name)
  if (!dir.exists(series_dir)) {
    dir.create(series_dir, recursive = TRUE, showWarnings = FALSE)
  }

  rds_name <- sprintf("GSEA_Capsule_[%s]_[%s].rds", custom_series_name, gsea_env$geneset$name)
  rds_path <- file.path(series_dir, rds_name)

  if (file.exists(rds_path) && !force) {
    message(sprintf("✅ 命中缓存！检测到已存在的 GSEA 胶囊: %s", rds_name))
    return(readRDS(rds_path))
  }

  # 缓存未命中，启动监控（如果启用）
  if (!is.null(perf_monitor)) {
    perf_monitor$start()

    # 设置退出时自动停止监控（防止遗漏）
    on.exit({
      if (!is.null(perf_monitor)) {
        perf_monitor$stop()
      }
    }, add = TRUE)
  }

  # 4. 准备任务列表（原有代码，保持不变）
  registry <- gsea_env$contrast_registry
  de_store <- gsea_env$de_store

  tasks <- list()
  for (i in 1:nrow(registry)) {
    row <- registry[i, ]
    cid <- row$contrast_id

    tasks[[cid]] <- list(
      contrast_id = cid,
      left = row$left_group,
      right = row$right_group,
      reverse = FALSE
    )

    if (bidirectional) {
      rev_cid <- paste(row$right_group, row$left_group, sep = "_vs_")
      tasks[[rev_cid]] <- list(
        contrast_id = rev_cid,
        left = row$right_group,
        right = row$left_group,
        reverse = TRUE
      )
    }
  }

  message(sprintf("🚀 准备就绪：共 %d 个对比任务待计算...", length(tasks)))

  # 5. 并行计算设置（原有代码）
  total_cores <- parallel::detectCores(logical = TRUE)
  use_cores <- min(total_cores, max(1, workers))
  message(sprintf("🖥️ 硬件调度: 使用 %d 核进行并行计算...", use_cores))

  options(future.globals.maxSize = 32000 * 1024^2)
  future::plan(future::multisession, workers = use_cores)

  # 记录开始时间（用于计算总耗时）
  global_start <- Sys.time()

  # 6. 核心计算循环（原有代码，保留所有逻辑）
  term2gene <- gsea_env$geneset$term2gene
  meta_dict <- gsea_env$geneset$meta_dict

  res_list <- future.apply::future_lapply(names(tasks), function(task_name) {
    # ... [保留原有所有计算逻辑] ...
    # 注意：内部使用.prepare_rank_vector等辅助函数
    # 最终返回list(name=, status=, data=, genelist=)
  }, future.seed = TRUE)

  names(res_list) <- names(tasks)
  future::plan(future::sequential)

  # 计算总耗时
  global_end <- Sys.time()
  total_duration <- as.numeric(difftime(global_end, global_start, units = "secs"))

  # 7. 停止监控并获取摘要（如启用）
  gsea_running_time_monitor <- NULL
  if (!is.null(perf_monitor)) {
    perf_monitor$stop()
    monitor_summary <- perf_monitor$get_summary()

    if (!is.null(monitor_summary)) {
      gsea_running_time_monitor <- list(
        duration_sec = monitor_summary$duration_sec,
        mem_peak_mb = monitor_summary$mem_peak_mb,
        cpu_max_pct = monitor_summary$cpu_max_pct,
        monitor_file = monitor_summary$monitor_file
      )

      # 打印监控摘要
      message("\n", rep("=", 50))
      message("           📊 性能监控摘要")
      message(rep("=", 50))
      message(sprintf("总耗时: %.2f 秒", gsea_running_time_monitor$duration_sec))
      message(sprintf("峰值内存: %.2f MB", gsea_running_time_monitor$mem_peak_mb))
      message(sprintf("最大CPU占用: %.2f%%", gsea_running_time_monitor$cpu_max_pct))
      message(sprintf("详细数据: %s", basename(gsea_running_time_monitor$monitor_file)))
      message(rep("=", 50))
    } else {
      # 监控失败，标记NA
      gsea_running_time_monitor <- list(
        duration_sec = total_duration,  # 至少记录总时间
        mem_peak_mb = NA_real_,
        cpu_max_pct = NA_real_,
        monitor_file = NA_character_
      )
    }

    # 从on.exit中移除（已手动停止）
    on.exit(NULL)
  }

  # 8. 封装结果对象（修改metadata部分）
  final_obj <- create_gsea_res(
    metadata = list(
      run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      workers_used = use_cores,
      parameters = list(
        bidirectional = bidirectional,
        minGSSize = minGSSize,
        enable_monitoring = enable_monitoring
      ),
      project_info = list(
        custom_series_name = custom_series_name,
        output_dir = normalizePath(output_dir, mustWork = FALSE),
        series_dir = normalizePath(series_dir, mustWork = FALSE),
        rds_path = normalizePath(rds_path, mustWork = FALSE)
      ),
      # 新增：性能监控数据（长命名，易查看）
      gsea_running_time_monitor = gsea_running_time_monitor
    ),
    backend_info = gsea_env$backend_info,
    contrast_registry = gsea_env$contrast_registry,
    de_store = gsea_env$de_store,
    expr_bundle = gsea_env$expr_bundle,
    geneset_info = gsea_env$geneset,
    results = res_list
  )

  # 9. 保存结果
  saveRDS(final_obj, rds_path)
  message(sprintf("\n✅ 计算完成！结果已保存至: %s", rds_path))

  if (isTRUE(enable_monitoring) && !is.null(gsea_running_time_monitor$monitor_file)) {
    message("📊 性能监控数据已保存，可用于基准测试分析")
  }

  return(final_obj)
}



#' @title 准备排序向量
#' @description 从 DE 表中提取 stat 列，进行清洗和排序。
#' @param de_table 差异分析表 (必须包含 gene_symbol 和 stat 列)
#' @param flip 是否翻转符号 (用于反向对比)
#' @return 排序后的命名向量
#' @keywords internal
.prepare_rank_vector <- function(de_table, flip = FALSE) {

  # 1. 清洗
  df <- de_table %>%
    dplyr::filter(!is.na(gene_symbol) & gene_symbol != "") %>%
    dplyr::mutate(gene_symbol = toupper(gene_symbol)) # 统一转大写

  # 2. 去重 (按 stat 绝对值最大的保留)
  df <- df %>%
    dplyr::arrange(dplyr::desc(abs(stat))) %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE)

  # 3. 提取向量
  vals <- df$stat
  names(vals) <- df$gene_symbol

  # 4. 翻转 (如果是反向对比)
  if (flip) {
    vals <- -vals
  }

  # 5. 降序排序
  vals <- sort(vals, decreasing = TRUE)

  return(vals)
}
