#==============================================================================
# 文件: R/utils-benchmark.R
# 描述: GSEA基准测试数据对齐与可视化工具
#==============================================================================

#' @title 对齐GSEA性能数据与系统监控数据
#' @description 将batch_calc_gsea生成的时间戳与独立监控脚本的CSV数据对齐，
#'   生成可用于发表文章的时间序列数据框。
#' @param gsea_res GseaRes对象（必须包含metadata$gsea_benchmark）
#' @param monitor_csv 监控脚本生成的CSV文件路径
#' @return 对齐后的数据框（包含relative_sec等时间轴统一字段）
#' @export
#' @examples
#' \dontrun{
#' res <- readRDS("GSEA_Capsule_[...].rds")
#' aligned <- align_benchmark_data(res, "system_monitor.csv")
#' plot(aligned$relative_sec, aligned$rss_mb, type = "l")
#' }
align_benchmark_data <- function(gsea_res, monitor_csv) {

  if (!inherits(gsea_res, "GseaRes")) {
    stop("输入必须是GseaRes对象")
  }

  if (!file.exists(monitor_csv)) {
    stop("监控CSV不存在: ", monitor_csv)
  }

  # 提取GSEA时间窗口
  bench <- gsea_res$metadata$gsea_benchmark
  start_ms <- bench$start_ms
  end_ms <- bench$end_ms

  # 读取Python生成的CSV（格式更简洁）
  monitor <- read.csv(monitor_csv, stringsAsFactors = FALSE)

  # Python和R的时间戳都是Unix毫秒，直接对齐
  aligned <- monitor[monitor$timestamp_ms >= (start_ms - 5000) &
                       monitor$timestamp_ms <= (end_ms + 5000), ]

  if (nrow(aligned) == 0) {
    warning("未找到时间重叠数据")
    return(NULL)
  }

  # 转换相对时间（秒）
  aligned$relative_sec <- (aligned$timestamp_ms - start_ms) / 1000

  # 阶段标记
  total_duration <- (end_ms - start_ms) / 1000
  aligned$phase <- ifelse(
    aligned$relative_sec < 0, "pre_start",
    ifelse(aligned$relative_sec > total_duration, "post_end",
           ifelse(aligned$relative_sec < total_duration * 0.2, "startup",
                  ifelse(aligned$relative_sec < total_duration * 0.8, "core_compute", "cleanup")
           )
    )
  )

  attr(aligned, "gsea_info") <- bench
  class(aligned) <- c("GseaBenchmarkAligned", "data.frame")

  message(sprintf("✅ 对齐完成: %d 采样点，跨度 %.1f 秒",
                  nrow(aligned), max(aligned$relative_sec) - min(aligned$relative_sec)))

  invisible(aligned)
}


#' @title 绘制GSEA基准测试内存曲线
#' @description 为发表文章生成标准内存使用曲线图
#' @param aligned_data align_benchmark_data返回的对齐数据框
#' @param highlight_phases 是否按阶段着色，默认TRUE
#' @return ggplot对象
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_vline labs theme_minimal scale_color_manual
plot_gsea_memory <- function(aligned_data, highlight_phases = TRUE) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("需要ggplot2包: install.packages('ggplot2')")
  }

  info <- attr(aligned_data, "gsea_info")

  p <- ggplot2::ggplot(aligned_data, ggplot2::aes(x = relative_sec, y = rss_mb)) +
    ggplot2::geom_line(linewidth = 0.8, color = "steelblue") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "darkgreen") +
    ggplot2::geom_vline(xintercept = info$duration_sec, linetype = "dashed", color = "darkred") +
    ggplot2::labs(
      title = sprintf("GSEA Memory Usage Profile (%d cores)", info$workers),
      subtitle = sprintf("Total duration: %.1f seconds", info$duration_sec),
      x = "Time (seconds from GSEA start)",
      y = "Memory (MB)",
      caption = sprintf("Peak: %.0f MB | Average: %.0f MB",
                        max(aligned_data$rss_mb, na.rm = TRUE),
                        mean(aligned_data$rss_mb[aligned_data$phase %in% c("startup", "core_compute", "cleanup")], na.rm = TRUE))
    ) +
    ggplot2::theme_minimal()

  if (highlight_phases) {
    # 添加阶段背景色（高级用法，可选）
    p <- p + ggplot2::aes(color = phase) +
      ggplot2::scale_color_manual(values = c(
        pre_start = "gray70",
        startup = "green3",
        core_compute = "steelblue",
        cleanup = "orange",
        post_end = "gray70"
      ), name = "Phase")
  }

  p
}


#' @title 打印GSEA基准测试对齐结果
#' @export
print.GseaBenchmarkAligned <- function(x, ...) {
  info <- attr(x, "gsea_info")

  cat("========================================\n")
  cat("   GSEA 基准测试对齐结果\n")
  cat("========================================\n")
  cat(sprintf("GSEA计算时间: %s\n", info$start_time))
  cat(sprintf("并行核心数:   %d\n", info$workers))
  cat(sprintf("GSEA总耗时:   %.2f 秒\n", info$duration_sec))
  cat("----------------------------------------\n")
  cat(sprintf("监控采样点:   %d\n", nrow(x)))
  cat(sprintf("监控时间跨度: %.1f 秒 (%.1f 至 %.1f)\n",
              max(x$relative_sec) - min(x$relative_sec),
              min(x$relative_sec), max(x$relative_sec)))
  cat(sprintf("峰值内存:     %.0f MB\n", max(x$rss_mb, na.rm = TRUE)))
  cat(sprintf("平均内存:     %.0f MB\n", mean(x$rss_mb[x$phase %in% c("startup", "core_compute", "cleanup")], na.rm = TRUE)))
  cat("========================================\n")

  # 各阶段统计
  cat("\n各阶段内存统计:\n")
  stage_stats <- tapply(x$rss_mb, x$phase, function(v)
    sprintf("均值 %.0f MB, 峰值 %.0f MB", mean(v, na.rm = TRUE), max(v, na.rm = TRUE))
  )
  for (phase in names(stage_stats)) {
    if (!is.na(stage_stats[phase])) {
      cat(sprintf("  [%s]: %s\n", phase, stage_stats[phase]))
    }
  }

  invisible(x)
}


#' @title 生成扩展性分析报告（多核对比）
#' @description 对比不同核心数下的GSEA性能（用于发表文章）
#' @param gsea_res_list 多个GseaRes对象的列表（命名如 list("1_core" = res1, "4_core" = res4)）
#' @param monitor_csv_list 对应的监控CSV路径列表（与gsea_res_list同名）
#' @return 对比数据框
#' @export
analyze_scalability <- function(gsea_res_list, monitor_csv_list = NULL) {

  if (is.null(names(gsea_res_list))) {
    names(gsea_res_list) <- paste0("run_", seq_along(gsea_res_list))
  }

  results <- lapply(names(gsea_res_list), function(name) {
    res <- gsea_res_list[[name]]

    if (is.null(res$metadata$gsea_benchmark)) {
      warning(name, " 缺少benchmark数据，跳过")
      return(NULL)
    }

    bench <- res$metadata$gsea_benchmark

    # 基础信息
    out <- data.frame(
      run_name = name,
      workers = bench$workers,
      duration_sec = bench$duration_sec,
      stringsAsFactors = FALSE
    )

    # 如果有监控数据，添加内存统计
    if (!is.null(monitor_csv_list) && name %in% names(monitor_csv_list)) {
      csv_path <- monitor_csv_list[[name]]
      if (file.exists(csv_path)) {
        monitor <- read.csv(csv_path, stringsAsFactors = FALSE)
        # 筛选GSEA时间段
        aligned <- monitor[monitor$timestamp_ms >= bench$start_ms &
                             monitor$timestamp_ms <= bench$end_ms, ]
        if (nrow(aligned) > 0) {
          out$mem_peak_mb <- max(aligned$rss_mb, na.rm = TRUE)
          out$mem_avg_mb <- mean(aligned$rss_mb, na.rm = TRUE)
        }
      }
    }

    out
  })

  do.call(rbind, results)
}
