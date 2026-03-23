#==============================================================================
# 文件: R/utils-performance.R
# 描述: GSEA运行性能监控工具（Windows/Linux/Mac跨平台）
# 作者: 独大立 (GSEAlens)
#==============================================================================

#------------------------------------------------------------------------------
# 监控守护进程（后台进程函数，通过callr调用）
#------------------------------------------------------------------------------

#' @title 监控守护进程主函数
#' @description 在独立R进程中运行，定期采样系统资源使用情况
#' @param parent_pid 父进程ID（用于检测父进程存活）
#' @param output_file CSV输出文件路径
#' @param interval_sec 采样间隔（秒）
#' @param max_duration_sec 最大监控时长（安全限制）
#' @keywords internal
monitor_daemon_main <- function(parent_pid, output_file, interval_sec = 0.5,
                                max_duration_sec = 36000) {

  # 加载必要的包（在子进程中）
  if (!requireNamespace("ps", quietly = TRUE)) {
    stop("监控进程：未安装ps包")
  }

  # 获取父进程句柄（用于存活检测）
  parent_handle <- tryCatch(
    ps::ps_handle(parent_pid),
    error = function(e) NULL
  )

  if (is.null(parent_handle)) {
    warning("监控进程：无法获取父进程句柄，启动失败")
    return(invisible(NULL))
  }

  # 获取自身进程句柄
  self_handle <- ps::ps_handle()

  # 初始化CSV文件（写入表头）
  header <- data.frame(
    timestamp = character(),
    rss_mb = numeric(),
    vms_mb = numeric(),
    cpu_percent_total = numeric(),
    cpu_percent_normalized = numeric(),
    workers = integer(),
    stringsAsFactors = FALSE
  )
  utils::write.csv(header, output_file, row.names = FALSE, quote = TRUE)

  # 主监控循环
  start_time <- Sys.time()
  sample_count <- 0

  while (TRUE) {
    # 检查父进程是否仍然存活
    if (!ps::ps_is_running(parent_handle)) {
      message("监控进程：父进程已终止，监控结束")
      break
    }

    # 检查超时
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    if (elapsed > max_duration_sec) {
      message("监控进程：达到最大监控时长，自动停止")
      break
    }

    # 采样系统状态
    tryCatch({
      mem_info <- ps::ps_memory_info(self_handle)
      cpu_info <- ps::ps_cpu_times(self_handle)

      # 计算CPU使用率（需要前后两次采样，首次为NA）
      current_time <- Sys.time()
      static_info <- ps::ps_info(self_handle)
      workers_count <- as.integer(static_info$n_threads %||% 1L)

      # 简化的CPU计算（基于ps的累积时间估算）
      # 注意：准确CPU需要连续采样差值，这里使用单次相对值
      cpu_total <- NA_real_
      cpu_normalized <- NA_real

      # 写入数据（追加模式）
      data_row <- data.frame(
        timestamp = format(current_time, "%Y-%m-%d %H:%M:%OS3"),
        rss_mb = round(mem_info$rss / 1024 / 1024, 2),  # 字节转MB
        vms_mb = round(mem_info$vms / 1024 / 1024, 2),
        cpu_percent_total = cpu_total,
        cpu_percent_normalized = cpu_normalized,
        workers = workers_count,
        stringsAsFactors = FALSE
      )

      # 使用write.table追加（比write.csv快，无表头）
      utils::write.table(
        data_row,
        file = output_file,
        append = TRUE,
        sep = ",",
        row.names = FALSE,
        col.names = FALSE,
        quote = TRUE
      )

      sample_count <- sample_count + 1

    }, error = function(e) {
      # 采样失败，记录警告但不中断
      warning("监控采样失败: ", e$message)
    })

    # 休眠指定间隔
    Sys.sleep(interval_sec)
  }

  # 返回统计信息
  invisible(list(
    samples = sample_count,
    duration = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
    output_file = output_file
  ))
}


#------------------------------------------------------------------------------
# 监控管理器（S3类）
#------------------------------------------------------------------------------

#' @title 创建性能监控管理器
#' @description 创建监控管理器实例，管理后台监控进程的生命周期
#' @param rds_path 关联的RDS文件路径（用于生成CSV路径）
#' @param workers 当前使用的worker数量
#' @param interval_sec 采样间隔
#' @return 监控管理器对象（S3类，继承自environment）
#' @keywords internal
new_performance_monitor <- function(rds_path, workers, interval_sec = 0.5) {

  # 创建环境作为对象容器（轻量级OOP）
  monitor <- new.env(parent = emptyenv())

  # 存储配置
  monitor$rds_path <- rds_path
  monitor$workers <- as.integer(workers)
  monitor$interval_sec <- interval_sec
  monitor$process <- NULL  # callr进程对象
  monitor$csv_path <- NULL

  # 生成CSV文件路径（与RDS同目录，包含workers数）
  monitor$csv_path <- .generate_monitor_path(rds_path, workers)

  # 绑定方法到环境
  monitor$start <- function() {
    if (!is.null(monitor$process) && monitor$process$is_alive()) {
      warning("监控进程已在运行")
      return(invisible(monitor))
    }

    # 检查ps包（严格模式）
    if (!requireNamespace("ps", quietly = TRUE)) {
      stop("性能监控已启用，但未安装ps包。请运行: install.packages('ps')")
    }

    if (!requireNamespace("callr", quietly = TRUE)) {
      stop("性能监控需要callr包。请运行: install.packages('callr')")
    }

    # 确保目录存在
    dir.create(dirname(monitor$csv_path), showWarnings = FALSE, recursive = TRUE)

    # 处理文件重名（追加_1, _2等）
    monitor$csv_path <- .resolve_monitor_path(monitor$csv_path)

    # 启动后台监控进程
    parent_pid <- Sys.getpid()

    tryCatch({
      monitor$process <- callr::r_bg(
        func = monitor_daemon_main,
        args = list(
          parent_pid = parent_pid,
          output_file = monitor$csv_path,
          interval_sec = monitor$interval_sec
        ),
        supervise = TRUE,  # 父进程崩溃时自动清理子进程
        cleanup = TRUE
      )

      # 等待进程启动（检查文件是否创建）
      max_wait <- 10  # 最多等5秒
      waited <- 0
      while (!file.exists(monitor$csv_path) && waited < max_wait) {
        Sys.sleep(0.5)
        waited <- waited + 1
      }

      if (!file.exists(monitor$csv_path)) {
        stop("监控进程启动失败：未检测到输出文件")
      }

      message("📊 性能监控已启动: ", basename(monitor$csv_path))

    }, error = function(e) {
      stop("无法启动性能监控: ", e$message)
    })

    invisible(monitor)
  }

  monitor$stop <- function() {
    if (is.null(monitor$process)) {
      return(invisible(monitor))
    }

    # 给监控进程时间完成最后写入（方案1：简单休眠）
    Sys.sleep(1)

    # 终止进程
    if (monitor$process$is_alive()) {
      monitor$process$kill()
    }

    # 验证输出文件
    if (file.exists(monitor$csv_path)) {
      file_size <- file.size(monitor$csv_path)
      message("📊 性能监控已停止: ", round(file_size / 1024, 1), " KB 数据")
    } else {
      warning("监控输出文件未找到")
    }

    invisible(monitor)
  }

  monitor$get_summary <- function() {
    # 读取CSV并计算摘要统计
    if (!file.exists(monitor$csv_path)) {
      return(NULL)
    }

    tryCatch({
      data <- utils::read.csv(monitor$csv_path, stringsAsFactors = FALSE)

      if (nrow(data) == 0) {
        return(NULL)
      }

      # 计算关键指标
      list(
        duration_sec = as.numeric(
          difftime(
            as.POSIXct(tail(data$timestamp, 1)),
            as.POSIXct(data$timestamp[1]),
            units = "secs"
          )
        ),
        mem_peak_mb = max(data$rss_mb, na.rm = TRUE),
        mem_avg_mb = mean(data$rss_mb, na.rm = TRUE),
        cpu_max_pct = max(data$cpu_percent_total, na.rm = TRUE),
        cpu_avg_pct = mean(data$cpu_percent_total, na.rm = TRUE),
        samples = nrow(data),
        monitor_file = monitor$csv_path
      )
    }, error = function(e) {
      warning("读取监控数据失败: ", e$message)
      NULL
    })
  }

  # 设置S3类属性
  class(monitor) <- c("GseaPerformanceMonitor", "environment")
  monitor
}


#------------------------------------------------------------------------------
# 辅助函数
#------------------------------------------------------------------------------

#' @title 生成监控文件路径
#' @keywords internal
.generate_monitor_path <- function(rds_path, workers) {
  # 解析RDS路径：GSEA_Capsule_[name]_[geneset].rds
  dir <- dirname(rds_path)
  base <- basename(rds_path)

  # 替换前缀
  monitor_base <- sub("^GSEA_Capsule_", "GSEA_Monitor_", base)

  # 移除.rds扩展名，添加_workers.csv
  monitor_base <- sub("\\.rds$", paste0("_", workers, ".csv"), monitor_base)

  file.path(dir, monitor_base)
}

#' @title 解决文件路径冲突（追加编号）
#' @keywords internal
.resolve_monitor_path <- function(csv_path) {
  if (!file.exists(csv_path)) {
    return(csv_path)
  }

  # 分解路径
  dir <- dirname(csv_path)
  base <- basename(csv_path)

  # 移除扩展名
  base_no_ext <- sub("\\.csv$", "", base)
  ext <- ".csv"

  # 查找可用编号
  counter <- 1
  while (TRUE) {
    new_path <- file.path(dir, paste0(base_no_ext, "_", counter, ext))
    if (!file.exists(new_path)) {
      return(new_path)
    }
    counter <- counter + 1
    if (counter > 100) {
      stop("无法生成监控文件路径：编号冲突过多")
    }
  }
}

#' @title 空值处理辅助函数
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x


#------------------------------------------------------------------------------
# 打印方法
#------------------------------------------------------------------------------

#' @export
print.GseaPerformanceMonitor <- function(x, ...) {
  cat("GSEA性能监控管理器\n")
  cat("  CSV路径: ", x$csv_path, "\n")
  cat("  Workers: ", x$workers, "\n")
  cat("  状态: ", if(!is.null(x$process) && x$process$is_alive()) "运行中" else "已停止", "\n")
  invisible(x)
}
