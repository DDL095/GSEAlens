#' @title 旧版本批量并行 GSEA 计算
#' @description 消费标准 GseaEnv 对象，执行并行 GSEA 计算。
#' @param gsea_env GseaEnv 对象 (由 setup_gsea_env 生成)
#' @param custom_series_name 字符串。分析系列名称，用于结果存档。
#' @param output_dir 字符串。结果输出路径，默认 "./GSEA_Output"。
#' @param workers 并行核心数。默认 4。
#' @param bidirectional 逻辑值。是否自动生成反向对比，默认 TRUE。
#' @param minGSSize 最小基因集大小，默认 10。
#' @param maxGSSize 最大基因集大小，默认 500。
#' @param pvalueCutoff P 值阈值，默认 1 (保留全量结果)。
#' @param force 逻辑值。是否强制重新计算，默认 FALSE。
#' @return GseaRes 对象，包含 metadata$gsea_benchmark 时间戳字段。
#' @noRd

batch_calc_gsea_old <- function(gsea_env,
                            custom_series_name = "Auto_Analysis",
                            output_dir = "./GSEA_Output",
                            workers = 4,
                            bidirectional = TRUE,
                            minGSSize = 10,
                            maxGSSize = 500,
                            pvalueCutoff = 1,
                            force = FALSE) {

  # 校验输入
  .check_gsea_env(gsea_env)

  # 记录高精度开始时间（用于外部监控对齐）
  start_time <- Sys.time()
  start_ms <- as.numeric(start_time) * 1000

  # 智能缓存拦截
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

  # 准备任务列表
  registry <- gsea_env$contrast_registry
  de_store <- gsea_env$de_store

  tasks <- list()

  for (i in 1:nrow(registry)) {
    row <- registry[i, ]
    cid <- row$contrast_id

    # 正向任务
    tasks[[cid]] <- list(
      contrast_id = cid,
      left = row$left_group,
      right = row$right_group,
      reverse = FALSE
    )

    # 反向任务
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

  # 并行计算设置
  total_cores <- parallel::detectCores(logical = TRUE)
  use_cores <- min(total_cores, max(1, workers))
  message(sprintf("🖥️ 硬件调度: 使用 %d 核进行并行计算...", use_cores))

  options(future.globals.maxSize = 32000 * 1024^2)
  future::plan(future::multisession, workers = use_cores)

  # 核心计算循环
  term2gene <- gsea_env$geneset$term2gene
  meta_dict <- gsea_env$geneset$meta_dict

  res_list <- future.apply::future_lapply(names(tasks), function(task_name) {

    task_info <- tasks[[task_name]]
    cid <- task_info$contrast_id

    # 获取 DE 表
    original_cid <- if (task_info$reverse) {
      paste(task_info$right, task_info$left, sep = "_vs_")
    } else {
      cid
    }

    de_table <- de_store[[original_cid]]

    if (is.null(de_table)) {
      warning(sprintf("无法找到 %s 的差异分析表，跳过。", original_cid))
      return(list(name = cid, status = "Failed", data = NULL, genelist = c()))
    }

    # 准备 Rank Vector
    genelist <- .prepare_rank_vector(de_table, flip = task_info$reverse)

    if (length(genelist) == 0) {
      return(list(name = cid, status = "Failed", data = NULL, genelist = c()))
    }

    # 运行 GSEA
    set.seed(123)
    gsea_res <- tryCatch({
      clusterProfiler::GSEA(
        geneList = genelist,
        TERM2GENE = term2gene,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        pvalueCutoff = pvalueCutoff,
        pAdjustMethod = "BH",
        verbose = FALSE,
        seed = 123,
        eps = 0
      )
    }, error = function(e) NULL)

    status <- if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) "Success" else "Failed/NoEnrich"

    # 注入元数据
    if (status == "Success" && !is.null(meta_dict)) {
      res_df <- as.data.frame(gsea_res@result)
      if ("Description" %in% colnames(res_df)) res_df <- dplyr::select(res_df, -Description)

      res_df <- res_df %>%
        dplyr::left_join(as.data.frame(meta_dict), by = "ID") %>%
        dplyr::mutate(
          Display_Collection = if("Combo_Name" %in% names(.)) Combo_Name else if("Collection" %in% names(.)) Collection else "Unknown",
          Display_Collection = as.factor(ifelse(is.na(Display_Collection), "Unknown", Display_Collection)),
          Pathway_Link = if("URL" %in% names(.)) {
            ifelse(is.na(URL) | URL == "", sprintf("<b>%s</b>", ID),
                   sprintf('<a href="%s" target="_blank">%s</a>', URL, ID))
          } else { sprintf("<b>%s</b>", ID) },
          Description = if("Description.y" %in% names(.)) Description.y else ID
        )

      rownames(res_df) <- res_df$ID
      gsea_res@result <- res_df
    }

    return(list(name = cid, status = status, data = gsea_res, genelist = genelist))

  }, future.seed = NULL)# use 123 as seed

  names(res_list) <- names(tasks)
  future::plan(future::sequential)

  # 记录结束时间
  end_time <- Sys.time()
  end_ms <- as.numeric(end_time) * 1000

  # 封装结果对象
  final_obj <- create_gsea_res(
    metadata = list(
      run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      workers_used = use_cores,
      parameters = list(bidirectional = bidirectional, minGSSize = minGSSize),
      project_info = list(
        custom_series_name = custom_series_name,
        output_dir = normalizePath(output_dir, mustWork = FALSE),
        series_dir = normalizePath(series_dir, mustWork = FALSE),
        rds_path = normalizePath(rds_path, mustWork = FALSE)
      ),
      # 新增：极简性能时间戳（用于外部监控对齐）
      gsea_benchmark = list(
        start_time = format(start_time, "%Y-%m-%d %H:%M:%OS3"),
        start_ms = start_ms,
        end_time = format(end_time, "%Y-%m-%d %H:%M:%OS3"),
        end_ms = end_ms,
        duration_sec = round((end_ms - start_ms) / 1000, 3),
        workers = use_cores
      )
    ),
    backend_info = gsea_env$backend_info,
    contrast_registry = gsea_env$contrast_registry,
    de_store = gsea_env$de_store,
    expr_bundle = gsea_env$expr_bundle,
    geneset_info = gsea_env$geneset,
    results = res_list
  )

  # 保存结果
  saveRDS(final_obj, rds_path)
  message(sprintf("\n✅ 计算完成！耗时 %.2f 秒，结果已保存至: %s",
                  final_obj$metadata$gsea_benchmark$duration_sec, rds_path))

  return(final_obj)
}


#' @title 准备排序向量_与旧版本配合
#' @description 从 DE 表中提取 stat 列，进行清洗和排序。
#' @param de_table 差异分析表 (必须包含 gene_symbol 和 stat 列)
#' @param flip 是否翻转符号 (用于反向对比)
#' @return 排序后的命名向量
#' @keywords internal

.prepare_rank_vector <- function(de_table, flip = FALSE) {

  # 清洗
  df <- de_table %>%
    dplyr::filter(!is.na(gene_symbol) & gene_symbol != "") %>%
    dplyr::mutate(gene_symbol = toupper(gene_symbol))

  # 去重 (按 stat 绝对值最大的保留)
  df <- df %>%
    dplyr::arrange(dplyr::desc(abs(stat))) %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE)

  # 提取向量
  vals <- df$stat
  names(vals) <- df$gene_symbol

  # 翻转 (如果是反向对比)
  if (flip) {
    vals <- -vals
  }

  # 降序排序
  vals <- sort(vals, decreasing = TRUE)

  return(vals)
}


#' @title 批量并行 GSEA 计算（优化版）
#' @description 消费标准 GseaEnv 对象，执行高效并行 GSEA 计算。
#' @param gsea_env GseaEnv 对象
#' @param custom_series_name 字符串。分析系列名称
#' @param output_dir 字符串。结果输出路径，默认 "./GSEA_Output"
#' @param workers 并行核心数。默认 4
#' @param bidirectional 逻辑值。是否自动生成反向对比，默认 TRUE
#' @param minGSSize 最小基因集大小，默认 10
#' @param maxGSSize 最大基因集大小，默认 500
#' @param pvalueCutoff P 值阈值，默认 1
#' @param force 逻辑值。是否强制重新计算，默认 FALSE
#' @param use_progress 逻辑值。是否显示进度条，默认 TRUE
#' @param chunk_size 整数。每个worker一次处理的任务数，默认 NULL（自动）
#' @return GseaRes 对象
#' @export
#' @title 批量并行 GSEA 计算（优化版 - 修复全局变量大小错误）
#' @description 消费标准 GseaEnv 对象，执行高效并行 GSEA 计算。
#' @export
#' @title 批量并行 GSEA 计算（带实时进度条）
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
                                      use_progress = F,  # 控制是否显示进度
                                      chunk_size = NULL) {


  # 第一步：参数校验和初始化（保持不变）

  .check_gsea_env(gsea_env)

  start_time <- Sys.time()
  start_ms <- as.numeric(start_time) * 1000

  series_dir <- file.path(output_dir, custom_series_name)
  if (!dir.exists(series_dir)) {
    dir.create(series_dir, recursive = TRUE, showWarnings = FALSE)
  }

  rds_name <- sprintf("GSEA_Capsule_[%s]_[%s].rds",
                      custom_series_name, gsea_env$geneset$name)
  rds_path <- file.path(series_dir, rds_name)

  if (file.exists(rds_path) && !force) {
    message(sprintf("✅ 命中缓存！检测到已存在的 GSEA 胶囊: %s", rds_name))
    return(readRDS(rds_path))
  }


  # 第二步：构建任务列表（保持不变）

  registry <- gsea_env$contrast_registry
  de_store <- gsea_env$de_store

  task_metadata <- list()
  for (i in seq_len(nrow(registry))) {
    row <- registry[i, ]
    cid <- row$contrast_id
    task_metadata[[cid]] <- list(
      task_id = cid,
      left_group = row$left_group,
      right_group = row$right_group,
      is_reversed = FALSE
    )
    if (bidirectional) {
      rev_cid <- paste(row$right_group, row$left_group, sep = "_vs_")
      task_metadata[[rev_cid]] <- list(
        task_id = rev_cid,
        left_group = row$right_group,
        right_group = row$left_group,
        is_reversed = TRUE
      )
    }
  }

  total_tasks <- length(task_metadata)
  message(sprintf("🚀 准备就绪：共 %d 个对比任务待计算...", total_tasks))


  # 第三步：并行设置（关键修复：progressr 处理器）


  # 设置全局变量大小限制（必须在 plan 之前）
  options(future.globals.maxSize = 192 * 1024^3)#设置无限大

  total_cores <- parallel::detectCores(logical = TRUE)
  use_cores <- min(total_cores, max(1, workers))
  message(sprintf("🖥️ 硬件调度: 使用 %d 核进行并行计算...", use_cores))

  if (is.null(chunk_size)) {
    chunk_size <- max(1, ceiling(total_tasks / (use_cores * 4)))
  }

  # 设置 future 计划
  future::plan(future::multisession, workers = use_cores)

  # 关键修复 1：启用 progressr 处理器（必须在 future_lapply 之前）
  if (use_progress) {
    if (!requireNamespace("progressr", quietly = TRUE)) {
      stop("请安装 progressr 包以显示进度条: install.packages('progressr')")
    }
    progressr::handlers(global = TRUE)  # 启用全局处理器
    progressr::handlers("progress")      # 使用 progress 包样式（或改为 "txtprogressbar"）
  }

  # 准备数据
  worker_term2gene <- gsea_env$geneset$term2gene
  worker_meta_dict <- gsea_env$geneset$meta_dict
  worker_de_list <- as.list(de_store)

  # 构建 chunks
  task_names <- names(task_metadata)
  task_chunks <- split(task_names, ceiling(seq_along(task_names) / chunk_size))

  message(sprintf("📦 分块策略：%d 个chunk，每个chunk最多 %d 个任务",
                  length(task_chunks), chunk_size))


  # 第四步：并行计算（关键修复：使用 progressr）


  # 创建进度条对象（progressr 风格）
  p <- progressr::progressor(steps = total_tasks, enable = use_progress)

  # 关键修复 2：在 worker 函数中调用 p() 来更新进度
  res_list <- future.apply::future_lapply(
    X = task_chunks,
    FUN = function(chunk_task_names,
                   task_metadata,
                   de_list,
                   term2gene,
                   meta_dict,
                   minGSSize,
                   maxGSSize,
                   pvalueCutoff,
                   progressor_fn) {  # 接收进度函数

      # 加载必要的包
      if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
          !requireNamespace("dplyr", quietly = TRUE)) {
        stop("Worker 缺少必要的包")
      }

      chunk_results <- list()

      for (task_name in chunk_task_names) {
        task_info <- task_metadata[[task_name]]
        task_id <- task_info$task_id

        # 获取原始对比ID
        if (task_info$is_reversed) {
          original_cid <- paste(task_info$right_group, task_info$left_group, sep = "_vs_")
        } else {
          original_cid <- task_id
        }

        de_table <- de_list[[original_cid]]

        if (is.null(de_table) || nrow(de_table) == 0) {
          chunk_results[[task_name]] <- list(
            name = task_id,
            status = "Failed",
            data = NULL,
            genelist = c()
          )
          # 关键：即使失败也要更新进度
          if (!is.null(progressor_fn)) progressor_fn()
          next
        }

        # 准备排序向量
        genelist <- tryCatch({
          .prepare_rank_vector_fast(de_table, flip = task_info$is_reversed)
        }, error = function(e) {
          c()
        })

        if (length(genelist) == 0) {
          chunk_results[[task_name]] <- list(
            name = task_id,
            status = "Failed",
            data = NULL,
            genelist = c()
          )
          if (!is.null(progressor_fn)) progressor_fn()
          next
        }

        # 执行 GSEA
        set.seed(123)
        gsea_res <- tryCatch({
          clusterProfiler::GSEA(
            geneList = genelist,
            TERM2GENE = term2gene,
            minGSSize = minGSSize,
            maxGSSize = maxGSSize,
            pvalueCutoff = pvalueCutoff,
            pAdjustMethod = "BH",
            verbose = FALSE,
            seed = 123,
            eps = 0
          )
        }, error = function(e) {
          NULL
        })

        status <- if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) "Success" else "Failed/NoEnrich"

        if (status == "Success" && !is.null(meta_dict)) {
          gsea_res@result <- .enrich_gsea_result(gsea_res@result, meta_dict)
        }

        chunk_results[[task_name]] <- list(
          name = task_id,
          status = status,
          data = gsea_res,
          genelist = genelist
        )

        # 关键：每完成一个任务，更新进度条
        if (!is.null(progressor_fn)) progressor_fn()
      }

      return(chunk_results)
    },
    task_metadata = task_metadata,
    de_list = worker_de_list,
    term2gene = worker_term2gene,
    meta_dict = worker_meta_dict,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    pvalueCutoff = pvalueCutoff,
    progressor_fn = if (use_progress) p else NULL,  # 传递进度函数
    future.seed = TRUE,
    future.scheduling = 1.0
  )


  # 第五步：清理与结果汇总（关键修复：关闭连接）


  # 关闭 future 计划（减少连接警告）
  future::plan(future::sequential)

  # 关键修复 3：主动关闭所有未使用的连接（消除警告）
  all_cons <- showConnections(all = TRUE)
  if (nrow(all_cons) > 1) {  # 第1行通常是 stdout
    for (con_idx in as.integer(rownames(all_cons))) {
      if (con_idx > 2) {  # 保留 stdin(0), stdout(1), stderr(2)
        try(suppressWarnings(close.connection(getConnection(con_idx))), silent = TRUE)
      }
    }
  }

  # 强制 GC 清理（可选，但有助于消除警告）
  invisible(gc(verbose = FALSE, full = TRUE))

  # 展平结果
  res_list_flat <- do.call(c, res_list)
  names(res_list_flat) <- names(task_metadata)

  end_time <- Sys.time()
  end_ms <- as.numeric(end_time) * 1000


  # 第六步：创建结果对象（保持不变）

  final_obj <- create_gsea_res(
    metadata = list(
      run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      workers_used = use_cores,
      chunk_size = chunk_size,
      parameters = list(
        bidirectional = bidirectional,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        pvalueCutoff = pvalueCutoff
      ),
      project_info = list(
        custom_series_name = custom_series_name,
        output_dir = normalizePath(output_dir, mustWork = FALSE),
        series_dir = normalizePath(series_dir, mustWork = FALSE),
        rds_path = normalizePath(rds_path, mustWork = FALSE)
      ),
      gsea_benchmark = list(
        start_time = format(start_time, "%Y-%m-%d %H:%M:%OS3"),
        start_ms = start_ms,
        end_time = format(end_time, "%Y-%m-%d %H:%M:%OS3"),
        end_ms = end_ms,
        duration_sec = round((end_ms - start_ms) / 1000, 3),
        workers = use_cores,
        total_tasks = length(task_metadata),
        successful_tasks = sum(sapply(res_list_flat, function(x) x$status == "Success"))
      )
    ),
    backend_info = gsea_env$backend_info,
    contrast_registry = gsea_env$contrast_registry,
    de_store = gsea_env$de_store,
    expr_bundle = gsea_env$expr_bundle,
    geneset_info = gsea_env$geneset,
    results = res_list_flat
  )

  saveRDS(final_obj, rds_path)

  success_count <- final_obj$metadata$gsea_benchmark$successful_tasks

  message(sprintf("\n✅ 计算完成！耗时 %.2f 秒",
                  final_obj$metadata$gsea_benchmark$duration_sec))
  message(sprintf("   成功分析: %d/%d 个对比", success_count, total_tasks))
  message(sprintf("   结果已保存至: %s", rds_path))

  return(final_obj)
}



#' @title 处理单个 GSEA Chunk（独立函数，避免闭包捕获大对象）
#' @description 这个函数独立于主函数，通过参数接收所有需要的数据，
#'   避免 future 自动捕获外部环境的巨大对象。
#' @param chunk_task_names 字符向量，当前 chunk 的任务名称
#' @param task_metadata 任务元数据列表
#' @param de_list DE 数据的 list 形式
#' @param term2gene 基因集映射数据框
#' @param meta_dict 元数据字典
#' @param minGSSize 最小基因集大小
#' @param maxGSSize 最大基因集大小
#' @param pvalueCutoff P值阈值
#' @return 当前 chunk 的结果列表
#' @keywords internal

.process_gsea_chunk <- function(chunk_task_names,
                                task_metadata,
                                de_list,
                                term2gene,
                                meta_dict,
                                minGSSize,
                                maxGSSize,
                                pvalueCutoff) {

  # 在 worker 进程内部加载必要的包（确保 worker 有这些包）
  if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE)) {
    stop("Worker 缺少必要的包: clusterProfiler 或 dplyr")
  }

  chunk_results <- list()

  for (task_name in chunk_task_names) {
    task_info <- task_metadata[[task_name]]
    task_id <- task_info$task_id

    # 获取原始对比ID（用于查询 de_list）
    if (task_info$is_reversed) {
      original_cid <- paste(task_info$right_group, task_info$left_group, sep = "_vs_")
    } else {
      original_cid <- task_id
    }

    # 查询DE表
    de_table <- de_list[[original_cid]]

    if (is.null(de_table) || nrow(de_table) == 0) {
      chunk_results[[task_name]] <- list(
        name = task_id,
        status = "Failed",
        data = NULL,
        genelist = c()
      )
      next
    }

    # 准备排序向量
    genelist <- tryCatch({
      .prepare_rank_vector_fast(de_table, flip = task_info$is_reversed)
    }, error = function(e) {
      warning(sprintf("排序向量准备失败: %s", e$message))
      c()
    })

    if (length(genelist) == 0) {
      chunk_results[[task_name]] <- list(
        name = task_id,
        status = "Failed",
        data = NULL,
        genelist = c()
      )
      next
    }

    # 执行 GSEA 计算
    set.seed(123)
    gsea_res <- tryCatch({
      clusterProfiler::GSEA(
        geneList = genelist,
        TERM2GENE = term2gene,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        pvalueCutoff = pvalueCutoff,
        pAdjustMethod = "BH",
        verbose = FALSE,
        seed = 123,
        eps = 0
      )
    }, error = function(e) {
      warning(sprintf("GSEA 计算失败 [%s]: %s", task_id, e$message))
      NULL
    })

    status <- if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) {
      "Success"
    } else {
      "Failed/NoEnrich"
    }

    # 注入元数据
    if (status == "Success" && !is.null(meta_dict)) {
      gsea_res@result <- .enrich_gsea_result(gsea_res@result, meta_dict)
    }

    chunk_results[[task_name]] <- list(
      name = task_id,
      status = status,
      data = gsea_res,
      genelist = genelist
    )
  }

  return(chunk_results)
}


#' @title 快速准备排序向量（优化版）
#' @description 使用 tidyverse 管道优化排序向量生成
#' @keywords internal

.prepare_rank_vector_fast <- function(de_table, flip = FALSE) {

  # 单管道操作，减少中间对象
  vals <- de_table %>%
    dplyr::filter(!is.na(gene_symbol), gene_symbol != "") %>%
    dplyr::mutate(
      gene_symbol = toupper(gene_symbol),
      abs_stat = abs(stat)
    ) %>%
    dplyr::arrange(dplyr::desc(abs_stat)) %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE) %>%
    {
      # 提取值和名称
      vec <- .$stat
      names(vec) <- .$gene_symbol
      vec
    }

  # 条件翻转
  if (flip) vals <- -vals

  # 排序
  sort(vals, decreasing = TRUE)
}

#' @title GSEA 结果元数据注入（优化版 - 修复 Description.y 错误）
#' @description 将 GSEA 结果与元数据字典合并，添加显示名称、链接等信息。
#'   修复了原版本中假设 Description.y 列存在导致的错误。
#' @param result_df GSEA 结果数据框（来自 clusterProfiler::GSEA 的 @result）
#' @param meta_dict 元数据字典，包含 ID, Collection, Combo_Name, URL, Description 等列
#' @return  enriched 的数据框
#' @keywords internal

.enrich_gsea_result_old <- function(result_df, meta_dict) {

  # 防御性检查：确保 meta_dict 是 data.frame/tibble
  if (is.null(meta_dict) || nrow(meta_dict) == 0) {
    warning("meta_dict 为空，返回原始结果")
    return(result_df)
  }
  # 保留原始 ID 用于后续恢复行名
  original_ids <- result_df$ID
  # 将 meta_dict 转换为 tibble 并检查列
  meta_tibble <- dplyr::as_tibble(meta_dict)

  # 检查 meta_dict 是否包含关键列
  required_cols <- c("ID")
  missing_required <- setdiff(required_cols, colnames(meta_tibble))
  if (length(missing_required) > 0) {
    stop(sprintf("meta_dict 缺少必需列: %s", paste(missing_required, collapse = ", ")))
  }

  # 检查可选列是否存在（用于后续逻辑）
  has_description_in_meta <- "Description" %in% colnames(meta_tibble)
  has_collection <- "Collection" %in% colnames(meta_tibble)
  has_combo_name <- "Combo_Name" %in% colnames(meta_tibble)
  has_url <- "URL" %in% colnames(meta_tibble)

  # 第一步：安全地处理 result_df 的 Description 列
  # 如果 result_df 有 Description，先备份到临时列
  result_processed <- result_df %>%
    dplyr::as_tibble()

  # 如果原始 result 有 Description，重命名为避免冲突
  if ("Description" %in% colnames(result_processed)) {
    result_processed <- result_processed %>%
      dplyr::rename(Original_Description = Description)
  } else {
    # 创建空的 Original_Description 列
    result_processed <- result_processed %>%
      dplyr::mutate(Original_Description = NA_character_)
  }

  # 第二步：执行 left join（meta_dict 中的 Description 会正常加入，不会重命名）
  joined_result <- result_processed %>%
    dplyr::left_join(
      meta_tibble,
      by = "ID",
      suffix = c("", "_meta")  # 如果 meta 中有同名列，加 _meta 后缀
    )

  # 第三步：安全地构建最终列（使用 coalesce 处理可能不存在的列）
  final_result <- joined_result %>%
    dplyr::mutate(
      # Display_Collection: 优先使用 Combo_Name，其次 Collection，最后 Unknown
      Display_Collection = dplyr::case_when(
        has_combo_name & !is.na(Combo_Name) ~ Combo_Name,
        has_collection & !is.na(Collection) ~ Collection,
        TRUE ~ "Unknown"
      ) %>% forcats::as_factor(),

      # Pathway_Link: 根据 URL 生成链接
      Pathway_Link = dplyr::case_when(
        !has_url | is.na(URL) | URL == "" ~ sprintf("<b>%s</b>", ID),
        TRUE ~ sprintf('<a href="%s" target="_blank">%s</a>', URL, ID)
      ),

      # Description: 优先级：meta_dict.Description > Original_Description > ID
      Description = dplyr::coalesce(
        if (has_description_in_meta) joined_result$Description else NULL,
        joined_result$Original_Description,
        joined_result$ID
      )
    ) %>%
    # 第四步：清理临时列
    dplyr::select(
      -dplyr::any_of(c("Original_Description")),  # 删除临时备份列
      #-dplyr::any_of(c("URL")),  # 删除 URL（已转为链接）
      -dplyr::ends_with("_meta")  # 删除可能的 _meta 后缀列（如果存在）
    )

  # 如果还有 Description_meta 列（当 meta_dict 有 Description 且 result 也有时），删除它
  if ("Description_meta" %in% colnames(final_result)) {
    final_result <- final_result %>%
      dplyr::select(-Description_meta)
  }

  # 🔧 关键修复：恢复行名
  final_result <- as.data.frame(final_result)
  if ("ID" %in% colnames(final_result)) {
    rownames(final_result) <- final_result$ID
  } else if (length(original_ids) == nrow(final_result)) {
    rownames(final_result) <- original_ids
  }

  return(final_result)
}



#' @title GSEA 结果元数据注入（安全版 - 保留行名）
#' @description 将 GSEA 结果与元数据字典合并，添加显示名称、链接等信息。
#'   此版本专门修复了行名丢失问题，确保 gseaResult 对象结构完整性。
#' @param result_df GSEA 结果数据框（来自 clusterProfiler::GSEA 的 @result）
#' @param meta_dict 元数据字典，包含 ID, Collection, Combo_Name, URL, Description 等列
#' @return enriched 的数据框，保持原始行名和行数
#' @keywords internal

.enrich_gsea_result <- function(result_df, meta_dict) {

  # --- 1. 输入验证与原始状态保存 ---

  if (!is.data.frame(result_df) || nrow(result_df) == 0) {
    warning("result_df 无效或为空，返回原样")
    return(result_df)
  }

  if (!is.data.frame(meta_dict) || nrow(meta_dict) == 0) {
    warning("meta_dict 无效或为空，返回原始结果")
    return(result_df)
  }

  if (!"ID" %in% colnames(result_df)) {
    stop("result_df 必须包含 'ID' 列用于匹配元数据")
  }

  if (!"ID" %in% colnames(meta_dict)) {
    stop("meta_dict 必须包含 'ID' 列用于匹配")
  }

  # 关键修复：立即保存原始行名，防止后续任何操作丢失
  original_rownames <- rownames(result_df)

  # 保存原始 Description（如果存在），用于后续合并
  has_original_desc <- "Description" %in% colnames(result_df)
  original_description <- if (has_original_desc) result_df$Description else NULL

  # --- 2. 准备元数据子集 ---

  # 确定 meta_dict 中可用的列，避免 select 时出错
  available_meta_cols <- colnames(meta_dict)
  cols_to_select <- c("ID")

  # 按优先级添加可选列
  optional_cols <- c("Combo_Name", "Collection", "URL", "Description")
  present_optional <- intersect(optional_cols, available_meta_cols)
  cols_to_select <- c(cols_to_select, present_optional)

  # 安全地选择列
  meta_subset <- meta_dict[, cols_to_select, drop = FALSE]

  # 检查重复 ID（会导致 join 后行数增加，破坏行名对应关系）
  if (any(duplicated(meta_subset$ID))) {
    warning("meta_dict 中存在重复的 ID，将保留首次出现项")
    meta_subset <- meta_subset[!duplicated(meta_subset$ID), ]
  }

  # --- 3. 执行数据合并（保持行顺序和行数） ---

  # 创建工作副本，移除原始的 Description 列（如果存在），避免与 meta_dict 冲突
  result_work <- result_df
  if (has_original_desc) {
    result_work$Description <- NULL
  }

  # 使用 base::merge 替代 dplyr::left_join 以更好控制行名行为
  # 或者使用 dplyr 但确保后续恢复行名
  merged_df <- dplyr::left_join(
    dplyr::as_tibble(result_work),  # 临时转为 tibble 进行 tidyverse 操作
    dplyr::as_tibble(meta_subset),
    by = "ID",
    suffix = c("", "_meta")  # 如果冲突，meta_dict 的列加 _meta 后缀
  )

  # 转回 data.frame 以便后续处理
  merged_df <- as.data.frame(merged_df)

  # --- 4. 构建展示列 ---

  # 4.1 Display_Collection: 优先级 Combo_Name > Collection > "Unknown"
  if ("Combo_Name" %in% colnames(merged_df)) {
    merged_df$Display_Collection <- dplyr::coalesce(merged_df$Combo_Name, merged_df$Collection)
  } else if ("Collection" %in% colnames(merged_df)) {
    merged_df$Display_Collection <- merged_df$Collection
  } else {
    merged_df$Display_Collection <- "Unknown"
  }

  # 处理 NA 并转换为因子
  na_collection <- is.na(merged_df$Display_Collection) | merged_df$Display_Collection == ""
  merged_df$Display_Collection[na_collection] <- "Unknown"
  merged_df$Display_Collection <- as.factor(merged_df$Display_Collection)

  # 4.2 Pathway_Link: 基于 URL 创建 HTML 链接
  if ("URL" %in% colnames(merged_df)) {
    merged_df$Pathway_Link <- ifelse(
      is.na(merged_df$URL) | merged_df$URL == "",
      sprintf("<b>%s</b>", merged_df$ID),
      sprintf('<a href="%s" target="_blank">%s</a>', merged_df$URL, merged_df$ID)
    )
  } else {
    merged_df$Pathway_Link <- sprintf("<b>%s</b>", merged_df$ID)
  }

  # 4.3 Description: 优先级 meta_dict.Description > original.Description > ID
  # 处理可能的 Description_meta 列（来自 suffix 冲突）
  if ("Description_meta" %in% colnames(merged_df)) {
    # 如果存在 Description_meta，说明 result_df 和 meta_dict 都有 Description
    meta_desc <- merged_df$Description_meta
    merged_df$Description_meta <- NULL  # 删除临时列
  } else if ("Description" %in% colnames(merged_df)) {
    # 只有 meta_dict 有 Description（result_work 中的已被删除）
    meta_desc <- merged_df$Description
  } else {
    meta_desc <- NULL
  }

  # 合并 Description 逻辑
  if (!is.null(meta_desc)) {
    # 有 meta_dict 的描述
    if (!is.null(original_description)) {
      merged_df$Description <- dplyr::coalesce(meta_desc, original_description)
    } else {
      merged_df$Description <- meta_desc
    }
  } else {
    # 没有 meta_dict 的描述，使用原始的
    if (!is.null(original_description)) {
      merged_df$Description <- original_description
    } else {
      merged_df$Description <- merged_df$ID
    }
  }

  # 最终后备：任何 NA 或空值都用 ID 填充
  na_desc_idx <- is.na(merged_df$Description) | merged_df$Description == ""
  merged_df$Description[na_desc_idx] <- merged_df$ID[na_desc_idx]

  # --- 5. 清理临时/中间列 ---

  # 删除从 meta_dict 引入的原始元数据列（保留我们构建的新列）
  cols_to_remove <- intersect(c("Combo_Name", "Collection", "URL"), colnames(merged_df))
  if (length(cols_to_remove) > 0) {
    merged_df <- merged_df[, !(colnames(merged_df) %in% cols_to_remove), drop = FALSE]
  }

  # --- 6. 关键修复：恢复原始行名 ---

  if (length(original_rownames) == nrow(merged_df)) {
    rownames(merged_df) <- original_rownames
  } else {
    # 如果行数不匹配（理论上不应发生，因为 left_join 保持左表行数）
    warning(sprintf(
      "行名长度不匹配: 原始 %d 行，处理后 %d 行。尝试使用 ID 列作为行名。",
      length(original_rownames),
      nrow(merged_df)
    ))
    # 后备方案：使用 ID 列
    if ("ID" %in% colnames(merged_df)) {
      rownames(merged_df) <- merged_df$ID
    }
  }

  return(merged_df)
}
