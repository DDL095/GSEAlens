#' @title Batch Parallel GSEA Calculation (Optimized Version)
#' @description Consumes a standard GseaEnv object and performs efficient parallel GSEA calculation.
#' @param gsea_env GseaEnv object
#' @param custom_series_name String. Analysis series name
#' @param output_dir String. Output directory for results, default "./GSEA_Output"
#' @param workers Number of parallel cores. Default 4
#' @param bidirectional Logical. Whether to automatically generate reverse contrasts, default TRUE
#' @param minGSSize Minimum gene set size, default 10
#' @param maxGSSize Maximum gene set size, default 500
#' @param pvalueCutoff P-value threshold, default 1
#' @param force Logical. Whether to force recalculation, default FALSE
#' @param use_progress Logical. Whether to show progress bar, default TRUE
#' @param chunk_size Integer. Number of tasks per worker per chunk, default NULL (auto)
#' @return GseaRes object
#' @export
#' @title Batch Parallel GSEA Calculation (Optimized Version - Fixed Global Variable Size Error)
#' @description Consumes a standard GseaEnv object and performs efficient parallel GSEA calculation.
#' @export
#' @title Batch Parallel GSEA Calculation (With Real-time Progress Bar)
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
                            use_progress = FALSE,  # 控制是否显示进度
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
    message(sprintf("Cache hit! Existing GSEA capsule detected: %s", rds_name))
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
  message(sprintf("Ready: %d contrast tasks pending calculation...", total_tasks))


  # 第三步：并行设置（关键修复：progressr 处理器）


  # 设置全局变量大小限制（必须在 plan 之前）
  options(future.globals.maxSize = 192 * 1024^3)#设置无限大

  total_cores <- parallel::detectCores(logical = TRUE)
  use_cores <- min(total_cores, max(1, workers))
  message(sprintf("Hardware scheduling: Using %d cores for parallel computation...", use_cores))

  if (is.null(chunk_size)) {
    chunk_size <- max(1, ceiling(total_tasks / (use_cores * 4)))
  }

  # 设置 future 计划
  future::plan(future::multisession, workers = use_cores)

  # 关键修复 1：启用 progressr 处理器（必须在 future_lapply 之前）
  if (use_progress) {
    if (!requireNamespace("progressr", quietly = TRUE)) {
      stop("Please install progressr package to show progress bar: install.packages('progressr')")
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

  message(sprintf("Chunk strategy: %d chunks, up to %d tasks per chunk",
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
        stop("Worker missing required packages")
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

  message(sprintf("\nCalculation complete! Time elapsed: %.2f seconds",
                  final_obj$metadata$gsea_benchmark$duration_sec))
  message(sprintf("   Successfully analyzed: %d/%d contrasts", success_count, total_tasks))
  message(sprintf("   Results saved to: %s", rds_path))

  return(final_obj)
}



#' @title Process Single GSEA Chunk (Independent Function to Avoid Closure Capturing Large Objects)
#' @description This function is independent from the main function and receives all required data
#'   via parameters, avoiding future automatically capturing huge objects from the external environment.
#' @param chunk_task_names Character vector, task names for current chunk
#' @param task_metadata Task metadata list
#' @param de_list DE data in list form
#' @param term2gene Gene set mapping dataframe
#' @param meta_dict Metadata dictionary
#' @param minGSSize Minimum gene set size
#' @param maxGSSize Maximum gene set size
#' @param pvalueCutoff P-value threshold
#' @return Result list for current chunk
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
    stop("Worker missing required packages: clusterProfiler or dplyr")
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
      warning(sprintf("Rank vector preparation failed: %s", e$message))
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
      warning(sprintf("GSEA calculation failed [%s]: %s", task_id, e$message))
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


#' @title Fast Rank Vector Preparation (Optimized Version)
#' @description Uses tidyverse pipes to optimize rank vector generation
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


#' @title GSEA Result Metadata Injection (Safe Version - Fixed for DDS Backend)
#' @description Merges GSEA results with metadata dictionary ensuring ALL columns are preserved.
#'   Specifically fixes the missing URL/Collection/Subcollection issue in DESeq2 backend.
#' @param result_df GSEA result dataframe from clusterProfiler::GSEA
#' @param meta_dict Metadata dictionary containing ID, Collection, Combo_Name, URL, Description, etc.
#' @return Enriched dataframe with 17 standard columns
#' @keywords internal

.enrich_gsea_result <- function(result_df, meta_dict) {

  # --- 0. 防御性检查 ---
  if (!is.data.frame(result_df) || nrow(result_df) == 0) {
    warning("result_df is invalid or empty")
    return(result_df)
  }

  if (!is.data.frame(meta_dict) || nrow(meta_dict) == 0) {
    warning("meta_dict is empty, cannot enrich results")
    return(result_df)
  }

  if (!"ID" %in% colnames(result_df) || !"ID" %in% colnames(meta_dict)) {
    stop("Both result_df and meta_dict must contain 'ID' column")
  }

  # 保存原始行名和ID
  original_rownames <- rownames(result_df)
  original_ids <- result_df$ID

  # --- 1. 标准化meta_dict列结构（确保所有必需列存在）---
  required_cols <- c("ID", "Description", "URL", "Collection", "Subcollection", "Combo_Name")

  # 如果meta_dict缺少某些列，创建空列占位
  for (col in required_cols) {
    if (!col %in% colnames(meta_dict)) {
      warning(sprintf("meta_dict missing column '%s', creating placeholder", col))
      meta_dict[[col]] <- NA_character_
    }
  }

  # 确保Subcollection不为NULL（DESeq2流程中可能出现）
  if (all(is.na(meta_dict$Subcollection)) || is.null(meta_dict$Subcollection)) {
    meta_dict$Subcollection <- ""
  }

  # 确保Combo_Name正确生成（如果缺失）
  if (all(is.na(meta_dict$Combo_Name))) {
    meta_dict$Combo_Name <- ifelse(
      is.na(meta_dict$Subcollection) | meta_dict$Subcollection == "",
      meta_dict$Collection,
      paste0(meta_dict$Collection, ":", meta_dict$Subcollection)
    )
  }

  # --- 2. 安全地准备result_df（移除可能与meta_dict冲突的列，但保留核心统计列）---
  # 定义GSEA核心列（统计结果，必须保留）
  core_stat_cols <- c("ID", "setSize", "enrichmentScore", "NES", "pvalue",
                      "p.adjust", "qvalue", "rank", "leading_edge", "core_enrichment")

  # 定义可能冲突的元数据列（这些应该从meta_dict获取，而非保留原值）
  conflict_cols <- c("Description", "URL", "Collection", "Subcollection", "Combo_Name")

  # 安全移除冲突列（如果存在）
  cols_to_keep <- setdiff(colnames(result_df), conflict_cols)
  result_work <- result_df[, cols_to_keep, drop = FALSE]

  # --- 3. 执行左连接（使用标准dplyr语法，确保列名不冲突）---
  # 选择meta_dict中需要的列，避免携带过多列
  meta_subset <- meta_dict[, required_cols, drop = FALSE]

  # 使用dplyr::left_join，自动处理后缀
  merged_df <- result_work %>%
    dplyr::left_join(
      meta_subset,
      by = "ID",
      suffix = c("", "_meta")  # 如果冲突，meta_dict的列加_meta后缀
    )

  # 检查合并结果
  if (nrow(merged_df) != nrow(result_work)) {
    warning(sprintf("Row count changed during merge: %d -> %d",
                    nrow(result_work), nrow(merged_df)))
  }

  # --- 4. 构建标准输出列（确保17列标准格式）---

  # 4.1 生成Display_Collection（从Combo_Name或Collection）
  merged_df$Display_Collection <- dplyr::coalesce(
    merged_df$Combo_Name,
    merged_df$Collection,
    "Unknown"
  )
  # 确保是因子类型（用于Shiny中的下拉筛选）
  merged_df$Display_Collection <- as.factor(merged_df$Display_Collection)

  # 4.2 生成Pathway_Link（HTML链接格式）
  merged_df$Pathway_Link <- ifelse(
    is.na(merged_df$URL) | merged_df$URL == "",
    sprintf("<b>%s</b>", merged_df$ID),
    sprintf('<a href="%s" target="_blank">%s</a>', merged_df$URL, merged_df$ID)
  )

  # 4.3 确保Description列存在（优先使用meta_dict的Description）
  if ("Description_meta" %in% colnames(merged_df)) {
    # 如果发生了列名冲突（Description.x/Description.y情况）
    # 使用meta_dict的Description（更完整的描述）
    merged_df$Description <- dplyr::coalesce(
      merged_df$Description_meta,
      merged_df$ID
    )
    # 删除临时列
    merged_df$Description_meta <- NULL
  } else if (!"Description" %in% colnames(merged_df)) {
    # 如果完全没有Description列，使用ID作为后备
    merged_df$Description <- merged_df$ID
  }

  # 处理NA值
  merged_df$Description[is.na(merged_df$Description)] <- merged_df$ID[is.na(merged_df$Description)]

  # --- 5. 清理临时列并恢复行名 ---

  # 删除所有_meta后缀的临时列（如果存在）
  meta_temp_cols <- grep("_meta$", colnames(merged_df), value = TRUE)
  if (length(meta_temp_cols) > 0) {
    merged_df <- merged_df[, !(colnames(merged_df) %in% meta_temp_cols), drop = FALSE]
  }

  # 确保行名正确（使用ID作为行名，这是clusterProfiler的标准）
  merged_df <- as.data.frame(merged_df)
  if ("ID" %in% colnames(merged_df)) {
    # 检查ID是否有重复
    if (any(duplicated(merged_df$ID))) {
      warning("Duplicate IDs found in result, using original rownames")
      rownames(merged_df) <- original_rownames
    } else {
      rownames(merged_df) <- merged_df$ID
    }
  } else {
    rownames(merged_df) <- original_rownames
  }

  # --- 6. 最终验证（确保17列标准）---
  standard_cols <- c("ID", "setSize", "enrichmentScore", "NES", "pvalue",
                     "p.adjust", "qvalue", "rank", "leading_edge", "core_enrichment",
                     "Description", "URL", "Collection", "Subcollection", "Combo_Name",
                     "Display_Collection", "Pathway_Link")

  missing_standard <- setdiff(standard_cols, colnames(merged_df))
  if (length(missing_standard) > 0) {
    warning(sprintf("Final result missing standard columns: %s",
                    paste(missing_standard, collapse = ", ")))
    # 为缺失的列创建空值
    for (col in missing_standard) {
      merged_df[[col]] <- NA_character_
    }
  }

  # 重新排序列以匹配标准顺序（可选，但有助于一致性）
  present_cols <- intersect(standard_cols, colnames(merged_df))
  extra_cols <- setdiff(colnames(merged_df), standard_cols)
  merged_df <- merged_df[, c(present_cols, extra_cols), drop = FALSE]

  message(sprintf("[.enrich_gsea_result] Successfully enriched: %d rows x %d columns",
                  nrow(merged_df), ncol(merged_df)))

  return(merged_df)
}
