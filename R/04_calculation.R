#' @title 批量并行 GSEA 计算
#' @description 消费标准 GseaEnv 对象，执行并行 GSEA 计算。
#' 自动从 de_store 中提取统计量进行排序，支持双向对比生成。
#' @param gsea_env GseaEnv 对象 (由 setup_gsea_env 生成)
#' @param custom_series_name 字符串。分析系列名称，用于结果存档。
#' @param output_dir 字符串。结果输出路径，默认 "./GSEA_Output"。
#' @param workers 并行核心数。默认 4。
#' @param bidirectional 逻辑值。是否自动生成反向对比，默认 TRUE。
#' @param minGSSize 最小基因集大小，默认 10。
#' @param maxGSSize 最大基因集大小，默认 500。
#' @param pvalueCutoff P 值阈值，默认 1 (保留全量结果)。
#' @param force 逻辑值。是否强制重新计算，默认 FALSE。
#' @return GseaRes 对象。
#' @export
batch_calc_gsea <- function(gsea_env,
                            custom_series_name = "Auto_Analysis",
                            output_dir = "./GSEA_Output",
                            workers = 4,
                            bidirectional = TRUE,
                            minGSSize = 10,
                            maxGSSize = 500,
                            pvalueCutoff = 1,
                            force = FALSE) {

  # 1. 校验输入
  .check_gsea_env(gsea_env)

  # 2. 智能缓存拦截
  series_dir <- file.path(output_dir, custom_series_name)
  if (!dir.exists(series_dir)) dir.create(series_dir, recursive = TRUE, showWarnings = FALSE)

  rds_name <- sprintf("GSEA_Capsule_[%s]_[%s].rds", custom_series_name, gsea_env$geneset$name)
  rds_path <- file.path(series_dir, rds_name)

  if (file.exists(rds_path) && !force) {
    message(sprintf("✅ 命中缓存！检测到已存在的 GSEA 胶囊: %s", rds_name))
    return(readRDS(rds_path))
  }

  # 3. 准备任务列表
  # 从 contrast_registry 获取任务
  registry <- gsea_env$contrast_registry
  de_store <- gsea_env$de_store

  # 构建任务映射表
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

  # 4. 并行计算设置
  total_cores <- parallel::detectCores(logical = TRUE)
  use_cores <- min(total_cores, max(1, workers))
  message(sprintf("🖥️ 硬件调度: 使用 %d 核进行并行计算...", use_cores))

  options(future.globals.maxSize = 32000 * 1024^2)
  future::plan(future::multisession, workers = use_cores)

  # 5. 核心计算循环
  term2gene <- gsea_env$geneset$term2gene
  meta_dict <- gsea_env$geneset$meta_dict

  res_list <- future.apply::future_lapply(names(tasks), function(task_name) {

    task_info <- tasks[[task_name]]
    cid <- task_info$contrast_id

    # 获取 DE 表
    # 如果是反向任务，需要找到原始的正向 ID 来获取 DE 表
    original_cid <- if (task_info$reverse) {
      paste(task_info$left, task_info$right, sep = "_vs_") # 反向任务的 original 是反向的
      # 注意：这里逻辑需要修正。
      # 如果任务是 B_vs_A (reverse=TRUE)，那么原始数据是 A_vs_B。
      # task_info$left 是 B, right 是 A。
      # 所以原始 ID 应该是 A_vs_B。
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
    # 使用统一的 .prepare_rank_vector 函数
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
          Description = if("Description.y" %in% names(.)) Description.y else ID # 优先使用字典里的描述
        )

      rownames(res_df) <- res_df$ID
      gsea_res@result <- res_df
    }

    return(list(name = cid, status = status, data = gsea_res, genelist = genelist))

  }, future.seed = TRUE)

  names(res_list) <- names(tasks)
  future::plan(future::sequential)

  # 6. 封装结果对象
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
      )
    ),
    backend_info = gsea_env$backend_info,
    contrast_registry = gsea_env$contrast_registry,
    de_store = gsea_env$de_store,
    expr_bundle = gsea_env$expr_bundle,
    geneset_info = gsea_env$geneset,
    results = res_list
  )

  # 7. 保存结果
  saveRDS(final_obj, rds_path)
  message(sprintf("\n✅ 计算完成！结果已保存至: %s", rds_path))

  return(final_obj)
}


# 内部辅助函数


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
