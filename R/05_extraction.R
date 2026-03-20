#' @title 提取特定对比组的 GSEA 结果
#' @description 从 GseaRes 胶囊中提取指定对比组的结果。
#' 支持按基因集子集 进行切片，并自动重新计算 FDR。
#' @param gsea_res GseaRes 对象 (由 batch_calc_gsea 生成)
#' @param contrast_id 字符串。对比组 ID (如 "A_vs_B")。
#' @param target_collection 字符串向量。欲提取的基因集集合 (如 "H", "C2:CP:KEGG_LEGACY")。
#' 默认 "ALL" 表示提取全部。
#' @return GseaTask 对象。
#' @export
extract_gsea_task <- function(gsea_res, contrast_id, target_collection = "ALL") {

  # 1. 校验输入
  .check_gsea_res(gsea_res)

  if (!(contrast_id %in% names(gsea_res$results))) {
    stop(sprintf("❌ 对比组 '%s' 不存在于结果中。可用: %s",
                 contrast_id, paste(names(gsea_res$results), collapse = ", ")))
  }

  task_info <- gsea_res$results[[contrast_id]]

  if (task_info$status != "Success") {
    stop(sprintf("❌ 对比组 '%s' 计算失败或无富集结果，无法提取。", contrast_id))
  }

  # 2. 获取元信息 (从 registry 中精确匹配，不再解析字符串)
  # 注意：contrast_id 可能是反向的 (B_vs_A)，registry 中可能只有正向的 (A_vs_B)
  # 但我们在 batch_calc_gsea 中已经把反向任务写入了 results，却没有写入 registry (这是设计上的权衡，registry 记录原始生物学对比)
  # 这里我们需要智能处理 left/right group

  # 优先从 results 里的 genelist 名称推断 (最准确)
  # 或者从 registry 查找

  reg <- gsea_res$contrast_registry
  reg_match <- reg[reg$contrast_id == contrast_id, ]

  left_group <- NA
  right_group <- NA

  if (nrow(reg_match) > 0) {
    left_group <- reg_match$left_group[1]
    right_group <- reg_match$right_group[1]
  } else {
    # 如果 registry 里没找到 (可能是反向生成的任务)，尝试从 ID 解析
    parts <- strsplit(contrast_id, "_vs_")[[1]]
    if (length(parts) == 2) {
      left_group <- parts[1]
      right_group <- parts[2]
    } else {
      left_group <- contrast_id
      right_group <- "Background"
    }
  }

  # 3. 提取 GSEA 结果对象
  gsea_obj <- task_info$data
  res_df <- as.data.frame(gsea_obj@result)

  # 4. 切片逻辑
  is_slice_mode <- length(target_collection) != 1 || toupper(target_collection[1]) != "ALL"

  if (is_slice_mode) {
    # 查找匹配的行
    # 支持 Collection (如 "H") 或 Combo_Name (如 "C2:CP:KEGG_LEGACY")
    match_idx <- which(res_df$Collection %in% target_collection |
                         res_df$Combo_Name %in% target_collection)

    if (length(match_idx) == 0) {
      stop(sprintf("❌ 切片失败！未在 '%s' 中找到属于 '%s' 的通路。",
                   contrast_id, paste(target_collection, collapse = ", ")))
    }

    res_df <- res_df[match_idx, , drop = FALSE]

    # 💥 关键：重算 FDR (消除全库背景的惩罚)
    res_df$p.adjust <- p.adjust(res_df$pvalue, method = "BH")
    if ("qvalue" %in% colnames(res_df)) res_df$qvalue <- res_df$p.adjust

    message(sprintf("✂️ [智能切片] 抽出 %d 条通路 | 已重算 FDR。", nrow(res_df)))

    # 更新 S4 对象
    gsea_obj@result <- res_df
  }

  # 5. 构建 GseaTask 对象
  # 整合 meta 信息，供下游绘图使用
  task_meta <- list(
    contrast_id = contrast_id,
    left_group = left_group,
    right_group = right_group,
    geneset_name = ifelse(is_slice_mode,
                          paste(target_collection, collapse = "_"),
                          gsea_res$geneset_info$name),
    meta_dict = gsea_res$geneset_info$meta_dict,
    expr_bundle = gsea_res$expr_bundle,
    project_info = gsea_res$metadata$project_info,
    backend_info = gsea_res$backend_info
  )

  obj <- create_gsea_task(gsea_res = gsea_obj, meta = task_meta)

  return(obj)
}


# 查看函数


#' @title 查看 GSEA 结果概览
#' @description 在控制台打印 GseaRes 对象的详细摘要。
#' @param gsea_res GseaRes 对象
#' @export
inspect_gsea_res <- function(gsea_res) {
  if (!inherits(gsea_res, "GseaRes")) stop("输入对象不是 GseaRes 类。")

  bi <- gsea_res$backend_info
  results <- gsea_res$results

  cat("\n", rep("=", 60), "\n", sep = "")
  cat("       📊 GSEAlens Result Summary\n")
  cat(rep("=", 60), "\n\n", sep = "")

  # 1. 后端信息
  cat("⚙️  [1] Backend\n")
  cat(sprintf("   • Type       : %s\n", bi$backend))
  cat(sprintf("   • Run Time   : %s\n", gsea_res$metadata$run_time))
  cat("\n")

  # 2. 计算状态统计
  cat("📉 [2] Calculation Status\n")
  status_count <- table(sapply(results, function(x) x$status))
  for (name in names(status_count)) {
    icon <- if (name == "Success") "✅" else "❌"
    cat(sprintf("   %s %s : %d\n", icon, name, status_count[[name]]))
  }
  cat("\n")

  # 3. 详细列表 (前5个)
  cat("⚖️  [3] Contrast Details (Top 5)\n")
  n_show <- min(5, length(results))
  for (i in 1:n_show) {
    task_name <- names(results)[i]
    task <- results[[i]]

    if (task$status == "Success") {
      n_sig <- sum(task$data@result$p.adjust < 0.05, na.rm = TRUE)
      cat(sprintf("   [%d] %-20s | FDR<0.05: %d\n", i, task_name, n_sig))
    } else {
      cat(sprintf("   [%d] %-20s | Failed\n", i, task_name))
    }
  }
  cat("\n")

  # 4. 下一步指引
  cat("🚀 [4] Next Step\n")
  cat("   提取特定结果进行可视化：\n")
  cat(sprintf('   > task <- extract_gsea_task(gsea_res, "%s")\n', names(results)[1]))
  cat('   > plot_directional_gsea(task, ...)\n')
  cat(rep("=", 60), "\n\n", sep = "")

  invisible(gsea_res)
}

#' @export
print.GseaRes <- function(x, ...) {
  inspect_gsea_res(x)
}
