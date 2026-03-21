#' @title 加载并智能归位 GSEA 计算胶囊
#' @description 安全载入计算胶囊。若文件脱离原始项目路径（比如被拷贝到了桌面），
#' 引擎将自动在当前工作目录复原标准文件夹，并将胶囊护送归位。
#' 完全兼容 GseaEnv 和 GseaRes 对象。
#' @param file_path 字符型。胶囊 rds 文件的绝对或相对路径。
#' @param auto_relocate 逻辑值。是否自动根据内置血统在当前目录下修复文件夹并归位？默认 TRUE。
#' @param inspect 逻辑值。加载后是否自动打印胶囊概览？默认 TRUE。
#' @return 返回解析后的 GseaRes 或 GseaEnv 胶囊对象。
#' @export
import_gsea_capsule <- function(file_path, auto_relocate = TRUE, inspect = TRUE) {

  # 1. 文件存在性检查
  if (!file.exists(file_path)) {
    stop(sprintf("❌ 文件不存在: %s", file_path))
  }

  message("🎁 正在唤醒计算胶囊...")
  capsule <- readRDS(file_path)

  # 2. 识别胶囊类型
  if (inherits(capsule, c("GseaEnv", "GseaEnvPro"))) {
    message("✅ 成功载入 [GseaEnv] 环境封装胶囊 (尚未进行并行计算)。")

    if (inspect && exists("inspect_gsea_env", mode = "function")) {
      inspect_gsea_env(capsule)
    }

    return(invisible(capsule))
  }

  # 3. 对于计算结果胶囊，执行智能归位逻辑
  if (!inherits(capsule, "GseaRes")) {
    warning("⚠️ 该文件似乎不是标准的 GSEAlens 胶囊！")
    return(invisible(capsule))
  }

  # 4. 检查是否有血统记忆（project_info）
  project_info <- capsule$metadata$project_info

  if (is.null(project_info)) {
    message("⚠️ 该胶囊为旧版本生成，缺乏项目血统记忆。直接载入（无法自动归位）。")
  } else {
    # 5. 获取当前路径与预期路径
    current_abs <- normalizePath(file_path, winslash = "/", mustWork = FALSE)
    expected_abs <- normalizePath(project_info$rds_path, winslash = "/", mustWork = FALSE)

    # 6. 路径不一致，且开启了自动归位
    if (current_abs != expected_abs && auto_relocate) {
      message(sprintf("🚨 [血统警报] 胶囊当前处于非标准路径"))
      message(sprintf("   📍 当前位置: %s", current_abs))
      message(sprintf("   🏠 原籍项目: [%s]", project_info$custom_series_name))

      # 7. 根据当前工作目录重建标准档案库
      local_series_dir <- file.path(getwd(), "GSEA_Output", project_info$custom_series_name)
      if (!dir.exists(local_series_dir)) {
        dir.create(local_series_dir, recursive = TRUE, showWarnings = FALSE)
        message(sprintf("   ✨ 已创建标准档案库目录: %s", local_series_dir))
      }

      # 8. 执行文件归位
      target_file <- file.path(local_series_dir, basename(file_path))

      if (!file.exists(target_file) || target_file == current_abs) {
        file.copy(from = file_path, to = target_file, overwrite = TRUE)
        message(sprintf("   🚀 已自动将胶囊遣返归位至标准档案库"))
        message(sprintf("   📦 新位置: %s", target_file))

        # 9. 更新胶囊内部的路径信息，防止下次加载时重复警告
        capsule$metadata$project_info$output_dir <- file.path(getwd(), "GSEA_Output")
        capsule$metadata$project_info$series_dir <- local_series_dir
        capsule$metadata$project_info$rds_path <- target_file

      } else {
        message("   ✅ 标准档案库中已有该文件备份（无需重复复制）。")
      }
    }
  }

  # 10. 自动调用探针函数打印概况
  if (inspect) {
    if (exists("inspect_gsea_res", mode = "function")) {
      inspect_gsea_res(capsule)
    } else {
      message("✅ 成功载入 [GseaRes] 计算完成结果胶囊！")
    }
  }

  return(invisible(capsule))
}



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
