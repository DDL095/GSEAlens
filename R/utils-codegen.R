#' @title GSEAlens 代码生成工具
#' @description 根据用户交互参数生成可复现的 R 代码，支持剪贴板复制和脚本导出。
#' @name utils-codegen
NULL

#' @title 生成 GSEA 绘图代码
#' @description 生成完整的、可直接运行的 R 脚本，用于复现当前的 GSEA 可视化。
#' @param contrast_id 字符串，对比组 ID（如 "IBAA_vs_PREA"）。
#' @param pathway_ids 字符向量，通路 ID 列表。
#' @param plot_subtype 整数，图像类型（1=仅富集图，2=富集+热图，3=完整带Rank）。
#' @param custom_colors 字符向量，自定义颜色代码。
#' @param layout_mode 字符串，排版模式（auto_grid/single_col/two_col/custom_grid）。
#' @param custom_nrow 整数，自定义行数（可选）。
#' @param custom_ncol 整数，自定义列数（可选）。
#' @param plot_width 数值，图像宽度（英寸）。
#' @param plot_height 数值，图像高度（英寸）。
#' @param plot_dpi 整数，图像分辨率。
#' @param rds_path 字符串，GseaRes 对象文件路径（用于代码中的注释）。
#' @return 字符串，格式化后的 R 代码。
#' @export
#' @examples
#' \dontrun{
#' code <- generate_plot_code(
#'   contrast_id = "IBAA_vs_PREA",
#'   pathway_ids = c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "KEGG_GLYCOLYSIS"),
#'   plot_subtype = 3,
#'   custom_colors = c("#E41A1C", "#377EB8")
#' )
#' cat(code)
#' }
generate_plot_code <- function(contrast_id,
                               pathway_ids,
                               plot_subtype = 3,
                               custom_colors = NULL,
                               layout_mode = "auto_grid",
                               custom_nrow = NULL,
                               custom_ncol = NULL,
                               plot_width = 12,
                               plot_height = 9,
                               plot_dpi = 300,
                               rds_path = "path/to/your/GSEA_Capsule.rds") {

  # 构建时间戳和元信息
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  # 处理颜色参数
  if (!is.null(custom_colors) && length(custom_colors) > 0) {
    # 确保颜色数量足够
    if (length(custom_colors) < length(pathway_ids)) {
      custom_colors <- rep(custom_colors, length.out = length(pathway_ids))
    }
    color_str <- paste0('c("', paste(custom_colors[seq_along(pathway_ids)], collapse = '", "'), '")')
  } else {
    color_str <- "NULL  # 使用默认颜色"
  }

  # 处理通路 ID 字符串
  if (length(pathway_ids) == 1) {
    pathway_str <- sprintf('"%s"', pathway_ids[1])
  } else {
    pathway_str <- paste0('c("', paste(pathway_ids, collapse = '", "'), '")')
  }

  # 构建布局代码片段
  layout_code <- ""
  if (length(pathway_ids) > 1) {
    layout_code <- sprintf("
# 计算排版布局 (%s 个通路)
layout_params <- switch('%s',
  'auto_grid' = {
    ncol <- ceiling(sqrt(%d))
    list(nrow = ceiling(%d / ncol), ncol = ncol)
  },
  'single_col' = list(nrow = %d, ncol = 1),
  'two_col' = list(nrow = ceiling(%d / 2), ncol = 2),
  'custom_grid' = list(nrow = %d, ncol = %d),
  list(nrow = 1, ncol = 1)
)

# 组合多图
p_combined <- patchwork::wrap_plots(
  plot_list,
  nrow = layout_params$nrow,
  ncol = layout_params$ncol,
  guides = 'collect'
) +
  patchwork::plot_annotation(
    title = 'GSEA Combined Plot: %s',
    subtitle = 'Contrast: %s',
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = 'bold'),
      plot.subtitle = ggplot2::element_text(size = 10, color = 'grey40')
    )
  )

p_final <- p_combined",
    length(pathway_ids),
    layout_mode,
    length(pathway_ids),
    length(pathway_ids),
    length(pathway_ids),
    length(pathway_ids),
    custom_nrow %||% 2,
    custom_ncol %||% 2,
    paste(pathway_ids, collapse = ", "),
    contrast_id
    )
  } else {
    layout_code <- "p_final <- p"
  }

  # 组合完整代码
  code <- sprintf('#!/usr/bin/env Rscript
# ============================================================================
# GSEAlens Auto-Generated Reproduction Code
# Generated: %s
# ============================================================================

# 1. 加载必要包 -----------------------------------------------------------
library(GSEAlens)
library(ggplot2)
library(patchwork)  # 用于组合多图

# 2. 加载 GSEA 结果胶囊 ----------------------------------------------------
# 请修改为你的实际文件路径
gsea_res <- readRDS("%s")

# 3. 提取特定对比组任务 ----------------------------------------------------
task <- extract_gsea_task(
  gsea_res = gsea_res,
  contrast_id = "%s",
  target_collection = "ALL"  # 如需切片，请修改此处
)

# 4. 定义绘图参数 ---------------------------------------------------------
target_pathways <- %s
custom_colors <- %s
plot_subtype <- %d  # 1=仅富集图, 2=富集+热图, 3=完整带Rank

# 5. 生成绘图 --------------------------------------------------------------
%s

# 单个通路绘图函数
create_gsea_plot <- function(pathway_id, color) {
  plot_directional_gsea(
    directional_gsea_obj = task,
    target_pathways = pathway_id,
    subPlot = plot_subtype,
    curveCol = color,
    main_title = pathway_id,
    add_pval = FALSE
  )
}

%s

# 6. 保存图像 --------------------------------------------------------------
output_file <- "GSEA_Combined_%s_%s.%s"

ggplot2::ggsave(
  filename = output_file,
  plot = p_final,
  width = %d,      # 英寸
  height = %d,     # 英寸
  dpi = %d,        # 分辨率
  bg = "white"     # 白色背景
)

message(sprintf("✅ 图像已保存至: %%s", output_file))
print(p_final)  # 在 R 会话中显示
',
timestamp,
rds_path,
contrast_id,
pathway_str,
color_str,
plot_subtype,
if (length(pathway_ids) > 1) sprintf('
# 生成单个图列表
plot_list <- lapply(seq_along(target_pathways), function(i) {
  create_gsea_plot(target_pathways[i], custom_colors[i])
})') else '
# 生成单个图
p <- create_gsea_plot(target_pathways[1], custom_colors[1])',
layout_code,
gsub("_vs_", "_", contrast_id),
format(Sys.time(), "%Y%m%d"),
"png",  # 默认格式，用户可自行修改
plot_width,
plot_height,
plot_dpi
  )

# 清理多余的空行
code <- gsub("\n{3,}", "\n\n", code)

return(code)
}

#' @title 生成 HTML 报告代码
#' @description 生成用于创建完整 HTML 报告的 R 代码。
#' @param contrast_id 字符串，对比组 ID。
#' @param p_adjust_cutoff 数值，FDR 过滤阈值。
#' @param top_plots 整数向量，c(正向数, 负向数)。
#' @param rds_path 字符串，胶囊文件路径。
#' @return 字符串，R 代码。
#' @export
generate_report_code <- function(contrast_id,
                                 p_adjust_cutoff = 1,
                                 top_plots = c(15, 15),
                                 rds_path = "path/to/your/GSEA_Capsule.rds") {

  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  sprintf('#!/usr/bin/env Rscript
# ============================================================================
# GSEAlens HTML Report Generation Code
# Generated: %s
# ============================================================================

library(GSEAlens)

# 加载数据
gsea_res <- readRDS("%s")

# 提取任务
task <- extract_gsea_task(
  gsea_res = gsea_res,
  contrast_id = "%s",
  target_collection = "ALL"
)

# 生成 HTML 报告
report_dir <- generate_gsea_html_report(
  res_obj = task,
  output_base_dir = NULL,  # 自动使用胶囊原地址
  p_adjust_cutoff = %f,
  top_plots = c(%d, %d),
  dpi = 200
)

message(sprintf("✅ 报告已生成至: %%s", report_dir))
',
timestamp,
rds_path,
contrast_id,
p_adjust_cutoff,
top_plots[1],
top_plots[2]
  )
}

#' @title 复制代码到剪贴板
#' @description 将生成的代码写入系统剪贴板，支持 Windows/Mac/Linux。
#' @param code 字符串，要复制的代码。
#' @param quiet 逻辑值，是否静默模式（不显示成功消息）。
#' @return 逻辑值，是否成功。
#' @export
#' @importFrom clipr write_clip
copy_code_to_clipboard <- function(code, quiet = FALSE) {
  if (!requireNamespace("clipr", quietly = TRUE)) {
    if (!quiet) {
      message("⚠️ 未安装 clipr 包，请手动安装: install.packages('clipr')")
      message("以下是代码内容，请手动复制：")
      message("========================================")
      cat(code)
      message("========================================")
    }
    return(FALSE)
  }

  tryCatch({
    clipr::write_clip(code)
    if (!quiet) {
      message("✅ 代码已复制到剪贴板！")
    }
    return(TRUE)
  }, error = function(e) {
    if (!quiet) {
      warning(sprintf("❌ 复制到剪贴板失败: %s", e$message))
      message("请手动复制上面的代码。")
    }
    return(FALSE)
  })
}

#' @title 保存代码到文件
#' @description 将生成的代码保存为 .R 脚本文件。
#' @param code 字符串，代码内容。
#' @param file_path 字符串，文件路径。如果为 NULL，使用交互式选择。
#' @return 字符串，保存的文件路径。
#' @export
save_code_to_file <- function(code, file_path = NULL) {
  if (is.null(file_path)) {
    # 生成默认文件名
    default_name <- sprintf("GSEAlens_code_%s.R", format(Sys.time(), "%Y%m%d_%H%M%S"))

    # 尝试使用交互式文件选择（如果在 Shiny 中可能不可用）
    if (interactive()) {
      file_path <- file.choose(new = TRUE)
      if (is.null(file_path) || file_path == "") {
        file_path <- default_name
      }
    } else {
      file_path <- default_name
    }
  }

  # 确保扩展名为 .R
  if (!grepl("\\.R$", file_path, ignore.case = TRUE)) {
    file_path <- paste0(file_path, ".R")
  }

  tryCatch({
    writeLines(code, con = file_path)
    message(sprintf("✅ 代码已保存至: %s", normalizePath(file_path)))
    return(normalizePath(file_path))
  }, error = function(e) {
    stop(sprintf("❌ 保存文件失败: %s", e$message))
  })
}

#' @title 生成 GseaRes 对象摘要代码
#' @description 生成用于查看和分析 GseaRes 对象的辅助代码。
#' @param rds_path 字符串，胶囊文件路径。
#' @return 字符串，R 代码。
#' @export
generate_summary_code <- function(rds_path = "path/to/your/GSEA_Capsule.rds") {
  sprintf('
# GSEAlens 数据对象分析代码
library(GSEAlens)

# 加载胶囊
res <- readRDS("%s")

# 查看概览
inspect_gsea_res(res)

# 查看所有对比组
names(res$results)

# 查看基因集信息
head(res$geneset_info$meta_dict)

# 提取特定对比组（示例）
task <- extract_gsea_task(res, contrast_id = names(res$results)[1])

# 查看该对比组的结果
head(as.data.frame(task$gsea_res@result))
',
          rds_path
  )
}
