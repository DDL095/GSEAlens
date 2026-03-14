# 请确保安装了以下依赖
library(dplyr)
library(tibble)
library(stringr)
library(limma)
library(clusterProfiler)
library(msigdbr)
library(BiocParallel)
library(GseaVis)
library(patchwork)
library(ggplot2)
library(DT)
library(htmlwidgets)
library(htmltools)


#' @title Build Custom Mixed MSigDB Dataframe
#' @description Helper function to fetch and row-bind multiple MSigDB collections while retaining all crucial metadata.
#' @param collection_list A list of character vectors, e.g., list(c("C2", "CP:KEGG"), c("C5", "GO:BP")).
#' @param species Character. Default "Homo sapiens".
#' @return A standardized dataframe ready for setup_gsea_env().
#' @export
build_custom_msigdb <- function(collection_list= list(c("C2", "CP:KEGG_LEGACY"),
                                                      c("C2", "CP:KEGG_MEDICUS"),
                                                      c("C5", "GO:BP"),
                                                      c("C5", "GO:CC"),
                                                      c("C5", "GO:MF")
),
species = "Homo sapiens") {
  message("🔧 Building custom mixed MSigDB database...")

  df_list <- lapply(collection_list, function(x) {
    cat <- x[1]
    subcat <- ifelse(length(x) > 1, x[2], "")
    message(sprintf("  -> Fetching %s : %s", cat, subcat))

    tmp <- msigdbr::msigdbr(species = species, category = cat, subcategory = subcat)
    # 强制保留核心列，并统一命名
    tmp %>% dplyr::select(
      gs_name,
      gene_symbol,
      gs_description,
      gs_url,
      gs_collection = gs_cat,
      gs_subcollection = gs_subcat
    )
  })

  final_df <- do.call(rbind, df_list)
  message(sprintf("✅ Built successfully! Total pathways: %d", length(unique(final_df$gs_name))))
  return(final_df)
}


#' @title Setup GSEA Environment (Polymorphic Engine)
#' @export
setup_gsea_env <- function(fit,
                           custom_df = NULL,
                           gmt_file = NULL,
                           custom_set_name = "Default_Set",
                           category = "H",
                           subcategory = NULL,
                           species = "Homo sapiens") {

  if (!inherits(fit, "MArrayLM")) stop("[Error] Input must be a limma 'fit' object.")

  # 1. 提取 Fit 中的比对组信息
  c_names <- colnames(fit)
  parsed <- lapply(c_names, function(x) {
    p <- strsplit(x, "\\s*-\\s*")[[1]]
    if(length(p) == 2) c(p[1], p[2]) else c(x, "Unknown")
  })
  parsed_df <- do.call(rbind, parsed)
  contrasts_df <- tibble::tibble(ID = seq_along(c_names), Contrast_Name = c_names, Num = parsed_df[,1], Den = parsed_df[,2])

  # 2. 核心：智能分发与元数据对齐
  term2gene <- NULL
  meta_dict <- NULL
  set_identifier <- custom_set_name

  if (!is.null(gmt_file)) {
    message("📥 Mode: External GMT File detected.")
    if (!file.exists(gmt_file)) stop(sprintf("[Error] GMT file not found: %s", gmt_file))
    if (custom_set_name == "Default_Set") stop("[Error] You must provide a meaningful 'custom_set_name' for external GMT.")

    lines <- readLines(gmt_file)
    gmt_list <- lapply(lines, function(l) strsplit(l, "\t")[[1]])

    meta_dict <- data.frame(
      ID = sapply(gmt_list, `[`, 1),
      Description = sapply(gmt_list, `[`, 2),
      URL = NA, Collection = "Custom_GMT", stringsAsFactors = FALSE
    )
    t2g_list <- lapply(gmt_list, function(x) data.frame(gs_name = x[1], gene_symbol = x[-(1:2)], stringsAsFactors = FALSE))
    term2gene <- do.call(rbind, t2g_list)

  } else if (!is.null(custom_df)) {
    message("🧬 Mode: Mixed Custom DataFrame detected.")
    req_cols <- c("gs_name", "gene_symbol")
    if (!all(req_cols %in% colnames(custom_df))) stop("[Error] custom_df is missing required columns.")

    term2gene <- custom_df %>% dplyr::select(gs_name, gene_symbol)

    # 提前补齐/统一缺失列 (在管道外部处理，避开 dplyr 解析器大坑)

    if (!"gs_description" %in% colnames(custom_df)) {
      custom_df$gs_description <- custom_df$gs_name
    }
    if (!"gs_url" %in% colnames(custom_df)) {
      custom_df$gs_url <- NA
    }
    if (!"gs_collection" %in% colnames(custom_df)) {
      custom_df$gs_collection <- "Mixed_Custom"
    }


    # 现在可以安全地进行 select 了，里面没有任何 if-else

    meta_dict <- custom_df %>%
      dplyr::select(
        ID = gs_name,
        Description = gs_description,
        URL = gs_url,
        Collection = gs_collection
      ) %>%
      dplyr::distinct(ID, .keep_all = TRUE)

    if (custom_set_name == "Default_Set") set_identifier <- "Mixed_Custom_Set"

  } else {
    # 模式 C：原生默认模式 (修复了空字符串报错问题)
    # 自动处理空子类 (比如 Hallmark 集合 H)
    if (!is.null(subcategory) && (is.na(subcategory) || nchar(trimws(subcategory)) == 0)) {
      subcategory <- NULL
    }

    if (is.null(subcategory)) {
      message(sprintf("🌐 Mode: Native MSigDB (Category: %s | No subcategory)", category))
      set_identifier <- category
    } else {
      message(sprintf("🌐 Mode: Native MSigDB (Category: %s | Subcategory: %s)", category, subcategory))
      set_identifier <- paste(category, subcategory, sep = "_")
    }

    df <- msigdbr::msigdbr(species = species, category = category, subcategory = subcategory)
    term2gene <- df %>% dplyr::select(gs_name, gene_symbol)
    meta_dict <- df %>% dplyr::select(ID = gs_name, Description = gs_description, URL = gs_url, Collection = gs_cat) %>% dplyr::distinct(ID, .keep_all = TRUE)
  }

  # 3. 组装坚如磐石的环境对象
  env_obj <- list(
    fit = fit,
    contrasts = contrasts_df,
    geneset = list(
      name = set_identifier,
      term2gene = term2gene,
      meta_dict = meta_dict
    )
  )
  class(env_obj) <- "GseaEnv"

  message(sprintf("✅ Environment Ready! Set Name: [%s], Total Pathways: %d", set_identifier, nrow(meta_dict)))
  return(env_obj)
}





#' @title Inspect GSEA Environment
#' @description Provides a beautifully formatted summary of the GseaEnv object, guiding the user on how to proceed to the next step.
#' @param gsea_env A GseaEnv object created by setup_gsea_env()
#' @export
inspect_gsea_env <- function(gsea_env) {

  if (!inherits(gsea_env, "GseaEnv")) {
    stop("[Error] Object is not a valid GseaEnv.")
  }

  # 提取基因集信息
  gs_name <- gsea_env$geneset$name
  total_pw <- nrow(gsea_env$geneset$meta_dict)
  total_genes <- length(unique(gsea_env$geneset$term2gene$gene_symbol))

  # 提取对比组信息
  c_df <- gsea_env$contrasts

  # 打印华丽的分割线和标题
  cat(rep("=", 65), "\n", sep = "")
  cat(sprintf("%-20s %s %-20s\n", "", "🧬 GSEA Environment Summary", ""))
  cat(rep("=", 65), "\n\n", sep = "")

  # Section 1: 基因集库信息
  cat("📦 [1] GENE SET DATABASE\n")
  cat(sprintf("  • Set Name       : %s\n", gs_name))
  cat(sprintf("  • Total Pathways : %s\n", format(total_pw, big.mark = ",")))
  cat(sprintf("  • Unique Genes   : %s\n\n", format(total_genes, big.mark = ",")))

  # Section 2: 核心指导 - 对比组清单
  cat("⚖️  [2] AVAILABLE CONTRASTS (From limma fit)\n")
  cat("  -> Use the 'ID' in run_directional_gsea(..., target_id = ID)\n\n")

  # 格式化表格打印
  header <- sprintf("  %-4s | %-20s | %-15s | %-15s", "ID", "Contrast_Name", "POSITIVE (Left)", "NEGATIVE (Right)")
  cat(header, "\n")
  cat("  ", rep("-", nchar(header)-2), "\n", sep="")

  for (i in 1:nrow(c_df)) {
    cat(sprintf("  %-4s | %-20s | %-15s | %-15s\n",
                c_df$ID[i], c_df$Contrast_Name[i], c_df$Num[i], c_df$Den[i]))
  }
  cat("\n")

  # Section 3: 下一步操作代码生成器
  cat("🚀 [3] NEXT STEP GUIDE (Copy & Paste)\n")

  # 默认拿第一个对比组做示例
  ex_id <- c_df$ID[1]
  ex_left <- c_df$Num[1]
  ex_right <- c_df$Den[1]

  cat(sprintf("  # Run default: POSITIVE NES = %s, NEGATIVE NES = %s\n", ex_left, ex_right))
  cat(sprintf("  res <- run_directional_gsea(gsea_env, target_id = %s)\n\n", ex_id))

  cat(sprintf("  # Force reversal: Make '%s' the POSITIVE (Left) group\n", ex_right))
  cat(sprintf("  res <- run_directional_gsea(gsea_env, target_id = %s, force_left_group = \"%s\")\n\n", ex_id, ex_right))

  cat(rep("=", 65), "\n", sep = "")

  return(invisible(gsea_env))
}

#' @export
print.GseaEnv <- function(x, ...) {
  inspect_gsea_env(x)
}

#' @title Run Directional GSEA
#' @export
run_directional_gsea <- function(gsea_env, target_id = 1, force_left_group = NULL,
                                 top_n = 40000, pvalueCutoff = 1, seed = 123) {

  if (!inherits(gsea_env, "GseaEnv")) stop("[Error] Input must be a GseaEnv object.")

  c_df <- gsea_env$contrasts
  if (!target_id %in% c_df$ID) stop("Invalid target_id.")

  t_row <- c_df %>% dplyr::filter(ID == target_id)
  coef_name <- t_row$Contrast_Name; left_grp <- t_row$Num; right_grp <- t_row$Den

  tt <- limma::topTable(gsea_env$fit, coef = coef_name, number = top_n) %>% as.data.frame()

  if (!is.null(force_left_group) && force_left_group == right_grp) {
    tt$t <- -tt$t
    left_grp <- t_row$Den; right_grp <- t_row$Num
  }
  genelist <- tt %>% dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    dplyr::mutate(SYMBOL = toupper(SYMBOL)) %>% dplyr::arrange(desc(abs(t))) %>%
    dplyr::distinct(SYMBOL, .keep_all = TRUE) %>% dplyr::select(SYMBOL, t) %>% tibble::deframe() %>% sort(decreasing = TRUE)

  message(sprintf("🚀 Running GSEA for %s vs %s...", left_grp, right_grp))
  set.seed(seed)
  res <- clusterProfiler::GSEA(geneList = genelist, TERM2GENE = gsea_env$geneset$term2gene, pvalueCutoff = pvalueCutoff, seed = seed)

  result_obj <- list(
    gsea_res = res,
    meta = list(
      left_group = left_grp, right_group = right_grp, contrast_name = coef_name,
      geneset_name = gsea_env$geneset$name, meta_dict = gsea_env$geneset$meta_dict
    )
  )
  class(result_obj) <- "DirectionalGSEA"
  return(result_obj)
}


#' @title Plot Directional GSEA (Bulletproof & Restores Classic Red/Blue Gradient)
#' @description 稳健的 GSEA 绘图函数，完美处理多通路合图、自动生成优雅的图例名称、
#'              并自带经典红蓝底部分布图。此版本通过 `subPlot` 参数控制是否显示 Rank List。
#' @param directional_gsea_obj 封装好的 DirectionalGSEA 对象 (需包含 gsea_res 和 meta)
#' @param target_pathways 需要绘制的通路 ID 向量 (如 c("HALLMARK_...", ...))
#' @param curveCol 自定义曲线颜色向量 (可选)
#' @param main_title 主标题名称
#' @param subPlot 控制 GseaVis::gseaNb 生成的子图数量。
#'                1: 仅GSEA富集图; 2: GSEA富集图 + 基因表达热图; 3: GSEA富集图 + 热图 + Rank List。
#' @param ... 传递给 GseaVis::gseaNb 的额外参数。
#' @return 一个 ggplot2 对象，表示 GSEA 结果图。
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_gradient2 geom_hline scale_x_continuous theme_bw theme element_blank element_text margin coord_cartesian labs geom_vline annotate
#' @importFrom grDevices colorRampPalette
#' @importFrom patchwork plot_annotation
#' @importFrom stringr str_to_title str_wrap
#' @importFrom GseaVis gseaNb
#' @export
plot_directional_gsea <- function(directional_gsea_obj, target_pathways, curveCol = NULL, main_title = "GSEA Plot", subPlot = 3, ...) {

  # 1. 解析对象与提取基础数据
  res <- directional_gsea_obj$gsea_res
  meta <- directional_gsea_obj$meta
  df <- as.data.frame(res)
  n_lines <- length(target_pathways)
  gList <- res@geneList

  # 2. 颜色池分配逻辑
  if (is.null(curveCol) || length(curveCol) < n_lines) {
    base_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D")
    curveCol <- if(n_lines <= length(base_colors)) base_colors[1:n_lines] else grDevices::colorRampPalette(base_colors)(n_lines)
  }

  curveCol_use <- curveCol[1:n_lines]
  curveCol_gsea <- curveCol_use
  if (n_lines == 1) curveCol_gsea <- c(curveCol_use[1], curveCol_use[1])

  # 3. 基础绘图 (调用 GseaVis::gseaNb) - 传入 subPlot 参数
  p_base <- GseaVis::gseaNb(
    object = res,
    geneSetID = target_pathways,
    subPlot = subPlot, # 💥 这里使用传入的 subPlot 参数
    addPval = FALSE,
    curveCol = curveCol_gsea,
    ...
  )


  # 4. 终极图例拦截与覆写逻辑
  if (n_lines > 1) {
    # a. 精准提取当前靶向通路的子集
    df_sub <- df[match(target_pathways, df$ID), ]

    # b. 锚定器与弹药库分离 (Separation of Concerns)
    raw_desc <- df_sub$Description  # 锚定器：底层 gseaNb 硬编码写在图上的名称
    name_id  <- df_sub$ID           # 弹药库：规整的标准 ID (永远有下划线)

    # c. 生成优雅短名称 (切除前缀 + 首字母大写 + 防爆换行)
    nice_labels <- sapply(name_id, function(x) {
      tit <- unlist(strsplit(x, split = "_"))
      if (length(tit) > 1) {
        formatted_text <- paste(stringr::str_to_title(tit[2:length(tit)]), collapse = " ")
      } else {
        formatted_text <- paste(stringr::str_to_title(tit), collapse = " ")
      }
      # 终极防线：自动在 45 个字符左右换行，防止 C5 通路图例挤爆画布
      stringr::str_wrap(formatted_text, width = 45)
    })

    # d. 执行 ggplot2 图层级强制覆写
    names(curveCol_use) <- raw_desc  # 必须使用 raw_desc 命名，确保映射准确

    override_scale <- ggplot2::scale_color_manual(
      name = "Term Name",
      values = curveCol_use,
      breaks = raw_desc,      # 锁定原始丑陋文本
      labels = nice_labels    # 替换为我们加工好的优雅文本
    )

    # 将新的 scale 注入到图形结构中 (应对 aplot 复合对象)
    # GseaVis::gseaNb 返回的是 an 'aplot' object, 其 ggplot2 对象存储在 plotlist 中
    if (inherits(p_base, "aplot")) {
      p_base$plotlist[[1]] <- p_base$plotlist[[1]] + override_scale
      # 只有当 subPlot 包含第二个子图时，才尝试修改其 scale
      if (subPlot >= 2 && !is.null(p_base$plotlist[[2]])) {
        p_base$plotlist[[2]] <- p_base$plotlist[[2]] + override_scale
      }
    } else { # 兼容 GseaVis 返回 ggplot 对象列表的情况
      p_base[[1]] <- p_base[[1]] + override_scale
      if (subPlot >= 2 && !is.null(p_base[[2]])) {
        p_base[[2]] <- p_base[[2]] + override_scale
      }
    }
  }


  # 5. 原生的经典红蓝底部分布图 (仅当 subPlot = 3 时生成和替换)
  if (subPlot == 3) { # 💥 条件判断：只有 subPlot = 3 时才生成 Rank List
    df_rank <- data.frame(x = 1:length(gList), y = as.numeric(gList))
    prank_classic <- ggplot2::ggplot(df_rank, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_col(ggplot2::aes(fill = y), width = 1, color = NA, show.legend = FALSE) +
      ggplot2::scale_fill_gradient2(low = "#08519C", mid = "white", high = "#A50F15", midpoint = 0) +
      ggplot2::geom_hline(yintercept = 0, size = 0.5, color = "black", linetype = "dashed") +
      ggplot2::scale_x_continuous(breaks = seq(0, length(gList), 5000)) +
      ggplot2::theme_bw(base_size = 14) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_text(colour = "black"), plot.margin = ggplot2::margin(t = -0.1, r = 0.2, b = 0.2, l = 0.2, unit = "cm")) +
      ggplot2::coord_cartesian(expand = 0) +
      ggplot2::labs(x = "Rank in Ordered Dataset", y = "Ranked List")

    z_cross <- sum(gList > 0); m_rank <- length(gList)
    anno_layers <- list(
      ggplot2::geom_vline(xintercept = z_cross, linetype = "dashed", color = "grey50"),
      ggplot2::annotate("text", x = z_cross, y = 0, label = paste0("Zero cross at ", z_cross), vjust = 1.5, hjust = -0.05, size = 3, color = "grey30"),
      ggplot2::annotate("text", x = m_rank * 0.01, y = max(gList) * 0.85, label = sprintf("'%s' (pos)", meta$left_group), color = "#A50F15", hjust = 0, size = 4, fontface = "italic"),
      ggplot2::annotate("text", x = m_rank * 0.99, y = min(gList) * 0.85, label = sprintf("'%s' (neg)", meta$right_group), color = "#08519C", hjust = 1, size = 4, fontface = "italic")
    )

    # 替换底部图层 (仅当 subPlot = 3 时)
    if (inherits(p_base, "aplot")) p_base$plotlist[[3]] <- prank_classic + anno_layers else p_base[[3]] <- prank_classic + anno_layers
  }

  # 6. 添加大标题与副标题
  p_final <- p_base + patchwork::plot_annotation(
    title = main_title, subtitle = sprintf("← %s | %s →", meta$left_group, meta$right_group),
    theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5, color = "grey40"))
  )

  return(p_final)
}







#' @title Generate GSEA HTML Report Bundle with Performance Control
#' @description Creates an independent bundle with main report, optionally generating detailed dashboards (Plots + Heatmaps) for top pathways to save time.
#' @param res_obj A completed DirectionalGSEA object from your pipeline.
#' @param output_base_dir The parent directory where the results bundle will be saved
#' @param p_adjust_cutoff Significance threshold for filtering pathways (Default: 1, keeps all pathways)
#' @param top_plots A numeric vector c(pos_N, neg_N) specifying how many top positive/negative pathways to plot. Set to c(0,0) to skip all plots. Set to c(Inf, Inf) for all.
#' @param save_pdf Logical, if TRUE, saves detailed plots as PDF as well as PNG.
#' @param dge_list Optional. An edgeR DGEList object to generate annotated expression heatmaps.
#' @export
#' @title 生成交互式 GSEA HTML 报告 (Pro终极优化版 - 扁平化目录 + 完整参数表)
#' @description 去除冗余文件夹嵌套，将 p.adjust/qvalue 引入表格，并强制 Description 在最后一列。
#' @export


#' @title Generate GSEA HTML Report Bundle (Pro 融合修复版)
#' @description 恢复原有的热图逻辑、HTML子网页、高亮红蓝表格，并严格采用 200 DPI PNG。
#' @export
generate_gsea_html_report <- function(res_obj, output_base_dir = "GSEA_Results",
                                      p_adjust_cutoff = 1, top_plots = c(15, 15),
                                      save_pdf = FALSE, dge_list = NULL, dpi = 200) {
  require(dplyr)
  if (!inherits(res_obj, "DirectionalGSEA")) stop("[Error] Invalid result object.")
  meta <- res_obj$meta
  res <- res_obj$gsea_res
  df <- as.data.frame(res)

  # 扁平化架构：不再套娃，直接使用母舰传进来的当前组别专属文件夹
  bundle_dir <- output_base_dir
  if (!dir.exists(bundle_dir)) dir.create(bundle_dir, recursive = TRUE)

  details_dir <- file.path(bundle_dir, "Details")
  if (!dir.exists(details_dir)) dir.create(details_dir, recursive = TRUE)

  # HTML 名与当前对比组文件夹同名
  main_html_path <- file.path(bundle_dir, paste0(basename(bundle_dir), ".html"))


  # 🌟 核心修改：构建 df_clean，确保 Pathway_Link 和 Description 各司其职
  df_clean <- df %>%
    dplyr::filter(p.adjust <= p_adjust_cutoff) %>%
    # 左连接 meta_dict，引入 long_description_for_html 和 URL
    dplyr::left_join(res_obj$meta$meta_dict, by = "ID") %>%
    dplyr::mutate(
      Enriched_In = factor(ifelse(NES > 0, meta$left_group, meta$right_group), levels = c(meta$left_group, meta$right_group)),
      Collection = as.factor(ifelse(is.na(Collection), "Unknown", Collection)),

      # 🌟 Pathway_Link 应该展示【短 ID】，并链接到 URL！
      Pathway_Link = ifelse(is.na(URL) | URL == "", sprintf("<b>%s</b>", ID), sprintf('<a href="%s" target="_blank" style="color: #0056b3;">%s</a>', URL, ID)),

      # 🌟 新增一个专门的 Description 列，用于显示【长文本】，纯文本不带链接！
      Description = long_description_for_html # 直接从连接来的 long_description_for_html 列赋值
    ) %>%
    dplyr::arrange(desc(abs(NES))) %>%
    dplyr::mutate(Rank = dplyr::row_number())

  if (nrow(df_clean) == 0) return(invisible(NULL))

  # 提取热图表达量 (完全沿用你的逻辑)
  has_expr <- FALSE; expr_mat <- NULL; cpm_mat <- NULL; sample_meta_sub <- NULL; n_left <- 0
  if (!is.null(dge_list)) {
    expr_mat <- edgeR::cpm(dge_list, log = TRUE)
    cpm_mat <- edgeR::cpm(dge_list, log = FALSE)
    sample_meta <- dge_list$samples
    left_samples <- rownames(sample_meta)[sample_meta$group == meta$left_group]
    right_samples <- rownames(sample_meta)[sample_meta$group == meta$right_group]
    target_samples <- c(left_samples, right_samples)

    if (length(target_samples) > 0) {
      expr_mat <- expr_mat[, target_samples, drop = FALSE]
      cpm_mat <- cpm_mat[, target_samples, drop = FALSE]
      sample_meta_sub <- sample_meta[target_samples, , drop = FALSE]
      n_left <- length(left_samples)
      has_expr <- TRUE
    }
  }

  detail_links <- character(nrow(df_clean))
  pos_plot_ids <- df_clean %>% dplyr::filter(NES > 0) %>% dplyr::arrange(desc(NES)) %>% head(top_plots[1]) %>% dplyr::pull(ID)
  neg_plot_ids <- df_clean %>% dplyr::filter(NES < 0) %>% dplyr::arrange(NES) %>% head(top_plots[2]) %>% dplyr::pull(ID)
  target_plot_ids <- c(pos_plot_ids, neg_plot_ids)

  # 遍历生成子图与详情页
  for (i in 1:nrow(df_clean)) {
    pw_id <- df_clean$ID[i]
    if (!(pw_id %in% target_plot_ids)) {
      detail_links[i] <- '<button class="btn btn-sm btn-outline-secondary" disabled style="padding: 2px 10px;">Skipped</button>'
      next
    }

    safe_pw_name <- gsub("[^A-Za-z0-9_.-]", "_", pw_id)
    detail_filename <- sprintf("Detail_Rank%03d_%s.html", i, safe_pw_name)
    gsea_png_name <- sprintf("GSEA_Rank%03d_%s.png", i, safe_pw_name)
    heat_png_name <- sprintf("Heatmap_Rank%03d_%s.png", i, safe_pw_name)

    # 1. 存 GSEA 图 (严格执行 200 DPI)
    p_gsea <- plot_directional_gsea(res_obj, pw_id, main_title = pw_id)
    ggplot2::ggsave(file.path(details_dir, gsea_png_name), plot = p_gsea, width = 8, height = 6, dpi = dpi, bg = "white")

    # 2. 基因统计表
    all_genes <- res@geneSets[[pw_id]]
    match_list_idx <- which(toupper(names(res@geneList)) %in% toupper(all_genes))
    valid_genes_list <- names(res@geneList)[match_list_idx]
    core_genes <- unlist(strsplit(as.character(df_clean$core_enrichment[i]), "/"))
    gene_table <- data.frame(
      Gene = valid_genes_list,
      Rank_Metric = unname(res@geneList[valid_genes_list]),
      Is_Core = ifelse(toupper(valid_genes_list) %in% toupper(core_genes), "✅ YES", "-"),
      stringsAsFactors = FALSE
    )
    gene_table <- gene_table[order(gene_table$Rank_Metric, decreasing = TRUE), ]
    colnames(gene_table) <- c("Gene_Name", "Rank_Metric", "Is_Core_Enrichment")
    html_table <- knitr::kable(gene_table, format = "html", row.names = FALSE, digits = 3, table.attr = 'class="table table-striped table-sm"')

    # 3. 绘制热图 (保留你神级的强制 -1 到 1 截断截断色彩)
    heat_html_tag <- "<p class='text-muted'>No expression data provided.</p>"
    if (has_expr) {
      expr_genes <- rownames(expr_mat)
      plot_genes <- expr_genes[which(toupper(expr_genes) %in% toupper(all_genes))]
      if (length(plot_genes) >= 2) {
        gene_metrics <- sapply(toupper(plot_genes), function(x) {
          idx <- match(x, toupper(names(res@geneList)))
          if (is.na(idx)) 0 else res@geneList[idx]
        })
        plot_genes <- plot_genes[order(gene_metrics, decreasing = TRUE)]
        plot_mat <- expr_mat[plot_genes, , drop = FALSE]
        plot_mat <- plot_mat[apply(plot_mat, 1, var) > 1e-6, , drop = FALSE]

        if (nrow(plot_mat) >= 2) {
          plot_cpm_numbers <- round(cpm_mat[rownames(plot_mat), , drop = FALSE])
          z_mat <- t(scale(t(plot_mat)))
          z_mat[is.na(z_mat)] <- 0
          z_mat[z_mat > 1] <- 1
          z_mat[z_mat < -1] <- -1

          ann_col <- data.frame(Group = sample_meta_sub$group, row.names = rownames(sample_meta_sub))
          grp_col <- c("#E41A1C", "#377EB8")
          names(grp_col) <- c(meta$left_group, meta$right_group)
          color_palette <- colorRampPalette(c("#67a9cf", "#f7f7f7", "#ef8a62"))(100)

          pheatmap::pheatmap(
            z_mat, scale = "none", cluster_cols = FALSE, cluster_rows = FALSE, gaps_col = n_left,
            color = color_palette, breaks = seq(-1, 1, length.out = 101),
            annotation_col = ann_col, annotation_colors = list(Group = grp_col),
            display_numbers = plot_cpm_numbers, number_color = "black", fontsize_number = 8,
            cellwidth = 35, cellheight = 16, filename = file.path(details_dir, heat_png_name),
            main = sprintf("Row-Scaled Z-Score [-1, 1]\nEnriched in: %s", as.character(df_clean$Enriched_In[i]))
          )
          heat_html_tag <- sprintf("<img src='%s' class='img-fluid'>", heat_png_name)
        }
      }
    }

    # 4. 生成子网页
    html_content <- sprintf('
      <!DOCTYPE html><html><head><meta charset="utf-8"><title>Detail: %s</title>
      <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet"></head>
      <body class="container mt-4 mb-5"><h2 class="text-primary">%s</h2>
      <p class="lead"><strong>Rank:</strong> %d | <strong>NES:</strong> %.3f | <strong>P-adj:</strong> %.4e</p><hr>
      <div class="row mb-5"><div class="col-md-5"><h4 class="mb-3">Enrichment Plot</h4><img src="%s" class="img-fluid border shadow-sm"></div>
      <div class="col-md-7"><h4 class="mb-3">Expression Heatmap</h4><div style="overflow:auto; max-height:800px;">%s</div></div></div>
      <div class="row"><div class="col-md-12"><h4 class="mb-3">Gene Statistics</h4><div style="max-height: 400px; overflow-y: auto;">%s</div></div></div>
      </body></html>
    ', pw_id, pw_id, i, df_clean$NES[i], df_clean$p.adjust[i], gsea_png_name, heat_html_tag, html_table)

    writeLines(html_content, con = file.path(details_dir, detail_filename))
    detail_links[i] <- sprintf('<a href="./Details/%s" target="_blank" class="btn btn-sm btn-success" style="padding: 2px 10px;">📊 Dashboard</a>', detail_filename)
  }

  df_clean$Detail_Page <- detail_links

  # 5. 生成主报表 (加入 pvalue，强制 Description 在最后一列)
  display_df <- df_clean %>%
    dplyr::select(Rank, Detail_Page, Pathway = Pathway_Link, Collection, Enriched_In, Size = setSize, NES, dplyr::any_of(c("pvalue", "p.adjust")), Description) # 🎯 这里的 Pathway 和 Description 列都是咱们精心构建的了！

  dt_table <- DT::datatable(
    display_df, rownames = FALSE, escape = FALSE, filter = "top",
    caption = htmltools::tags$caption(style = 'caption-side: top; text-align: center; font-size: 150%; font-weight: bold;', sprintf("GSEA Report: %s vs %s", meta$left_group, meta$right_group)),
    extensions = c('Buttons', 'Scroller'),
    options = list(dom = 'Bfrtip', deferRender = TRUE, scrollY = 600, scroller = TRUE, pageLength = -1, buttons = c('copy', 'csv', 'excel'), autoWidth = TRUE)
  ) %>%
    DT::formatRound(columns = c('NES'), digits = 3) %>%
    DT::formatSignif(columns = intersect(colnames(display_df), c('pvalue', 'p.adjust')), digits = 4) %>%
    DT::formatStyle('Enriched_In', backgroundColor = DT::styleEqual(c(meta$left_group, meta$right_group), c('#fee0d2', '#deebf7'))) %>%
    DT::formatStyle('NES', color = DT::styleInterval(0, c('blue', 'red')), fontWeight = 'bold')

  # 保存时分离 lib，防止卡顿
  old_wd <- getwd()
  setwd(bundle_dir)
  tryCatch({
    htmlwidgets::saveWidget(dt_table, file = basename(main_html_path), selfcontained = FALSE, libdir = "lib", title = sprintf("GSEA Report: %s vs %s [%s]", meta$left_group, meta$right_group, meta$geneset_name)#sprintf("🧬 %s vs %s", meta$left_group, meta$right_group)



    )
  }, finally = { setwd(old_wd) })

  return(invisible(bundle_dir))
}




gsea_genelist <- function(){




}



#' 打印并返回 msigdb_relation 数据
#'
#' @return 返回 msigdb_relation 数据框
#' @export
print_msigdb_relation <- function() {
  # 在 R 包内部，直接使用包内同名数据即可
  df_relation |>
    print() |>
    invisible() # 隐式返回，允许赋值但不在控制台刷屏两遍
}


#' @title Plot Directional GSEA (Ribbon Version)
#' @description 稳健的 GSEA 绘图函数，完美处理多通路合图、自动生成优雅的图例名称、
#'              并采用平滑的带状填充 Rank List，优化 SVG/PDF 性能。
#'              此版本使用连续颜色渐变，并通过 `subPlot` 参数控制是否显示 Rank List。
#' @param directional_gsea_obj 封装好的 DirectionalGSEA 对象 (需包含 gsea_res 和 meta)
#' @param target_pathways 需要绘制的通路 ID 向量 (如 c("HALLMARK_...", ...))
#' @param curveCol 自定义曲线颜色向量 (可选)
#' @param main_title 主标题名称
#' @param n_smooth_points 用于平滑 Rank List 的点数，点数越多越平滑，文件越大。
#' @param subPlot 控制 GseaVis::gseaNb 生成的子图数量。
#'                1: 仅GSEA富集图; 2: GSEA富集图 + 基因表达热图; 3: GSEA富集图 + 热图 + Rank List。
#' @param ... 传递给 GseaVis::gseaNb 的额外参数。
#' @return 一个 ggplot2 对象，表示 GSEA 结果图。
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line scale_fill_gradient2 geom_hline scale_x_continuous theme_bw theme element_blank element_text margin coord_cartesian labs geom_vline annotate
#' @importFrom grDevices colorRampPalette
#' @importFrom patchwork plot_annotation
#' @importFrom stringr str_to_title str_wrap
#' @importFrom GseaVis gseaNb
#' @export
plot_directional_gsea_ribbon <- function(directional_gsea_obj, target_pathways, curveCol = NULL, main_title = "GSEA Plot", n_smooth_points = 400, subPlot = 3, ...) {

  # 1. 解析对象与提取基础数据
  res <- directional_gsea_obj$gsea_res
  meta <- directional_gsea_obj$meta
  df <- as.data.frame(res)
  n_lines <- length(target_pathways)
  gList <- res@geneList

  # 2. 颜色池分配逻辑
  if (is.null(curveCol) || length(curveCol) < n_lines) {
    base_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D")
    curveCol <- if(n_lines <= length(base_colors)) base_colors[1:n_lines] else grDevices::colorRampPalette(base_colors)(n_lines)
  }

  curveCol_use <- curveCol[1:n_lines]
  curveCol_gsea <- curveCol_use
  if (n_lines == 1) curveCol_gsea <- c(curveCol_use[1], curveCol_use[1])

  # 3. 基础绘图 (调用 GseaVis::gseaNb) - 传入 subPlot 参数
  p_base <- GseaVis::gseaNb(
    object = res,
    geneSetID = target_pathways,
    subPlot = subPlot, # 💥 这里使用传入的 subPlot 参数
    addPval = FALSE,
    curveCol = curveCol_gsea,
    ...
  )


  # 4. 终极图例拦截与覆写逻辑 (保持不变)
  if (n_lines > 1) {
    # a. 精准提取当前靶向通路的子集
    df_sub <- df[match(target_pathways, df$ID), ]

    # b. 锚定器与弹药库分离 (Separation of Concerns)
    raw_desc <- df_sub$Description  # 锚定器：底层 gseaNb 硬编码写在图上的名称
    name_id  <- df_sub$ID           # 弹药库：规整的标准 ID (永远有下划线)

    # c. 生成优雅短名称 (切除前缀 + 首字母大写 + 防爆换行)
    nice_labels <- sapply(name_id, function(x) {
      tit <- unlist(strsplit(x, split = "_"))
      if (length(tit) > 1) {
        formatted_text <- paste(stringr::str_to_title(tit[2:length(tit)]), collapse = " ")
      } else {
        formatted_text <- paste(stringr::str_to_title(tit), collapse = " ")
      }
      # 终极防线：自动在 45 个字符左右换行，防止 C5 通路图例挤爆画布
      stringr::str_wrap(formatted_text, width = 45)
    })

    # d. 执行 ggplot2 图层级强制覆写
    names(curveCol_use) <- raw_desc  # 必须使用 raw_desc 命名，确保映射准确

    override_scale <- ggplot2::scale_color_manual(
      name = "Term Name",
      values = curveCol_use,
      breaks = raw_desc,      # 锁定原始丑陋文本
      labels = nice_labels    # 替换为我们加工好的优雅文本
    )

    # 将新的 scale 注入到图形结构中 (应对 aplot 复合对象)
    if (inherits(p_base, "aplot")) {
      p_base$plotlist[[1]] <- p_base$plotlist[[1]] + override_scale
      # 只有当 subPlot 包含第二个子图时，才尝试修改其 scale
      if (subPlot >= 2 && !is.null(p_base$plotlist[[2]])) {
        p_base$plotlist[[2]] <- p_base$plotlist[[2]] + override_scale
      }
    } else { # 兼容 GseaVis 返回 ggplot 对象列表的情况
      p_base[[1]] <- p_base[[1]] + override_scale
      if (subPlot >= 2 && !is.null(p_base[[2]])) {
        p_base[[2]] <- p_base[[2]] + override_scale
      }
    }
  }


  # 5. 【核心优化】平滑的带状填充 Rank List 图 (仅当 subPlot = 3 时生成和替换)
  if (subPlot == 3) { # 💥 条件判断：只有 subPlot = 3 时才生成 Rank List
    # 获取原始的排序基因列表数据
    df_rank_original <- data.frame(x = 1:length(gList), y = as.numeric(gList))

    # 创建一个平滑样条函数 (spline function)
    s_fun <- splinefun(x = df_rank_original$x, y = df_rank_original$y)

    # 生成新的 x 值，这些值在原始数据的范围内均匀分布
    x_smooth <- seq(min(df_rank_original$x), max(df_rank_original$x), length.out = n_smooth_points)

    # 从样条函数中获取对应的平滑 y 值
    y_smooth <- s_fun(x_smooth)

    # 构建用于绘图的平滑数据框
    df_rank_smoothed <- data.frame(x = x_smooth, y = y_smooth)

    # 绘制平滑的 Rank List (使用连续颜色渐变)
    prank_classic <- ggplot2::ggplot(df_rank_smoothed, ggplot2::aes(x = x, y = y)) + # 使用平滑后的数据
      # 使用 geom_ribbon 绘制填充区域，从 y=0 到曲线 y 值，并通过 fill 映射 y 值颜色
      ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = y, fill = y), show.legend = FALSE) +
      # 额外添加 geom_line 绘制平滑曲线的黑色轮廓，增强视觉清晰度
      ggplot2::geom_line(color = "black", size = 0.3) +

      ggplot2::scale_fill_gradient2(low = "#08519C", mid = "white", high = "#A50F15", midpoint = 0) + # 💥 回到连续渐变
      ggplot2::geom_hline(yintercept = 0, size = 0.5, color = "black", linetype = "dashed") +
      ggplot2::scale_x_continuous(breaks = seq(0, length(gList), 5000)) + # X 轴刻度保持与原始 gList 长度一致
      ggplot2::theme_bw(base_size = 14) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_text(colour = "black"), plot.margin = ggplot2::margin(t = -0.1, r = 0.2, b = 0.2, l = 0.2, unit = "cm")) +
      ggplot2::coord_cartesian(expand = 0) +
      ggplot2::labs(x = "Rank in Ordered Dataset", y = "Ranked List")

    # 计算 Zero Cross 点和总基因数 (这些仍然基于原始 gList 的索引和长度)
    z_cross <- sum(gList > 0)
    m_rank <- length(gList)

    # 确定平滑曲线的最大最小值，用于注释的定位，使其更精确
    max_y_smooth <- max(y_smooth)
    min_y_smooth <- min(y_smooth)

    anno_layers <- list(
      ggplot2::geom_vline(xintercept = z_cross, linetype = "dashed", color = "grey50"),
      ggplot2::annotate("text", x = z_cross, y = 0, label = paste0("Zero cross at ", z_cross), vjust = 1.5, hjust = -0.05, size = 3, color = "grey30"),
      # 调整注释的 y 坐标，使其相对于平滑曲线的 extrema 更合理
      ggplot2::annotate("text", x = m_rank * 0.01, y = max_y_smooth * 0.85, label = sprintf("'%s' (pos)", meta$left_group), color = "#A50F15", hjust = 0, size = 4, fontface = "italic"),
      ggplot2::annotate("text", x = m_rank * 0.99, y = min_y_smooth * 0.85, label = sprintf("'%s' (neg)", meta$right_group), color = "#08519C", hjust = 1, size = 4, fontface = "italic")
    )

    # 替换底部图层 (仅当 subPlot = 3 时)
    if (inherits(p_base, "aplot")) p_base$plotlist[[3]] <- prank_classic + anno_layers else p_base[[3]] <- prank_classic + anno_layers
  }

  # 6. 添加大标题与副标题 (保持不变)
  p_final <- p_base + patchwork::plot_annotation(
    title = main_title, subtitle = sprintf("← %s | %s →", meta$left_group, meta$right_group),
    theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5, color = "grey40"))
  )

  return(p_final)
}








#' @title 交互式构建 GSEA 基因集与超级标签 (Pro 引擎)
#' @description 基于 msigdbr 提供交互式菜单，选择基因集并生成智能语义化标签。支持未知未来子集的动态截断。
#' @param species 物种，默认 "Homo sapiens"
#' @param auto_select 可选。如果不想交互，可传入所选序号向量，例如 c(17, 25)
#' @return 返回包含 TERM2GENE、TERM2NAME、SuperTag 的列表对象
#' @export
build_gsea_pathways_pro <- function(species = "Homo sapiens", auto_select = NULL) {

  suppressPackageStartupMessages(require(msigdbr))
  suppressPackageStartupMessages(require(dplyr))

  # 动态拉取库
  avail_colls <- msigdbr::msigdbr_collections() %>%
    dplyr::arrange(gs_collection, gs_subcollection)

  # 白名单字典 (25 种核心集合)
  dict <- data.frame(
    gs_collection = c(
      "C1", "C2", "C2", "C2", "C2", "C2", "C2", "C2", "C2",
      "C3", "C3", "C3", "C3", "C4", "C4", "C4",
      "C5", "C5", "C5", "C5", "C6", "C7", "C7", "C8", "H"
    ),
    gs_subcollection = c(
      "", "CGP", "CP", "CP:BIOCARTA", "CP:KEGG_LEGACY", "CP:KEGG_MEDICUS",
      "CP:PID", "CP:REACTOME", "CP:WIKIPATHWAYS", "MIR:MIRDB", "MIR:MIR_LEGACY",
      "TFT:GTRD", "TFT:TFT_LEGACY", "3CA", "CGN", "CM",
      "GO:BP", "GO:CC", "GO:MF", "HPO", "", "IMMUNESIGDB", "VAX", "", ""
    ),
    short_tag = c(
      "C1Pos", "CGP", "CP", "BioC", "KeggL", "KeggM",
      "PID", "Reac", "Wiki", "MirDB", "MirL",
      "TFTG", "TFTL", "3CA", "CGN", "CM",
      "GoBP", "GoCC", "GoMF", "HPO", "C6Onc", "ImmS", "VAX", "C8Cell", "Hal"
    ),
    description = c(
      "C1 位置基因集", "C2 化学和遗传微扰", "C2 经典通路(总)", "C2 BioCarta", "C2 KEGG(旧)", "C2 KEGG(医学)",
      "C2 PID 通路", "C2 Reactome", "C2 WikiPathways", "C3 miRNA 靶标", "C3 miRNA(旧)",
      "C3 TF 靶标(GTRD)", "C3 TF(旧)", "C4 3CA 癌症", "C4 癌症邻居", "C4 癌症模块",
      "C5 GO 生物过程(BP)", "C5 GO 细胞组分(CC)", "C5 GO 分子功能(MF)", "C5 人类表型(HPO)", "C6 肿瘤特征", "C7 免疫特征", "C7 疫苗特征", "C8 细胞类型", "H 癌症标志(Hallmark)"
    ),
    stringsAsFactors = FALSE
  )

  avail_colls$gs_subcollection_clean <- gsub(".*:NO_SUB", "", avail_colls$gs_subcollection)
  menu_df <- dplyr::left_join(avail_colls, dict,
                              by = c("gs_collection", "gs_subcollection_clean" = "gs_subcollection"))

  # 🔮 动态截断机制 (应对未来 C9 等新集合)
  missing_idx <- is.na(menu_df$short_tag)
  if (any(missing_idx)) {
    menu_df$short_tag[missing_idx] <- paste0(
      menu_df$gs_collection[missing_idx], "_",
      substr(gsub("[^A-Za-z]", "", menu_df$gs_subcollection[missing_idx]), 1, 4)
    )
    menu_df$description[missing_idx] <- paste0("新集合: ", menu_df$gs_subcollection[missing_idx])
  }

  # 交互界面
  if (is.null(auto_select)) {
    message("\n", rep("=", 60))
    message(sprintf("🌟 欢迎使用 DudaliRnaseq PRO 基因集向导 (%s)", species))
    message(rep("=", 60))
    for (i in 1:nrow(menu_df)) {
      cat(sprintf("[%2d] %-6s | %-16s | %s\n",
                  i, menu_df$short_tag[i],
                  paste0(menu_df$gs_collection[i], ifelse(menu_df$gs_subcollection_clean[i]=="", "", paste0(":", menu_df$gs_subcollection_clean[i]))),
                  menu_df$description[i]))
    }
    user_input <- readline(prompt = "\n👉 请输入编号 (用逗号分隔, 如 17,25): ")
    selected_idx <- as.integer(unlist(strsplit(user_input, "[, ]+")))
    selected_idx <- selected_idx[!is.na(selected_idx) & selected_idx >= 1 & selected_idx <= nrow(menu_df)]
    if (length(selected_idx) == 0) stop("❌ 无效的输入！")
  } else {
    selected_idx <- auto_select
  }

  selected_rows <- menu_df[selected_idx, ]
  super_tag <- ifelse(nrow(selected_rows) <= 4,
                      paste(selected_rows$short_tag, collapse = "_"),
                      sprintf("Mix%d_%s", nrow(selected_rows), selected_rows$short_tag[1]))

  message(sprintf("\n✅ 已选定 %d 个集合。智能批次 Tag: [%s]", nrow(selected_rows), super_tag))

  # 提取数据 (前置拉取，节省后续多核重复读取的时间)
  pathway_list <- lapply(1:nrow(selected_rows), function(i) {
    c_coll <- selected_rows$gs_collection[i]
    c_sub <- selected_rows$gs_subcollection[i]
    if (grepl("NO_SUB", c_sub) || c_sub == "") {
      msigdbr::msigdbr(species = species, category = c_coll)
    } else {
      msigdbr::msigdbr(species = species, category = c_coll, subcategory = c_sub)
    }
  })

  all_pathways <- dplyr::bind_rows(pathway_list)
  TERM2GENE <- all_pathways %>% dplyr::select(gs_name, gene_symbol)
  TERM2NAME <- all_pathways %>% dplyr::select(ID = gs_name, Description = gs_description, URL = gs_url, Collection = gs_cat) %>% dplyr::distinct(ID, .keep_all = TRUE)

  return(list(
    TERM2GENE = TERM2GENE,
    meta_dict = TERM2NAME,
    SuperTag = super_tag,
    collections_used = selected_rows
  ))
}











#' @title 组装计算胶囊 (Pro 引擎)
#' @description 将差异结果(fit)、基因集与表达矩阵完美焊死在一起，实现一次打包，终身复现。
#' @param fit limma 分析得到的 MArrayLM 对象
#' @param pathway_obj build_gsea_pathways_pro() 返回的基因集对象
#' @param expr_data 你的 DGEList 或者标准化后的表达矩阵 (用于画热图，可为 NULL)
#' @export
setup_gsea_env_pro <- function(fit, pathway_obj, expr_data = NULL) {

  if (!inherits(fit, "MArrayLM")) stop("❌ fit 必须是 limma 的对象！")

  # 1. 自动捕获对比组（产生 12 组的源头！）
  c_names <- colnames(fit)
  parsed <- lapply(c_names, function(x) {
    p <- strsplit(x, "\\s*-\\s*")[[1]]
    if(length(p) == 2) c(p[1], p[2]) else c(x, "Unknown")
  })
  parsed_df <- do.call(rbind, parsed)
  contrasts_df <- tibble::tibble(ID = seq_along(c_names), Contrast_Name = c_names, Num = parsed_df[,1], Den = parsed_df[,2])

  # 2. 组装终极数据胶囊
  env_obj <- list(
    fit = fit,
    contrasts = contrasts_df,
    geneset = list(
      name = pathway_obj$SuperTag,       # 超级 Tag 被正式接管
      term2gene = pathway_obj$TERM2GENE, # 提前准备好的数据
      meta_dict = pathway_obj$meta_dict,
      used_collections = pathway_obj$collections_used
    ),
    expr_data = expr_data                # 🌟 表达矩阵被永远封存进胶囊中
  )

  class(env_obj) <- "GseaEnvPro"

  message(sprintf("✅ 胶囊封装完毕！对比组数: %d | Tag: [%s]", nrow(contrasts_df), pathway_obj$SuperTag))
  return(env_obj)
}



#' @title 检查计算胶囊 (Pro 引擎)
#' @description 极其美观地打印当前工作流所装载的所有元数据和对比任务。
#' @param env_pro setup_gsea_env_pro 创建的对象
#' @export
inspect_gsea_env_pro <- function(env_pro) {
  if (!inherits(env_pro, "GseaEnvPro")) stop("❌ 需要传入 GseaEnvPro 胶囊！")

  cat(rep("=", 65), "\n", sep = "")
  cat(sprintf("%-18s %s %-18s\n", "", "📦 GSEA PRO CAPSULE 封装报告", ""))
  cat(rep("=", 65), "\n\n", sep = "")

  cat("🧬 [1] 基因集字典 (Super Tag):", env_pro$geneset$name, "\n")
  cat(sprintf("  • 包含 %d 条通路映射关系。\n", nrow(env_pro$geneset$meta_dict)))

  cat("\n⚖️  [2] 从 fit 自动捕获的对比组 (将裂变运算):\n")
  for (i in 1:nrow(env_pro$contrasts)) {
    cat(sprintf("  [%d] %s  (左: %s | 右: %s)\n",
                env_pro$contrasts$ID[i], env_pro$contrasts$Contrast_Name[i],
                env_pro$contrasts$Num[i], env_pro$contrasts$Den[i]))
  }

  cat("\n🌡️  [3] 表达矩阵冻结状态:\n")
  if (is.null(env_pro$expr_data)) {
    cat("  ⚠️ 未传入 expr_data，后续 HTML 报告将只出折线图，不出热图。\n")
  } else {
    cat("  ✅ 表达矩阵已成功冻结入囊，随时可调取绘制热图！\n")
  }

  cat(rep("-", 65), "\n")
  cat("💡 下一步操作指引:\n")
  cat("   calc_res <- batch_parallel_gsea_pro(env_pro)\n")
  cat(rep("=", 65), "\n")
}



#' @title 并行计算 GSEA 核心引擎 (Pro 版本)
#' @description 自动感知硬件多核算力，对封装在胶囊中的所有对比组（支持双向真实计算，无翻转作弊）进行严格的 GSEA 置换检验。
#' @param gsea_env 由 setup_gsea_env_pro() 创建的计算胶囊
#' @param bidirectional 逻辑值，是否进行双向对比（自动生成 A_vs_B 和 B_vs_A），默认 TRUE
#' @param minGSSize 基因集最小包含基因数，默认 10
#' @param maxGSSize 基因集最大包含基因数，默认 500
#' @param pvalueCutoff GSEA 的 P 值阈值（推荐 1，保留全貌用于后续过滤），默认 1
#' @return 返回 DudaliGseaResPro 对象，包含元数据和所有计算结果
#' @export
#' @title 极速批量计算 GSEA (Pro - 修正版，回归原生稳定逻辑)
#' @description
#' 基于 topTable + SYMBOL + t 值 构建 geneList，
#' 支持真实双向计算，并将结果封装为 DudaliGseaResPro。
#'
#' @param gsea_env 由 setup_gsea_env_pro() 生成的 DudaliGseaEnvPro 对象
#' @param bidirectional 是否进行双向计算，默认 TRUE
#' @param workers 并行核心数，默认 NULL（自动检测并保留 4 核）
#' @param top_n topTable 提取基因数，默认 Inf
#' @param pvalueCutoff GSEA 的 pvalueCutoff，默认 1
#' @param minGSSize 最小基因集大小，默认 10
#' @param maxGSSize 最大基因集大小，默认 500
#' @param seed 随机种子，默认 123
#' @param future_maxsize_gb future 允许传输对象上限（GB），默认 32
#' @return DudaliGseaResPro 对象
#' @export
#' @title 并行计算 GSEA 核心引擎 (Pro 完美复刻版)
#' @description 结合原生代码最稳健的 t 值去重与大写转换逻辑，支持硬件核数自定与双向翻转计算。
#' @param gsea_env 由 setup_gsea_env_pro() 创建的计算胶囊
#' @param bidirectional 逻辑值，是否进行双向对比，默认 TRUE
#' @param workers 整数，手动指定使用的 CPU 核心数。如果为 NULL，则自动探测并保留 4 核。
#' @param minGSSize 基因集最小包含基因数，默认 10
#' @param maxGSSize 基因集最大包含基因数，默认 500
#' @param pvalueCutoff GSEA 的 P 值阈值（强制默认 1，保留全貌用于后续过滤）
#' @return 返回 DudaliGseaResPro 对象
#' @export
#' @title Batch Calculate GSEA Pro (Smart Caching Edition)
#' @export
batch_calc_gsea_pro <- function(gsea_env, custom_series_name = "Auto_Analysis", output_dir = "./GSEA_Output",
                                force = FALSE, bidirectional = TRUE, workers = 20,
                                minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1) {

  if (!inherits(gsea_env, "GseaEnvPro")) stop("❌ 请传入标准 GseaEnvPro 胶囊对象！")

  # ✅ 正确写法：从 gsea_env 里提取标签，用于给 RDS 命名
  super_tag <- gsea_env$geneset$name

  # 🔴 智能缓存拦截系统 (Smart Checkpoint)

  series_dir <- file.path(output_dir, custom_series_name)
  if (!dir.exists(series_dir)) dir.create(series_dir, recursive = TRUE, showWarnings = FALSE)

  rds_name <- sprintf("GSEA_Capsule_[%s]_[%s].rds", custom_series_name, super_tag)
  rds_path <- file.path(series_dir, rds_name)

  if (file.exists(rds_path) && force == FALSE) {
    message(sprintf("✅ 命中缓存！检测到已存在的 GSEA 胶囊: %s", rds_name))
    message("   触发极速载入模式，直接跳过计算环节...")
    return(readRDS(rds_path))
  }


  # 🟡 苦力计算系统

  total_cores <- parallel::detectCores(logical = TRUE)
  use_cores <- if (is.null(workers)) max(1, total_cores - 4) else min(total_cores, max(1, workers))
  message(sprintf("🖥️ 硬件侦测: 发现 %d 个逻辑核心。调度 %d 核执行并行计算！", total_cores, use_cores))

  future::plan(future::multisession, workers = use_cores)

  fit <- gsea_env$fit
  contrasts <- gsea_env$contrasts
  tasks <- list()

  message("🔍 正在提取 topTable 并进行去重与大小写清洗...")
  for (i in 1:nrow(contrasts)) {
    c_name <- contrasts$Contrast_Name[i]
    num <- contrasts$Num[i]; den <- contrasts$Den[i]

    tt <- limma::topTable(fit, coef = c_name, number = Inf) %>% as.data.frame()
    if (!"SYMBOL" %in% colnames(tt)) tt$SYMBOL <- rownames(tt)

    genelist_base <- tt %>%
      dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") %>%
      dplyr::mutate(SYMBOL = toupper(SYMBOL)) %>%
      dplyr::arrange(dplyr::desc(abs(t))) %>%
      dplyr::distinct(SYMBOL, .keep_all = TRUE) %>%
      dplyr::select(SYMBOL, t) %>%
      tibble::deframe()

    tasks[[paste0(num, "_vs_", den)]] <- sort(genelist_base, decreasing = TRUE)
    if (bidirectional) {
      tasks[[paste0(den, "_vs_", num)]] <- sort(-genelist_base, decreasing = TRUE)
    }
  }

  message(sprintf("🚀 洗练完毕，生成 %d 个独立对比 genelist。发射置换计算...", length(tasks)))

  term2gene <- gsea_env$geneset$term2gene
  res_list <- future.apply::future_lapply(names(tasks), function(task_name) {
    genelist <- tasks[[task_name]]
    set.seed(123)
    gsea_res <- tryCatch({
      clusterProfiler::GSEA(
        geneList = genelist, TERM2GENE = term2gene, minGSSize = minGSSize,
        maxGSSize = maxGSSize, pvalueCutoff = pvalueCutoff, pAdjustMethod = "BH", verbose = FALSE, seed = 123
      )
    }, error = function(e) NULL)

    status <- if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) "Success" else "Failed/NoEnrich"
    return(list(name = task_name, status = status, data = gsea_res, genelist = genelist))
  }, future.seed = TRUE)

  names(res_list) <- names(tasks)

  final_obj <- list(
    metadata = list(
      run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      workers_used = use_cores,
      parameters = list(bidirectional = bidirectional, minGSSize = minGSSize)
    ), # 👈 这里的逗号还在吗？
    geneset_info = gsea_env$geneset,  # 👈 这是我们新加的，注意尾部要有逗号！
    results = res_list,
    expr_data = gsea_env$expr_data
  )
  class(final_obj) <- "DudaliGseaResPro"
  future::plan(future::sequential)


  # 🟢 资产落袋与教学输出系统 (Archiving)

  saveRDS(final_obj, rds_path)

  example_task <- names(tasks)[1]
  cat(sprintf("\n=============================================================\n"))
  cat(sprintf("✅ [数据资产已存档] 计算胶囊已安全保存，防崩溃保护已启动！\n"))
  cat(sprintf("📁 文件路径: %s\n\n", rds_path))
  cat(sprintf("💡 未来如需写文章画单图，请复制以下 R 代码直接调取：\n"))
  cat(sprintf("-------------------------------------------------------------\n"))
  cat(sprintf("my_capsule <- readRDS(\"%s\")\n", rds_path))
  cat(sprintf("my_task <- DudaliRnaseq::extract_gsea_task_pro(my_capsule, \"%s\")\n", example_task))
  cat(sprintf("DudaliRnaseq::plot_directional_gsea(my_task, target_pathways = c(\"ID_1\", \"ID_2\"))\n"))
  cat(sprintf("=============================================================\n\n"))

  return(final_obj)
}







#' @title 查看 GSEA 计算胶囊结果 (Pro 引擎)
#' @description 格式化打印 DudaliGseaResPro 对象的元数据和各个对比组的富集状态，替代容易导致 RStudio 崩溃的 View()。
#' @param res_pro 由 batch_calc_gsea_pro() 生成的 DudaliGseaResPro 对象
#' @export
inspect_gsea_res_pro <- function(res_pro) {

  if (!inherits(res_pro, "DudaliGseaResPro")) {
    stop("❌ 必须传入 DudaliGseaResPro 对象！")
  }

  cat("\n", rep("=", 65), "\n", sep = "")
  cat(sprintf("%-18s %s %-18s\n", "", "🌟 GSEA Pro 运行报告 (Capsule)", ""))
  cat(rep("=", 65), "\n\n", sep = "")

  # 1. 基础元数据展示
  cat(sprintf("📅 运行时间 : %s\n", res_pro$metadata$run_time))
  cat(sprintf("🖥️ 调度核数 : %s Threads 并行\n", res_pro$metadata$workers_used))
  cat(sprintf("🧬 基因集合 : [%s]\n", res_pro$pathway_info$super_tag))

  expr_status <- if (!is.null(res_pro$expr_data)) "✅ 已封装 (支持自动渲染热图)" else "⚠️ 未封装 (仅出折线图)"
  cat(sprintf("🌡️ 表达矩阵 : %s\n", expr_status))
  cat("\n", rep("-", 65), "\n", sep = "")

  # 2. 组别运行状态展示
  cat(sprintf("📊 组别运行状态 (共 %d 组):\n\n", length(res_pro$results)))

  for (i in seq_along(res_pro$results)) {
    task_name <- names(res_pro$results)[i]
    task_info <- res_pro$results[[i]]

    if (task_info$status == "Success") {
      # 智能统计通路情况
      gsea_obj <- task_info$data
      total_count <- nrow(gsea_obj@result)
      sig_count <- sum(gsea_obj@result$p.adjust < 0.05, na.rm = TRUE)

      cat(sprintf("  [%02d] ✅ %-18s | 成功 (总富集: %4d 条, p.adj < 0.05: %3d 条)\n",
                  i, task_name, total_count, sig_count))
    } else {
      cat(sprintf("  [%02d] ❌ %-18s | 失败 (无显著富集结果或抛出异常)\n",
                  i, task_name))
    }
  }

  cat("\n", rep("=", 65), "\n", sep = "")
  cat("💡 下一步操作指引: \n")
  cat("  1. 批量出图生成网页: batch_plot_gsea_pro(obj)\n")
  cat("  2. 单独提取某组画图: my_data <- extract_gsea_task_pro(obj, '组别名')\n")
  cat(rep("=", 65), "\n\n", sep = "")
}














#' @title GSEA 模块化结果渲染器 (Pro 引擎 - 并行高控版)
#' @description 完全解耦的并行绘图流程。自动调度多核调用 generate_gsea_html_report()。
#' @param calc_batch_obj DudaliGseaResPro 计算结果胶囊
#' @param custom_series_name 自定义整体项目名称
#' @param output_dir 输出的总文件夹路径
#' @param top_plots 正负 NES top N 截断 c(pos_n, neg_n)
#' @param workers 渲染逻辑核数
#' @export
#' @title Batch Plot GSEA Pro
#' @export
batch_plot_gsea_pro <- function(calc_batch_obj, custom_series_name = "GSEA_Analysis",
                                output_dir = "./GSEA_Output", top_plots = c(15, 15),
                                workers = 20, dpi = 200) {


  if (!inherits(calc_batch_obj, "DudaliGseaResPro")) stop("❌ 必须传入 DudaliGseaResPro 胶囊！")
  super_tag <- calc_batch_obj$geneset_info$name
  master_dict <- calc_batch_obj$geneset_info$meta_dict
  expr_data <- calc_batch_obj$expr_data

  if (length(top_plots) == 1) top_plots <- c(top_plots, top_plots)

  # 将 HTML 报告文件夹统一放置在项目名称的子目录中，与 RDS 团聚
  series_dir <- file.path(output_dir, custom_series_name)
  if (!dir.exists(series_dir)) dir.create(series_dir, recursive = TRUE)

  success_tasks <- Filter(function(x) x$status == "Success", calc_batch_obj$results)
  if (length(success_tasks) == 0) {
    message("⚠️ 没有富集成功的对比组。")
    return(invisible(NULL))
  }

  total_cores <- parallel::detectCores(logical = TRUE)
  use_cores <- if (is.null(workers)) max(1, total_cores - 4) else min(total_cores, max(1, workers))
  message(sprintf("🎨 启动并行渲染引擎: 调度 %d 核对 %d 个对比组执行【全级联渲染】...", use_cores, length(success_tasks)))

  future::plan(future::multisession, workers = use_cores)

  future.apply::future_lapply(success_tasks, function(task) {
    # (在 future_lapply 内部)
    comp_name <- task$name
    parts <- strsplit(comp_name, "_vs_")[[1]]

    # 【完美旁路注入】：去原装字典里把长名字找回来！
    current_ids <- task$data@result$ID
    matched_idx <- match(current_ids, master_dict$ID)

    # 🌟 核心修改：将原始的 Description 列重命名为 long_description_for_html，
    # 🌟 这样在下游 generate_gsea_html_report 中就可以直接使用这个明确的列名。
    meta_dict <- data.frame(
      ID = current_ids,
      long_description_for_html = master_dict$Description[matched_idx], # 🎯 这里：长描述的正名
      URL = master_dict$URL[matched_idx],
      Collection = master_dict$Collection[matched_idx],
      stringsAsFactors = FALSE
    )

    res_obj <- list(
      gsea_res = task$data,
      meta = list(
        left_group = parts[1],
        right_group = parts[2],
        geneset_name = super_tag,
        meta_dict = meta_dict # 包含 long_description_for_html 的 meta_dict 已经打包进去了！
      )
    )
    class(res_obj) <- "DirectionalGSEA"

    # 存入专属子文件夹
    # 🌟 文件夹命名优化响应：
    # 结合对比组和基因集缩写！
    # 最终生成的文件夹名会是极其直观的："IBAA_vs_LTA_[Hal]" 或 "IBAA_vs_LTA_[C2KEGG_C5GO]"
    sub_folder_name <- sprintf("%s_[%s]", comp_name, super_tag)
    bundle_dir <- file.path(series_dir, sub_folder_name)

    tryCatch({
      generate_gsea_html_report(
        res_obj = res_obj, output_base_dir = bundle_dir, p_adjust_cutoff = 1,
        top_plots = top_plots, save_pdf = FALSE, dge_list = expr_data, dpi = dpi
      )
    }, error = function(e) warning(sprintf("任务 %s 渲染失败: %s", comp_name, e$message)))

    return(comp_name)
  }, future.seed = TRUE)

  future::plan(future::sequential)
  message(sprintf("\n🎉 带有热图与高亮表格的终极 HTML 报告已全部生成完毕！\n📂 点击直达: %s/", normalizePath(series_dir, winslash = "/", mustWork = FALSE)))
}

































#' @title Pro 数据提取器 (文章精细绘图专用 - 缝合外壳版)
#' @description 直接输入对比组名称，提取并包装成原生 DirectionalGSEA 对象，完美兼容底层 plot_directional_gsea 画图。
#' @param calc_res DudaliGseaResPro 计算胶囊
#' @param task_name 想提取的组别名 (如 "LTA_vs_POSTA")
#' @export
#' @title Extract Single GSEA Task Pro for Plotting
#' @description 从 DudaliGseaResPro 胶囊中提取单个对比的 GSEA 结果，并格式化为 DirectionalGSEA 对象，
#'              供 plot_directional_gsea 或 generate_gsea_html_report 单独使用。
#' @param gsea_capsule DudaliGseaResPro 对象。
#' @param task_name 想要提取的对比任务名称 (e.g., "Treatment_vs_Control")。
#' @return DirectionalGSEA 对象。
#' @export
extract_gsea_task_pro <- function(gsea_capsule, task_name) {

  if (!inherits(gsea_capsule, "DudaliGseaResPro")) stop("❌ 请传入标准 DudaliGseaResPro 胶囊对象！")
  if (!(task_name %in% names(gsea_capsule$results))) stop(sprintf("❌ 任务 '%s' 不存在于胶囊中！", task_name))

  # 1. 提取目标任务
  task <- gsea_capsule$results[[task_name]]

  # 2. 从任务名称中解析出左右组别
  parts <- strsplit(task_name, "_vs_")[[1]]
  left_grp <- parts[1]
  right_grp <- parts[2]

  # 3. 从总胶囊中获取原始的通路元数据字典 (这里面包含了长描述、URL、Collection等！)
  master_geneset_meta <- gsea_capsule$geneset_info$meta_dict

  # 4. 获取当前 GSEA 结果的通路 ID
  current_gsea_ids <- task$data@result$ID

  # 5. 将当前 GSEA 结果的 ID 匹配回总的通路元数据字典
  matched_indices <- match(current_gsea_ids, master_geneset_meta$ID)

  # 6. 构建 DirectionalGSEA 对象所需的 meta_dict
  meta_dict_for_res_obj <- data.frame(
    ID = current_gsea_ids,
    # 🌟 从 master_geneset_meta 中提取我们需要的长描述、URL、Collection！
    long_description_for_html = master_geneset_meta$Description[matched_indices], # 沿用我们之前确定的名字
    URL = master_geneset_meta$URL[matched_indices],
    Collection = master_geneset_meta$Collection[matched_indices],
    stringsAsFactors = FALSE
  )

  # 7. 组装 DirectionalGSEA 对象
  res_obj <- list(
    gsea_res = task$data,
    meta = list(
      left_group = left_grp,
      right_group = right_grp,
      # 🌟 从总胶囊中获取正确的 geneset 名称！
      geneset_name = gsea_capsule$geneset_info$name,
      meta_dict = meta_dict_for_res_obj # 使用构建好的 meta_dict
    )
  )
  class(res_obj) <- "DirectionalGSEA"

  return(res_obj)
}







#' @title GSEA 全流程一键终极总揽引擎 (All-in-One)
#' @description 从环境胶囊一键直达最终的高清 PNG 与自治 HTML 网页，并自动保存复现用的 RDS 文件。
#' @param gsea_env setup_gsea_env_pro 创建的环境胶囊
#' @param custom_series_name 自定义分析批次名称（如 "2026_03_09_最终汇报"）
#' @param save_rds 逻辑值，是否自动在根目录留存包含所有信息的数据胶囊
#' @export
#' @title 终极管线：全自动 GSEA 计算与并行渲染总控
#' @export
batch_parallel_gsea_pro <- function(
    gsea_env,
    custom_series_name = "Auto_Analysis",
    output_dir = "./GSEA_Output",
    force = FALSE,               # 👈 控制是否无视缓存强制重算
    bidirectional = TRUE,
    top_plots = c(15, 15),
    workers = 20,
    dpi = 200,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1
) {

  # 步骤 1：呼叫智能计算引擎 (它会自己决定是瞬间载入 RDS 还是老实计算)
  calc_res <- DudaliRnaseq::batch_calc_gsea_pro(
    gsea_env = gsea_env,
    custom_series_name = custom_series_name,
    output_dir = output_dir,
    force = force,
    bidirectional = bidirectional,
    workers = workers,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    pvalueCutoff = pvalueCutoff
  )

  # 步骤 2：呼叫渲染引擎 (传入计算所得的胶囊对象)
  DudaliRnaseq::batch_plot_gsea_pro(
    calc_batch_obj = calc_res,
    custom_series_name = custom_series_name,
    output_dir = output_dir,
    top_plots = top_plots,
    workers = workers,
    dpi = dpi
  )

  # 隐式返回胶囊，方便在控制台继续折腾
  return(invisible(calc_res))
}





#' @title Batch Parallel GSEA Computation (Only Math)
#' @description 仅进行并行的 GSEA 计算，返回计算结果列表，不出图（速度极快）
#' @export
batch_calc_gsea <- function(gsea_env, bidirectional = TRUE, workers = 25, series_name = NULL, ...) {

  con_df <- gsea_env$contrasts
  if (is.null(con_df)) con_df <- gsea_env$contrasts_info
  if (is.null(con_df)) stop("❌ 找不到对比表格！")

  # 优先使用用户手动传入的 series_name
  if (!is.null(series_name)) {
    set_name <- series_name
  } else {
    set_name <- gsea_env$set_name
    if (is.null(set_name)) set_name <- "GSEA_Result"
  }

  tasks <- list()
  for (i in 1:nrow(con_df)) {
    id <- con_df[[1]][i]
    left <- con_df[[3]][i]
    right <- con_df[[4]][i]

    tasks[[length(tasks) + 1]] <- list(id = id, force_left = left, name = sprintf("[%s]_%s_vs_%s", set_name, left, right))
    if (bidirectional) {
      tasks[[length(tasks) + 1]] <- list(id = id, force_left = right, name = sprintf("[%s]_%s_vs_%s", set_name, right, left))
    }
  }

  message(sprintf("🧮 [阶段 1] 开始并行计算 GSEA (共 %d 个任务, %d 线程)...", length(tasks), workers))
  message(sprintf("   🏷️ 当前分析批次名称 (Series Name): %s", set_name))

  # 🔓 解除 future 包的 500MB 数据传输限制（设定为 32GB）
  options(future.globals.maxSize = 32000 * 1024^2)
  future::plan(future::multisession, workers = workers)

  calc_results <- future.apply::future_lapply(tasks, function(task) {
    suppressPackageStartupMessages({
      require(dplyr, quietly = TRUE)
      require(magrittr, quietly = TRUE)
      require(clusterProfiler, quietly = TRUE)
      require(DudaliRnaseq, quietly = TRUE)
    })

    tryCatch({
      res <- run_directional_gsea(gsea_env = gsea_env, target_id = task$id, force_left_group = task$force_left, ...)
      return(list(name = task$name, status = "Success", data = res))
    }, error = function(e) {
      return(list(name = task$name, status = "Failed", error = e$message, data = NULL))
    })
  }, future.seed = TRUE, future.packages = c("dplyr", "magrittr", "clusterProfiler", "DudaliRnaseq"))

  future::plan(future::sequential)
  message("✅ GSEA 计算阶段完成！")
  return(calc_results)
}

#' @title Batch Parallel GSEA HTML Report Generation
#' @description 接收 batch_calc_gsea 的结果，多核并行生成 HTML 报告和折线图
#' @export
batch_plot_gsea <- function(calc_results,
                            output_root = "./GSEA_Final_Reports",
                            series_name = NULL,
                            workers = 4,
                            p_adjust_cutoff = 1,
                            top_plots = c(20, 20),
                            save_pdf = FALSE,
                            dge_list = NULL,
                            ...) {

  valid_tasks <- Filter(function(x) x$status == "Success" && !is.null(x$data), calc_results)
  message(sprintf("📊 [阶段 2] 开始并行生成 HTML 报告 (共 %d 个任务, %d 线程)...", length(valid_tasks), workers))

  # 🔓 再次解除出图阶段的 500MB 数据传输限制（设定为 8GB）
  options(future.globals.maxSize = 8000 * 1024^2)
  future::plan(future::multisession, workers = workers)

  plot_results <- future.apply::future_lapply(valid_tasks, function(task) {
    suppressPackageStartupMessages({
      require(dplyr, quietly = TRUE)
      require(magrittr, quietly = TRUE)
      require(DudaliRnaseq, quietly = TRUE)
    })

    tryCatch({
      df <- as.data.frame(task$data$gsea_res)
      if (nrow(df) > 0) {

        current_name <- task$name
        if (!is.null(series_name)) {
          current_name <- sub("^\\[.*?\\]_", sprintf("[%s]_", series_name), current_name)
        }

        out_dir <- file.path(output_root, current_name)
        if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

        generate_gsea_html_report(
          res_obj = task$data,
          output_base_dir = out_dir,
          p_adjust_cutoff = p_adjust_cutoff,
          top_plots = top_plots,
          save_pdf = save_pdf,
          dge_list = dge_list,
          ...
        )

        return(list(name = current_name, status = "✅ 图表生成成功"))
      } else {
        return(list(name = task$name, status = "⚠️ 无显著通路，跳过出图"))
      }
    }, error = function(e) {
      return(list(name = task$name, status = sprintf("❌ 出图失败: %s", e$message)))
    })
  }, future.seed = TRUE, future.packages = c("dplyr", "magrittr", "DudaliRnaseq"))

  future::plan(future::sequential)

  message("\n=======================================================")
  for (res in plot_results) cat(sprintf("%s : %s\n", res$status, res$name))
  message("=======================================================")
  message(sprintf("📂 报告总目录: %s", normalizePath(output_root, mustWork = FALSE)))
}

#' @title All-in-One: Batch Parallel GSEA & Report
#' @description 一键全自动：先并行计算，后并行出图，并允许精细控制报告生成参数
#' @export
batch_parallel_gsea <- function(gsea_env,
                                output_root = "./GSEA_Final_Reports",
                                series_name = NULL,
                                bidirectional = TRUE,
                                workers = 4,
                                p_adjust_cutoff = 1,
                                top_plots = c(10, 10),
                                save_pdf = FALSE,
                                dge_list = NULL,
                                ...) {
  calc_res <- batch_calc_gsea(
    gsea_env = gsea_env,
    bidirectional = bidirectional,
    workers = workers,
    series_name = series_name,
    ...
  )

  batch_plot_gsea(
    calc_results = calc_res,
    output_root = output_root,
    series_name = series_name,
    workers = workers,
    p_adjust_cutoff = p_adjust_cutoff,
    top_plots = top_plots,
    save_pdf = save_pdf,
    dge_list = dge_list,
    ...
  )

  invisible(calc_res)
}
