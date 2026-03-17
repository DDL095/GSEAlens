#' @title Plot Directional GSEA (Bulletproof & Restores Classic Red/Blue Gradient)
#' @description 极客级 GSEA 绘图引擎，完美兼容单通路/多通路合图。
#'              自动生成优雅图例、拦截并重绘经典红蓝基因分布带。统计学标注默认关闭。
#' @param directional_gsea_obj 封装好的 DirectionalGSEA 对象 (提取器直接返回的结果)
#' @param target_pathways 需要绘制的通路 ID 向量 (支持单个或多个，如 c("ID_1", "ID_2"))
#' @param curveCol 自定义曲线颜色向量 (可选，如果不传会自动分配高级调色盘)
#' @param main_title 主标题名称
#' @param subPlot 控制 GseaVis::gseaNb 生成的子图数量 (1: 仅富集图; 2: 富集+热图; 3: 完整带Rank)
#' @param add_pval 是否在图上添加统计学标注 (P-value/FDR)。默认 FALSE (保持画面纯净)。
#' @param ... 传递给 GseaVis::gseaNb 的额外参数 (例如 pvalX, pvalY 控制文本位置)
#' @return ggplot2 复合对象
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_gradient2 geom_hline scale_x_continuous theme_bw theme element_blank element_text margin coord_cartesian labs geom_vline annotate scale_color_manual
#' @importFrom grDevices colorRampPalette
#' @importFrom patchwork plot_annotation
#' @importFrom stringr str_to_title str_wrap
#' @export
plot_directional_gsea <- function(directional_gsea_obj, target_pathways, curveCol = NULL,
                                  main_title = "GSEA Plot", subPlot = 3, add_pval = FALSE, ...) {

  # 1. 解析对象与提取基础数据
  res <- directional_gsea_obj$gsea_res
  meta <- directional_gsea_obj$meta
  df <- as.data.frame(res)
  n_lines <- length(target_pathways)
  gList <- res@geneList

  # 2. 高级颜色池分配逻辑 (单线纯色，多线渐变或高对比)
  if (is.null(curveCol) || length(curveCol) < n_lines) {
    base_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#1B9E77", "#D95F02", "#7570B3", "#E7298A")
    curveCol <- if(n_lines <= length(base_colors)) base_colors[1:n_lines] else grDevices::colorRampPalette(base_colors)(n_lines)
  }

  curveCol_use <- curveCol[1:n_lines]
  curveCol_gsea <- curveCol_use
  if (n_lines == 1) curveCol_gsea <- c(curveCol_use[1], curveCol_use[1])

  # 3. 基础绘图 (原生调用 gseaNb) - 动态传入 add_pval
  p_base <- GseaVis::gseaNb(
    object = res,
    geneSetID = target_pathways,
    subPlot = subPlot,
    addPval = add_pval, # 💡 开放接口，默认 FALSE
    curveCol = curveCol_gsea,
    ...
  )

  # 4. 终极图例拦截与覆写逻辑 (专治多通路名字难看)
  if (n_lines > 1) {
    df_sub <- df[match(target_pathways, df$ID), ]
    raw_desc <- df_sub$Description  # 底层硬编码的描述
    name_id  <- df_sub$ID           # 我们的原生 ID

    # 智能截断与标题化 (您的这段逻辑非常完美，保留！)
    nice_labels <- sapply(name_id, function(x) {
      tit <- unlist(strsplit(x, split = "_"))
      if (length(tit) > 1) {
        formatted_text <- paste(stringr::str_to_title(tit[2:length(tit)]), collapse = " ")
      } else {
        formatted_text <- paste(stringr::str_to_title(tit), collapse = " ")
      }
      stringr::str_wrap(formatted_text, width = 45) # 防爆换行
    })

    names(curveCol_use) <- raw_desc
    override_scale <- ggplot2::scale_color_manual(
      name = "Term Name",
      values = curveCol_use,
      breaks = raw_desc,
      labels = nice_labels
    )

    # 兼容 aplot 与 ggplot 列表结构强制注入 scale
    if (inherits(p_base, "aplot")) {
      p_base$plotlist[[1]] <- p_base$plotlist[[1]] + override_scale
      if (subPlot >= 2 && !is.null(p_base$plotlist[[2]])) {
        p_base$plotlist[[2]] <- p_base$plotlist[[2]] + override_scale
      }
    } else {
      p_base[[1]] <- p_base[[1]] + override_scale
      if (subPlot >= 2 && !is.null(p_base[[2]])) {
        p_base[[2]] <- p_base[[2]] + override_scale
      }
    }
  }

  # 5. 原生经典红蓝底部分布图 (强力替换底层自带的刻板条码)
  if (subPlot == 3) {
    df_rank <- data.frame(x = 1:length(gList), y = as.numeric(gList))
    prank_classic <- ggplot2::ggplot(df_rank, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_col(ggplot2::aes(fill = y), width = 1, color = NA, show.legend = FALSE) +
      ggplot2::scale_fill_gradient2(low = "#08519C", mid = "white", high = "#A50F15", midpoint = 0) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.5, color = "black", linetype = "dashed") + # 修复了 size 警告
      ggplot2::scale_x_continuous(breaks = seq(0, length(gList), 5000)) +
      ggplot2::theme_bw(base_size = 14) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_text(colour = "black"),
                     plot.margin = ggplot2::margin(t = -0.1, r = 0.2, b = 0.2, l = 0.2, unit = "cm")) +
      ggplot2::coord_cartesian(expand = 0) +
      ggplot2::labs(x = "Rank in Ordered Dataset", y = "Ranked List")

    z_cross <- sum(gList > 0); m_rank <- length(gList)
    anno_layers <- list(
      ggplot2::geom_vline(xintercept = z_cross, linetype = "dashed", color = "grey50"),
      ggplot2::annotate("text", x = z_cross, y = 0, label = paste0("Zero cross at ", z_cross), vjust = 1.5, hjust = -0.05, size = 3, color = "grey30"),
      ggplot2::annotate("text", x = m_rank * 0.01, y = max(gList) * 0.85, label = sprintf("'%s' (pos)", meta$left_group), color = "#A50F15", hjust = 0, size = 4, fontface = "italic"),
      ggplot2::annotate("text", x = m_rank * 0.99, y = min(gList) * 0.85, label = sprintf("'%s' (neg)", meta$right_group), color = "#08519C", hjust = 1, size = 4, fontface = "italic")
    )

    if (inherits(p_base, "aplot")) {
      p_base$plotlist[[3]] <- prank_classic + anno_layers
    } else {
      p_base[[3]] <- prank_classic + anno_layers
    }
  }

  # 6. 添加结构化大标题
  p_final <- p_base + patchwork::plot_annotation(
    title = main_title, subtitle = sprintf("← %s | %s →", meta$left_group, meta$right_group),
    theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
                           plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5, color = "grey40"))
  )

  return(p_final)
}


#' @title Generate GSEA HTML Report Bundle (终极智能寻路版)
#' @description 自动规避列名污染，智能抓取对象内的表达矩阵，输出带双P值展示的舒适网页报表与热图。支持根据计算胶囊原地址自动归位报表。
#' @param res_obj extract_gsea_task_pro 返回的 DirectionalGSEA 对象
#' @param output_base_dir 输出文件夹路径。若留空(NULL)，将自动跟随计算胶囊的项目原址建立专属 HTML 文件夹！
#' @param p_adjust_cutoff FDR 过滤阈值，默认 1
#' @param top_plots 绘制详细子网页的通路数量，格式为 c(正向数量, 负向数量)
#' @param save_pdf 预留参数
#' @param dge_list edgeR DGEList 对象。若为空，将全网自动在 res_obj$meta 内雷达搜寻表达矩阵！
#' @param dpi GSEA 富集图的分辨率，默认 200
#' @export
generate_gsea_html_report <- function(res_obj, output_base_dir = NULL,
                                      p_adjust_cutoff = 1, top_plots = c(15, 15),
                                      save_pdf = FALSE, dge_list = NULL, dpi = 200) {
  require(dplyr)
  if (!inherits(res_obj, "DirectionalGSEA")) stop("❌ [Error] 传入的对象不是标准的 DirectionalGSEA。")

  meta <- res_obj$meta
  res <- res_obj$gsea_res
  df <- as.data.frame(res)

  # 🌟 1. 自动寻路逻辑：追溯原发地
  if (is.null(output_base_dir)) {
    if (!is.null(meta$project_info) && !is.null(meta$project_info$series_dir)) {
      output_base_dir <- file.path(meta$project_info$series_dir, sprintf("HTML_Report_%s_vs_%s", meta$left_group, meta$right_group))
    } else {
      output_base_dir <- "GSEA_Results" # 备用回退文件夹
    }
  }

  bundle_dir <- output_base_dir
  if (!dir.exists(bundle_dir)) dir.create(bundle_dir, recursive = TRUE)

  details_dir <- file.path(bundle_dir, "Details")
  if (!dir.exists(details_dir)) dir.create(details_dir, recursive = TRUE)

  main_html_path <- file.path(bundle_dir, paste0(basename(bundle_dir), ".html"))


  # 🛡️ 2. 核心清洗区：智能 left_join 避免列名污染
  dict_cols <- colnames(meta$meta_dict)
  df_cols <- colnames(df)
  cols_to_add <- setdiff(dict_cols, df_cols) # 取差集
  safe_dict <- meta$meta_dict %>% dplyr::select(ID, dplyr::any_of(cols_to_add))

  df_clean <- df %>%
    dplyr::filter(p.adjust <= p_adjust_cutoff) %>%
    dplyr::left_join(safe_dict, by = "ID") %>%
    dplyr::mutate(
      Enriched_In = factor(ifelse(NES > 0, meta$left_group, meta$right_group), levels = c(meta$left_group, meta$right_group)),
      Display_Collection = if("Combo_Name" %in% names(.)) Combo_Name else if("Collection" %in% names(.)) Collection else "Unknown",
      Display_Collection = as.factor(ifelse(is.na(Display_Collection), "Unknown", Display_Collection)),
      Pathway_Link = if("URL" %in% names(.)) {
        ifelse(is.na(URL) | URL == "", sprintf("<b>%s</b>", ID), sprintf('<a href="%s" target="_blank" style="color: #0056b3; text-decoration: none;">%s</a>', URL, ID))
      } else { sprintf("<b>%s</b>", ID) },
      Description = if("long_description_for_html" %in% names(.)) long_description_for_html else ID
    ) %>%
    dplyr::arrange(desc(abs(NES))) %>%
    dplyr::mutate(Rank = dplyr::row_number())

  if (nrow(df_clean) == 0) {
    message("⚠️ 警告：当前截断值下没有显著通路，跳过报告生成！")
    return(invisible(NULL))
  }

  # 📡 3. 表达矩阵智能雷达：探测热图数据
  message("🔍 [智能雷达] 正在探测表达矩阵用于绘制热图...")
  if (is.null(dge_list)) {
    if (!is.null(meta$expr_data)) {
      dge_list <- meta$expr_data
      message("   ✅ 成功从 res_obj$meta 提取 expr_data")
    } else {
      message("   ⚠️ 未找到内置的 expr_data，将跳过热图绘制。")
    }
  }

  has_expr <- FALSE; expr_mat <- NULL; cpm_mat <- NULL; sample_meta_sub <- NULL; n_left <- 0

  if (!is.null(dge_list)) {
    sample_meta <- dge_list$samples
    if (!"group" %in% colnames(sample_meta)) {
      message("   ❌ [错误] 表达矩阵中缺少 'group' 样本列，无法区分分组！")
    } else {
      left_samples <- rownames(sample_meta)[sample_meta$group == meta$left_group]
      right_samples <- rownames(sample_meta)[sample_meta$group == meta$right_group]
      target_samples <- c(left_samples, right_samples)

      if (length(target_samples) == 0) {
        message(sprintf("   ❌ [错误] 无法在矩阵中找到对比组 '%s' 或 '%s'！", meta$left_group, meta$right_group))
      } else {
        message(sprintf("   ✅ [热图就绪] 匹配到样本: 左组(%s) %d个 | 右组(%s) %d个", meta$left_group, length(left_samples), meta$right_group, length(right_samples)))
        expr_mat <- edgeR::cpm(dge_list, log = TRUE)[, target_samples, drop = FALSE]
        cpm_mat <- edgeR::cpm(dge_list, log = FALSE)[, target_samples, drop = FALSE]
        sample_meta_sub <- sample_meta[target_samples, , drop = FALSE]
        n_left <- length(left_samples)
        has_expr <- TRUE
      }
    }
  }

  detail_links <- character(nrow(df_clean))
  pos_plot_ids <- df_clean %>% dplyr::filter(NES > 0) %>% dplyr::arrange(desc(NES)) %>% head(top_plots[1]) %>% dplyr::pull(ID)
  neg_plot_ids <- df_clean %>% dplyr::filter(NES < 0) %>% dplyr::arrange(NES) %>% head(top_plots[2]) %>% dplyr::pull(ID)
  target_plot_ids <- c(pos_plot_ids, neg_plot_ids)

  # 🖼️ 4. 遍历生成子图与详情页
  for (i in seq_len(nrow(df_clean))) {
    pw_id <- df_clean$ID[i]
    if (!(pw_id %in% target_plot_ids)) {
      detail_links[i] <- '<button class="btn btn-sm btn-outline-secondary" disabled style="padding: 2px 10px;">Skipped</button>'
      next
    }

    safe_pw_name <- gsub("[^A-Za-z0-9_.-]", "_", pw_id)
    detail_filename <- sprintf("Detail_Rank%03d_%s.html", i, safe_pw_name)
    gsea_png_name <- sprintf("GSEA_Rank%03d_%s.png", i, safe_pw_name)
    heat_png_name <- sprintf("Heatmap_Rank%03d_%s.png", i, safe_pw_name)

    known_prefixes <- c("KEGG_", "REACTOME_", "WP_", "BIOCARTA_", "GO_", "HP_", "HALLMARK_")
    clean_title <- pw_id
    for (pref in known_prefixes) { if (startsWith(pw_id, pref)) { clean_title <- sub(pref, "", pw_id); break } }
    clean_title <- stringr::str_to_title(gsub("_", " ", clean_title))

    # 4.1 绘制 GSEA 主图
    p_gsea <- plot_directional_gsea(res_obj, target_pathways = pw_id, main_title = clean_title, subPlot = 3, add_pval = FALSE)
    ggplot2::ggsave(file.path(details_dir, gsea_png_name), plot = p_gsea, width = 8, height = 6, dpi = dpi, bg = "white")

    # 4.2 基因统计表
    all_genes <- res@geneSets[[pw_id]]
    match_list_idx <- which(toupper(names(res@geneList)) %in% toupper(all_genes))
    valid_genes_list <- names(res@geneList)[match_list_idx]
    core_genes <- unlist(strsplit(as.character(df_clean$core_enrichment[i]), "/"))
    gene_table <- data.frame(
      Gene = valid_genes_list, Rank_Metric = unname(res@geneList[valid_genes_list]),
      Is_Core = ifelse(toupper(valid_genes_list) %in% toupper(core_genes), "✅ YES", "-"), stringsAsFactors = FALSE
    )
    gene_table <- gene_table[order(gene_table$Rank_Metric, decreasing = TRUE), ]
    colnames(gene_table) <- c("Gene_Name", "Rank_Metric", "Is_Core_Enrichment")
    html_table <- knitr::kable(gene_table, format = "html", row.names = FALSE, digits = 3, table.attr = 'class="table table-striped table-sm"')

    # 4.3 绘制热图
    heat_html_tag <- "<p class='text-muted' style='margin-top:20px;'>No expression data matched for heatmap.</p>"
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
          z_mat[is.na(z_mat)] <- 0; z_mat[z_mat > 1] <- 1; z_mat[z_mat < -1] <- -1

          ann_col <- data.frame(Group = sample_meta_sub$group, row.names = rownames(sample_meta_sub))
          grp_col <- c("#E41A1C", "#377EB8"); names(grp_col) <- c(meta$left_group, meta$right_group)
          color_palette <- grDevices::colorRampPalette(c("#67a9cf", "#f7f7f7", "#ef8a62"))(100)

          pheatmap::pheatmap(
            z_mat, scale = "none", cluster_cols = FALSE, cluster_rows = FALSE, gaps_col = n_left,
            color = color_palette, breaks = seq(-1, 1, length.out = 101),
            annotation_col = ann_col, annotation_colors = list(Group = grp_col),
            display_numbers = plot_cpm_numbers, number_color = "black", fontsize_number = 8,
            cellwidth = 35, cellheight = 16, filename = file.path(details_dir, heat_png_name),
            main = sprintf("Row-Scaled Z-Score [-1, 1]\nEnriched in: %s", as.character(df_clean$Enriched_In[i]))
          )
          heat_html_tag <- sprintf("<img src='%s' class='img-fluid shadow-sm border'>", heat_png_name)
        }
      }
    }

    # 🌟 4.4 拼接 Dashboard HTML (引入 P-value 和 P-adj 双展示)
    html_content <- sprintf('
      <!DOCTYPE html><html><head><meta charset="utf-8"><title>Detail: %s</title>
      <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet"></head>
      <body class="container mt-4 mb-5"><h2 class="text-primary">%s</h2>
      <h6 class="text-muted mb-4">Original ID: %s</h6>
      <p class="lead"><strong>Rank:</strong> %d | <strong>NES:</strong> %.3f | <strong>P-value:</strong> %.4e | <strong>FDR (P-adj):</strong> <span class="badge bg-danger">%.4e</span></p><hr>
      <div class="row mb-5"><div class="col-md-5"><h4 class="mb-3">Enrichment Plot</h4><img src="%s" class="img-fluid border shadow-sm"></div>
      <div class="col-md-7"><h4 class="mb-3">Expression Heatmap</h4><div style="overflow:auto; max-height:800px;">%s</div></div></div>
      <div class="row"><div class="col-md-12"><h4 class="mb-3">Gene Statistics</h4><div style="max-height: 400px; overflow-y: auto;">%s</div></div></div>
      </body></html>
    ', clean_title, clean_title, pw_id, i, df_clean$NES[i], df_clean$pvalue[i], df_clean$p.adjust[i], gsea_png_name, heat_html_tag, html_table)

    writeLines(html_content, con = file.path(details_dir, detail_filename))
    detail_links[i] <- sprintf('<a href="./Details/%s" target="_blank" class="btn btn-sm btn-success" style="padding: 2px 10px; text-decoration: none;">📊 Dashboard</a>', detail_filename)
  }

  df_clean$Detail_Page <- detail_links

  # 🌐 5. 生成主交互式数据表 (加上 Pvalue 和 P.adjust)
  display_df <- df_clean %>%
    dplyr::select(Rank, Detail_Page, Pathway = Pathway_Link, Collection = Display_Collection, Enriched_In, Size = setSize, NES, pvalue, p.adjust, Description)

  dt_table <- DT::datatable(
    display_df, rownames = FALSE, escape = FALSE, filter = "top",
    caption = htmltools::tags$caption(style = 'caption-side: top; text-align: center; font-size: 150%; font-weight: bold;', sprintf("GSEA Report: %s vs %s", meta$left_group, meta$right_group)),
    extensions = c('Buttons', 'Scroller'),
    options = list(dom = 'Bfrtip', deferRender = TRUE, scrollY = 600, scroller = TRUE, pageLength = -1, buttons = c('copy', 'csv', 'excel'), autoWidth = TRUE)
  ) %>%
    DT::formatRound(columns = c('NES'), digits = 3) %>%
    DT::formatSignif(columns = c('pvalue', 'p.adjust'), digits = 4) %>%
    DT::formatStyle('Enriched_In', backgroundColor = DT::styleEqual(c(meta$left_group, meta$right_group), c('#fee0d2', '#deebf7'))) %>%
    DT::formatStyle('NES', color = DT::styleInterval(0, c('blue', 'red')), fontWeight = 'bold')

  old_wd <- getwd()
  setwd(bundle_dir)
  tryCatch({
    htmlwidgets::saveWidget(dt_table, file = basename(main_html_path), selfcontained = FALSE, libdir = "lib", title = sprintf("GSEA Report: %s vs %s [%s]", meta$left_group, meta$right_group, meta$geneset_name))
  }, finally = { setwd(old_wd) })

  message(sprintf("✅ 完美！HTML 报告已自动寻址并生成至: %s", bundle_dir))
  return(invisible(bundle_dir))
}




