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





