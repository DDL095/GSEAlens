#' @title GSEAlens 全局导入
#' @description 集中管理所有外部依赖导入，避免重复定义
#' @keywords internal
#' @name utils-imports
NULL


# 基础 R 包导入（不需要 @importFrom，因为已在 Depends/Imports 中）


#' @importFrom grDevices dev.off png colorRampPalette
#' @importFrom graphics plot text par lines points abline legend
#' @importFrom methods is new slot slotNames
NULL

# stats 包函数（显式导入常用函数）
#' @importFrom stats p.adjust terms var na.omit setNames
#' @importFrom stats sd median quantile reorder aggregate
NULL

# utils 包函数
#' @importFrom utils head tail packageVersion data
NULL


# Tidyverse 核心


#' @importFrom dplyr filter arrange mutate select left_join right_join
#' @importFrom dplyr distinct desc group_by summarise ungroup bind_rows
#' @importFrom dplyr rename pull slice slice_head slice_tail if_else coalesce
#' @importFrom tidyr pivot_longer pivot_wider separate unite drop_na
#' @importFrom tibble tibble as_tibble enframe deframe rownames_to_column
#' @importFrom stringr str_split str_trim str_detect str_replace str_wrap
#' @importFrom stringr str_replace_all str_to_title str_sub str_c
NULL


# 可视化


#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_col geom_boxplot
#' @importFrom ggplot2 geom_jitter geom_hline geom_vline geom_text annotate
#' @importFrom ggplot2 geom_tile geom_segment geom_rect
#' @importFrom ggplot2 scale_color_manual scale_fill_manual scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous scale_size_manual scale_alpha_manual
#' @importFrom ggplot2 labs theme element_text element_blank element_rect
#' @importFrom ggplot2 theme_bw theme_minimal coord_cartesian ggsave
#' @importFrom ggplot2 facet_wrap facet_grid
NULL


# Shiny 生态系统


#' @importFrom shiny shinyApp fluidPage titlePanel sidebarLayout sidebarPanel
#' @importFrom shiny mainPanel tabsetPanel tabPanel br tags HTML icon
#' @importFrom shiny selectInput textInput actionButton observeEvent reactiveVal
#' @importFrom shiny reactive renderPlot renderText showModal modalDialog
#' @importFrom shiny updateSelectInput updateSelectizeInput invalidateLater
#' @importFrom shiny isolate moduleServer NS tagList div span hr h4
#' @importFrom shiny req validate need showNotification
#' @importFrom shinycssloaders withSpinner
NULL

#' @importFrom DT datatable renderDataTable dataTableOutput formatStyle
#' @importFrom DT styleInterval styleEqual
NULL

#' @importFrom plotly plot_ly renderPlotly plotlyOutput layout add_annotations
#' @importFrom plotly event_data ggplotly style
NULL


# 生物信息学核心


#' @importFrom limma topTable eBayes contrasts.fit makeContrasts
NULL

#' @importFrom edgeR cpm rpkm DGEList calcNormFactors
NULL

#' @importFrom DESeq2 DESeq results resultsNames counts
#' @importFrom DESeq2 vst estimateSizeFactors
NULL

#' @importFrom SummarizedExperiment assay colData rowData
NULL

#' @importFrom S4Vectors metadata metadata<-
NULL


# GSEA 相关


# 找到这一行（约第98行）并删除：
# @importFrom clusterProfiler enrichplot  <-- 删除这一行

# 改为：
#' @importFrom clusterProfiler GSEA
NULL

# 如果您确实需要使用 enrichplot 包的功能，请单独导入：
#' @importFrom enrichplot dotplot cnetplot heatplot ridgeplot gseaplot upsetplot
NULL

# 如果需要使用基础 barplot，从 graphics 导入：
#' @importFrom graphics barplot
NULL

# GseaVis 不导入具体函数，使用 :: 调用避免冲突


# 热图与可视化


#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation draw rowAnnotation
NULL

#' @importFrom circlize colorRamp2
NULL

#' @importFrom grid gpar unit grid.text
NULL


# 并行计算与工具


#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
NULL

#' @importFrom htmltools tags HTML
NULL

#' @importFrom htmlwidgets saveWidget JS
NULL

#' @importFrom jsonlite fromJSON toJSON
NULL

#' @importFrom clipr write_clip
NULL

#' @importFrom msigdbr msigdbr msigdbr_collections
NULL

#' @importFrom BiocParallel bplapply bpparam
NULL
