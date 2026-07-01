#' @title GSEAlens Global Imports

#' @description Centralized management of all external dependency imports to avoid redundant definitions

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

#' @importFrom utils head tail packageVersion data combn read.csv

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





#' @importFrom clusterProfiler GSEA

NULL



#' @importFrom enrichplot dotplot cnetplot heatplot ridgeplot gseaplot upsetplot

NULL



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



#' @importFrom rlang .data

#' @importFrom rlang `%||%`

NULL





# ---- Native R utilities (base packages) -------------------------------------

# `tools` ships with base R but Bioconductor requires explicit importFrom

# for any `::` access in package code. See utils-accessors.R ::file_ext use.

#' @importFrom tools file_ext

NULL





# ---- RNG encapsulation ------------------------------------------------------

# withr::with_seed / local_seed replace set.seed() (reviewer comment #4).

# Declared in Imports (not Suggests) because R/ code calls them directly.

#' @importFrom withr with_seed local_seed

NULL





# ---- GSEA scoring helper ----------------------------------------------------

# enrichit::gseaScores is the GSEA running-score core used by vis-gsea-core.R.

#' @importFrom enrichit gseaScores

NULL





# ---- Network / graph --------------------------------------------------------

# igraph powers the Pathway Relation module and Hubgene Network layouts.

#' @importFrom igraph graph_from_data_frame layout_with_fr layout_with_kk

#' @importFrom igraph layout_in_circle V vcount ecount

NULL





# ---- Interactive network rendering ------------------------------------------

# visNetwork powers the Hubgene Network shiny module.

#' @importFrom visNetwork visNetwork visNetworkOutput renderVisNetwork

#' @importFrom visNetwork visPhysics visLayout visInteraction visOptions visEvents

NULL





# ---- Shiny UI helpers -------------------------------------------------------

# shinyjs powers enable/disable of the edge-detail control in Pathway Relation.

#' @importFrom shinyjs enable disable

NULL





# ---- Progress reporting -----------------------------------------------------

# progressr powers the progressor inside batch_calc_gsea.

#' @importFrom progressr handlers progressor

NULL





utils::globalVariables(c(

  ".", "gene_symbol", "stat", "abs_stat", "ID", "NES", "Collection",

  "Combo_Name", "Display_Collection", "Enriched_In", "Expression",

  "Group", "Sample", "Subcollection", "URL", "Description", "Description.y",

  "Detail_Page", "Pathway_Link", "Rank", "setSize", "pvalue", "gs_collection",

  "gs_description", "gs_name", "gs_subcollection", "gs_url", "phase",

  "read.csv", "relative_sec", "rss_mb", "x", "y"

))

