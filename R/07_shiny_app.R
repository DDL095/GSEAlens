#' @title GSEAlens 交互式分析工作站 (新版架构)
#' @description 基于 GseaRes 对象的模块化 Shiny 应用，支持 Limma 和 DESeq2 双后端。
#' @param gsea_res GseaRes 对象 (由 batch_calc_gsea() 生成)
#' @param lang 语言设置，"zh" 或 "en"
#' @return Shiny 应用对象
#' @export
#' @importFrom shiny shinyApp fluidPage titlePanel sidebarLayout sidebarPanel mainPanel
#' @importFrom shiny tabsetPanel tabPanel br tags HTML icon
#' @importFrom shinycssloaders withSpinner
launch_gsea_app <- function(gsea_res, lang = "zh") {

  # 设置语言
  if (!is.null(lang)) {
    set_gsealens_lang(lang)
  }

  # 校验输入
  if (!inherits(gsea_res, "GseaRes")) {
    stop("❌ 请传入标准的 GseaRes 对象 (由 batch_calc_gsea() 生成)")
  }

  # 动态生成 UI
  ui <- shiny::fluidPage(
    shiny::tags$head(
      shiny::tags$style(shiny::HTML("
        .modal-dialog { max-width: 90vw !important; width: 90vw !important; }
        .modal-body { max-height: 85vh; overflow-y: auto; background-color: #fafafa; }
        .white-box {
          background-color: white; padding: 15px; border-radius: 8px;
          box-shadow: 0 2px 8px rgba(0,0,0,0.1); margin-bottom: 15px;
        }
        .description-cell {
          white-space: normal !important; word-wrap: break-word !important;
          max-width: 300px; font-size: 12px; line-height: 1.4;
        }
        .dt-center { text-align: center !important; }
        .well { background-color: #f8f9fa; border: 1px solid #e9ecef; border-radius: 8px; padding: 15px; }
      "))
    ),
    shiny::titlePanel(shiny::HTML("🧬 GSEAlens PRO 4.0: 双后端全息分析工作站")),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        width = 3,
        mod_data_prep_ui("data_prep")
      ),
      shiny::mainPanel(
        width = 9,
        shiny::tabsetPanel(
          id = "main_tabs",
          shiny::tabPanel(
            title = shiny::HTML("📊 主工作台"),
            shiny::br(),
            mod_master_table_ui("master_table"),
            shiny::hr(),
            mod_multi_plot_ui("multi_plot")
          ),
          shiny::tabPanel(
            title = shiny::HTML("🌋 全息四重联动"),
            shiny::br(),
            mod_quadrant_ui("quadrant")
          )
        )
      )
    ),
    # 隐藏的弹窗占位
    mod_pathway_modal_ui("pathway_modal")
  )

  # Server 逻辑
  # Server 逻辑
  # Server 逻辑
  server <- function(input, output, session) {

    # 1. 数据预处理模块（返回列表）
    data_prep_list <- mod_data_prep_server("data_prep", gsea_res)

    # 2. 主表格模块（传入 data reactive）
    table_result <- mod_master_table_server("master_table", data_prep_list$data)

    # 3. 联合绘图模块（传入 data reactive）
    mod_multi_plot_server("multi_plot", data_prep_list$data, table_result)

    # 4. 四重联动模块（传入完整列表）
    mod_quadrant_server("quadrant", data_prep_list, gsea_res)

    # 5. 详情弹窗模块（传入 data reactive）
    mod_pathway_modal_server("pathway_modal", data_prep_list$data, table_result$show_modal, gsea_res)
  }

  message("🚀 GSEAlens PRO 4.0 已启动")
  message(sprintf("   后端类型: %s", gsea_res$backend_info$backend))
  message(sprintf("   对比组数: %d", nrow(gsea_res$contrast_registry)))

  shiny::shinyApp(ui, server)
}
