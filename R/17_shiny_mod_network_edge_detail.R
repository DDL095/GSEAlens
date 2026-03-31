#' @title Network Edge Detail Modal Module
#' @description Modal for displaying detailed edge information including shared genes
#'   and multiple similarity metrics (Jaccard, Overlap, Dice).
#' @keywords internal

#' @title Network Edge Detail Modal UI
#' @keywords internal
mod_network_edge_detail_ui <- function(id) {
  ns <- shiny::NS(id)
  NULL  # Modal dynamically generated in server
}


#' @title Network Edge Detail Modal Server
#' @param id Module ID
#' @param parent_session Parent module session for receiving edge click events
#' @keywords internal
mod_network_edge_detail_server <- function(id, parent_session) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ============================================================
    # 监听来自主模块的边点击事件
    # ============================================================

    shiny::observeEvent(parent_session$input[["network_edge_clicked"]], {
      message_data <- parent_session$input[["network_edge_clicked"]]

      if (is.null(message_data) || !is.list(message_data)) return()

      # 显示 Modal
      show_edge_detail_modal(message_data)
    })

    # ============================================================
    # Modal 显示函数
    # ============================================================

    show_edge_detail_modal <- function(edge_info) {

      if (!is.list(edge_info) || is.null(edge_info$from) || is.null(edge_info$to)) {
        return()
      }

      from_pw <- edge_info$from
      to_pw <- edge_info$to
      shared_count <- edge_info$shared
      jaccard <- edge_info$jaccard
      shared_genes <- edge_info$shared_genes

      # 计算其他相似度指标
      if (!is.null(edge_info$pathway_a_genes) && !is.null(edge_info$pathway_b_genes)) {
        genes_a <- toupper(edge_info$pathway_a_genes)
        genes_b <- toupper(edge_info$pathway_b_genes)

        overlap_coef <- shared_count / min(length(genes_a), length(genes_b))
        dice_coef <- (2 * shared_count) / (length(genes_a) + length(genes_b))
      } else {
        overlap_coef <- NA
        dice_coef <- NA
      }

      # 格式化共享基因
      if (is.null(shared_genes) || length(shared_genes) == 0) {
        gene_display <- "No shared genes found"
      } else {
        gene_tags <- sapply(shared_genes, function(g) {
          sprintf('<span style="display:inline-block; background:#e3f2fd; padding:4px 8px; margin:2px; border-radius:4px; font-size:12px;">%s</span>', g)
        })
        gene_display <- paste(gene_tags, collapse = "")
      }

      # Modal 内容
      modal_content <- shiny::div(
        style = "padding: 10px;",

        # 相似度指标卡片
        shiny::div(
          style = "background: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 15px;",
          shiny::h4("Similarity Metrics", style = "margin-top: 0;"),
          shiny::fluidRow(
            shiny::column(3,
                          shiny::div(style = "text-align: center;",
                                     shiny::h3(sprintf("%.3f", jaccard), style = "color: #007bff; margin: 0;"),
                                     shiny::small("Jaccard Index"))
            ),
            shiny::column(3,
                          shiny::div(style = "text-align: center;",
                                     shiny::h3(sprintf("%.3f", overlap_coef), style = "color: #28a745; margin: 0;"),
                                     shiny::small("Overlap Coef"))
            ),
            shiny::column(3,
                          shiny::div(style = "text-align: center;",
                                     shiny::h3(sprintf("%.3f", dice_coef), style = "color: #dc3545; margin: 0;"),
                                     shiny::small("Dice Coeff"))
            ),
            shiny::column(3,
                          shiny::div(style = "text-align: center;",
                                     shiny::h3(sprintf("%d", shared_count), style = "color: #6c757d; margin: 0;"),
                                     shiny::small("Shared Genes"))
            )
          )
        ),

        # 共享基因列表
        shiny::div(
          style = "background: #fff; padding: 15px; border-radius: 8px; border: 1px solid #dee2e6;",
          shiny::h4(sprintf("Shared Core Genes (%d)", shared_count), style = "margin-top: 0; color: #333;"),
          shiny::HTML(gene_display)
        ),

        # 导出按钮
        shiny::div(
          style = "margin-top: 15px; text-align: right;",
          shiny::downloadButton(ns("export_shared_genes"),
                                label = "Export Gene List (CSV)",
                                class = "btn-sm btn-primary")
        )
      )

      # 显示 Modal
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML(sprintf(
          '<span style="font-size: 18px;">Edge Detail</span><br>
           <small style="color: #666;">%s &harr; %s</small>',
          from_pw, to_pw
        )),
        size = "l",
        easyClose = TRUE,
        footer = shiny::modalButton("Close"),
        modal_content
      ))
    }

    # ============================================================
    # 导出功能
    # ============================================================

    output$export_shared_genes <- shiny::downloadHandler(
      filename = function() {
        "shared_genes.csv"
      },
      content = function(file) {
        # 获取当前 Modal 的数据
        message_data <- parent_session$input[["network_edge_clicked"]]

        if (!is.null(message_data) && !is.null(message_data$shared_genes)) {
          genes_df <- data.frame(
            Gene = message_data$shared_genes,
            stringsAsFactors = FALSE
          )
          readr::write_csv(genes_df, file)
        }
      }
    )

  })
}
