# ============================================================
# 文件: R/mod_ai_abs_page.R
# 功能: AI 提示词生成模块（支持自定义模板）
# ============================================================

#' @title AI Interpretation Page UI
#' @param id Module ID
#' @keywords internal

mod_ai_abs_page_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::fluidRow(
      # ===== 左侧控制栏 =====
      shiny::column(
        width = 3,
        shiny::div(
          class = "well",
          style = "padding: 15px;",

          shiny::h4("AI Prompt Generator"),
          shiny::hr(),

          # 当前对比组信息
          shiny::h5("Current Comparison"),
          shiny::uiOutput(ns("current_contrast_info")),

          shiny::hr(),

          # 选中通路统计
          shiny::h5("Selected Pathways"),
          shiny::uiOutput(ns("selection_stats")),

          shiny::hr(),

          # ===== 自定义模板区域 =====
          shiny::h5("Custom Template"),
          shiny::div(
            style = "background: #f0f8ff; padding: 10px; border-radius: 5px; margin-bottom: 10px;",
            shiny::checkboxInput(
              ns("use_custom_template"),
              label = "Use custom template",
              value = FALSE
            )
          ),

          shiny::conditionalPanel(
            condition = sprintf("input['%s'] == true", ns("use_custom_template")),
            shiny::div(
              style = "background: #fff; padding: 10px; border: 1px solid #ddd; border-radius: 5px;",
              shiny::p(shiny::strong("占位符说明:")),
              shiny::tags$small(
                "{left_group} - 左组名称", shiny::br(),
                "{right_group} - 右组名称", shiny::br(),
                "{comparison} - 对比组描述", shiny::br(),
                "{total} - 总通路数", shiny::br(),
                "{high} - 高置信度数量", shiny::br(),
                "{mod} - 中置信度数量", shiny::br(),
                "{low} - 低置信度数量", shiny::br(),
                "{table} - 通路数据表格", shiny::br(),
                shiny::hr(),
                "示例: {left_group} 在通路 X 中呈现激活趋势"
              ),
              shiny::hr(),
              shiny::actionButton(
                ns("import_template"),
                label = "Import Template",
                icon = shiny::icon("upload"),
                class = "btn-outline-primary btn-sm",
                style = "width: 100%;"
              ),
              shiny::actionButton(
                ns("export_template"),
                label = "Export Current",
                icon = shiny::icon("download"),
                class = "btn-outline-secondary btn-sm",
                style = "width: 100%; margin-top: 5px;"
              ),
              shiny::hr(),
              shiny::textAreaInput(
                ns("custom_template_text"),
                label = "Edit Template:",
                value = "",
                width = "100%",
                height = "300px"
              )
            )
          ),

          shiny::hr(),

          # 生成按钮
          shiny::actionButton(
            ns("generate_prompt"),
            label = "Generate Prompt",
            icon = shiny::icon("wand-magic-sparkles"),
            class = "btn-primary",
            style = "width: 100%; font-weight: bold;"
          ),

          shiny::hr(),

          # 复制按钮
          shiny::actionButton(
            ns("copy_to_clipboard"),
            label = "Copy to Clipboard",
            icon = shiny::icon("copy"),
            class = "btn-success",
            style = "width: 100%;"
          )
        )
      ),

      # ===== 右侧提示词输出 =====
      shiny::column(
        width = 9,
        shiny::div(
          class = "white-box",
          style = "min-height: 800px; padding: 15px;",

          shiny::h4("AI Prompt Output"),
          shiny::hr(),

          # 提示词输出区域
          shiny::div(
            style = "background: #f8f9fa; padding: 15px; border-radius: 5px;",
            shiny::textAreaInput(
              ns("prompt_output"),
              label = NULL,
              value = "Click 'Generate Prompt' button to create AI interpretation prompt.",
              width = "100%",
              height = "700px"
            )
          )
        )
      )
    )
  )
}


#' @title AI Interpretation Page Server
#' @param id Module ID
#' @param gsea_res GseaRes object
#' @param data_prep_list Data preprocessing list
#' @param table_controller Table controller
#' @keywords internal

mod_ai_abs_page_server <- function(id, gsea_res, data_prep_list, table_controller) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ============================================================
    # 默认提示词模板（中文优化版）
    # ============================================================

    default_template <- 'Role Setting: GSEA Enrichment Direction Interpretation Expert

language: English
description: You are an expert in interpreting GSEA results. You have removed the sign from the Normalized Enrichment Score (NES) and focus solely on the absolute value (|NES|) for analysis. Your core task is to objectively describe the enrichment trend of gene sets in the predefined {left_group} group or {right_group} group based on the degree of enrichment, allowing the use of terms such as "upregulated" and "activated", while ensuring all analyses are based on absolute values.

### Gene Set Names and Their Analysis
In the analysis, in addition to using absolute values, please pay attention to interpreting the names of gene sets. For example, the gene set name "FLORIO_NEOCORTEX_BASAL_RADIAL_GLIA_DN" should be understood as this gene set being downregulated in a certain group. Therefore, analysis should combine the meaning of the gene set name with its absolute value.

Core Principles [Must Be Strictly Followed]
- The absolute value of the gene set (|NES|) directly indicates its enrichment intensity in the ranked gene list; the larger the absolute value, the more pronounced the enrichment trend of that gene set.
- When using terms such as "upregulated" or "activated" in the analysis, it must be stated that the sign has been removed and only the absolute value is being used for description.

Expression Rules (Golden Rule):
- The use of expressions such as "inhibited" or "decreased" is strictly prohibited.
- The use of terms such as "activated" or "upregulated" is permitted, based on the absolute value (|NES|) of each gene set to describe the enrichment trend, combined with analysis of the meaning of the gene set name.

Mandatory Objective Description Templates:
- "This gene set shows an activation trend in the [{left_group} group] (or: the members of this gene set have overall higher expression levels in the [{left_group} group], with |NES| = {abs_nes})."
- "This gene set shows an activation trend in the [{right_group} group] (or: the members of this gene set have overall higher expression levels in the [{right_group} group], with |NES| = {abs_nes})."

Constraints and Correction Mechanisms:
- In any interpretation, if language that is hypothetical or inferential is used, it must be stopped immediately and strictly traced back to the Golden Rule for expression correction.
- All descriptions must be based on absolute values and ensure accurate description to avoid ambiguity.

Workflow
1. Confirm Background (Prerequisites):
   - Group A = {left_group} group
   - Group B = {right_group} group
2. Assess Significance: Use the FDR q-val (typically < 0.25) indicator to determine whether the enrichment direction is statistically significant, and judge the strength of enrichment based on the absolute value of |NES|.
3. Classification and Description:
   - Strictly classify and describe significant enrichment results (e.g., FDR < 0.25) as follows:
     - Enriched in {left_group}: Classified as "gene sets activated in the {left_group} group".
     - Enriched in {right_group}: Classified as "gene sets activated in the {right_group} group".
4. Cautious Interpretation: After completing the objective directional description, analyze the biological significance of the two enrichment directions and provide constructive insights and explanations,and provide a possible functional interpretation of gene pathways in combination with specific leadingedge genes.

Initialization
As your GSEA enrichment direction interpretation expert, I am ready. Please analyze the following GSEA results.

### Enrichment Analysis Data

**Comparison Group**: {comparison}

**Statistical Summary**: Total {total} pathways | High confidence: {high} | Medium confidence: {mod} | Low confidence: {low}

| # | Pathway ID | |NES| | Enrichment Direction | FDR | Description | Leading Edge Genes |
|:--:|:-----------|:----:|:-----------:|:------:|:--------------|:-------------------|
{table}

**Field Descriptions**:
- **|NES|**: Absolute Normalized Enrichment Score. |NES| >= 1.5 indicates significant enrichment; |NES| >= 2.0 indicates strong enrichment
- **Enrichment Direction**: Enrichment in {left_group} indicates activation in the {left_group} group; enrichment in {right_group} indicates activation in the {right_group} group
- **FDR**: Multiple testing corrected P-value. FDR < 0.25 is the MSigDB standard threshold
- **Leading Edge Genes**: Core contributing genes, key to understanding regulatory mechanisms'

    # 初始化自定义模板为空
    shiny::observe({
      shiny::updateTextAreaInput(
        session,
        "custom_template_text",
        value = default_template
      )
    })

    # ============================================================
    # 辅助函数
    # ============================================================

    #' 提取 Leading Edge 基因并格式化为逗号分隔
    extract_leading_genes <- function(core_str, max_genes = 15) {
      if (is.null(core_str) || is.na(core_str) || core_str == "") {
        return(list(genes = "N/A", count = 0))
      }
      genes <- unlist(strsplit(as.character(core_str), "/"))
      genes <- trimws(genes)
      genes <- genes[genes != ""]
      count <- length(genes)
      if (count > max_genes) {
        display <- paste(paste(genes[1:max_genes], collapse = ", "),
                         paste0("(+", count - max_genes, " more)"),
                         sep = " ")
      } else {
        display <- paste(genes, collapse = ", ")
      }
      return(list(genes = display, count = count))
    }

    #' 判断置信度
    get_confidence <- function(nes, fdr) {
      abs_nes <- abs(nes)
      if (abs_nes >= 1.5 && fdr < 0.05) {
        return("High confidence")
      } else if (abs_nes >= 1.0 && fdr < 0.25) {
        return("Medium confidence")
      } else {
        return("Low confidence")
      }
    }

    # ============================================================
    # 当前对比组信息
    # ============================================================

    output$current_contrast_info <- shiny::renderUI({
      data_list <- data_prep_list$data()
      shiny::req(data_list)

      left <- data_list$left_group
      right <- data_list$right_group

      shiny::tagList(
        shiny::strong("Left: ", style = "color: #E41A1C;"),
        shiny::code(left), shiny::br(),
        shiny::strong("Right: ", style = "color: #377EB8;"),
        shiny::code(right)
      )
    })

    # ============================================================
    # 选中通路统计
    # ============================================================

    output$selection_stats <- shiny::renderUI({
      shiny::req(table_controller)
      selected_ids <- table_controller$selected_pathways()
      n_selected <- length(selected_ids)

      if (n_selected == 0) {
        shiny::div(
          style = "background: #fff3cd; padding: 10px; border-radius: 5px; color: #856404;",
          shiny::icon("exclamation-triangle"),
          " No pathways selected",
          shiny::br(),
          shiny::tags$small("Check 'Joint Plot' in main table")
        )
      } else {
        shiny::div(
          style = "background: #d4edda; padding: 10px; border-radius: 5px; color: #155724;",
          shiny::icon("check-circle"),
          sprintf(" %d pathway(s) selected", n_selected)
        )
      }
    })

    # ============================================================
    # 导入模板
    # ============================================================

    shiny::observeEvent(input$import_template, {
      shiny::showModal(shiny::modalDialog(
        title = "Import Custom Template",
        size = "l",
        easyClose = TRUE,
        shiny::textAreaInput(
          ns("import_text"),
          label = "Paste your template here:",
          value = "",
          width = "100%",
          height = "400px"
        ),
        footer = shiny::tagList(
          shiny::actionButton(ns("confirm_import"), "Confirm Import", class = "btn-primary"),
          shiny::modalButton("Cancel")
        )
      ))
    })

    shiny::observeEvent(input$confirm_import, {
      imported_text <- input$import_text
      if (!is.null(imported_text) && imported_text != "") {
        shiny::updateTextAreaInput(
          session,
          "custom_template_text",
          value = imported_text
        )
        shiny::showNotification("Template imported successfully!", type = "message", duration = 3)
      }
      shiny::removeModal()
    })

    # ============================================================
    # 导出模板
    # ============================================================

    shiny::observeEvent(input$export_template, {
      current_template <- input$custom_template_text
      if (is.null(current_template) || current_template == "") {
        current_template <- default_template
      }

      # 保存到文件
      output_path <- file.path(getwd(), "gsea_prompt_template.txt")
      writeLines(current_template, output_path, useBytes = TRUE)

      shiny::showNotification(
        sprintf("Template exported to: %s", output_path),
        type = "message", duration = 5
      )
    })

    # ============================================================
    # 生成提示词
    # ============================================================

    shiny::observeEvent(input$generate_prompt, {
      shiny::req(data_prep_list$data())

      data_list <- data_prep_list$data()

      # 获取选中的 pathways
      selected_ids <- table_controller$selected_pathways()

      # 如果没有选中，使用 Top 15
      if (length(selected_ids) == 0) {
        selected_ids <- data_list$df$ID[1:min(15, nrow(data_list$df))]
        shiny::showNotification(
          "No pathways selected - using Top 15 by |NES|",
          type = "warning", duration = 3
        )
      }

      # 筛选数据
      df_subset <- data_list$df[data_list$df$ID %in% selected_ids, ]
      df_subset <- df_subset[order(-abs(df_subset$NES)), ]

      # 限制数量
      max_pathways <- 20
      if (nrow(df_subset) > max_pathways) {
        df_subset <- df_subset[1:max_pathways, ]
      }

      # 基本信息
      left_group <- data_list$left_group
      right_group <- data_list$right_group
      comparison <- paste0(left_group, " vs ", right_group)

      # 统计
      n_total <- nrow(df_subset)
      n_high <- sum(sapply(seq_len(n_total), function(i) {
        get_confidence(df_subset$NES[i], df_subset$p.adjust[i]) == "高置信度"
      }))
      n_mod <- sum(sapply(seq_len(n_total), function(i) {
        get_confidence(df_subset$NES[i], df_subset$p.adjust[i]) == "中置信度"
      }))
      n_low <- n_total - n_high - n_mod

      # ============================================================
      # 构建通路表格行
      # ============================================================

      table_rows <- c()
      for (i in seq_len(nrow(df_subset))) {
        row <- df_subset[i, ]

        # 提取 Leading Edge 基因
        leading <- extract_leading_genes(row$core_enrichment)

        # 确定富集方向
        enriched_in <- ifelse(row$NES > 0, left_group, right_group)

        # 置信度
        conf <- get_confidence(row$NES, row$p.adjust)

        # 截断 Description
        desc <- if (!is.null(row$Description) && !is.na(row$Description)) {
          substr(as.character(row$Description), 1, 9999)
        } else {
          row$ID
        }

        # 构建表格行
        table_rows <- c(table_rows, sprintf(
          "| %d | `%s` | %.2f | %s | %.2e | %s | %s [%s] |",
          i,
          row$ID,
          abs(row$NES),
          enriched_in,
          row$p.adjust,
          substr(desc, 1, 9999),
          leading$genes,
          conf
        ))
      }

      table_body <- paste(table_rows, collapse = "\n")

      # ============================================================
      # 获取模板（自定义或默认）
      # ============================================================

      if (isTRUE(input$use_custom_template) &&
          !is.null(input$custom_template_text) &&
          input$custom_template_text != "") {
        template_text <- input$custom_template_text
      } else {
        template_text <- default_template
      }

      # ============================================================
      # 替换占位符
      # ============================================================

      prompt_text <- template_text

      # 基础占位符替换
      prompt_text <- gsub("\\{left_group\\}", left_group, prompt_text, fixed = FALSE)
      prompt_text <- gsub("\\{right_group\\}", right_group, prompt_text, fixed = FALSE)
      prompt_text <- gsub("\\{comparison\\}", comparison, prompt_text, fixed = FALSE)
      prompt_text <- gsub("\\{total\\}", as.character(n_total), prompt_text, fixed = FALSE)
      prompt_text <- gsub("\\{high\\}", as.character(n_high), prompt_text, fixed = FALSE)
      prompt_text <- gsub("\\{mod\\}", as.character(n_mod), prompt_text, fixed = FALSE)
      prompt_text <- gsub("\\{low\\}", as.character(n_low), prompt_text, fixed = FALSE)
      prompt_text <- gsub("\\{table\\}", table_body, prompt_text, fixed = FALSE)

      # 特殊占位符：abs_nes（用于描述模板中）
      # 这里需要逐行处理，可以扩展更多占位符
      prompt_text <- gsub("\\{abs_nes\\}", sprintf("%.2f", abs(df_subset$NES[1])), prompt_text, fixed = FALSE)

      # 更新输出
      shiny::updateTextAreaInput(
        session,
        "prompt_output",
        value = prompt_text
      )

      shiny::showNotification(
        sprintf("Prompt generated: %d pathways", n_total),
        type = "message", duration = 3
      )
    })

    # ============================================================
    # 复制到剪贴板
    # ============================================================

    shiny::observeEvent(input$copy_to_clipboard, {
      prompt <- input$prompt_output

      if (is.null(prompt) || prompt == "") {
        shiny::showNotification("No prompt to copy!", type = "warning", duration = 3)
        return()
      }

      # 使用 clipr 或系统剪贴板
      if (requireNamespace("clipr", quietly = TRUE)) {
        clipr::write_clip(prompt)
        shiny::showNotification("Copied!", type = "message", duration = 2)
      } else {
        writeLines(prompt, con <- file("clipboard", "w"))
        close(con)
        shiny::showNotification("Copied!", type = "message", duration = 2)
      }
    })

  })
}
