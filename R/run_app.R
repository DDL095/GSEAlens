#' @title GSEAlens PRO 模块化重构版 (Shiny Modules架构)
#' @description 修复core_genes作用域、Description换行、基因名大小写匹配、高亮置顶等问题
#' @importFrom shiny moduleServer NS tagList tags HTML actionButton fluidRow column
#' @importFrom DT datatable renderDataTable dataTableOutput formatStyle styleInterval
#' @importFrom plotly plot_ly renderPlotly plotlyOutput layout add_annotations event_data
#' @importFrom dplyr filter arrange mutate select left_join case_when
#' @importFrom stringr str_wrap



# ==================== 模块1: 数据预处理引擎 (修订版 - 移除HTML生成) ====================
#' @title 数据预处理UI模块
#' @description 集成对比组选择、联合绘图控制、表达量类型选择及数据切片功能
#' @param id 模块命名空间ID
#' @return Shiny UI元素列表
modDataPrepUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    # 对比组选择（核心控制）
    shiny::selectInput(
      inputId = ns("selected_contrast"),
      label = "⚖️ 选择对比组 (Contrast):",
      choices = NULL  # 在server中动态填充
    ),

    shiny::hr(),

    # 联合绘图控制台
    shiny::h4("🎨 联合绘图控制"),
    shiny::selectInput(
      inputId = ns("plot_subtype"),
      label = "GSEAvis 样式:",
      choices = c("1: 仅经典富集" = "1", "2: 富集+热图带" = "2", "3: 完整带Rank" = "3"),
      selected = "3"
    ),
    shiny::textInput(
      inputId = ns("custom_colors"),
      label = "多通路自定义颜色 (逗号分隔):",
      value = "#E41A1C, #377EB8, #4DAF4A, #984EA3",
      placeholder = "例如: #FF0000, #00FF00, #0000FF"
    ),

    shiny::hr(),

    # 表达量类型选择器
    shiny::h4("📊 表达量度量"),
    shiny::selectInput(
      inputId = ns("expression_type"),
      label = "选择表达数据类型:",
      choices = c(
        "log2(CPM)" = "log2cpm",
        "CPM (原始)" = "cpm",
        "log2(FPKM)" = "log2fpkm",
        "FPKM (原始)" = "fpkm",
        "log2(TPM)" = "log2tpm",
        "TPM (原始)" = "tpm"
      ),
      selected = "log2cpm"
    ),

    shiny::hr(),

    # 数据切片与排序
    shiny::h4("🎯 数据切片与排序"),
    shiny::selectizeInput(
      inputId = ns("selected_collections"),
      label = "选择基因集亚组 (支持多选):",
      choices = NULL,
      multiple = TRUE,
      options = list(
        plugins = list("remove_button"),
        placeholder = "选择亚组或保留ALL"
      )
    ),
    shiny::selectInput(
      inputId = ns("sort_by"),
      label = "全局排序策略:",
      choices = c(
        "按 NES (降序)" = "nes_desc",
        "按 NES (升序)" = "nes_asc",
        "按 NES 绝对值 (降序)" = "abs_nes_desc",
        "按 P-value (升序)" = "pval_asc",
        "按 FDR (升序)" = "fdr_asc"
      ),
      selected = "abs_nes_desc"
    ),

    shiny::actionButton(
      inputId = ns("run_btn"),
      label = "🚀 确认配置 / 更新工作台",
      class = "btn-success",
      style = "width: 100%; font-weight: bold; margin-top: 15px;"
    )
  )
}

#' @title 数据预处理Server模块 (修订版)
#' @description 处理数据流，动态匹配limma系数，准备联合绘图参数
#' @param id 模块ID
#' @param res_capsule GSEA结果胶囊对象
#' @return 响应式数据列表，包含处理后的数据框及配置参数（不再包含HTML按钮列）

#' @title 数据预处理Server模块 (修订版 - 支持双向对比匹配)
#' @description 处理数据流，动态匹配limma系数（支持正向/反向对比），准备联合绘图参数
#' @param id 模块ID
#' @param res_capsule GSEA结果胶囊对象
#' @return 响应式数据列表，包含处理后的数据框、配置参数及系数方向标志

modDataPrepServer <- function(id, res_capsule) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # 存储处理结果
    result_data <- shiny::reactiveVal(NULL)
    # 标志位：确保只自动初始化一次
    has_initialized <- shiny::reactiveVal(FALSE)

    # 动态填充对比组选择器
    shiny::observe({
      shiny::req(res_capsule$results)
      contrast_choices <- names(res_capsule$results)
      if (length(contrast_choices) == 0) return()

      shiny::updateSelectInput(
        session = session,
        inputId = "selected_contrast",
        choices = contrast_choices,
        selected = contrast_choices[1]
      )
    })

    # 核心数据处理函数（参数化版本 - 支持双向系数匹配）
    process_data_core <- function(selected_contrast_val,
                                  selected_collections_val,
                                  sort_by_val = "abs_nes_desc",
                                  expression_type_val = "log2cpm",
                                  plot_subtype_val = "3",
                                  custom_colors_val = "#E41A1C, #377EB8, #4DAF4A, #984EA3",
                                  is_auto_init = FALSE) {

      if (is.null(selected_contrast_val)) {
        message("[DEBUG] process_data_core: selected_contrast_val is NULL")
        return(NULL)
      }

      task_info <- res_capsule$results[[selected_contrast_val]]
      if (is.null(task_info$data)) {
        message("[DEBUG] process_data_core: task_info$data is NULL")
        return(NULL)
      }

      # 解析对比组名称
      contrast_parts <- strsplit(selected_contrast_val, "_vs_")[[1]]
      left_group <- contrast_parts[1]
      right_group <- if (length(contrast_parts) >= 2) contrast_parts[2] else "Background"

      # ========== 关键修复：动态匹配limma系数名称（支持双向对比） ==========
      target_coef_pattern <- gsub("_vs_", " - ", selected_contrast_val)  # 正向：A_vs_B -> A - B
      reverse_coef_pattern <- paste(right_group, "-", left_group)        # 反向：A_vs_B -> B - A
      coef_name <- NULL
      coef_reversed <- FALSE  # 新增标志：是否使用反向系数（需要反转logFC）

      if (!is.null(res_capsule$limma_fit)) {
        fit_coefficients <- res_capsule$limma_fit$coefficients
        available_coefs <- colnames(fit_coefficients)

        # 1. 尝试正向精确匹配（A_vs_B -> A - B）
        exact_match <- which(available_coefs == target_coef_pattern)
        if (length(exact_match) > 0) {
          coef_name <- available_coefs[exact_match[1]]
          message(sprintf("[DEBUG] Matched coefficient (forward): %s", coef_name))
        } else {
          # 2. 尝试正向模糊匹配（去除空格后比较）
          normalized_target <- gsub("\\s+", "", toupper(target_coef_pattern))
          normalized_coefs <- gsub("\\s+", "", toupper(available_coefs))
          fuzzy_match <- which(normalized_coefs == normalized_target)

          if (length(fuzzy_match) > 0) {
            coef_name <- available_coefs[fuzzy_match[1]]
            message(sprintf("[DEBUG] Matched coefficient (forward, fuzzy): %s", coef_name))
          } else {
            # 3. 关键修复：尝试反向精确匹配（A_vs_B -> B - A）
            exact_match_rev <- which(available_coefs == reverse_coef_pattern)
            if (length(exact_match_rev) > 0) {
              coef_name <- available_coefs[exact_match_rev[1]]
              coef_reversed <- TRUE
              message(sprintf("[DEBUG] Matched coefficient (reversed): %s, will flip logFC", coef_name))
            } else {
              # 4. 尝试反向模糊匹配
              normalized_target_rev <- gsub("\\s+", "", toupper(reverse_coef_pattern))
              fuzzy_match_rev <- which(normalized_coefs == normalized_target_rev)
              if (length(fuzzy_match_rev) > 0) {
                coef_name <- available_coefs[fuzzy_match_rev[1]]
                coef_reversed <- TRUE
                message(sprintf("[DEBUG] Matched coefficient (reversed, fuzzy): %s, will flip logFC", coef_name))
              } else if (length(available_coefs) > 0) {
                # 最终回退：使用第一个可用系数（并警告）
                coef_name <- available_coefs[1]
                warning(sprintf("无法匹配对比组 [%s] 的系数，回退到第一个可用系数 [%s]",
                                selected_contrast_val, coef_name))
              }
            }
          }
        }
      }

      # 数据切片与清洗
      raw_df <- as.data.frame(task_info$data)
      meta_dict <- res_capsule$geneset_info$meta_dict

      # 关键修复：兼容URL列缺失的情况
      if (!"URL" %in% colnames(raw_df)) {
        raw_df$URL <- ""
        message("[DEBUG] URL column missing, created empty column")
      }

      if (!is.null(meta_dict) && "Description" %in% colnames(meta_dict)) {
        desc_map <- stats::setNames(meta_dict$Description, meta_dict$ID)
        raw_df$Brief_Description <- desc_map[raw_df$ID]
        raw_df$Brief_Description[is.na(raw_df$Brief_Description)] <- raw_df$ID[is.na(raw_df$Brief_Description)]
      } else {
        raw_df$Brief_Description <- raw_df$ID
      }

      raw_df$Brief_Description <- vapply(
        raw_df$Brief_Description,
        function(desc_text) stringr::str_wrap(as.character(desc_text), width = 200),
        character(1)
      )

      # 预留blank1-10
      for (blank_idx in 1:10) {
        raw_df[[sprintf("blank%d", blank_idx)]] <- NA_character_
      }

      # 子集过滤逻辑
      effective_collections <- selected_collections_val
      if (is.null(effective_collections) || length(effective_collections) == 0) {
        effective_collections <- "ALL"
        message("[DEBUG] No collections selected, defaulting to ALL")
      }

      if (!("ALL" %in% effective_collections)) {
        selected_chars <- as.character(effective_collections)
        display_col <- if ("Combo_Name" %in% colnames(raw_df)) "Combo_Name" else "Collection"

        message(sprintf("[DEBUG] Filtering by collections: %s (display_col: %s)",
                        paste(selected_chars, collapse=", "), display_col))

        collection_mask <- Reduce("|", lapply(selected_chars, function(sel_pattern) {
          if (display_col %in% colnames(raw_df)) {
            col_vals <- as.character(raw_df[[display_col]])
            # 关键修复：使用startsWith避免正则转义问题
            if (grepl(":", sel_pattern, fixed = TRUE)) {
              startsWith(col_vals, sel_pattern)
            } else {
              startsWith(col_vals, sel_pattern)
            }
          } else {
            rep(FALSE, nrow(raw_df))
          }
        }))

        if (any(collection_mask)) {
          raw_df <- raw_df[collection_mask, ]
          message(sprintf("[DEBUG] Filtered to %d rows", nrow(raw_df)))
        } else {
          warning(sprintf("Collection filter returned 0 rows for: %s",
                          paste(selected_chars, collapse=", ")))
          return(NULL)
        }
      }

      if (nrow(raw_df) == 0) {
        message("[DEBUG] raw_df has 0 rows after filtering")
        return(NULL)
      }

      # FDR重算
      raw_df$p.adjust <- stats::p.adjust(raw_df$pvalue, method = "BH")

      # 排序逻辑
      sort_config <- switch(sort_by_val,
                            "nes_desc" = list(col = "NES", desc = TRUE),
                            "nes_asc" = list(col = "NES", desc = FALSE),
                            "abs_nes_desc" = list(col = "abs_NES", desc = TRUE),
                            "pval_asc" = list(col = "pvalue", desc = FALSE),
                            "fdr_asc" = list(col = "p.adjust", desc = FALSE),
                            list(col = "abs_NES", desc = TRUE))

      if (sort_config$col == "abs_NES") {
        raw_df$abs_NES <- abs(raw_df$NES)
      }

      sort_order <- order(raw_df[[sort_config$col]], decreasing = sort_config$desc)
      raw_df <- raw_df[sort_order, ]
      raw_df$Rank <- seq_len(nrow(raw_df))

      # 构建显示列 - 移除HTML生成，仅保留数据字段
      raw_df$Safe_ID <- gsub("'", "\\\\'", raw_df$ID)

      raw_df$Pathway_Link <- ifelse(
        !is.na(raw_df$URL) & raw_df$URL != "",
        sprintf('<a href="%s" target="_blank">%s</a>', raw_df$URL, raw_df$ID),
        sprintf("<b>%s</b>", raw_df$ID)
      )

      raw_df$Enriched_In <- ifelse(raw_df$NES > 0, left_group, right_group)

      # 返回纯净数据结构（新增 coef_reversed 标志）
      list(
        df = raw_df,
        gsea_res = task_info$data,
        left_group = left_group,
        right_group = right_group,
        contrast_name = selected_contrast_val,
        coef_name = coef_name,
        coef_reversed = coef_reversed,  # 新增：系数方向标志
        expression_type = expression_type_val,
        plot_subtype = as.numeric(plot_subtype_val),
        custom_colors = custom_colors_val,
        is_preview = is_auto_init
      )
    }

    # 动态填充基因集子集选择器
    shiny::observe({
      shiny::req(input$selected_contrast)
      task_info <- res_capsule$results[[input$selected_contrast]]
      if (is.null(task_info$data)) return()

      gsea_df <- as.data.frame(task_info$data)
      available_collections <- unique(c(
        if ("Collection" %in% colnames(gsea_df)) gsea_df$Collection,
        if ("Combo_Name" %in% colnames(gsea_df)) gsea_df$Combo_Name
      ))
      available_collections <- available_collections[!is.na(available_collections)]
      available_collections <- setdiff(available_collections, "ALL")

      if (length(available_collections) == 0) return()

      shiny::updateSelectizeInput(
        session = session,
        inputId = "selected_collections",
        choices = c("ALL", sort(available_collections)),
        selected = if (is.null(input$selected_collections)) character(0) else input$selected_collections
      )
    })

    # ==================== 自动初始化逻辑（保留） ====================
    init_observer <- shiny::observeEvent(input$selected_contrast, {
      # 如果已初始化，销毁自身并返回
      if (has_initialized()) {
        init_observer$destroy()
        return()
      }

      # 关键修复：使用 observe + invalidateLater 实现延迟，但移除 once = TRUE
      # 改用局部标志位确保单次执行
      delay_done <- shiny::reactiveVal(FALSE)

      shiny::observe({
        # 延迟100ms确保UI就绪
        shiny::invalidateLater(100, session)

        # 双重检查防止重复执行
        if (delay_done() || has_initialized()) return()

        delay_done(TRUE)

        # 在 isolate 中执行初始化逻辑
        shiny::isolate({
          selected_contrast <- input$selected_contrast
          if (is.null(selected_contrast)) return()

          task_info <- res_capsule$results[[selected_contrast]]
          if (is.null(task_info$data)) return()

          gsea_df <- as.data.frame(task_info$data)

          # 计算可用集合
          available_collections <- unique(c(
            if ("Collection" %in% colnames(gsea_df)) gsea_df$Collection,
            if ("Combo_Name" %in% colnames(gsea_df)) gsea_df$Combo_Name
          ))
          available_collections <- available_collections[!is.na(available_collections)]
          available_collections <- setdiff(available_collections, "ALL")

          if (length(available_collections) == 0) return()

          # 计算最小集合（使用startsWith避免正则问题）
          collection_counts <- sapply(available_collections, function(col_name) {
            if ("Combo_Name" %in% colnames(gsea_df)) {
              sum(startsWith(as.character(gsea_df$Combo_Name), col_name))
            } else {
              sum(startsWith(as.character(gsea_df$Collection), col_name))
            }
          })

          min_idx <- which.min(collection_counts)
          if (length(min_idx) == 0 || is.na(min_idx)) return()

          min_collection <- available_collections[min_idx]
          min_count <- collection_counts[min_idx]

          message(sprintf("[AUTO-INIT] Selected min collection: %s (%d pathways)",
                          min_collection, min_count))

          # 更新UI选择
          shiny::updateSelectizeInput(
            session = session,
            inputId = "selected_collections",
            choices = c("ALL", sort(available_collections)),
            selected = min_collection
          )

          # 执行数据处理
          auto_data <- process_data_core(
            selected_contrast_val = selected_contrast,
            selected_collections_val = min_collection,
            sort_by_val = "abs_nes_desc",
            expression_type_val = "log2cpm",
            plot_subtype_val = "3",
            custom_colors_val = "#E41A1C, #377EB8, #4DAF4A, #984EA3",
            is_auto_init = TRUE
          )

          if (!is.null(auto_data)) {
            result_data(auto_data)
            has_initialized(TRUE)
            message("[AUTO-INIT] Data successfully loaded")

            tryCatch({
              shiny::showNotification(
                paste("已自动加载最小基因集:", min_collection, "(", min_count, "条通路)"),
                type = "default",
                duration = 5
              )
            }, error = function(e) NULL)
          } else {
            message("[AUTO-INIT] Failed to load data")
          }
        }) # end isolate
      }) # end observe (延迟执行)

    }, ignoreNULL = TRUE, ignoreInit = FALSE)

    # 手动确认按钮
    shiny::observeEvent(input$run_btn, {
      message("[MANUAL] Run button clicked")
      shiny::req(input$selected_contrast)

      manual_data <- process_data_core(
        selected_contrast_val = input$selected_contrast,
        selected_collections_val = input$selected_collections,
        sort_by_val = input$sort_by,
        expression_type_val = input$expression_type,
        plot_subtype_val = input$plot_subtype,
        custom_colors_val = input$custom_colors,
        is_auto_init = FALSE
      )

      if (is.null(manual_data)) {
        shiny::showNotification("所选配置无数据或过滤条件过严", type = "error", duration = 5)
        return()
      }

      result_data(manual_data)
      shiny::showNotification("工作台已更新", type = "message", duration = 3)
    })

    # 返回结果
    return(shiny::reactive({
      result_data()
    }))
  })
}






# ==================== 模块2: 主工作台表格 (修订版 - 接管HTML生成) ====================
#' @title 主工作台表格UI模块
#' @description 提供可调整列宽、支持联合展示选择的数据表格界面
#' @param id 模块命名空间ID
#' @return Shiny UI元素
modMasterTableUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::div(
      class = "master-table-container",
      style = "width: 100%; overflow-x: auto;",
      # 添加操作提示
      shiny::tags$div(
        style = "margin-bottom: 10px; color: #666; font-size: 12px;",
        shiny::HTML("💡 <b>操作提示：</b>拖拽列标题调整顺序 | 拖拽列边框调整宽度 | 勾选\"联合展示\"列选择通路作图 | 点击\"Dashboard\"查看详情")
      ),
      DT::dataTableOutput(outputId = ns("table"))
    ),
    # 隐藏元素用于存储联合展示选择状态
    shiny::tags$div(
      style = "display: none;",
      shiny::textInput(
        inputId = ns("joint_selection_store"),
        label = NULL,
        value = ""
      )
    )
  )
}

#' @title 主工作台表格Server模块 (修订版)
#' @description 渲染交互式DT表格，动态生成HTML按钮，统一命名空间事件监听
#' @param id 模块ID
#' @param data_prep_module 来自模块1的数据流
#' @return 响应式向量，包含当前选中的通路ID（用于联合绘图）

modMasterTableServer <- function(id, data_prep_module) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # 维护联合展示选择状态（持久化存储）
    # 维护联合展示选择状态（持久化存储）
    joint_selected_ids <- shiny::reactiveVal(character(0))
    has_joint_interaction <- shiny::reactiveVal(FALSE)

    # Safe_ID 与原始 ID 转换
    to_safe_id <- function(x) {
      gsub("'", "\\\\'", x, fixed = TRUE)
    }

    to_original_id <- function(x) {
      gsub("\\\\'", "'", x, fixed = TRUE)
    }

    # 供外部模块调用：删除指定通路
    remove_pathways <- function(ids) {
      ids <- ids[!is.na(ids) & nzchar(ids)]
      if (length(ids) == 0) return(invisible(NULL))

      has_joint_interaction(TRUE)
      safe_ids <- to_safe_id(ids)
      new_selection <- setdiff(joint_selected_ids(), safe_ids)
      joint_selected_ids(new_selection)

      message(sprintf(
        "[JOINT-SELECT] Removed from external selector: %s (remain: %d)",
        paste(ids, collapse = ", "),
        length(new_selection)
      ))

      invisible(new_selection)
    }

    # 供外部模块调用：清空所有通路
    clear_selection <- function() {
      has_joint_interaction(TRUE)
      joint_selected_ids(character(0))
      message("[JOINT-SELECT] Cleared all from external selector")
      invisible(character(0))
    }

    # 关键修复：监听本模块生成的复选框事件（命名空间：master_table-joint_plot_toggle）
    shiny::observeEvent(input$joint_plot_toggle, {
      has_joint_interaction(TRUE)
      toggle_info <- input$joint_plot_toggle

      if (!is.null(toggle_info) && is.list(toggle_info)) {
        pathway_id <- toggle_info$id
        is_checked <- toggle_info$checked

        current_selection <- joint_selected_ids()

        if (isTRUE(is_checked)) {
          # 添加到选择集（去重）
          if (!(pathway_id %in% current_selection)) {
            new_selection <- c(current_selection, pathway_id)
            joint_selected_ids(new_selection)
            message(sprintf("[JOINT-SELECT] Added: %s (total: %d)", pathway_id, length(new_selection)))
          }
        } else {
          # 从选择集移除
          new_selection <- setdiff(current_selection, pathway_id)
          joint_selected_ids(new_selection)
          message(sprintf("[JOINT-SELECT] Removed: %s (total: %d)", pathway_id, length(new_selection)))
        }
      }
    })

    # 全选/清空功能（可选API，供外部调用）
    shiny::observeEvent(input$select_all_visible, {
      has_joint_interaction(TRUE)
      data_list <- data_prep_module()
      shiny::req(data_list)

      all_ids <- data_list$df$Safe_ID
      joint_selected_ids(all_ids)
      message(sprintf("[JOINT-SELECT] Select all: %d items", length(all_ids)))

      # 触发表格重绘以更新复选框状态（通过JS）
      session$sendCustomMessage(type = "checkAllBoxes", message = list(ids = all_ids))
    })

    shiny::observeEvent(input$clear_all_selection, {
      has_joint_interaction(TRUE)
      joint_selected_ids(character(0))
      message("[JOINT-SELECT] Cleared all")
      session$sendCustomMessage(type = "uncheckAllBoxes", message = list())
    })

    # 主表格渲染（关键修复：动态添加HTML列，统一命名空间）
    output$table <- DT::renderDataTable({
      # 数据验证
      shiny::validate(shiny::need(
        data_prep_module(),
        "等待数据加载...\n请在左侧选择对比组并点击\"确认配置\""
      ))

      data_list <- data_prep_module()
      shiny::validate(shiny::need(
        !is.null(data_list$df) && nrow(data_list$df) > 0,
        "当前筛选条件下无数据，请调整基因集亚组选择"
      ))

      display_df <- data_list$df

      message(sprintf("[TABLE-RENDER] Rendering %d rows", nrow(display_df)))

      # ========== 关键修复：动态生成HTML交互列（使用master_table命名空间） ==========

      # 1. 生成联合展示复选框（使用ns()确保命名空间正确）
      display_df$Select_for_Plot <- sprintf(
        '<input type="checkbox" class="joint-plot-checkbox" data-id="%s" id="%s_%s" %s onclick="Shiny.setInputValue(&#39;%s&#39;, {id: &#39;%s&#39;, checked: this.checked}, {priority: &#39;event&#39;});"/>',
        display_df$Safe_ID,
        ns("chk"),           # 使用本模块命名空间
        display_df$Safe_ID,
        ifelse(display_df$Safe_ID %in% joint_selected_ids(), 'checked="checked"', ''),
        ns("joint_plot_toggle"),  # 事件名：master_table-joint_plot_toggle
        display_df$Safe_ID
      )

      # 2. 生成详情按钮（使用ns()确保命名空间正确）
      display_df$Detail_Btn <- sprintf(
        '<button class="btn btn-sm btn-success pathway-detail-btn" onClick="Shiny.setInputValue(&#39;%s&#39;, &#39;%s&#39;, {priority: &#39;event&#39;})">🔍 Dashboard</button>',
        ns("show_modal"),    # 事件名：master_table-show_modal
        display_df$Safe_ID
      )

      # 确定展示列顺序（Select_for_Plot和Detail_Btn现在由本模块生成）
      base_cols <- c("Rank", "Select_for_Plot", "Detail_Btn", "Pathway_Link",
                     "Enriched_In", "setSize", "NES", "pvalue", "p.adjust",
                     "Brief_Description")
      blank_cols <- paste0("blank", 1:10)
      id_col <- "ID"

      show_cols <- c(base_cols, blank_cols, id_col)
      show_cols <- intersect(show_cols, colnames(display_df))

      message(sprintf("[TABLE-RENDER] Columns: %s", paste(show_cols, collapse=", ")))

      # 数值格式化（保留原始值用于排序）
      display_df$NES_formatted <- round(display_df$NES, 3)
      display_df$pvalue_formatted <- signif(display_df$pvalue, 3)
      display_df$p.adjust_formatted <- signif(display_df$p.adjust, 3)

      # 构建最终显示数据框（按show_cols顺序）
      display_df_ordered <- display_df[, show_cols, drop = FALSE]

      # 替换为格式化版本（但保留原始列用于排序）
      if ("NES" %in% colnames(display_df_ordered)) {
        display_df_ordered$NES <- display_df$NES_formatted
      }
      if ("pvalue" %in% colnames(display_df_ordered)) {
        display_df_ordered$pvalue <- display_df$pvalue_formatted
      }
      if ("p.adjust" %in% colnames(display_df_ordered)) {
        display_df_ordered$p.adjust <- display_df$p.adjust_formatted
      }

      # 关键修复：安全计算列索引（避免integer(0)导致DT报错）
      # 使用match替代which，找不到时返回NA而非integer(0)
      get_col_idx <- function(col_name) {
        idx <- match(col_name, show_cols)
        if (is.na(idx)) return(NULL)
        return(idx - 1)  # DT使用0-based索引
      }

      rank_idx <- get_col_idx("Rank")
      select_idx <- get_col_idx("Select_for_Plot")
      detail_idx <- get_col_idx("Detail_Btn")
      pathway_idx <- get_col_idx("Pathway_Link")
      enriched_idx <- get_col_idx("Enriched_In")
      desc_idx <- get_col_idx("Brief_Description")
      blank_idx <- which(grepl("^blank", show_cols)) - 1
      id_idx <- get_col_idx("ID")

      # 构建columnDefs（过滤掉NULL）
      col_defs <- list()

      if (!is.null(rank_idx)) {
        col_defs[[length(col_defs) + 1]] <- list(width = '50px', targets = rank_idx)
      }

      if (!is.null(select_idx)) {
        col_defs[[length(col_defs) + 1]] <- list(
          width = '80px',
          targets = select_idx,
          className = 'dt-center',
          orderable = TRUE,
          render = htmlwidgets::JS("function(data, type, row) {
            if (type === 'sort') {
              var match = data.match(/checked\\s*=/i);
              return match ? '1' : '0';
            }
            return data;
          }")
        )
      }

      if (!is.null(detail_idx)) {
        col_defs[[length(col_defs) + 1]] <- list(
          width = '100px',
          targets = detail_idx,
          orderable = FALSE
        )
      }

      if (!is.null(pathway_idx)) {
        col_defs[[length(col_defs) + 1]] <- list(
          width = '200px',
          targets = pathway_idx
        )
      }

      if (!is.null(enriched_idx)) {
        col_defs[[length(col_defs) + 1]] <- list(
          width = '100px',
          targets = enriched_idx
        )
      }

      # 数值列（setSize, NES, pvalue, p.adjust）
      num_cols <- c("setSize", "NES", "pvalue", "p.adjust")
      num_idx <- sapply(num_cols, get_col_idx)
      num_idx <- num_idx[!sapply(num_idx, is.null)]
      if (length(num_idx) > 0) {
        col_defs[[length(col_defs) + 1]] <- list(
          width = '80px',
          targets = unlist(num_idx)
        )
      }

      if (!is.null(desc_idx)) {
        col_defs[[length(col_defs) + 1]] <- list(
          width = '400px',
          targets = desc_idx,
          className = 'dt-left description-cell',
          render = htmlwidgets::JS("function(data, type, row) {
            if (type === 'display' && data.length > 100) {
              return '<div style=\"white-space: normal; word-wrap: break-word; max-height: 60px; overflow: hidden; text-overflow: ellipsis;\" title=\"' + data + '\">' + data.substring(0, 100) + '...</div>';
            }
            return '<div style=\"white-space: normal; word-wrap: break-word;\">' + data + '</div>';
          }")
        )
      }

      if (length(blank_idx) > 0) {
        col_defs[[length(col_defs) + 1]] <- list(
          visible = FALSE,
          targets = blank_idx
        )
      }

      if (!is.null(id_idx)) {
        col_defs[[length(col_defs) + 1]] <- list(
          visible = FALSE,
          targets = id_idx
        )
      }

      # 构建DT表格（关键修复：添加rowCallback确保复选框事件绑定）
      dt <- DT::datatable(
        data = display_df_ordered,
        escape = FALSE,
        selection = "multiple",
        rownames = FALSE,
        extensions = c('Buttons', 'Scroller', 'ColReorder'),
        callback = htmlwidgets::JS(sprintf("
          table.on('column-reorder', function(e, settings, details) {
            console.log('Column reordered:', details);
          });

          // 关键修复：确保每次重绘后重新绑定复选框事件（使用正确命名空间）
          table.on('draw', function() {
            var checkboxes = document.querySelectorAll('.joint-plot-checkbox');
            checkboxes.forEach(function(chk) {
              // 移除旧事件避免重复绑定
              var newChk = chk.cloneNode(true);
              chk.parentNode.replaceChild(newChk, chk);

              // 绑定新事件（使用master_table命名空间）
              newChk.addEventListener('change', function() {
                var eventData = {
                  id: this.getAttribute('data-id'),
                  checked: this.checked
                };
                Shiny.setInputValue('%s', eventData, {priority: 'event'});
              });
            });
          });
        ", ns("joint_plot_toggle"))),
        options = list(
          scrollX = TRUE,
          scrollY = "65vh",
          scroller = TRUE,
          deferRender = TRUE,
          pageLength = -1,
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel', 'colvis'),
          colReorder = TRUE,
          autoWidth = TRUE,
          ordering = TRUE,
          columnDefs = col_defs,
          order = list(list(0, 'asc'))  # 默认按第一列（Rank）升序
        )
      ) %>%
        DT::formatStyle(
          columns = "Enriched_In",
          backgroundColor = DT::styleEqual(
            unique(display_df$Enriched_In),
            c('#fee0d2', '#deebf7', '#c6dbef', '#fdd0a2')[1:length(unique(display_df$Enriched_In))]
          ),
          fontWeight = 'bold'
        ) %>%
        DT::formatStyle(
          columns = "NES",
          color = DT::styleInterval(0, c('#1052bd', '#cc212f')),
          fontWeight = 'bold'
        ) %>%
        DT::formatStyle(
          columns = "p.adjust",
          backgroundColor = DT::styleInterval(c(0.01, 0.05), c('#fc9272', '#fee0d2', 'transparent')),
          fontWeight = DT::styleInterval(0.05, c('bold', 'normal'))
        )

      dt
    }, server = TRUE)

    # 返回选中的通路ID（优先使用联合展示选择，回退到行选择）
    selected_pathways <- shiny::reactive({
      joint_ids <- joint_selected_ids()

      # 只要用户已经使用过“联合展示”复选框/外部删除器，
      # 就只认 joint_selected_ids，即使它当前为空也不回退到行选择
      if (isTRUE(has_joint_interaction())) {
        original_ids <- to_original_id(joint_ids)
        message(sprintf("[SELECTION] Returning %d joint-selected pathways", length(original_ids)))
        return(original_ids)
      }

      # 回退到传统行选择（兼容保留）
      data_list <- data_prep_module()
      if (is.null(data_list) || is.null(data_list$df)) {
        return(character(0))
      }

      if (is.null(input$table_rows_selected) || length(input$table_rows_selected) == 0) {
        return(character(0))
      }

      selected_rows <- input$table_rows_selected
      selected_ids <- data_list$df$ID[selected_rows]
      message(sprintf("[SELECTION] Returning %d row-selected pathways", length(selected_ids)))
      return(selected_ids)
    })

    # 新增：返回详情按钮触发事件（供父模块监听弹窗）
    # 通过 reactive 暴露 show_modal 事件
    show_modal_trigger <- shiny::reactive({
      input$show_modal
    })

    # 返回联合展示选择状态（供模块5使用）和弹窗触发器
    return(list(
      selected_pathways = selected_pathways,
      show_modal = show_modal_trigger,
      remove_pathways = remove_pathways,
      clear_selection = clear_selection
    ))
  })
}

   # ==================== 模块====================

# ==================== 模块3: 全息四重联动 (修订版 - 支持双向对比方向校正) ====================
#' @title 全息四重联动UI模块
#' @description 包含通路火山图、基因Rank分布、差异表达火山图和表达箱线图
#' @param id 模块命名空间ID
#' @return Shiny UI元素列表
modQuadrantUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      shiny::column(
        width = 6,
        shiny::div(
          class = "white-box",
          shiny::h4("1. 宏观: 通路火山图 (点击持久标签)"),
          plotly::plotlyOutput(outputId = ns("volcano_pathway"), height = "450px")
        )
      ),
      shiny::column(
        width = 6,
        shiny::div(
          class = "white-box",
          shiny::h4("2. 微观: 基因 Rank 分布 (无描边)"),
          plotly::plotlyOutput(outputId = ns("volcano_gene"), height = "450px")
        )
      )
    ),
    shiny::fluidRow(
      shiny::column(
        width = 6,
        shiny::div(
          class = "white-box",
          shiny::h4("3. 差异表达火山图 (高亮置顶)"),
          plotly::plotlyOutput(outputId = ns("limma_volcano"), height = "450px")
        )
      ),
      shiny::column(
        width = 6,
        shiny::div(
          class = "white-box",
          shiny::h4("4. 全量表达箱线图 (大小写不敏感)"),
          plotly::plotlyOutput(outputId = ns("gene_expr_box"), height = "450px")
        )
      )
    )
  )
}

modQuadrantServer <- function(id, data_prep_module, res_capsule) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    selected_labels <- shiny::reactiveVal(character(0))

    current_coef_name <- shiny::reactive({
      data_list <- data_prep_module()
      if (is.null(data_list) || is.null(data_list$coef_name)) return(NULL)
      return(data_list$coef_name)
    })

    current_expression_type <- shiny::reactive({
      data_list <- data_prep_module()
      if (is.null(data_list) || is.null(data_list$expression_type)) {
        return("log2cpm")
      }
      return(data_list$expression_type)
    })

    # 辅助函数：优化对比组显示格式（IBAA_vs_PREA → IBAA vs PREA）
    format_contrast_display <- function(contrast_name) {
      if (is.null(contrast_name)) return("")
      # 先处理 _vs_，再处理剩余下划线
      display_name <- gsub("_vs_", " vs ", contrast_name)
      display_name <- gsub("_", " ", display_name)
      return(display_name)
    }

    # 1. 宏观通路火山图
    output$volcano_pathway <- plotly::renderPlotly({
      data_list <- data_prep_module()
      shiny::req(data_list)

      pathway_df <- data_list$df
      pathway_df$Color <- dplyr::case_when(
        pathway_df$NES > 0 ~ "#cc212f",
        pathway_df$NES < 0 ~ "#1052bd",
        TRUE ~ "#C0C0C0"
      )
      pathway_df$Size <- 10

      # 标题格式优化
      contrast_display <- format_contrast_display(data_list$contrast_name)
      plot_title <- sprintf("Contrast: %s", contrast_display)

      p <- plotly::plot_ly(
        data = pathway_df,
        x = ~NES,
        y = ~-log10(p.adjust),
        type = "scatter",
        mode = "markers",
        marker = list(
          color = ~Color,
          size = ~Size,
          opacity = 0.8,
          line = list(color = "white", width = 1)
        ),
        text = ~sprintf("%s<br>NES: %.3f<br>FDR: %.2e", ID, NES, p.adjust),
        hoverinfo = "text",
        key = ~ID,
        source = ns("pathway_volcano")
      ) %>%
        plotly::layout(
          title = list(text = plot_title, font = list(size = 12)),
          xaxis = list(title = "NES", zeroline = FALSE),
          yaxis = list(title = "-log10 (FDR)", zeroline = FALSE),
          dragmode = "select",
          showlegend = FALSE
        )

      labs <- selected_labels()
      if (length(labs) > 0) {
        lab_df <- pathway_df[pathway_df$ID %in% labs, ]
        if (nrow(lab_df) > 0) {
          p <- p %>% plotly::add_annotations(
            x = lab_df$NES,
            y = -log10(lab_df$p.adjust),
            text = lab_df$ID,
            showarrow = TRUE,
            arrowhead = 2,
            ax = 20,
            ay = -30,
            font = list(size = 10, color = "black"),
            bgcolor = "rgba(255,255,255,0.9)",
            bordercolor = "#333",
            borderwidth = 1
          )
        }
      }
      p
    })

    shiny::observeEvent(plotly::event_data("plotly_click", source = ns("pathway_volcano")), {
      click_data <- plotly::event_data("plotly_click", source = ns("pathway_volcano"))
      if (!is.null(click_data$key)) {
        current_labels <- selected_labels()
        if (click_data$key %in% current_labels) {
          selected_labels(setdiff(current_labels, click_data$key))
        } else {
          selected_labels(c(current_labels, click_data$key))
        }
      }
    })

    # 2. 微观基因Rank分布
    output$volcano_gene <- plotly::renderPlotly({
      data_list <- data_prep_module()
      shiny::req(data_list)

      gene_list <- data_list$gsea_res@geneList

      rank_df <- data.frame(
        Rank = seq_along(gene_list),
        Metric = as.numeric(gene_list),
        Gene = names(gene_list),
        Gene_Upper = toupper(names(gene_list)),
        stringsAsFactors = FALSE
      )

      rank_df$Point_Color <- "#CFD8DC"
      rank_df$Point_Size <- 4
      rank_df$Point_Alpha <- 0.3

      pathway_click <- plotly::event_data("plotly_click", source = ns("pathway_volcano"))
      if (!is.null(pathway_click) && !is.null(pathway_click$key)) {
        pathway_genes_upper <- toupper(data_list$gsea_res@geneSets[[pathway_click$key]])
        match_idx <- which(rank_df$Gene_Upper %in% pathway_genes_upper)

        if (length(match_idx) > 0) {
          rank_df$Point_Color[match_idx] <- "#FF9800"
          rank_df$Point_Size[match_idx] <- 12
          rank_df$Point_Alpha[match_idx] <- 1
        }
      }

      # 标题显示当前通路和对比组
      contrast_display <- format_contrast_display(data_list$contrast_name)
      plot_title <- sprintf("Rank Dist. | %s", contrast_display)
      if (!is.null(pathway_click$key)) {
        plot_title <- sprintf("Pathway: %s | %s", pathway_click$key, contrast_display)
      }

      plotly::plot_ly(
        data = rank_df,
        x = ~Rank,
        y = ~Metric,
        type = "scattergl",
        mode = "markers",
        marker = list(
          color = ~Point_Color,
          size = ~Point_Size,
          opacity = ~Point_Alpha,
          line = list(width = 0)
        ),
        text = ~Gene,
        hoverinfo = "text",
        key = ~Gene_Upper
      ) %>%
        plotly::layout(
          title = list(text = plot_title, font = list(size = 11)),
          xaxis = list(title = "Gene Rank"),
          yaxis = list(title = "Ranking Metric"),
          showlegend = FALSE
        )
    })

    # 3. 差异表达火山图（关键修复：支持双向对比方向校正）
    output$limma_volcano <- plotly::renderPlotly({
      shiny::req(res_capsule$limma_fit)
      data_list <- data_prep_module()
      shiny::req(data_list)

      fit_object <- res_capsule$limma_fit
      coef_name <- current_coef_name()

      if (is.null(coef_name)) {
        return(plotly::plot_ly() %>%
                 plotly::layout(title = "未找到匹配的对比系数"))
      }

      tryCatch({
        top_table <- limma::topTable(
          fit = fit_object,
          coef = coef_name,
          number = Inf,
          sort.by = "none"
        )

        # ========== 关键修复：根据coef_reversed标志反转logFC方向 ==========
        if (!is.null(data_list$coef_reversed) && isTRUE(data_list$coef_reversed)) {
          top_table$logFC <- -top_table$logFC
          message(sprintf("[VOLCANO] Reversing logFC for %s using reversed coefficient [%s]",
                          data_list$contrast_name, coef_name))
        } else {
          message(sprintf("[VOLCANO] Using coefficient [%s] for %s (forward direction)",
                          coef_name, data_list$contrast_name))
        }
      },
      error = function(error_obj) {
        warning(sprintf("提取系数 [%s] 失败: %s", coef_name, error_obj$message))
        return(plotly::plot_ly() %>%
                 plotly::layout(title = sprintf("系数提取失败: %s", coef_name)))
      })

      top_table$Gene_Symbol <- rownames(top_table)
      if ("SYMBOL" %in% colnames(top_table)) {
        top_table$Gene_Symbol <- top_table$SYMBOL
      }
      top_table$Gene_Upper <- toupper(top_table$Gene_Symbol)

      pval_col <- if ("adj.P.Val" %in% colnames(top_table)) "adj.P.Val" else "P.Value"
      top_table$negLogP <- -log10(top_table$P.Value)

      top_table$Base_Color <- dplyr::case_when(
        top_table$logFC > 0.5 & top_table[[pval_col]] < 0.05 ~ "#cc212f",
        top_table$logFC < -0.5 & top_table[[pval_col]] < 0.05 ~ "#1052bd",
        TRUE ~ "#C0C0C0"
      )
      top_table$Point_Size <- 10
      top_table$Is_Highlight <- FALSE

      pathway_click <- plotly::event_data("plotly_click", source = ns("pathway_volcano"))
      if (!is.null(pathway_click) && !is.null(pathway_click$key)) {
        pathway_genes_upper <- toupper(data_list$gsea_res@geneSets[[pathway_click$key]])
        highlight_idx <- which(top_table$Gene_Upper %in% pathway_genes_upper)

        if (length(highlight_idx) > 0) {
          top_table$Base_Color[highlight_idx] <- "#FF9800"
          top_table$Point_Size[highlight_idx] <- 16
          top_table$Is_Highlight[highlight_idx] <- TRUE
        }
      }

      top_table <- top_table[order(top_table$Is_Highlight), ]

      # 标题优化（显示格式化的对比组名）
      contrast_display <- format_contrast_display(data_list$contrast_name)
      plot_title <- sprintf("DEG: %s (n=%d)", contrast_display, nrow(top_table))

      plotly::plot_ly(
        data = top_table,
        x = ~logFC,
        y = ~negLogP,
        type = "scattergl",
        mode = "markers",
        marker = list(
          color = ~Base_Color,
          size = ~Point_Size,
          opacity = 0.8,
          line = list(color = "white", width = 1)
        ),
        text = ~sprintf("%s<br>logFC: %.2f<br>FDR: %.2e", Gene_Symbol, logFC, .data[[pval_col]]),
        hoverinfo = "text",
        key = ~Gene_Upper,
        source = ns("gene_volcano")
      ) %>%
        plotly::layout(
          title = list(text = plot_title, font = list(size = 12)),
          xaxis = list(title = "log2 Fold Change", zeroline = FALSE),
          yaxis = list(title = "-log10 P-value", zeroline = FALSE),
          showlegend = FALSE
        )
    })

    # 4. 全量表达箱线图
    output$gene_expr_box <- plotly::renderPlotly({
      shiny::req(res_capsule$expr_data)
      data_list <- data_prep_module()

      expression_type <- current_expression_type()

      gene_click <- plotly::event_data("plotly_click", source = ns("gene_volcano"))

      if (is.null(gene_click) || is.null(gene_click$key)) {
        return(plotly::plot_ly() %>%
                 plotly::layout(
                   title = list(
                     text = "👈 请在左侧火山图点击基因",
                     font = list(size = 14)
                   )
                 ))
      }

      target_gene_upper <- gene_click$key
      dge_list <- res_capsule$expr_data

      # 表达量计算
      expression_matrix <- switch(
        expression_type,
        "log2cpm" = edgeR::cpm(dge_list, log = TRUE),
        "cpm" = edgeR::cpm(dge_list, log = FALSE),
        "log2fpkm" = {
          if (!is.null(dge_list$genes) && "Length" %in% colnames(dge_list$genes)) {
            log2(edgeR::rpkm(dge_list, gene.length = dge_list$genes$Length) + 1)
          } else {
            warning("未找到基因长度信息，回退到log2(CPM)")
            edgeR::cpm(dge_list, log = TRUE)
          }
        },
        "fpkm" = {
          if (!is.null(dge_list$genes) && "Length" %in% colnames(dge_list$genes)) {
            edgeR::rpkm(dge_list, gene.length = dge_list$genes$Length)
          } else {
            warning("未找到基因长度信息，回退到CPM")
            edgeR::cpm(dge_list, log = FALSE)
          }
        },
        "log2tpm" = {
          if (!is.null(dge_list$genes) && "Length" %in% colnames(dge_list$genes)) {
            cpm_vals <- edgeR::cpm(dge_list, log = FALSE)
            kb_lengths <- dge_list$genes$Length / 1000
            rpk_vals <- cpm_vals / kb_lengths
            tpm_vals <- t(t(rpk_vals) * 1e6 / colSums(rpk_vals))
            log2(tpm_vals + 1)
          } else {
            warning("未找到基因长度信息，回退到log2(CPM)")
            edgeR::cpm(dge_list, log = TRUE)
          }
        },
        "tpm" = {
          if (!is.null(dge_list$genes) && "Length" %in% colnames(dge_list$genes)) {
            cpm_vals <- edgeR::cpm(dge_list, log = FALSE)
            kb_lengths <- dge_list$genes$Length / 1000
            rpk_vals <- cpm_vals / kb_lengths
            t(t(rpk_vals) * 1e6 / colSums(rpk_vals))
          } else {
            warning("未找到基因长度信息，回退到CPM")
            edgeR::cpm(dge_list, log = FALSE)
          }
        },
        edgeR::cpm(dge_list, log = TRUE)
      )

      gene_names_upper <- toupper(rownames(expression_matrix))
      match_idx <- which(gene_names_upper == target_gene_upper)

      if (length(match_idx) == 0) {
        return(plotly::plot_ly() %>%
                 plotly::layout(
                   title = list(
                     text = sprintf("基因 '%s' 未找到", target_gene_upper),
                     font = list(size = 14)
                   )
                 ))
      }

      actual_gene <- rownames(expression_matrix)[match_idx[1]]
      gene_expression <- expression_matrix[actual_gene, ]

      sample_metadata <- dge_list$samples

      # 构建绘图数据框 - 保留所有组别
      plot_data <- data.frame(
        Sample_ID = names(gene_expression),
        Expression_Value = as.numeric(gene_expression),
        Group = as.character(sample_metadata$group),
        stringsAsFactors = FALSE
      )

      # 动态Y轴标签
      y_axis_label <- switch(
        expression_type,
        "log2cpm" = "log2(CPM)",
        "cpm" = "CPM",
        "log2fpkm" = "log2(FPKM)",
        "fpkm" = "FPKM",
        "log2tpm" = "log2(TPM)",
        "tpm" = "TPM",
        "Expression"
      )

      # 标题优化
      contrast_display <- format_contrast_display(data_list$contrast_name)
      plot_title <- sprintf("%s | %s | Ref: %s", actual_gene, y_axis_label, contrast_display)

      # 颜色映射扩展以支持更多组别（最多8种颜色）
      all_groups <- unique(plot_data$Group)
      group_colors <- c("#1052bd", "#cc212f", "#4DAF4A", "#984EA3",
                        "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
      color_mapping <- stats::setNames(
        group_colors[seq_along(all_groups)],
        all_groups
      )

      p <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(x = Group, y = Expression_Value, fill = Group)
      ) +
        ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        ggplot2::geom_jitter(
          width = 0.2,
          size = 3,
          ggplot2::aes(text = Sample_ID),
          alpha = 0.8
        ) +
        ggplot2::scale_fill_manual(values = color_mapping) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::labs(
          title = plot_title,
          y = y_axis_label,
          x = ""
        ) +
        ggplot2::theme(
          legend.position = "none",
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          plot.title = ggplot2::element_text(size = 10, face = "bold")
        )

      plotly::ggplotly(p, tooltip = c("text", "y"))
    })
  })
}
# ==================== 模块4: Dashboard弹窗（修复版 + ComplexHeatmap重构） ====================
modPathwayModalUI <- function(id) {
  ns <- shiny::NS(id)
  # 模态框通过showModal动态生成，此处仅占位
  NULL
}

modPathwayModalServer <- function(id, data_prep_module, trigger_event, res_capsule) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ========== 关键修复1：使用reactiveVal存储当前通路数据 ==========
    current_pathway_data <- shiny::reactiveVal(NULL)

    # 渲染GSEA图 (保留原有逻辑)
    output$modal_gsea_plot <- shiny::renderPlot({
      pdata <- current_pathway_data()
      shiny::req(pdata)

      tryCatch({
        print(plot_directional_gsea(
          directional_gsea_obj = list(
            gsea_res = pdata$data_list$gsea_res,
            meta = list(
              left_group = pdata$data_list$left_group,
              right_group = pdata$data_list$right_group
            )
          ),
          target_pathways = pdata$pathway_id,
          subPlot = 3,
          add_pval = TRUE
        ))
      }, error = function(error_object) {
        graphics::plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        graphics::text(1, 1, paste("绘图失败:", error_object$message), col = "red", cex = 1.2)
      })
    })

    # ================= 核心重构：渲染热图 (使用 ComplexHeatmap，保持原有美学) =================
    output$modal_heatmap <- shiny::renderPlot({
      pdata <- current_pathway_data()
      shiny::req(pdata)

      data_list <- pdata$data_list
      pathway_genes <- pdata$pathway_genes
      core_genes <- pdata$core_genes

      shiny::req(res_capsule$expr_data)
      dge_list <- res_capsule$expr_data
      sample_meta <- dge_list$samples

      # 恢复原有逻辑：只显示当前对比组的两组样本
      target_samples <- rownames(sample_meta)[
        sample_meta$group %in% c(data_list$left_group, data_list$right_group)
      ]

      if (length(target_samples) == 0) return()

      expr_matrix <- edgeR::cpm(dge_list, log = TRUE)[, target_samples, drop = FALSE]
      cpm_matrix <- edgeR::cpm(dge_list, log = FALSE)[, target_samples, drop = FALSE]

      expr_genes <- rownames(expr_matrix)
      match_idx <- which(toupper(expr_genes) %in% toupper(pathway_genes))

      if (length(match_idx) < 2) return()

      plot_genes <- expr_genes[match_idx]

      # 恢复原有排序逻辑
      gene_metrics <- sapply(plot_genes, function(gene_symbol) {
        idx <- match(toupper(gene_symbol), toupper(names(data_list$gsea_res@geneList)))
        if (is.na(idx)) return(0)
        return(data_list$gsea_res@geneList[idx])
      })

      is_leading <- toupper(plot_genes) %in% toupper(core_genes)
      sort_order <- order(is_leading, gene_metrics, decreasing = TRUE)
      plot_genes <- plot_genes[sort_order]

      plot_matrix <- expr_matrix[plot_genes, , drop = FALSE]
      plot_matrix <- plot_matrix[apply(plot_matrix, 1, stats::var) > 1e-6, , drop = FALSE]

      if (nrow(plot_matrix) < 2) return()

      # Z-score 标准化与截断（保持原有参数）
      z_matrix <- t(scale(t(plot_matrix)))
      z_matrix[is.na(z_matrix)] <- 0
      z_matrix[z_matrix > 1.5] <- 1.5
      z_matrix[z_matrix < -1.5] <- -1.5

      # 提取需要展示的 CPM 数字
      display_numbers <- round(cpm_matrix[rownames(z_matrix), , drop = FALSE])

      # 构造 ComplexHeatmap 注释与参数（保持原有美学）
      group_split <- factor(sample_meta[target_samples, "group"],
                            levels = c(data_list$left_group, data_list$right_group))

      # 恢复原有颜色映射 (淡蓝-白-亮粉橙)
      col_fun <- circlize::colorRamp2(c(-1.5, 0, 1.5), c("#67a9cf", "#f7f7f7", "#ef8a62"))

      # 列注释颜色（红蓝对比）
      group_colors <- c("#da130f", "#447fba")
      names(group_colors) <- c(data_list$left_group, data_list$right_group)

      top_ann <- ComplexHeatmap::HeatmapAnnotation(
        Group = sample_meta[target_samples, "group"],
        col = list(Group = group_colors),
        annotation_name_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        simple_anno_size = grid::unit(0.6, "cm")
      )

      # 行注释 (Leading Edge)
      leading_status <- ifelse(toupper(rownames(z_matrix)) %in% toupper(core_genes), "YES", "NO")
      leading_colors <- c("YES" = "#FF9800", "NO" = "#E0E0E0")

      right_ann <- ComplexHeatmap::rowAnnotation(
        LeadingEdge = leading_status,
        col = list(LeadingEdge = leading_colors),
        annotation_name_gp = grid::gpar(fontsize = 12, fontface = "bold"),
        simple_anno_size = grid::unit(0.4, "cm")
      )

      # 恢复原有单元格渲染：固定黑色字体，大字号
      cell_fun <- function(j, i, x, y, width, height, fill) {
        val <- display_numbers[i, j]
        # 保持原有设置：固定黑色字体
        grid::grid.text(val, x, y, gp = grid::gpar(fontsize = 13, col = "black", fontface = "bold"))
      }

      # 绘制 Heatmap 对象 (保持原有参数)
      ht <- ComplexHeatmap::Heatmap(
        z_matrix,
        name = "Z-Score",
        col = col_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_split = group_split,
        cluster_column_slices = FALSE,
        top_annotation = top_ann,
        right_annotation = right_ann,
        cell_fun = cell_fun,
        row_names_gp = grid::gpar(fontsize = 15, fontface = ifelse(leading_status == "YES", "bold", "plain")),
        column_names_gp = grid::gpar(fontsize = 15, fontface = "bold"),
        rect_gp = grid::gpar(col = "white", lwd = 1), # 白色网格线
        show_heatmap_legend = TRUE,
        width = NULL,
        height = NULL
      )

      ComplexHeatmap::draw(ht, merge_legend = TRUE)

    }, height = function() {
      # 动态高度计算：保持原有逻辑
      pdata <- current_pathway_data()
      if (is.null(pdata)) return(400)

      n_genes <- length(intersect(toupper(rownames(res_capsule$expr_data$counts)), toupper(pdata$pathway_genes)))
      return(max(400, n_genes * 28 + 150))
    })

    # 渲染基因统计表 (保留原有逻辑，添加Rank_in_List)
    output$modal_gene_table <- DT::renderDataTable({
      pdata <- current_pathway_data()
      shiny::req(pdata)

      gene_list <- pdata$data_list$gsea_res@geneList
      all_genes_in_set <- pdata$pathway_genes
      core_genes <- pdata$core_genes

      # 构建基因排名表（添加Rank_in_List）
      gene_ranks <- data.frame(
        Gene_Symbol = names(gene_list),
        Rank_Metric = round(as.numeric(gene_list), 3),
        stringsAsFactors = FALSE
      )

      # 计算在全基因组列表中的排名
      gene_ranks$Rank_in_List <- seq_len(nrow(gene_ranks))

      gene_table <- gene_ranks[toupper(gene_ranks$Gene_Symbol) %in% toupper(all_genes_in_set), ]
      gene_table$Is_Core <- ifelse(
        toupper(gene_table$Gene_Symbol) %in% toupper(core_genes),
        "✅ YES", "—"
      )

      # 排序：Core基因在前，其次按Metric降序
      gene_table <- gene_table[order(gene_table$Is_Core == "✅ YES", gene_table$Rank_Metric, decreasing = TRUE), ]
      rownames(gene_table) <- NULL

      # 调整列顺序：Gene, Rank_in_List, Rank_Metric, Is_Core
      display_cols <- c("Gene_Symbol", "Rank_in_List", "Rank_Metric", "Is_Core")

      DT::datatable(
        gene_table[, display_cols],
        rownames = FALSE,
        escape = FALSE,
        extensions = c('Scroller'),
        options = list(
          pageLength = -1,
          scrollX = TRUE,
          scrollY = "40vh",
          scroller = TRUE,
          deferRender = TRUE,
          dom = 'frtip',
          order = list(list(3, 'desc'), list(2, 'desc'))  # 先按Is_Core排序，再按Metric
        ),
        colnames = c("Gene", "Rank in List", "Metric", "Leading Edge")
      ) %>%
        DT::formatStyle(
          columns = "Is_Core",
          backgroundColor = DT::styleEqual(c("✅ YES", "—"), c("#FF9800", "transparent")),
          fontWeight = DT::styleEqual("✅ YES", "bold"),
          color = DT::styleEqual("✅ YES", "white")
        )
    })

    # ========== 关键修复3：observeEvent只负责更新数据和显示modal ==========
    shiny::observeEvent(trigger_event(), {
      pathway_id <- trigger_event()
      shiny::req(pathway_id)

      data_list <- data_prep_module()
      shiny::req(data_list)

      res_df <- as.data.frame(data_list$gsea_res, stringsAsFactors = FALSE)
      core_str <- res_df$core_enrichment[res_df$ID == pathway_id]

      core_genes <- if (length(core_str) > 0 && !is.na(core_str[1])) {
        unlist(strsplit(as.character(core_str[1]), "/"))
      } else {
        character(0)
      }

      pathway_genes <- data_list$gsea_res@geneSets[[pathway_id]]

      current_pathway_data(list(
        pathway_id = pathway_id,
        core_genes = core_genes,
        pathway_genes = pathway_genes,
        data_list = data_list
      ))

      # 标题格式优化：IBAA_vs_PREA -> IBAA vs PREA
      contrast_display <- gsub("_vs_", " vs ", data_list$contrast_name)

      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML(sprintf("<b style='color:#0056b3;'>Pathway Dashboard: %s</b><br><small>Contrast: %s</small>", pathway_id, contrast_display)),
        size = "l",
        easyClose = TRUE,
        shiny::fluidRow(
          shiny::column(
            width = 5,
            shiny::div(
              class = "white-box",
              shiny::h4("经典 GSEA 富集轮廓"),
              shiny::plotOutput(outputId = ns("modal_gsea_plot"), height = "400px")
            )
          ),
          shiny::column(
            width = 7,
            shiny::div(
              class = "white-box",
              shiny::h4("核心基因表达热图 (CPM)"),
              # 关键修复：外层 div 限制最大高度为 650px 并支持内部滚动；
              # plotOutput height 设为 "auto"，由 server 端真实基因数计算出的高度直接撑开
              shiny::div(
                style = "height: 650px; overflow-y: auto; overflow-x: auto;",
                shiny::plotOutput(outputId = ns("modal_heatmap"), height = "auto", width = "100%")
              )
            )
          )
        ),
        shiny::hr(),
        shiny::div(
          class = "white-box",
          shiny::h4("Leading Edge 基因统计表 (连续滚动)"),
          DT::dataTableOutput(outputId = ns("modal_gene_table"))
        ),
        footer = shiny::modalButton("关闭")
      ))
    })
  })
}




# ==================== 模块5: 联合绘图控制台 (修订版) ====================
#' @title 联合绘图控制台UI模块
#' @description 展示多通路联合GSEA图（控制参数已移至模块1）
#' @param id 模块命名空间ID
#' @return Shiny UI元素
modMultiPlotUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::div(
      style = "margin-top: 20px; padding: 15px; background-color: #f8f9fa; border-radius: 8px;",
      shiny::h4("🎨 联合通路绘图", style = "color: #333; margin-bottom: 15px;"),
      shiny::uiOutput(outputId = ns("selection_info")),

      shiny::selectizeInput(
        inputId = ns("pathway_selector"),
        label = "当前已选择通路（可直接点 × 删除）:",
        choices = character(0),
        multiple = TRUE,
        options = list(
          plugins = list("remove_button"),
          placeholder = "当前勾选的通路会显示在这里，可直接删除",
          persist = FALSE,
          create = FALSE
        )
      ),

      shiny::actionButton(
        inputId = ns("clear_selected_pathways"),
        label = "🗑️ 清空当前选择通路",
        class = "btn-warning",
        style = "margin-bottom: 15px;"
      ),

      shiny::plotOutput(outputId = ns("multi_plot"), height = "600px")
    )
  )
}

#' @title 联合绘图控制台Server模块
#' @description 根据模块1的配置和模块2的选中通路执行联合绘图
#' @param id 模块ID
#' @param data_prep_module 来自模块1的数据流（包含plot_subtype和custom_colors）
#' @param selected_ids 来自模块2的选中通路ID向量（reactive）
modMultiPlotServer <- function(id, data_prep_module, selection_controller) {
  shiny::moduleServer(id, function(input, output, session) {

    selected_ids <- selection_controller$selected_pathways
    updating_selector <- shiny::reactiveVal(FALSE)

    # 同步当前已选通路到 selectizeInput
    shiny::observe({
      sel <- selected_ids()

      updating_selector(TRUE)
      on.exit(updating_selector(FALSE), add = TRUE)

      shiny::updateSelectizeInput(
        session = session,
        inputId = "pathway_selector",
        choices = sel,
        selected = sel,
        server = TRUE
      )
    })

    # 用户在 selectize 中手动删除通路
    shiny::observeEvent(input$pathway_selector, {
      if (isTRUE(updating_selector())) return()

      current_sel <- selected_ids()
      editor_sel <- input$pathway_selector
      if (is.null(editor_sel)) editor_sel <- character(0)

      removed_ids <- setdiff(current_sel, editor_sel)

      if (length(removed_ids) > 0) {
        selection_controller$remove_pathways(removed_ids)
        shiny::showNotification(
          paste("已移除通路:", paste(removed_ids, collapse = ", ")),
          type = "message",
          duration = 3
        )
      }
    }, ignoreInit = TRUE)

    # 一键清空
    shiny::observeEvent(input$clear_selected_pathways, {
      selection_controller$clear_selection()
      shiny::showNotification("已清空当前联合绘图通路", type = "message", duration = 3)
    })

    # 显示当前选中状态
    output$selection_info <- shiny::renderUI({
      sel <- selected_ids()
      count <- length(sel)

      if (count == 0) {
        shiny::tags$div(
          style = "color: #856404; background-color: #fff3cd; padding: 10px; border-radius: 4px; margin-bottom: 10px;",
          "⚠️ 请在工作台表格中勾选“联合展示”列以添加通路；添加后可在下方手动删除"
        )
      } else {
        shiny::tags$div(
          style = "color: #155724; background-color: #d4edda; padding: 10px; border-radius: 4px; margin-bottom: 10px;",
          sprintf("✅ 已选择 %d 条通路进行联合绘图，可在下方删除单条或全部清空", count)
        )
      }
    })

    # 联合绘图主输出
    output$multi_plot <- shiny::renderPlot({
      sel_ids <- selected_ids()
      shiny::req(sel_ids, length(sel_ids) > 0)

      data_list <- data_prep_module()
      shiny::req(data_list)

      plot_subtype <- data_list$plot_subtype
      custom_colors <- data_list$custom_colors

      color_vector <- trimws(strsplit(custom_colors, ",")[[1]])
      if (length(color_vector) < length(sel_ids)) {
        color_vector <- rep(color_vector, length.out = length(sel_ids))
      }

      contrast_display <- gsub("_vs_", " vs ", data_list$contrast_name)

      tryCatch({
        if (exists("plot_directional_gsea", mode = "function")) {
          p <- plot_directional_gsea(
            directional_gsea_obj = list(
              gsea_res = data_list$gsea_res,
              meta = list(
                left_group = data_list$left_group,
                right_group = data_list$right_group
              )
            ),
            target_pathways = sel_ids,
            subPlot = plot_subtype,
            curveCol = color_vector,
            main_title = sprintf("联合展示: %d 条通路 | %s", length(sel_ids), contrast_display)
          )
          print(p)
        } else {
          graphics::par(mfrow = c(2, 2))
          for (i in seq_len(min(4, length(sel_ids)))) {
            graphics::plot(
              1:10, stats::rnorm(10), type = "l", col = color_vector[i],
              main = sel_ids[i], xlab = "Rank", ylab = "ES"
            )
          }
        }
      }, error = function(error_obj) {
        graphics::par(mar = c(2, 2, 2, 2))
        graphics::plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        graphics::text(
          1, 1,
          labels = sprintf("绘图错误: %s", error_obj$message),
          col = "red", cex = 1.2
        )
      })
    })
  })
}


# ==================== 主应用组装 ====================
#' @export
launch_gsea_app <- function(res_capsule) {

  if (!inherits(res_capsule, "GseaResPro")) {
    stop("传入的对象不是标准的 GseaResPro 计算胶囊！")
  }

  ui <- shiny::fluidPage(
    shiny::tags$head(
      shiny::tags$style(shiny::HTML("
        .modal-dialog { max-width: 85vw !important; width: 85vw !important; }
        .modal-body { min-height: 80vh; overflow-y: auto; background-color: #fafafa; }
        .white-box {
          background-color: white;
          padding: 15px;
          border-radius: 8px;
          box-shadow: 0 2px 8px rgba(0,0,0,0.1);
          margin-bottom: 15px;
        }
        .master-table-container { width: 100%; overflow-x: auto; }
        /* Description自动换行样式 */
        .description-cell {
          white-space: normal !important;
          word-wrap: break-word !important;
          max-width: 300px;
          font-size: 12px;
          line-height: 1.4;
        }
        /* 复选框列居中 */
        .dt-center { text-align: center !important; }
        /* 侧边栏样式优化 */
        .well {
          background-color: #f8f9fa;
          border: 1px solid #e9ecef;
          border-radius: 8px;
          padding: 15px;
        }
      "))
    ),
    shiny::titlePanel("🧬 GSEAlens PRO 3.0: 模块化全息分析工作站"),

    shiny::sidebarLayout(
      shiny::sidebarPanel(
        width = 3,
        # 模块1：数据预处理（已集成联合绘图控制）
        modDataPrepUI("data_prep")
      ),
      shiny::mainPanel(
        width = 9,
        shiny::tabsetPanel(
          shiny::tabPanel(
            title = "📊 主工作台",
            shiny::br(),
            # 模块2：主工作台表格
            modMasterTableUI("master_table"),
            shiny::hr(),
            # 模块5：联合绘图控制台（显示在表格下方）
            modMultiPlotUI("multi_plot")
          ),
          shiny::tabPanel(
            title = "🌋 全息四重联动",
            shiny::br(),
            # 模块3：四重联动图表
            modQuadrantUI("quadrant")
          )
        )
      )
    ),
    # 模块4：Dashboard弹窗占位（内容由server动态控制）
    modPathwayModalUI("pathway_modal")
  )

  server <- function(input, output, session) {

    # 1. 数据预处理模块（全局共享，包含联合绘图配置）
    data_prep <- modDataPrepServer("data_prep", res_capsule)

    # 2. 主表格模块（返回列表：包含 selected_pathways 和 show_modal）
    master_table_result <- modMasterTableServer("master_table", data_prep)

    # 提取 reactive（这样命名保持清晰）
    selected_pathways <- master_table_result$selected_pathways
    show_modal_trigger <- master_table_result$show_modal

    # 3. 联合绘图模块（传入具体的 reactive，而非整个列表）
    modMultiPlotServer("multi_plot", data_prep, master_table_result)

    # 4. 四重联动模块
    modQuadrantServer("quadrant", data_prep, res_capsule)

    # 5. Dashboard弹窗模块（直接使用提取的 reactive）
    modPathwayModalServer("pathway_modal", data_prep, show_modal_trigger, res_capsule)

    # 可选：保留对 joint_plot_toggle 的监听（用于调试）
    shiny::observeEvent(input[["master_table-joint_plot_toggle"]], {
      message("[DEBUG] Joint plot toggle event received")
      invisible(NULL)
    })
  }

  message("🚀 GSEAlens PRO 3.0 模块化版本已启动")
  message("架构组成: DataPrep | MasterTable | QuadrantPlot | PathwayModal | MultiPlot")
  message("版本特性:")
  message("  ✓ 动态对比组切换（IBAA vs PREA格式显示）")
  message("  ✓ 多表达量类型支持（CPM/FPKM/TPM及其log2）")
  message("  ✓ 联合绘图控制集成于左侧栏")
  message("  ✓ 工作台列宽可调（ColReorder）")
  message("  ✓ 联合展示复选框排序筛选")
  shiny::shinyApp(ui, server)
}
