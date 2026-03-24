#' @title 数据预处理模块 UI（联合画布控制增强版）
#' @keywords internal

mod_data_prep_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::selectInput(
      ns("selected_contrast"),
      label = "选择对比组 (Contrast)",
      choices = NULL
    ),
    shiny::h4("🎯 数据切片"),
    shiny::selectizeInput(
      ns("selected_collections"),
      label = "选择基因集亚组:",
      choices = NULL,
      multiple = TRUE,
      options = list(plugins = list("remove_button"), placeholder = "选择亚组或保留ALL")
    ),
    shiny::selectInput(
      ns("sort_by"),
      label = "全局排序策略:",
      choices = c(
        "按 NES 绝对值 (降序)" = "abs_nes_desc",
        "按 NES (降序)" = "nes_desc",
        "按 P-value (升序)" = "pval_asc"
      ),
      selected = "abs_nes_desc"
    ),
    shiny::actionButton(
      ns("run_btn"),
      label = "🚀 确认对比组与富集基因集",
      class = "btn-success",
      style = "width: 100%; font-weight: bold; margin-top: 15px; font-size: 16px;"
    ),

    shiny::helpText(
      style = "margin-top: 10px; text-align: center; color: #28a745;",
      "💡 修改设置后必须点击此按钮"
    ),

    shiny::hr(),

    # 🔧 新增：联合GSEA画布控制（完全移入侧边栏）
    shiny::h4("🖼️ 联合GSEA画布"),
    shiny::selectizeInput(
      ns("joint_contrasts"),
      label = "选择对比组（多选，支持排列）:",
      choices = NULL,
      multiple = TRUE,
      options = list(
        plugins = list("remove_button"),
        placeholder = '选择对比组...',
        maxItems = 999
      )
    ),
    shiny::numericInput(
      ns("joint_ncol"),
      "每行数量 (ncol):",
      value = 2,
      min = 1,
      max = 10,
      step = 1
    ),
    shiny::actionButton(
      ns("joint_generate"),
      "🎨 生成/更新多通路画布",
      class = "btn-success",
      style = "width: 100%; font-weight: bold; margin-top: 10px;"
    ) ,
    shiny::hr(),

    shiny::h4("🎯 差异表达基因标记"),
    shiny::div(
      style = "background-color: #e8f4f8; padding: 10px; border-radius: 5px; border-left: 4px solid #17a2b8;",
      shiny::selectizeInput(
        ns("pending_genes"),
        label = "选择感兴趣基因（支持多选）:",
        choices = NULL,
        multiple = TRUE,
        options = list(
          plugins = list('remove_button'),
          placeholder = '输入基因名（如 TP53）',
          maxItems = 50,
          closeAfterSelect = FALSE,
          selectOnTab = TRUE
        )
      ),
      shiny::actionButton(
        ns("apply_genes_btn"),
        label = "🎯 确认应用基因标记",
        class = "btn-info",
        style = "width: 100%; margin-top: 10px; font-weight: bold;"
      ),
      shiny::helpText(
        style = "margin-top: 8px; color: #666;",
        "选择基因后点击确认，避免实时刷新"
      )
    ),

    shiny::div(
      style = "background-color: #d4edda; padding: 10px; border-radius: 5px; margin-top: 10px;",
      shiny::HTML(
        "<strong style='color: #155724;'>💾 已确认基因标记：</strong><br>
    <small style='color: #666;'>切换对比组后仍然保留</small>"
      ),
      shiny::uiOutput(ns("confirmed_genes_display"))
    )
    ,
    shiny::uiOutput(ns("applied_genes_display")),
    shiny::hr(),

    shiny::h4("🔄 选择组别显示顺序"),
    shiny::div(
      style = "background-color: #fff3cd; padding: 10px; border-radius: 5px; border-left: 4px solid #ffc107;",
      shiny::selectInput(
        ns("boxplot_order_pending"),
        label = "组别排序选项:",
        choices = c("默认顺序" = "default"),
        selected = "default"
      ),
      shiny::actionButton(
        ns("apply_order_btn"),
        label = "✓ 确认排序（刷新箱线图）",
        class = "btn-warning",
        style = "width: 100%; margin-top: 10px; font-weight: bold;"
      ),
      # ✅ 新增：实时排序状态显示
      shiny::uiOutput(ns("order_status_display")),
      shiny::helpText(
        style = "margin-top: 8px; color: #666;",
        "选择顺序后点击确认，箱线图将自动排序"
      )
    ),
    shiny::hr(),

    shiny::h4("🎨 联合绘图控制"),
    shiny::selectInput(
      ns("plot_subtype"),
      label = "GSEAvis 样式:",
      choices = c("1: 仅经典富集" = "1", "2: 富集+热图带" = "2", "3: 完整带Rank" = "3"),
      selected = "3"
    ),
    shiny::textInput(
      ns("custom_colors"),
      label = "多通路自定义颜色:",
      value = "#E41A1C, #377EB8, #4DAF4A, #984EA3",
      placeholder = "例如: #FF0000, #00FF00"
    ),
    shiny::hr(),

    # 🔧 修复：动态表达数据类型（根据后端更新）
    shiny::h4("📊 表达量度量"),
    shiny::selectInput(
      ns("expression_type"),
      label = "选择表达数据类型:",
      choices = NULL,  # 动态填充
      selected = NULL
    )
  )
}

#' @title 数据预处理模块 Server（联合画布参数传递版）
#' @keywords internal

mod_data_prep_server <- function(id, gsea_res) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    result_data <- shiny::reactiveVal(NULL)
    has_initialized <- shiny::reactiveVal(FALSE)

    applied_genes <- shiny::reactiveVal(character(0))
    applied_boxplot_order <- shiny::reactiveVal("default")

    # 🔧 新增：联合画布参数缓存
    canvas_contrasts_val <- shiny::reactiveVal(character(0))
    canvas_ncol_val <- shiny::reactiveVal(3)

    current_contrast_cache <- shiny::reactiveVal(NULL)
    pending_genes_internal <- shiny::reactiveVal(character(0))

    # 🔧 修复：根据后端类型动态设置表达数据类型选项
    shiny::observe({
      backend <- gsea_res$backend_info$backend

      if (backend == "limma_voom") {
        # Limma-voom流程：支持CPM、log2CPM，不支持VST
        choices <- c(
          "log2(CPM)" = "logcpm",
          "CPM (原始)" = "cpm",
          "log2(FPKM)" = "logfpkm",
          "FPKM (原始)" = "fpkm"
        )
        selected <- "logcpm"
        message("📊 检测到Limma-voom后端，加载对应表达类型")
      } else if (backend == "deseq2") {
        # DESeq2流程：支持VST、CPM等
        choices <- c(
          "log2(CPM)" = "logcpm",
          "CPM (原始)" = "cpm",
          "VST (方差稳定变换)" = "vst",
          "log2(Normalized Counts + 1)" = "lognorm"
        )
        selected <- "logcpm"
        message("📊 检测到DESeq2后端，加载对应表达类型")
      } else {
        choices <- c("log2(CPM)" = "logcpm", "CPM (原始)" = "cpm")
        selected <- "logcpm"
      }

      shiny::updateSelectInput(
        session,
        "expression_type",
        choices = choices,
        selected = selected
      )
    })

    # 初始化对比组（包含排列：正向+反向）
    shiny::observe({
      registry <- gsea_res$contrast_registry
      all_contrasts <- list()

      # 构建排列：包含正向和反向
      for (i in 1:nrow(registry)) {
        row <- registry[i, ]
        # 正向
        all_contrasts[[row$contrast_id]] <- paste(row$left_group, "vs", row$right_group)
        # 反向（排列）
        rev_id <- paste(row$right_group, row$left_group, sep = "_vs_")
        rev_display <- paste(row$right_group, "vs", row$left_group, "(反向)")
        all_contrasts[[rev_id]] <- rev_display
      }

      choices <- setNames(names(all_contrasts), unlist(all_contrasts))

      # 更新单选（主工作台）
      shiny::updateSelectInput(
        session,
        "selected_contrast",
        choices = choices,
        selected = registry$contrast_id[1]
      )

      # 🔧 更新多选（联合画布）- 默认全选所有排列
      shiny::updateSelectizeInput(
        session,
        "joint_contrasts",
        choices = choices,
        selected = names(all_contrasts)[1:min(4, length(all_contrasts))]
      )
    })

    # 基因列表更新（根治版：已确认基因独立维护）
    shiny::observeEvent(input$selected_contrast, {
      contrast_id <- input$selected_contrast
      shiny::req(contrast_id)

      if (!is.null(current_contrast_cache()) && current_contrast_cache() == contrast_id) {
        return()
      }

      message(sprintf("🔄 切换对比组: %s", contrast_id))
      current_contrast_cache(contrast_id)

      tryCatch({
        de_df <- get_de_table(gsea_res, contrast_id)
        if (!is.null(de_df) && "gene_symbol" %in% colnames(de_df)) {
          gene_choices <- de_df$gene_symbol
          gene_choices <- gene_choices[!is.na(gene_choices)]
          gene_choices <- sort(unique(gene_choices))

          # ✅ 关键：切换对比组时不清空 UI 勾选；把"已确认基因"回填到 pending_genes
          current_applied <- applied_genes()  # 已确认基因（永不清空）

          # 由于当前对比组的 gene_choices 可能不包含全部 applied 基因，所以取交集（大小写不敏感）
          applied_upper <- toupper(current_applied)
          choices_upper <- toupper(gene_choices)

          selected_upper <- applied_upper[applied_upper %in% choices_upper]
          pending_selected <- gene_choices[choices_upper %in% selected_upper]

          pending_selected <- unique(pending_selected)

          shiny::updateSelectizeInput(
            session,
            "pending_genes",
            choices = gene_choices,
            selected = pending_selected,     # <- 不清空：允许用户在框里删改
            server = length(gene_choices) > 1000
          )

          # pending 先跟 applied 同步，让用户删改后点确认才会真正刷新 applied
          pending_genes_internal(pending_selected)

          message(sprintf(
            "✅ [对比切换] pending 已回填 applied（当前 choices 内 %d 个；applied 总数 %d 个）",
            length(pending_selected), length(current_applied)
          ))

          message(sprintf("✅ 已确认的基因标记保持不变（%d个）", length(current_applied)))
        }
      }, error = function(e) {
        message("更新基因列表失败: ", e$message)
      })
    }, ignoreInit = FALSE)

    # 生成排列组合
    generate_all_permutations <- function(groups) {
      n <- length(groups)
      if (n <= 1) return(list(default = "默认顺序"))

      perms <- list()
      permute <- function(arr, l, r) {
        if (l == r) {
          perm_str <- paste(arr, collapse = ",")
          label_str <- paste(arr, collapse = "→")
          perms[[perm_str]] <<- label_str
          return()
        }
        for (i in l:r) {
          tmp <- arr[l]; arr[l] <- arr[i]; arr[i] <- tmp
          permute(arr, l + 1, r)
          tmp <- arr[l]; arr[l] <- arr[i]; arr[i] <- tmp
        }
      }
      permute(groups, 1, n)
      return(perms)
    }

    generate_limited_perms <- function(groups, max_perms = 100) {
      perms <- list()
      n <- length(groups)
      count <- 0

      perms[[paste(groups, collapse = ",")]] <- paste(groups, collapse = "→")
      count <- count + 1

      perms[[paste(rev(groups), collapse = ",")]] <- paste(rev(groups), collapse = "→")
      count <- count + 1

      for (i in 2:min(n, max_perms/2)) {
        shifted <- groups[c(i:n, 1:(i-1))]
        perms[[paste(shifted, collapse = ",")]] <- paste(shifted, collapse = "→")
        count <- count + 1
        if (count >= max_perms) break
      }

      if (count < max_perms) {
        set.seed(123)
        for (i in 1:(max_perms - count)) {
          shuffled <- sample(groups)
          perm_str <- paste(shuffled, collapse = ",")
          if (!(perm_str %in% names(perms))) {
            perms[[perm_str]] <- paste(shuffled, collapse = "→")
          }
        }
      }

      return(perms)
    }

    # 动态更新Boxplot排序选项
    shiny::observe({
      sample_meta <- tryCatch(get_sample_meta(gsea_res), error = function(e) NULL)
      if (is.null(sample_meta) || !"group" %in% colnames(sample_meta)) return()

      all_groups <- levels(sample_meta$group)
      if (length(all_groups) <= 1) return()

      total_perms <- factorial(length(all_groups))

      if (total_perms > 100) {
        perms <- generate_limited_perms(all_groups, max_perms = 100)
        perms[["limited"]] <- "--- 仅显示前100种 ---"
      } else {
        perms <- generate_all_permutations(all_groups)
      }

      perms <- c(list(default = "默认顺序"), perms)

      shiny::updateSelectInput(
        session,
        "boxplot_order_pending",
        choices = perms,
        selected = "default"
      )
    })

    # 基因确认按钮
    shiny::observeEvent(input$pending_genes, {
      genes <- input$pending_genes
      if (is.null(genes)) genes <- character(0)
      pending_genes_internal(genes)
    }, ignoreNULL = FALSE, ignoreInit = TRUE)

    shiny::observeEvent(input$apply_genes_btn, {
      genes_to_apply <- pending_genes_internal()
      applied_genes(genes_to_apply)
      message(sprintf("🎯 已确认基因标记: %d个", length(genes_to_apply)))
      shiny::showNotification(sprintf("已标记 %d 个基因", length(genes_to_apply)), type = "message", duration = 3)
    })

    # 排序确认按钮（关键修改）
    shiny::observeEvent(input$apply_order_btn, {
      order_to_apply <- input$boxplot_order_pending
      applied_boxplot_order(order_to_apply)  # ✅ 更新 reactiveVal

      # 强制刷新四重联动模块（这是关键！）
      message(sprintf("🔄 [排序确认] 用户选择了: %s | 四重联动模块已触发更新",
                      if(order_to_apply == "default") "默认" else order_to_apply))

      shiny::showNotification(
        sprintf("✅ 排序已应用: %s | 箱线图已刷新",
                if(order_to_apply == "default") "默认顺序" else gsub(",", " → ", order_to_apply)),
        type = "message",
        duration = 2
      )
    }, ignoreInit = TRUE)  # 👈 关键：添加 ignoreInit = TRUE
    # 在 mod_data_prep_server 的关键处添加：
    #shiny::observe({
    #  applied <- applied_genes()
    #  pending <- pending_genes_internal()
    #  contrast <- current_contrast_cache()
    #
    #  message(sprintf(
    #    "📊 [基因管理状态] 对比组=%s | 已确认=%d个 | 待选=%d个",
    #    contrast, length(applied), length(pending)
    #  ))
    #})

    # 显示实时排序状态
    output$order_status_display <- shiny::renderUI({
      pending <- input$boxplot_order_pending
      confirmed <- applied_boxplot_order()

      status_text <- if (pending == "default") {
        "📌 默认顺序"
      } else {
        sprintf("📌 待应用: %s", gsub(",", " → ", pending))
      }

      confirmed_text <- if (confirmed != "default" && !is.na(confirmed) && confirmed != "") {
        sprintf("<br>✅ 已应用: %s", gsub(",", " → ", confirmed))
      } else {
        ""
      }

      shiny::div(
        style = "background: #e7f3ff; padding: 8px; border-radius: 4px; font-size: 12px; color: #004085; margin-top: 10px;",
        shiny::HTML(paste0(status_text, confirmed_text))
      )
    })

    output$confirmed_genes_display <- shiny::renderUI({
      genes <- applied_genes()
      if (length(genes) == 0) {
        shiny::span("（无）", style = "color: #999;")
      } else {
        shiny::div(
          style = "margin-top: 5px; padding: 8px; background: #f1f8f4; border-radius: 3px;",
          shiny::strong(length(genes), "个基因：", style = "color: #28a745;"),
          shiny::br(),
          shiny::span(paste(genes, collapse = ", "), style = "font-size: 12px; color: #555;")
        )
      }
    })

    # 🔧 新增：监听画布参数变化
    shiny::observeEvent(input$canvas_contrasts, {
      canvas_contrasts_val(input$canvas_contrasts)
    }, ignoreNULL = FALSE)

    shiny::observeEvent(input$canvas_ncol, {
      canvas_ncol_val(input$canvas_ncol)
    })



    # 核心数据处理
    process_data_core <- function(contrast_id, collections, sort_by, expr_type, plot_subtype, colors, is_auto = FALSE) {
      if (is.null(contrast_id)) return(NULL)

      task_obj <- tryCatch({
        extract_gsea_task(gsea_res, contrast_id, target_collection = collections)
      }, error = function(e) {
        message(sprintf("提取任务失败: %s", e$message))
        return(NULL)
      })

      if (is.null(task_obj)) return(NULL)

      meta <- task_obj$meta
      gsea_res_obj <- task_obj$gsea_res
      df <- as.data.frame(gsea_res_obj@result)

      if (nrow(df) == 0) return(NULL)

      df$abs_NES <- abs(df$NES)
      sort_config <- switch(sort_by,
                            "nes_desc" = list(col = "NES", desc = TRUE),
                            "abs_nes_desc" = list(col = "abs_NES", desc = TRUE),
                            "pval_asc" = list(col = "pvalue", desc = FALSE),
                            list(col = "abs_NES", desc = TRUE)
      )

      sort_order <- order(df[[sort_config$col]], decreasing = sort_config$desc)
      df <- df[sort_order, ]
      df$Rank <- seq_len(nrow(df))
      df$Safe_ID <- gsub("'", "\\\\'", df$ID)

      reg_row <- gsea_res$contrast_registry[gsea_res$contrast_registry$contrast_id == contrast_id, ]

      if (nrow(reg_row) == 0) {
        parts <- strsplit(contrast_id, "_vs_")[[1]]
        if (length(parts) == 2) {
          left <- parts[1]; right <- parts[2]
        } else {
          left <- contrast_id; right <- "Background"
        }
      } else {
        left <- reg_row$left_group[1]
        right <- reg_row$right_group[1]
      }

      list(
        df = df,
        gsea_res = gsea_res_obj,
        task_obj = task_obj,
        left_group = left,
        right_group = right,
        contrast_id = contrast_id,
        contrast_name = sprintf("%s vs %s", left, right),
        expression_type = expr_type,
        plot_subtype = as.numeric(plot_subtype),
        custom_colors = colors,
        is_preview = is_auto,
        backend = gsea_res$backend_info$backend,
        # 🔧 新增：传递画布参数
        joint_contrasts = canvas_contrasts_val(),
        joint_ncol = canvas_ncol_val()
      )
    }

    # 自动初始化
    shiny::observeEvent(input$selected_contrast, {
      if (has_initialized()) return()

      shiny::invalidateLater(200, session)

      shiny::isolate({
        if (is.null(input$selected_contrast)) return()

        task_temp <- tryCatch(
          extract_gsea_task(gsea_res, input$selected_contrast, "ALL"),
          error = function(e) NULL
        )
        if (is.null(task_temp)) return()

        df_temp <- as.data.frame(task_temp$gsea_res@result)
        available <- unique(c(df_temp$Collection, df_temp$Combo_Name))
        available <- setdiff(available, c(NA, "Unknown"))

        if (length(available) == 0) return()

        counts <- sapply(available, function(x) sum(startsWith(as.character(df_temp$Combo_Name), x)))
        min_col <- available[which.min(counts)[1]]

        shiny::updateSelectizeInput(
          session,
          "selected_collections",
          choices = c("ALL", sort(available)),
          selected = min_col
        )

        auto_data <- process_data_core(
          input$selected_contrast,
          min_col,
          "abs_nes_desc",
          input$expression_type %||% "logcpm",
          "3",
          "#E41A1C, #377EB8, #4DAF4A, #984EA3",
          TRUE
        )

        if (!is.null(auto_data)) {
          result_data(auto_data)
          has_initialized(TRUE)
          shiny::showNotification(sprintf("已自动加载: %s", min_col), type = "default", duration = 3)
        }
      })
    }, ignoreNULL = TRUE, ignoreInit = TRUE)

    # 手动确认按钮
    shiny::observeEvent(input$run_btn, {
      shiny::req(input$selected_contrast)

      cols <- input$selected_collections
      if (is.null(cols) || length(cols) == 0) cols <- "ALL"

      manual_data <- process_data_core(
        input$selected_contrast,
        cols,
        input$sort_by,
        input$expression_type,
        input$plot_subtype,
        input$custom_colors,
        FALSE
      )

      if (is.null(manual_data)) {
        shiny::showNotification("所选配置无数据", type = "error", duration = 5)
        return()
      }

      result_data(manual_data)
      shiny::showNotification("✅ 工作台已更新", type = "message", duration = 3)
    })

    # 🔧 新增：画布生成触发器（事件反应式）
    joint_generate_event <- shiny::eventReactive(input$joint_generate, {
      list(
        contrasts = input$joint_contrasts,
        ncol = input$joint_ncol,
        timestamp = Sys.time()
      )
    }, ignoreNULL = TRUE)

    return(list(
      data = shiny::reactive({ result_data() }),
      highlight_genes = applied_genes,
      boxplot_order = applied_boxplot_order,  # reactiveVal
      # 🔧 新增：联合画布控制项（统一命名）
      joint_contrasts = shiny::reactive({ input$joint_contrasts }),
      joint_ncol = shiny::reactive({ input$joint_ncol }),
      joint_generate = joint_generate_event  # 事件反应式
    ))
  })
}
