# ==============================================================================

# Section: Data Preprocessing Module UI

# ==============================================================================



#' @title Data Preprocessing Module UI

#' @description

#' User interface for the data preprocessing module, providing controls for:

#' \itemize{

#'   \item Contrast selection and gene set subgroup filtering

#'   \item Joint GSEA canvas for multi-pathway visualization

#'   \item Differential expression gene marker selection

#'   \item Group display order customization

#'   \item Expression metrics configuration

#' }

#' @param id Module namespace ID

#' @keywords internal



mod_data_prep_ui <- function(id) {

  ns <- shiny::NS(id)

  shiny::tagList(

    shiny::selectInput(

      ns("selected_contrast"),

      label = "Select Contrast",

      choices = NULL

    ),

    shiny::h4("Data Filtering"),

    shiny::selectizeInput(

      ns("selected_collections"),

      label = "Select Gene Set Subgroup:",

      choices = NULL,

      multiple = TRUE,

      options = list(plugins = list("remove_button"), placeholder = "Select subgroup or keep ALL")

    ),

    shiny::selectInput(

      ns("sort_by"),

      label = "Global Sorting Strategy:",

      choices = c(

        "By NES Absolute Value (Desc)" = "abs_nes_desc",

        "By NES (Desc)" = "nes_desc",

        "By P-value (Asc)" = "pval_asc"

      ),

      selected = "abs_nes_desc"

    ),

    shiny::actionButton(

      ns("run_btn"),

      label = "Confirm Contrast and Gene Sets",

      class = "btn-success",

      style = "width: 100%; font-weight: bold; margin-top: 15px; font-size: 16px;"

    ),

    shiny::helpText(

      style = "margin-top: 10px; text-align: center; color: #28a745;",

      "Click this button after modifying settings"

    ),

    shiny::hr(),

    shiny::h4("Joint GSEA Canvas"),

    shiny::selectizeInput(

      ns("joint_contrasts"),

      label = "Select Contrasts (Multi-select, supports permutations):",

      choices = NULL,

      multiple = TRUE,

      options = list(

        plugins = list("remove_button"),

        placeholder = "Select contrasts...",

        maxItems = 999

      )

    ),

    shiny::numericInput(

      ns("joint_ncol"),

      "Number per row (ncol):",

      value = 2,

      min = 1,

      max = 10,

      step = 1

    ),

    shiny::actionButton(

      ns("joint_generate"),

      "Generate/Update Multi-Pathway Canvas",

      class = "btn-success",

      style = "width: 100%; font-weight: bold; margin-top: 10px;"

    ),

    shiny::hr(),

    shiny::h4("Differential Expression Gene Markers"),

    shiny::div(

      style = "background-color: #e8f4f8; padding: 10px; border-radius: 5px; border-left: 4px solid #17a2b8;",

      shiny::selectizeInput(

        ns("pending_genes"),

        label = "Select Genes of Interest (Multi-select supported):",

        choices = NULL,

        multiple = TRUE,

        options = list(

          plugins = list("remove_button"),

          placeholder = "Enter gene names (e.g., TP53)",

          maxItems = 999,

          closeAfterSelect = FALSE,

          selectOnTab = TRUE

        )

      ),

      shiny::actionButton(

        ns("apply_genes_btn"),

        label = "Confirm and Apply Gene Markers",

        class = "btn-info",

        style = "width: 100%; margin-top: 10px; font-weight: bold;"

      ),

      shiny::actionButton(

        ns("clear_all_genes_btn"),

        label = "Clear All",

        class = "btn-warning",

        style = "width: 100%; margin-top: 10px; font-weight: bold;"

      ),

      shiny::div(

        style = "background-color: #d4edda; padding: 10px; border-radius: 5px; margin-top: 10px;",

        shiny::HTML(

          "<strong style='color: #155724;'>Confirmed Gene Markers:</strong><br>

    <small style='color: #666;'>Retained after switching contrasts</small>"

        ),

        shiny::uiOutput(ns("confirmed_genes_display"))

      ),

      shiny::uiOutput(ns("applied_genes_display")),

      shiny::hr(),

      shiny::h4("Group Display Order"),

      shiny::div(

        style = "background-color: #fff3cd; padding: 10px; border-radius: 5px; border-left: 4px solid #ffc107;",

        shiny::selectInput(

          ns("boxplot_order_pending"),

          label = "Group Ordering Options:",

          choices = c("Default Order" = "default"),

          selected = "default"

        ),

        shiny::actionButton(

          ns("apply_order_btn"),

          label = "Confirm Order (Refresh Boxplot)",

          class = "btn-warning",

          style = "width: 100%; margin-top: 10px; font-weight: bold;"

        ),

        shiny::uiOutput(ns("order_status_display")),

        shiny::helpText(

          style = "margin-top: 8px; color: #666;",

          "Click confirm after selecting order; boxplot will auto-sort"

        )

      ),

      shiny::hr(),

      shiny::h4("Joint Plotting Control"),

      shiny::selectInput(

        ns("plot_subtype"),

        label = "GSEAvis Style:",

        choices = c("1: Classic Enrichment Only" = "1", "2: Enrichment + Heatmap Strip" = "2", "3: Full Strip with Rank" = "3"),

        selected = "3"

      ),

      shiny::textInput(

        ns("custom_colors"),

        label = "Multi-Pathway Custom Colors:",

        value = "#E41A1C, #377EB8, #4DAF4A, #984EA3",

        placeholder = "e.g., #FF0000, #00FF00"

      ),

      shiny::hr(),

      shiny::h4("Expression Metrics"),

      shiny::selectInput(

        ns("expression_type"),

        label = "Select Expression Data Type:",

        choices = NULL,

        selected = NULL

      )

    )

  )

}



# ==============================================================================

# Section: Data Preprocessing Module Server

# ==============================================================================



#' @title Data Preprocessing Module Server

#' @description

#' Server logic for the data preprocessing module. Handles:

#' \itemize{

#'   \item Contrast and gene set selection

#'   \item Data processing and sorting

#'   \item Gene marker management (add, confirm, clear)

#'   \item Boxplot ordering configuration

#'   \item Joint canvas parameter passing

#'   \item Expression data type selection based on backend

#' }

#' @param id Module namespace ID

#' @param gsea_res GSEA results object containing contrast registry and backend info

#' @return A list of reactive functions and values:

#' \itemize{

#'   \item \code{data}: Reactive data object with processed GSEA results

#'   \item \code{highlight_genes}: Reactive character vector of confirmed gene markers

#'   \item \code{boxplot_order}: Reactive value for boxplot group ordering

#'   \item \code{joint_contrasts}: Reactive vector of contrasts for joint canvas

#'   \item \code{joint_ncol}: Reactive number of columns for joint canvas

#'   \item \code{joint_generate}: Event reactive for canvas regeneration trigger

#'   \item \code{update_pending_genes}: Function to update pending genes from modal

#'   \item \code{get_pending_genes}: Function to retrieve current pending genes

#' }

#' @keywords internal



mod_data_prep_server <- function(id, gsea_res) {

  shiny::moduleServer(id, function(input, output, session) {

    ns <- session$ns



    # === 新增：获取物种信息并更新 label ===

    species <- gsea_res$geneset_info$species %||% "HS"



    shiny::observe({

      shiny::updateSelectizeInput(

        session,

        "selected_collections",

        label = sprintf("Select Gene Set Subgroup [%s]:", species)

      )

    })

    # === 新增结束 ===



    result_data <- shiny::reactiveVal(NULL)

    has_initialized <- shiny::reactiveVal(FALSE)



    applied_genes <- shiny::reactiveVal(character(0))

    applied_boxplot_order <- shiny::reactiveVal("default")

    pending_genes_internal <- shiny::reactiveVal(character(0))



    canvas_contrasts_val <- shiny::reactiveVal(character(0))

    canvas_ncol_val <- shiny::reactiveVal(3)



    current_contrast_cache <- shiny::reactiveVal(NULL)

    pending_genes_internal <- shiny::reactiveVal(character(0))



    # 根据后端类型动态设置表达数据类型选项

    shiny::observe({

      backend <- gsea_res$backend_info$backend



      if (backend == "limma_voom") {

        choices <- c(

          "log2(CPM)" = "logcpm",

          "CPM (raw)" = "cpm",

          "log2(FPKM)" = "logfpkm",

          "FPKM (raw)" = "fpkm"

        )

        selected <- "cpm"

        message("Detected Limma-voom backend, loading corresponding expression types")

      } else if (backend == "deseq2") {

        choices <- c(

          "log2(CPM)" = "logcpm",

          "CPM (raw)" = "cpm",

          "VST (Variance Stabilizing Transform)" = "vst",

          "log2(Normalized Counts + 1)" = "lognorm"

        )

        selected <- "cpm"

        message("Detected DESeq2 backend, loading corresponding expression types")

      } else {

        choices <- c("log2(CPM)" = "logcpm", "CPM (raw)" = "cpm")

        selected <- "cpm"

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

      for (i in seq_len(nrow(registry))) {

        row <- registry[i, ]

        # 正向

        all_contrasts[[row$contrast_id]] <- paste(row$left_group, "vs", row$right_group)

        # 反向（排列）

        rev_id <- paste(row$right_group, row$left_group, sep = "_vs_")

        rev_display <- paste(row$right_group, "vs", row$left_group, "(reversed)")

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



      # 更新多选（联合画布）- 默认全选所有排列

      shiny::updateSelectizeInput(

        session,

        "joint_contrasts",

        choices = choices,

        selected = names(all_contrasts)[seq_len(min(4, length(all_contrasts)))]

      )

    })





    # 基因列表更新（根治版：已确认基因独立维护）

    shiny::observeEvent(input$selected_contrast,

      {

        contrast_id <- input$selected_contrast

        shiny::req(contrast_id)



        if (!is.null(current_contrast_cache()) && current_contrast_cache() == contrast_id) {

          return()

        }



        message(sprintf("Switching contrast: %s", contrast_id))

        current_contrast_cache(contrast_id)



        tryCatch(

          {

            de_df <- get_de_table(gsea_res, contrast_id)

            if (!is.null(de_df) && "gene_symbol" %in% colnames(de_df)) {

              gene_choices <- de_df$gene_symbol

              gene_choices <- gene_choices[!is.na(gene_choices)]

              gene_choices <- sort(unique(gene_choices))



              current_applied <- applied_genes()



              applied_upper <- toupper(current_applied)

              choices_upper <- toupper(gene_choices)



              selected_upper <- applied_upper[applied_upper %in% choices_upper]

              pending_selected <- gene_choices[choices_upper %in% selected_upper]



              pending_selected <- unique(pending_selected)



              shiny::updateSelectizeInput(

                session,

                "pending_genes",

                choices = gene_choices,

                selected = pending_selected,

                server = length(gene_choices) > 1000

              )



              pending_genes_internal(pending_selected)



              message(sprintf(

                "Pending refilled from applied (%d in current choices; %d total in applied)",

                length(pending_selected), length(current_applied)

              ))



              message(sprintf("Confirmed gene markers remain unchanged (%d genes)", length(current_applied)))

            }

          },

          error = function(e) {

            message("Failed to update gene list: ", e$message)

          }

        )

      },

      ignoreInit = FALSE

    )



    # 生成排列组合

    generate_all_permutations <- function(groups) {

      n <- length(groups)

      if (n <= 1) {

        return(list(default = "Default Order"))

      }



      # Use an environment as a mutable accumulator to avoid the `<<-` super-assignment.

      # Semantics are equivalent to the previous recursive <<- pattern, but the

      # scope of the mutation is explicit and does not leak into parent frames.

      state <- new.env(parent = emptyenv())

      state$perms <- list()



      permute <- function(arr, l, r) {

        if (l == r) {

          perm_str <- paste(arr, collapse = ",")

          label_str <- paste(arr, collapse = "->")

          state$perms[[perm_str]] <- label_str

          return()

        }

        for (i in l:r) {

          tmp <- arr[l]

          arr[l] <- arr[i]

          arr[i] <- tmp

          permute(arr, l + 1, r)

          tmp <- arr[l]

          arr[l] <- arr[i]

          arr[i] <- tmp

        }

      }

      permute(groups, 1, n)

      return(as.list(state$perms))

    }



    generate_limited_perms <- function(groups, max_perms = 100, seed = 123) {

      perms <- list()

      n <- length(groups)

      count <- 0



      perms[[paste(groups, collapse = ",")]] <- paste(groups, collapse = "->")

      count <- count + 1



      perms[[paste(rev(groups), collapse = ",")]] <- paste(rev(groups), collapse = "->")

      count <- count + 1



      for (i in 2:min(n, max_perms / 2)) {

        shifted <- groups[c(i:n, seq_len(i - 1))]

        perms[[paste(shifted, collapse = ",")]] <- paste(shifted, collapse = "->")

        count <- count + 1

        if (count >= max_perms) break

      }



      if (count < max_perms) {

        # Local RNG scope: seed is required for reproducible random permutation

        # sampling when max_perms exceeds deterministic (shift-based) permutations.

        # The seed is exposed to users via the `seed` parameter (default 123).

        # withr::local_seed restores the caller's RNG state on function exit.

        withr::local_seed(seed)

        for (i in seq_len(max_perms - count)) {

          shuffled <- sample(groups)

          perm_str <- paste(shuffled, collapse = ",")

          if (!(perm_str %in% names(perms))) {

            perms[[perm_str]] <- paste(shuffled, collapse = "->")

          }

        }

      }



      return(perms)

    }



    # 动态更新Boxplot排序选项

    shiny::observe({

      sample_meta <- tryCatch(get_sample_meta(gsea_res), error = function(e) NULL)

      if (is.null(sample_meta) || !"group" %in% colnames(sample_meta)) {

        return()

      }



      all_groups <- levels(sample_meta$group)

      if (length(all_groups) <= 1) {

        return()

      }



      total_perms <- factorial(length(all_groups))



      if (total_perms > 100) {

        perms <- generate_limited_perms(all_groups, max_perms = 100)

        perms[["limited"]] <- "--- Showing first 100 only ---"

      } else {

        perms <- generate_all_permutations(all_groups)

      }



      perms <- c(list(default = "Default Order"), perms)



      shiny::updateSelectInput(

        session,

        "boxplot_order_pending",

        choices = perms,

        selected = "default"

      )

    })



    # 基因确认按钮

    # 关键修复：input$pending_genes observer应该追加而非覆盖

    # 当用户在输入框手动输入时，只更新pending_genes_internal

    # 不再无条件覆盖Modal添加的基因

    shiny::observeEvent(input$pending_genes,

      {

        genes <- input$pending_genes

        if (is.null(genes)) genes <- character(0)



        # 清理基因名

        genes_clean <- toupper(trimws(genes))

        genes_clean <- genes_clean[genes_clean != ""]



        # 获取当前pending

        current_pending <- pending_genes_internal()



        # 判断：是追加还是覆盖？

        # 如果用户输入的基因不在当前pending中，认为是手动追加

        # 如果用户输入的基因与当前pending完全不同（通过手动清空后输入），认为是覆盖

        if (length(genes_clean) == 0) {

          # 用户清空了输入框，保持pending不变

          pending_genes_internal(current_pending)

        } else if (length(current_pending) == 0) {

          # 当前pending为空，直接设置

          pending_genes_internal(genes_clean)

        } else {

          # 当前pending有值，比较差异

          # 如果输入的内容比pending多很多，认为是手动输入（覆盖）

          # 否则认为是追加

          genes_not_in_pending <- setdiff(genes_clean, current_pending)

          pending_not_in_genes <- setdiff(current_pending, genes_clean)



          # 如果大部分pending都在输入中，且输入中有新基因，认为是追加

          if (length(pending_not_in_genes) < length(current_pending) * 0.5 && length(genes_not_in_pending) > 0) {

            # 追加模式：保留原有 + 添加新的

            pending_genes_internal(unique(c(current_pending, genes_clean)))

          } else {

            # 覆盖模式：直接替换

            pending_genes_internal(genes_clean)

          }

        }



        message(sprintf("[PendingInput] User input updated: %d genes", length(genes_clean)))

      },

      ignoreNULL = FALSE,

      ignoreInit = TRUE

    )



    shiny::observeEvent(input$apply_genes_btn, {

      genes_to_apply <- pending_genes_internal()



      # 合并：已有的 + 新输入的

      current_applied <- applied_genes()

      new_applied <- union(toupper(current_applied), toupper(genes_to_apply))

      applied_genes(new_applied)



      message(sprintf(

        "Applied genes: %d (merged: %d existing + %d new)",

        length(new_applied), length(current_applied), length(genes_to_apply)

      ))

      shiny::showNotification(sprintf("Applied genes: %d (merged)", length(new_applied)),

        type = "message", duration = 3

      )

    })



    # 【新增】Clear All 按钮

    shiny::observeEvent(input$clear_all_genes_btn, {

      applied_genes(character(0))

      pending_genes_internal(character(0))

      shiny::updateSelectizeInput(

        session,

        "pending_genes",

        selected = character(0)

      )

      message("Cleared all applied genes")

      shiny::showNotification("Cleared all gene markers", type = "message", duration = 2)

    })



    # 排序确认按钮

    shiny::observeEvent(input$apply_order_btn,

      {

        order_to_apply <- input$boxplot_order_pending

        applied_boxplot_order(order_to_apply)



        message(sprintf(

          "Order Confirm: User selected: %s | Quad-link module triggered update",

          if (order_to_apply == "default") "default" else order_to_apply

        ))



        shiny::showNotification(

          sprintf(

            "Order applied: %s | Boxplot refreshed",

            if (order_to_apply == "default") "Default order" else gsub(",", " -> ", order_to_apply)

          ),

          type = "message",

          duration = 2

        )

      },

      ignoreInit = TRUE

    )



    # 显示实时排序状态

    output$order_status_display <- shiny::renderUI({

      pending <- input$boxplot_order_pending

      confirmed <- applied_boxplot_order()



      status_text <- if (pending == "default") {

        "Default Order"

      } else {

        sprintf("Pending: %s", gsub(",", " -> ", pending))

      }



      confirmed_text <- if (confirmed != "default" && !is.na(confirmed) && confirmed != "") {

        sprintf("<br>Applied: %s", gsub(",", " -> ", confirmed))

      } else {

        ""

      }



      shiny::div(

        style = "background: #e7f3ff; padding: 8px; border-radius: 4px; font-size: 12px; color: #004085; margin-top: 10px;",

        shiny::HTML(paste0(status_text, confirmed_text))

      )

    })



    # ═══════════════════════════════════════════════════════════════════════════

    # 修复：显示原始大小写的基因名

    # ═══════════════════════════════════════════════════════════════════════════

    output$confirmed_genes_display <- shiny::renderUI({

      genes <- applied_genes()



      if (length(genes) == 0) {

        shiny::span("(None)", style = "color: #999;")

      } else {

        # 按需重建 symbol_map（获取原始大小写）

        current_contrast <- input$selected_contrast

        if (!is.null(current_contrast) && current_contrast != "") {

          symbol_map <- .rebuild_symbol_map(gsea_res, current_contrast)

          display_genes <- .get_display_symbols(genes, symbol_map)

        } else {

          display_genes <- genes

        }



        shiny::div(

          style = "margin-top: 5px; padding: 8px; background: #f1f8f4; border-radius: 3px;",

          shiny::strong(length(genes), " genes:", style = "color: #28a745;"),

          shiny::br(),

          shiny::span(paste(display_genes, collapse = ", "), style = "font-size: 12px; color: #555;")

        )

      }

    })



    # 监听画布参数变化

    shiny::observeEvent(input$canvas_contrasts,

      {

        canvas_contrasts_val(input$canvas_contrasts)

      },

      ignoreNULL = FALSE

    )



    shiny::observeEvent(input$canvas_ncol, {

      canvas_ncol_val(input$canvas_ncol)

    })





    # 核心数据处理

    process_data_core <- function(contrast_id, collections, sort_by, expr_type, plot_subtype, colors, is_auto = FALSE) {

      if (is.null(contrast_id)) {

        return(NULL)

      }



      task_obj <- tryCatch(

        {

          extract_gsea_task(gsea_res, contrast_id, target_collection = collections)

        },

        error = function(e) {

          message(sprintf("Task extraction failed: %s", e$message))

          return(NULL)

        }

      )



      if (is.null(task_obj)) {

        return(NULL)

      }



      meta <- task_obj$meta

      gsea_res_obj <- task_obj$gsea_res

      df <- as.data.frame(gsea_res_obj@result)



      if (nrow(df) == 0) {

        return(NULL)

      }



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

          left <- parts[1]

          right <- parts[2]

        } else {

          left <- contrast_id

          right <- "Background"

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

        joint_contrasts = canvas_contrasts_val(),

        joint_ncol = canvas_ncol_val(),

        collections = collections

      )

    }



    # 自动初始化

    shiny::observeEvent(input$selected_contrast,

      {

        if (has_initialized()) {

          return()

        }



        shiny::invalidateLater(200, session)



        shiny::isolate({

          if (is.null(input$selected_contrast)) {

            return()

          }



          task_temp <- tryCatch(

            extract_gsea_task(gsea_res, input$selected_contrast, "ALL"),

            error = function(e) NULL

          )

          if (is.null(task_temp)) {

            return()

          }



          df_temp <- as.data.frame(task_temp$gsea_res@result)

          available <- unique(c(df_temp$Collection, df_temp$Combo_Name))

          available <- setdiff(available, c(NA, "Unknown"))



          if (length(available) == 0) {

            return()

          }



          # Sort collections: Hallmark (H/MH) first, then by main code

          # (C1/M1, C2/M2, ...) in natural numeric order; within the same

          # main code the main collection precedes its sub-collections

          # (e.g. C2 before C2:CGP, C2 before C2:CP:REACTOME).

          avail_chr <- as.character(available)

          main_code <- sub(":.*$", "", avail_chr)

          coll_num <- ifelse(

            main_code %in% c("H", "MH"), 0L,

            suppressWarnings(as.integer(gsub("[^0-9]", "", main_code)))

          )

          coll_num[is.na(coll_num)] <- 99L

          has_sub <- grepl(":", avail_chr)

          sub_part <- ifelse(has_sub, sub("^[^:]+:", "", avail_chr), "")

          ord <- order(coll_num, has_sub, sub_part)

          sorted_available <- avail_chr[ord]



          # Default to Hallmark when present; otherwise fall back to the

          # collection with the fewest entries (original heuristic).

          hallmark <- intersect(c("H", "MH"), sorted_available)

          if (length(hallmark) > 0) {

            min_col <- hallmark[1]

          } else {

            counts <- vapply(sorted_available, function(x) {

              sum(startsWith(as.character(df_temp$Combo_Name), x))

            }, numeric(1))

            min_col <- sorted_available[which.min(counts)[1]]

          }



          shiny::updateSelectizeInput(

            session,

            "selected_collections",

            choices = c("ALL", sorted_available),

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

            shiny::showNotification(sprintf("Auto-loaded: %s", min_col), type = "default", duration = 3)

          }

        })

      },

      ignoreNULL = TRUE,

      ignoreInit = TRUE

    )



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

        shiny::showNotification("No data for selected configuration", type = "error", duration = 5)

        return()

      }



      result_data(manual_data)

      shiny::showNotification("Workspace updated", type = "message", duration = 3)

    })



    # 画布生成触发器

    joint_generate_event <- shiny::eventReactive(input$joint_generate,

      {

        list(

          contrasts = input$joint_contrasts,

          ncol = input$joint_ncol,

          timestamp = Sys.time()

        )

      },

      ignoreNULL = TRUE

    )



    # Modal传递基因接口（追加模式）- 完全重写

    update_pending_genes <- function(new_genes_from_modal) {

      # 防御性检查

      if (is.null(new_genes_from_modal) || length(new_genes_from_modal) == 0) {

        message("[UpdatePending] No new genes from modal, skipping")

        return(invisible(NULL))

      }



      # 清理基因名：去除空白 + 转大写

      new_genes_clean <- toupper(trimws(as.character(new_genes_from_modal)))

      new_genes_clean <- new_genes_clean[!is.na(new_genes_clean)]

      new_genes_clean <- new_genes_clean[new_genes_clean != ""]

      new_genes_clean <- unique(new_genes_clean)



      # 直接获取当前pending_genes_internal的值

      current_pending <- tryCatch(

        {

          shiny::isolate(pending_genes_internal())

        },

        error = function(e) {

          message(sprintf("[UpdatePending] Error: %s", e$message))

          character(0)

        }

      )



      # 追加合并（不是替换）

      updated_pending <- unique(c(current_pending, new_genes_clean))



      # 更新reactive值

      pending_genes_internal(updated_pending)



      message(sprintf(

        "[UpdatePending] Merged: %d existing + %d new = %d total (pending)",

        length(current_pending), length(new_genes_clean), length(updated_pending)

      ))



      # 获取当前contrast的基因choices

      gene_choices_upper <- character(0)



      tryCatch(

        {

          contrast_id <- shiny::isolate(input$selected_contrast)

          if (!is.null(contrast_id) && contrast_id != "") {

            de_df <- get_de_table(gsea_res, contrast_id)

            if (!is.null(de_df) && "gene_symbol" %in% colnames(de_df)) {

              gene_choices <- de_df$gene_symbol

              gene_choices <- gene_choices[!is.na(gene_choices)]

              gene_choices <- trimws(gene_choices)

              gene_choices <- gene_choices[gene_choices != ""]

              gene_choices_upper <- toupper(sort(unique(gene_choices)))

            }

          }

        },

        error = function(e) {

          message(sprintf("[UpdatePending] Error getting gene choices: %s", e$message))

        }

      )



      # 如果没有choices，使用pending作为choices

      if (length(gene_choices_upper) == 0) {

        gene_choices_upper <- updated_pending

      }



      # 刷新selectizeInput

      tryCatch(

        {

          shiny::updateSelectizeInput(

            session,

            "pending_genes",

            choices = gene_choices_upper,

            selected = updated_pending,

            options = list(

              plugins = list("remove_button"),

              placeholder = "Enter gene names (e.g., TP53)",

              maxItems = 999,

              closeAfterSelect = FALSE,

              selectOnTab = TRUE

            ),

            server = FALSE

          )

          message("[UpdatePending] UI refreshed")

        },

        error = function(e) {

          message(sprintf("[UpdatePending] UI update error: %s", e$message))

        }

      )



      invisible(NULL)

    }





    return(list(

      data = shiny::reactive({

        result_data()

      }),

      highlight_genes = applied_genes,

      boxplot_order = applied_boxplot_order,

      joint_contrasts = shiny::reactive({

        input$joint_contrasts

      }),

      joint_ncol = shiny::reactive({

        input$joint_ncol

      }),

      joint_generate = joint_generate_event,

      update_pending_genes = update_pending_genes,

      get_pending_genes = function() {

        shiny::isolate(pending_genes_internal())

      }

    ))

  })

}

