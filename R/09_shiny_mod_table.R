#' @title Master Workspace Table UI
#' @description Interactive data table module for displaying GSEA pathway results with checkbox selection and modal integration.
#' @keywords internal

mod_master_table_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::tags$head(
      shiny::tags$script(shiny::HTML(sprintf("
        Shiny.addCustomMessageHandler('%s', function(message) {
          var ids = message.ids;
          var ns = message.ns;
          var checkboxes = document.querySelectorAll('.joint-plot-checkbox');
          checkboxes.forEach(function(cb) {
            var id = cb.getAttribute('data-id');
            var shouldBeChecked = ids.includes(id);
            if (cb.checked !== shouldBeChecked) {
              cb.checked = shouldBeChecked;
            }
          });
        });
      ", ns("updateCheckbox"))))
    ),

    shiny::div(
      class = "master-table-container",
      style = "width: 100%; overflow-x: auto;",
      shiny::tags$div(
        style = "margin-bottom: 10px; color: #666; font-size: 12px;",
        shiny::HTML("Drag column headers to reorder | Check 'Joint Plot' column to select pathways | Click 'Dashboard' for details")
      ),
      DT::dataTableOutput(ns("table"))
    ),
    shiny::tags$div(
      style = "display: none;",
      shiny::textInput(ns("joint_selection_store"), label = NULL, value = "")
    )
  )
}


#' @title Master Workspace Table Server
#' @description Provides data display, CSV annotation merging, and column display control with decoupled checkbox state.
#' @param id Module ID
#' @param data_prep Reactive data from the data preprocessing module
#' @return List containing selected_pathways and show_modal
#' @keywords internal

mod_master_table_server <- function(id, data_prep) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    joint_selected <- shiny::reactiveVal(character(0))
    has_interaction <- shiny::reactiveVal(FALSE)

    to_safe <- function(x) gsub("'", "\\\\'", x, fixed = TRUE)
    to_original <- function(x) gsub("\\\\'", "'", x, fixed = TRUE)

    # 动态加载CSV注释文件
    pathway_annotations <- shiny::reactive({
      possible_paths <- c(
        file.path("inst", "extdata", "pathway_annotations.csv"),
        file.path("extdata", "pathway_annotations.csv"),
        system.file("extdata", "pathway_annotations.csv", package = "GSEAlens")
      )

      csv_path <- NULL
      for (path in possible_paths) {
        if (file.exists(path)) {
          csv_path <- path
          break
        }
      }

      if (is.null(csv_path)) {
        message("pathway_annotations.csv not found, skipping annotation loading")
        return(NULL)
      }

      tryCatch({
        # 使用read.csv，支持UTF-8和中文
        anno_df <- read.csv(csv_path, stringsAsFactors = FALSE,
                            check.names = FALSE,
                            encoding = "UTF-8")

        if (nrow(anno_df) == 0) {
          message("CSV file is empty")
          return(NULL)
        }

        if (!"ID" %in% colnames(anno_df)) {
          warning("CSV file missing ID column, cannot merge")
          return(NULL)
        }

        message(sprintf("Successfully loaded annotation file: %d rows x %d columns", nrow(anno_df), ncol(anno_df)))
        return(anno_df)
      }, error = function(e) {
        warning(sprintf("Failed to read CSV annotation file: %s", e$message))
        return(NULL)
      })
    })

    shiny::observeEvent(input$joint_plot_toggle, {
      has_interaction(TRUE)
      toggle <- input$joint_plot_toggle
      if (is.null(toggle) || !is.list(toggle)) return()

      current <- joint_selected()
      if (isTRUE(toggle$checked)) {
        if (!(toggle$id %in% current)) joint_selected(c(current, toggle$id))
      } else {
        joint_selected(setdiff(current, toggle$id))
      }
    })


    # p值格式化函数
    format_pvalue <- function(x, threshold = 0.001) {
      ifelse(x < threshold,
             format(x, digits = 2, scientific = TRUE),
             format(round(x, 3), nsmall = 3))
    }


    # 表格渲染 - 使用 isolate() 阻止 checkbox 状态触发刷新
    output$table <- DT::renderDataTable({
      data_list <- data_prep()
      shiny::validate(shiny::need(data_list, "Waiting for data to load..."))
      shiny::validate(shiny::need(nrow(data_list$df) > 0, "No data under current filter criteria"))

      df <- data_list$df
      left_grp <- data_list$left_group
      right_grp <- data_list$right_group

      # 计算 Enriched_In 列
      df$Enriched_In <- ifelse(df$NES > 0, left_grp, right_grp)
      df$Enriched_In <- factor(df$Enriched_In, levels = c(left_grp, right_grp))

      # 合并 Description
      if (!is.null(data_list$task_obj$meta$meta_dict)) {
        meta_dict <- data_list$task_obj$meta$meta_dict
        if ("Brief_Description" %in% colnames(meta_dict)) {
          desc_map <- setNames(meta_dict$Brief_Description, meta_dict$ID)
          df$Description <- desc_map[df$ID]
        } else if ("Description" %in% colnames(meta_dict)) {
          desc_map <- setNames(meta_dict$Description, meta_dict$ID)
          df$Description <- desc_map[df$ID]
        } else {
          df$Description <- df$ID
        }
      } else {
        df$Description <- df$ID
      }
      df$Description[is.na(df$Description)] <- df$ID[is.na(df$Description)]

      df$NES_display <- round(df$NES, 2)
      df$pvalue_display <- format_pvalue(df$pvalue, threshold = 0.001)
      df$padj_display <- format_pvalue(df$p.adjust, threshold = 0.001)

      # 这样勾选复选框时，表格不会重新渲染，保持滚动位置和搜索状态
      current_selection <- isolate(joint_selected())

      # 交互列
      df$Select_for_Plot <- sprintf(
        '<input type="checkbox" class="joint-plot-checkbox" data-id="%s" %s onclick="Shiny.setInputValue(&#39;%s&#39;, {id: &#39;%s&#39;, checked: this.checked}, {priority: &#39;event&#39;});"/>',
        df$Safe_ID,
        ifelse(df$Safe_ID %in% current_selection, 'checked="checked"', ''),
        ns("joint_plot_toggle"),
        df$Safe_ID
      )

      df$Detail_Btn <- sprintf(
        '<button class="btn btn-sm btn-success" onclick="Shiny.setInputValue(&#39;%s&#39;, &#39;%s&#39;, {priority: &#39;event&#39;})">Dashboard</button>',
        ns("show_modal"),
        df$Safe_ID
      )

      # 合并CSV注释数据
      anno_df <- pathway_annotations()
      if (!is.null(anno_df)) {
        cols_to_keep <- sapply(anno_df, function(col) {
          !all(is.na(col) | col == "" | col == "NA")
        })
        cols_to_keep["ID"] <- TRUE

        if (sum(cols_to_keep) > 1) {
          anno_df_filtered <- anno_df[, cols_to_keep, drop = FALSE]
          df <- dplyr::left_join(df, anno_df_filtered, by = "ID")
          message(sprintf("Merged annotation columns: %s",
                          paste(setdiff(colnames(anno_df_filtered), "ID"), collapse = ", ")))
        }
      }

      # 构建显示列
      base_cols <- c("Rank", "Select_for_Plot", "Detail_Btn", "ID", "Enriched_In",
                     "NES_display", "pvalue_display", "padj_display", "setSize", "Description")

      if (!is.null(anno_df)) {
        csv_cols <- setdiff(colnames(anno_df), "ID")
        csv_cols_non_empty <- sapply(csv_cols, function(col) {
          if (!col %in% colnames(df)) return(FALSE)
          vals <- df[[col]]
          !all(is.na(vals) | vals == "" | vals == "NA")
        })
        csv_cols_to_show <- csv_cols[csv_cols_non_empty]
        base_cols <- c(base_cols, csv_cols_to_show)
      }

      # 预留BLANK1-10列
      blank_cols <- paste0("BLANK", 1:10)
      existing_blank <- blank_cols[blank_cols %in% colnames(df)]
      if (length(existing_blank) > 0) {
        blank_non_empty <- sapply(existing_blank, function(col) {
          !all(is.na(df[[col]]) | df[[col]] == "")
        })
        base_cols <- c(base_cols, existing_blank[blank_non_empty])
      }

      display_cols <- intersect(base_cols, colnames(df))
      dt_data <- df[, display_cols]

      col_name_map <- c(
        "Rank" = "Rank",
        "Select_for_Plot" = "Joint Plot",
        "Detail_Btn" = "Dashboard",
        "ID" = "Pathway",
        "Enriched_In" = "Enriched In",
        "NES_display" = "NES",
        "pvalue_display" = "P-value",
        "padj_display" = "FDR",
        "setSize" = "Size",
        "Description" = "Description"
      )

      # 为CSV列和BLANK列保留原始名称
      for (col in display_cols) {
        if (!col %in% names(col_name_map)) {
          col_name_map[col] <- col
        }
      }

      names(dt_data) <- col_name_map[display_cols]


      col_defs <- list()

      if ("Select_for_Plot" %in% display_cols) {
        idx <- which(display_cols == "Select_for_Plot") - 1
        col_defs[[length(col_defs) + 1]] <- list(
          width = '40px',
          targets = idx,
          className = 'dt-center',
          orderable = FALSE
        )
      }

      if ("Detail_Btn" %in% display_cols) {
        idx <- which(display_cols == "Detail_Btn") - 1
        col_defs[[length(col_defs) + 1]] <- list(
          width = '60px',
          targets = idx,
          className = 'dt-center',
          orderable = FALSE
        )
      }

      if ("Rank" %in% display_cols) {
        idx <- which(display_cols == "Rank") - 1
        col_defs[[length(col_defs) + 1]] <- list(
          width = '20px',
          targets = idx,
          className = 'dt-center'
        )
      }

      if ("ID" %in% display_cols) {
        idx <- which(display_cols == "ID") - 1
        col_defs[[length(col_defs) + 1]] <- list(
          width = '200px',
          targets = idx
        )
      }

      if ("Enriched_In" %in% display_cols) {
        idx <- which(display_cols == "Enriched_In") - 1
        col_defs[[length(col_defs) + 1]] <- list(
          width = '80px',
          targets = idx,
          className = 'dt-center'
        )
      }

      # 数值列设置较窄宽度
      numeric_cols <- c("NES_display", "pvalue_display", "padj_display", "setSize")
      for (numeric_col in numeric_cols) {
        if (numeric_col %in% display_cols) {
          idx <- which(display_cols == numeric_col) - 1
          col_defs[[length(col_defs) + 1]] <- list(
            width = '60px',
            targets = idx,
            className = 'dt-center'
          )
        }
      }

      if ("Description" %in% display_cols) {
        idx <- which(display_cols == "Description") - 1
        col_defs[[length(col_defs) + 1]] <- list(
          width = '800px',
          targets = idx,
          className = 'description-cell'
        )
      }

      csv_and_blank_cols <- setdiff(display_cols, c("Rank", "Select_for_Plot", "Detail_Btn", "ID",
                                                    "Enriched_In", "NES_display", "pvalue_display",
                                                    "padj_display", "setSize", "Description"))

      for (col in csv_and_blank_cols) {
        idx <- which(display_cols == col) - 1
        col_defs[[length(col_defs) + 1]] <- list(
          width = '800px',
          targets = idx,
          className = 'dt-csv-column'
        )
      }


      # 关键修复：优化 DT::datatable 配置
      dt <- DT::datatable(
        dt_data,
        escape = FALSE,
        selection = "multiple",
        rownames = FALSE,
        extensions = c('Buttons', 'Scroller', 'ColReorder', 'FixedColumns'),
        options = list(
          scrollX = TRUE,
          scrollY = "60vh",
          scroller = TRUE,
          deferRender = TRUE,
          pageLength = -1,
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel', 'colvis'),
          colReorder = TRUE,
          columnDefs = col_defs,
          order = list(list(0, 'asc')),
          autoWidth = TRUE,
          fixedColumns = list(
            leftColumns = 2
          ),
          scrollCollapse = TRUE,
          language = list(
            emptyTable = "No data available",
            zeroRecords = "No matching records found"
          )
        )
      )


      if ("Enriched In" %in% names(dt_data)) {
        dt <- dt %>% DT::formatStyle(
          columns = "Enriched In",
          backgroundColor = DT::styleEqual(c(left_grp, right_grp), c('#fee0d2', '#deebf7')),
          color = DT::styleEqual(c(left_grp, right_grp), c('#cc212f', '#1052bd')),
          fontWeight = 'bold'
        )
      }

      # 2. NES 列方向着色（使用十六进制色值：负值蓝色，正值红色）
      if ("NES" %in% names(dt_data)) {
        dt <- dt %>% DT::formatStyle(
          columns = "NES",
          color = DT::styleInterval(0, c('#1052bd', '#cc212f')),
          fontWeight = 'bold'
        )
      }

      # 3. FDR 列分段着色（使用styleInterval：根据阈值设置不同颜色）
      if ("FDR" %in% names(dt_data)) {
        dt <- dt %>% DT::formatStyle(
          columns = "FDR",
          backgroundColor = DT::styleInterval(c(0.01,0.05, 0.25), c('#fc9272',"#fdb9a2", '#fee0d2', 'transparent')),
          fontWeight = DT::styleInterval(0.05, c('bold', 'normal'))
        )
      }

      if ("P-value" %in% names(dt_data)) {
        dt <- dt %>% DT::formatStyle(
          columns = "P-value",
          backgroundColor = DT::styleInterval(c(0.01, 0.05), c('#fc9272', '#fee0d2', 'transparent')),
          fontWeight = DT::styleInterval(0.05, c('bold', 'normal'))
        )
      }


      dt
    }, server = TRUE)

    # 监听 joint_selected 变化，通过 JavaScript 同步复选框状态
    shiny::observe({
      sel <- joint_selected()
      # 发送自定义消息到前端，更新复选框状态
      session$sendCustomMessage(type = ns("updateCheckbox"), message = list(
        ids = sel,
        ns = ns("")
      ))
    })

    selected_pathways <- shiny::reactive({
      if (isTRUE(has_interaction())) return(to_original(joint_selected()))
      data_list <- data_prep()
      if (is.null(data_list) || is.null(input$table_rows_selected)) return(character(0))
      return(data_list$df$ID[input$table_rows_selected])
    })

    show_modal_trigger <- shiny::reactive({ input$show_modal })

    remove_pathways <- function(ids) {
      ids <- ids[!is.na(ids)]
      if (length(ids) == 0) return()
      has_interaction(TRUE)
      safe_ids <- to_safe(ids)
      new_sel <- setdiff(joint_selected(), safe_ids)
      joint_selected(new_sel)
    }

    clear_selection <- function() {
      has_interaction(TRUE)
      joint_selected(character(0))
    }

    return(list(
      selected_pathways = selected_pathways,
      show_modal = show_modal_trigger,
      remove_pathways = remove_pathways,
      clear_selection = clear_selection
    ))
  })
}
