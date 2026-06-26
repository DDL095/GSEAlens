# Section: Shiny Module - Master Table ----

#' @title Master Workspace Table UI
#' @description Interactive data table module for displaying GSEA pathway results
#'   with checkbox selection and modal integration.
#' @param id Module ID
#' @return Shiny UI tagList
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
      ", ns("updateCheckbox")))),
      # External selection update message handler
      shiny::tags$script(shiny::HTML(sprintf("
        Shiny.addCustomMessageHandler('%s', function(message) {
          var action = message.action;
          var ids = message.ids || [];
          var checkboxes = document.querySelectorAll('.joint-plot-checkbox');

          if (action === 'set') {
            checkboxes.forEach(function(cb) {
              var id = cb.getAttribute('data-id');
              cb.checked = ids.includes(id);
            });
          } else if (action === 'add') {
            checkboxes.forEach(function(cb) {
              var id = cb.getAttribute('data-id');
              if (ids.includes(id)) {
                cb.checked = true;
              }
            });
          } else if (action === 'remove') {
            checkboxes.forEach(function(cb) {
              var id = cb.getAttribute('data-id');
              if (ids.includes(id)) {
                cb.checked = false;
              }
            });
          } else if (action === 'clear') {
            checkboxes.forEach(function(cb) {
              cb.checked = false;
            });
          }
        });
      ", ns("externalUpdate"))))
    ),
    shiny::div(
      id = ns("table-wrapper"),
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
#' @description Provides data display, column display control with decoupled checkbox state.
#'   Supports merging addition_data columns into the main table.
#'   Includes external selection update capabilities for bidirectional sync.
#' @param id Module ID
#' @param data_prep Reactive data from the data preprocessing module
#' @param addition_data Optional. A data frame with pathway annotations to merge.
#' @return List containing selected_pathways, show_modal, and selection control functions
#' @keywords internal
#' @importFrom dplyr left_join
mod_master_table_server <- function(id, data_prep, addition_data = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    joint_selected <- shiny::reactiveVal(character(0))
    has_interaction <- shiny::reactiveVal(FALSE)


    to_safe <- function(x) gsub("'", "\\'", x, fixed = TRUE)
    to_original <- function(x) gsub("\\'", "'", x, fixed = TRUE)

    # Validate and prepare addition_data
    .validate_addition_data_internal <- function(add_data) {
      if (is.null(add_data)) {
        return(NULL)
      }
      if (!is.data.frame(add_data)) {
        warning("[MasterTable] addition_data must be a data.frame, ignoring")
        return(NULL)
      }
      if (!"ID" %in% colnames(add_data)) {
        warning("[MasterTable] addition_data must contain 'ID' column, ignoring")
        return(NULL)
      }
      add_data$ID <- as.character(add_data$ID)
      return(add_data)
    }

    addition_data_validated <- .validate_addition_data_internal(addition_data)


    # =========================================
    # Key fix: Joint plot toggle event - use identical() for boolean comparison
    # =========================================
    shiny::observeEvent(input$joint_plot_toggle, {
      toggle <- input$joint_plot_toggle
      if (is.null(toggle) || !is.list(toggle)) {
        return()
      }

      # Key fix: JavaScript passes strings "true"/"false", need to convert
      is_checked <- identical(toggle$checked, TRUE) || identical(toggle$checked, "true")

      has_interaction(TRUE) # Mark as interacted when event fires

      current <- joint_selected()

      if (is_checked) {
        # Add to selected list
        if (!(toggle$id %in% current)) {
          joint_selected(c(current, toggle$id))
          message(sprintf("[MasterTable] Added pathway: %s", toggle$id))
        }
      } else {
        # Remove from selected list
        new_sel <- setdiff(current, toggle$id)
        joint_selected(new_sel)
        message(sprintf("[MasterTable] Removed pathway: %s", toggle$id))
      }
    })

    # External selection update interface (for other modules to call)
    update_selection <- function(action = c("set", "add", "remove", "clear"), ids = character(0)) {
      action <- match.arg(action)
      current <- joint_selected()

      switch(action,
        "set" = {
          has_interaction(TRUE)
          joint_selected(unique(ids))
          session$sendCustomMessage(type = ns("externalUpdate"), message = list(
            action = "set",
            ids = to_safe(unique(ids))
          ))
          message("[MasterTable] Selection set to: ", paste(ids, collapse = ", "))
        },
        "add" = {
          has_interaction(TRUE)
          new_sel <- union(current, ids)
          joint_selected(unique(new_sel))
          session$sendCustomMessage(type = ns("externalUpdate"), message = list(
            action = "add",
            ids = to_safe(unique(ids))
          ))
          message("[MasterTable] Added to selection: ", paste(ids, collapse = ", "))
        },
        "remove" = {
          has_interaction(TRUE)
          new_sel <- setdiff(current, ids)
          joint_selected(new_sel)
          session$sendCustomMessage(type = ns("externalUpdate"), message = list(
            action = "remove",
            ids = to_safe(unique(ids))
          ))
          message("[MasterTable] Removed from selection: ", paste(ids, collapse = ", "))
        },
        "clear" = {
          has_interaction(TRUE)
          joint_selected(character(0))
          session$sendCustomMessage(type = ns("externalUpdate"), message = list(
            action = "clear",
            ids = character(0)
          ))
          message("[MasterTable] Selection cleared")
        }
      )
    }

    # =========================================
    # Table rendering
    # =========================================
    output$table <- DT::renderDataTable(
      {
        data_list <- data_prep()
        shiny::validate(shiny::need(data_list, "Waiting for data to load..."))
        shiny::validate(shiny::need(nrow(data_list$df) > 0, "No data under current filter criteria"))

        df <- data_list$df
        left_grp <- data_list$left_group
        right_grp <- data_list$right_group

        # Calculate Enriched_In
        df$Enriched_In <- ifelse(df$NES > 0, left_grp, right_grp)
        df$Enriched_In <- factor(df$Enriched_In, levels = c(left_grp, right_grp))

        # Merge Description
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

        # Merge addition_data columns
        if (!is.null(addition_data_validated)) {
          add_cols <- setdiff(colnames(addition_data_validated), "ID")
          df <- dplyr::left_join(df, addition_data_validated, by = "ID", suffix = c("", "_add"))
        }

        df$NES_display <- round(df$NES, 2)

        # Get current selection state
        current_selection <- isolate(joint_selected())

        # Interactive columns
        df$Select_for_Plot <- sprintf(
          '<input type="checkbox" class="joint-plot-checkbox" data-id="%s" %s onclick="Shiny.setInputValue(\'%s\', {id: \'%s\', checked: this.checked}, {priority: \'event\'});"/>',
          df$Safe_ID,
          ifelse(df$Safe_ID %in% current_selection, 'checked="checked"', ""),
          ns("joint_plot_toggle"),
          df$Safe_ID
        )

        df$Detail_Btn <- sprintf(
          '<button class="btn btn-sm btn-success" onclick="Shiny.setInputValue(\'%s\', \'%s\', {priority: \'event\'})">Dashboard</button>',
          ns("show_modal"),
          df$Safe_ID
        )

        # Display columns
        base_cols <- c(
          "Rank", "Select_for_Plot", "Detail_Btn", "ID", "Enriched_In",
          "NES_display", "pvalue", "p.adjust", "setSize", "Description"
        )

        if (!is.null(addition_data_validated)) {
          add_cols_display <- setdiff(colnames(addition_data_validated), "ID")
          add_cols_valid <- add_cols_display[vapply(add_cols_display, function(col) {
            !all(is.na(df[[col]]) | df[[col]] == "" | is.null(df[[col]]))
          }, logical(1))]
          base_cols <- c(base_cols, add_cols_valid)
        }

        display_cols <- intersect(base_cols, colnames(df))
        dt_data <- df[, display_cols, drop = FALSE]

        col_name_map <- c(
          "Rank" = "Rank",
          "Select_for_Plot" = "Joint Plot",
          "Detail_Btn" = "Dashboard",
          "ID" = "Pathway",
          "Enriched_In" = "Enriched In",
          "NES_display" = "NES",
          "pvalue" = "P-value",
          "p.adjust" = "FDR",
          "setSize" = "Size",
          "Description" = "Description"
        )

        for (col in display_cols) {
          if (!col %in% names(col_name_map)) {
            col_name_map[col] <- col
          }
        }

        names(dt_data) <- col_name_map[display_cols]

        # Column definitions
        col_defs <- list()

        if ("Select_for_Plot" %in% display_cols) {
          idx <- which(display_cols == "Select_for_Plot") - 1
          col_defs[[length(col_defs) + 1]] <- list(
            width = "40px", targets = idx, className = "dt-center", orderable = FALSE
          )
        }

        if ("Detail_Btn" %in% display_cols) {
          idx <- which(display_cols == "Detail_Btn") - 1
          col_defs[[length(col_defs) + 1]] <- list(
            width = "60px", targets = idx, className = "dt-center", orderable = FALSE
          )
        }

        if ("Rank" %in% display_cols) {
          idx <- which(display_cols == "Rank") - 1
          col_defs[[length(col_defs) + 1]] <- list(
            width = "20px", targets = idx, className = "dt-center"
          )
        }

        if ("ID" %in% display_cols) {
          idx <- which(display_cols == "ID") - 1
          col_defs[[length(col_defs) + 1]] <- list(
            width = "200px", targets = idx
          )
        }

        if ("Enriched_In" %in% display_cols) {
          idx <- which(display_cols == "Enriched_In") - 1
          col_defs[[length(col_defs) + 1]] <- list(
            width = "80px", targets = idx, className = "dt-center"
          )
        }

        numeric_cols <- c("NES_display", "pvalue", "p.adjust", "setSize")
        for (numeric_col in numeric_cols) {
          if (numeric_col %in% display_cols) {
            idx <- which(display_cols == numeric_col) - 1
            col_defs[[length(col_defs) + 1]] <- list(
              width = "60px", targets = idx, className = "dt-center"
            )
          }
        }

        if ("Description" %in% display_cols) {
          idx <- which(display_cols == "Description") - 1
          col_defs[[length(col_defs) + 1]] <- list(
            width = "800px", targets = idx, className = "description-cell"
          )
        }

        # JS render function: scientific notation for small values
        js_pval_render <- "
        function(data, type) {
          if (type === 'sort' || type === 'type') return data;
          if (data === null || isNaN(data)) return '-';
          if (Math.abs(data) < 0.001) return data.toExponential(2);
          return data.toFixed(4);
        }
      "

        # Find P-value and FDR column indices
        pval_idx <- which(display_cols == "pvalue") - 1
        padj_idx <- which(display_cols == "p.adjust") - 1

        # Add P-value and FDR render definitions
        col_defs[[length(col_defs) + 1]] <- list(
          targets = pval_idx,
          render = htmlwidgets::JS(js_pval_render)
        )
        col_defs[[length(col_defs) + 1]] <- list(
          targets = padj_idx,
          render = htmlwidgets::JS(js_pval_render)
        )


        # Render table
        dt <- DT::datatable(
          dt_data,
          escape = FALSE,
          selection = "none",
          rownames = FALSE,
          extensions = c("Buttons", "Scroller", "ColReorder", "FixedColumns"),
          options = list(
            scrollX = TRUE,
            scrollY = "70vh",
            scroller = TRUE,
            deferRender = TRUE,
            pageLength = -1,
            dom = "Bfrtip",
            buttons = c("copy", "csv", "excel", "colvis"),
            colReorder = TRUE,
            columnDefs = col_defs,
            order = list(list(0, "asc")),
            autoWidth = FALSE,
            fixedColumns = list(leftColumns = 2),
            scrollCollapse = TRUE,
            language = list(emptyTable = "No data available", zeroRecords = "No matching records found")
          )
        )

        if ("Enriched In" %in% names(dt_data)) {
          dt <- dt |> DT::formatStyle(
            columns = "Enriched In",
            backgroundColor = DT::styleEqual(c(left_grp, right_grp), c("#fee0d2", "#deebf7")),
            color = DT::styleEqual(c(left_grp, right_grp), c("#cc212f", "#1052bd")),
            fontWeight = "bold"
          )
        }

        if ("NES" %in% names(dt_data)) {
          dt <- dt |> DT::formatStyle(
            columns = "NES",
            color = DT::styleInterval(0, c("#1052bd", "#cc212f")),
            fontWeight = "bold"
          )
        }

        if ("FDR" %in% names(dt_data)) {
          dt <- dt |> DT::formatStyle(
            columns = "FDR",
            backgroundColor = DT::styleInterval(c(0.01, 0.05, 0.25), c("#fc9272", "#fdb9a2", "#fee0d2", "transparent")),
            fontWeight = DT::styleInterval(0.05, c("bold", "normal"))
          )
        }

        if ("P-value" %in% names(dt_data)) {
          dt <- dt |> DT::formatStyle(
            columns = "P-value",
            backgroundColor = DT::styleInterval(c(0.01, 0.05), c("#fc9272", "#fee0d2", "transparent")),
            fontWeight = DT::styleInterval(0.05, c("bold", "normal"))
          )
        }

        dt
      },
      server = TRUE
    )

    # Sync checkbox state
    shiny::observe({
      sel <- joint_selected()
      session$sendCustomMessage(type = ns("updateCheckbox"), message = list(
        ids = sel, ns = ns("")
      ))
    })

    # =========================================
    # Key fix: selected_pathways reactive - ensure checkbox state is returned
    # =========================================
    selected_pathways <- shiny::reactive({
      # Return checkbox interaction state first
      if (isTRUE(has_interaction())) {
        return(to_original(joint_selected()))
      }
      # Fall back to row click
      data_list <- data_prep()
      if (is.null(data_list) || is.null(input$table_rows_selected)) {
        return(character(0))
      }
      return(data_list$df$ID[input$table_rows_selected])
    })

    show_modal_trigger <- shiny::reactive({
      input$show_modal
    })

    remove_pathways <- function(ids) {
      ids <- ids[!is.na(ids)]
      if (length(ids) == 0) {
        return()
      }
      has_interaction(TRUE)
      safe_ids <- to_safe(ids)
      new_sel <- setdiff(joint_selected(), safe_ids)
      joint_selected(new_sel)
    }

    clear_selection <- function() {
      has_interaction(TRUE)
      joint_selected(character(0))
    }

    # Return all interfaces
    return(list(
      selected_pathways = selected_pathways,
      show_modal = show_modal_trigger,
      remove_pathways = remove_pathways,
      clear_selection = clear_selection,
      update_selection = update_selection
    ))
  })
}
