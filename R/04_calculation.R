#' @title Batch Parallel GSEA Calculation (Optimized Version)
#' @description Consumes a standard GseaEnv object and performs efficient parallel GSEA calculation.
#' @param gsea_env GseaEnv object
#' @param custom_series_name String. Analysis series name
#' @param output_dir String. Output directory for results, default "./GSEA_Output"
#' @param workers Number of parallel cores. Default 4
#' @param bidirectional Logical. Whether to automatically generate reverse contrasts, default TRUE
#' @param minGSSize Minimum gene set size, default 10
#' @param maxGSSize Maximum gene set size, default 500
#' @param pvalueCutoff P-value threshold, default 1
#' @param force Logical. Whether to force recalculation, default FALSE
#' @param use_progress Logical. Whether to show progress bar, default TRUE
#' @param chunk_size Integer. Number of tasks per worker per chunk, default NULL (auto)
#' @return GseaRes object
#' @export
batch_calc_gsea <- function(gsea_env,
                            custom_series_name = "Auto_Analysis",
                            output_dir = "./GSEA_Output",
                            workers = 4,
                            bidirectional = TRUE,
                            minGSSize = 10,
                            maxGSSize = 500,
                            pvalueCutoff = 1,
                            force = FALSE,
                            use_progress = FALSE,
                            chunk_size = NULL) {

  .check_gsea_env(gsea_env)

  start_time <- Sys.time()
  start_ms <- as.numeric(start_time) * 1000

  series_dir <- file.path(output_dir, custom_series_name)
  if (!dir.exists(series_dir)) {
    dir.create(series_dir, recursive = TRUE, showWarnings = FALSE)
  }

  rds_name <- sprintf("GSEA_Capsule_[%s]_[%s].rds",
                      custom_series_name, gsea_env$geneset$name)
  rds_path <- file.path(series_dir, rds_name)

  if (file.exists(rds_path) && !force) {
    message(sprintf("Cache hit! Existing GSEA capsule detected: %s", rds_name))
    return(readRDS(rds_path))
  }

  registry <- gsea_env$contrast_registry
  de_store <- gsea_env$de_store

  task_metadata <- list()
  for (i in seq_len(nrow(registry))) {
    row <- registry[i, ]
    cid <- row$contrast_id
    task_metadata[[cid]] <- list(
      task_id = cid,
      left_group = row$left_group,
      right_group = row$right_group,
      is_reversed = FALSE
    )
    if (bidirectional) {
      rev_cid <- paste(row$right_group, row$left_group, sep = "_vs_")
      task_metadata[[rev_cid]] <- list(
        task_id = rev_cid,
        left_group = row$right_group,
        right_group = row$left_group,
        is_reversed = TRUE
      )
    }
  }

  total_tasks <- length(task_metadata)
  message(sprintf("Ready: %d contrast tasks pending calculation...", total_tasks))

  options(future.globals.maxSize = 192 * 1024^3)

  total_cores <- parallel::detectCores(logical = TRUE)
  use_cores <- min(total_cores, max(1, workers))
  message(sprintf("Hardware scheduling: Using %d cores for parallel computation...", use_cores))

  if (is.null(chunk_size)) {
    chunk_size <- max(1, ceiling(total_tasks / (use_cores * 4)))
  }

  future::plan(future::multisession, workers = use_cores)

  if (use_progress) {
    if (!requireNamespace("progressr", quietly = TRUE)) {
      stop("Please install progressr package: install.packages('progressr')")
    }
    progressr::handlers(global = TRUE)
    progressr::handlers("progress")
  }

  worker_term2gene <- gsea_env$geneset$term2gene
  worker_meta_dict <- gsea_env$geneset$meta_dict
  worker_de_list <- as.list(de_store)

  task_names <- names(task_metadata)
  task_chunks <- split(task_names, ceiling(seq_along(task_names) / chunk_size))

  message(sprintf("Chunk strategy: %d chunks, up to %d tasks per chunk",
                  length(task_chunks), chunk_size))

  p <- progressr::progressor(steps = total_tasks, enable = use_progress)

  res_list <- future.apply::future_lapply(
    X = task_chunks,
    FUN = function(chunk_task_names,
                   task_metadata,
                   de_list,
                   term2gene,
                   meta_dict,
                   minGSSize,
                   maxGSSize,
                   pvalueCutoff,
                   progressor_fn) {

      if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
          !requireNamespace("dplyr", quietly = TRUE)) {
        stop("Worker missing required packages")
      }

      chunk_results <- list()

      for (task_name in chunk_task_names) {
        task_info <- task_metadata[[task_name]]
        task_id <- task_info$task_id

        if (task_info$is_reversed) {
          original_cid <- paste(task_info$right_group, task_info$left_group, sep = "_vs_")
        } else {
          original_cid <- task_id
        }

        de_table <- de_list[[original_cid]]

        if (is.null(de_table) || nrow(de_table) == 0) {
          chunk_results[[task_name]] <- list(
            name = task_id,
            status = "Failed",
            data = NULL,
            genelist = c()
          )
          if (!is.null(progressor_fn)) progressor_fn()
          next
        }

        genelist <- tryCatch({
          .prepare_rank_vector_fast(de_table, flip = task_info$is_reversed)
        }, error = function(e) {
          c()
        })

        if (length(genelist) == 0) {
          chunk_results[[task_name]] <- list(
            name = task_id,
            status = "Failed",
            data = NULL,
            genelist = c()
          )
          if (!is.null(progressor_fn)) progressor_fn()
          next
        }

        detect_case_format <- function(genes) {
          sample <- head(unique(genes), 100)
          sample_filtered <- sample[!grepl("^[0-9]", sample)]
          if (length(sample_filtered) == 0) return("mixed")
          n_title <- sum(grepl("^[A-Z][a-z]", sample_filtered))
          n_upper <- sum(grepl("^[A-Z]{2,}$", sample_filtered))
          if (n_title > n_upper) return("title_case")
          if (n_upper > n_title) return("upper_case")
          return("mixed")
        }

        geneset_species <- gsea_env$geneset$species %||% "HS"
        term2gene_format <- detect_case_format(term2gene$gene_symbol)
        de_format <- detect_case_format(names(genelist))

        if (term2gene_format == "upper_case") {
          if (de_format == "title_case") {
            names(genelist) <- toupper(names(genelist))
          }
        } else if (term2gene_format == "title_case") {
          if (de_format == "upper_case") {
            stop(paste0(
              "[Species Mismatch Error]\n",
              "Gene set species: Mouse (MM) - requires title case gene symbols\n",
              "DE gene format: Uppercase (e.g., GAPDH, IRF7) - appears to be human data\n\n",
              "Please use human gene sets (species='HS') for uppercase DE data,\n",
              "or provide mouse DE data with title case gene symbols (e.g., Gapdh, Irf7)."
            ))
          }
        }

        set.seed(123)
        gsea_res <- tryCatch({
          clusterProfiler::GSEA(
            geneList = genelist,
            TERM2GENE = term2gene,
            minGSSize = minGSSize,
            maxGSSize = maxGSSize,
            pvalueCutoff = pvalueCutoff,
            pAdjustMethod = "BH",
            verbose = FALSE,
            seed = 123,
            eps = 0
          )
        }, error = function(e) {
          NULL
        })

        status <- if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) "Success" else "Failed/NoEnrich"

        if (status == "Success" && !is.null(meta_dict)) {
          gsea_res@result <- .enrich_gsea_result(gsea_res@result, meta_dict)
        }

        chunk_results[[task_name]] <- list(
          name = task_id,
          status = status,
          data = gsea_res,
          genelist = genelist
        )

        if (!is.null(progressor_fn)) progressor_fn()
      }

      return(chunk_results)
    },
    task_metadata = task_metadata,
    de_list = worker_de_list,
    term2gene = worker_term2gene,
    meta_dict = worker_meta_dict,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    pvalueCutoff = pvalueCutoff,
    progressor_fn = if (use_progress) p else NULL,
    future.seed = TRUE,
    future.scheduling = 1.0
  )

  future::plan(future::sequential)

  all_cons <- showConnections(all = TRUE)
  if (nrow(all_cons) > 1) {
    for (con_idx in as.integer(rownames(all_cons))) {
      if (con_idx > 2) {
        try(suppressWarnings(close.connection(getConnection(con_idx))), silent = TRUE)
      }
    }
  }

  invisible(gc(verbose = FALSE, full = TRUE))

  res_list_flat <- do.call(c, res_list)
  names(res_list_flat) <- names(task_metadata)

  end_time <- Sys.time()
  end_ms <- as.numeric(end_time) * 1000

  final_obj <- create_gsea_res(
    metadata = list(
      run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      workers_used = use_cores,
      chunk_size = chunk_size,
      parameters = list(
        bidirectional = bidirectional,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        pvalueCutoff = pvalueCutoff
      ),
      project_info = list(
        custom_series_name = custom_series_name,
        output_dir = normalizePath(output_dir, mustWork = FALSE),
        series_dir = normalizePath(series_dir, mustWork = FALSE),
        rds_path = normalizePath(rds_path, mustWork = FALSE)
      ),
      gsea_benchmark = list(
        start_time = format(start_time, "%Y-%m-%d %H:%M:%OS3"),
        start_ms = start_ms,
        end_time = format(end_time, "%Y-%m-%d %H:%M:%OS3"),
        end_ms = end_ms,
        duration_sec = round((end_ms - start_ms) / 1000, 3),
        workers = use_cores,
        total_tasks = length(task_metadata),
        successful_tasks = sum(sapply(res_list_flat, function(x) x$status == "Success"))
      )
    ),
    backend_info = gsea_env$backend_info,
    contrast_registry = gsea_env$contrast_registry,
    de_store = gsea_env$de_store,
    expr_bundle = gsea_env$expr_bundle,
    geneset_info = gsea_env$geneset,
    results = res_list_flat
  )

  saveRDS(final_obj, rds_path)

  success_count <- final_obj$metadata$gsea_benchmark$successful_tasks

  message(sprintf("\nCalculation complete! Time elapsed: %.2f seconds",
                  final_obj$metadata$gsea_benchmark$duration_sec))
  message(sprintf("   Successfully analyzed: %d/%d contrasts", success_count, total_tasks))
  message(sprintf("   Results saved to: %s", rds_path))

  return(final_obj)
}


#' @title Process Single GSEA Chunk
#' @param chunk_task_names Character vector
#' @param task_metadata Task metadata list
#' @param de_list DE data in list form
#' @param term2gene Gene set mapping dataframe
#' @param meta_dict Metadata dictionary
#' @param minGSSize Minimum gene set size
#' @param maxGSSize Maximum gene set size
#' @param pvalueCutoff P-value threshold
#' @return Result list for current chunk
#' @keywords internal
.process_gsea_chunk <- function(chunk_task_names,
                                task_metadata,
                                de_list,
                                term2gene,
                                meta_dict,
                                minGSSize,
                                maxGSSize,
                                pvalueCutoff) {

  if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE)) {
    stop("Worker missing required packages: clusterProfiler or dplyr")
  }

  chunk_results <- list()

  for (task_name in chunk_task_names) {
    task_info <- task_metadata[[task_name]]
    task_id <- task_info$task_id

    if (task_info$is_reversed) {
      original_cid <- paste(task_info$right_group, task_info$left_group, sep = "_vs_")
    } else {
      original_cid <- task_id
    }

    de_table <- de_list[[original_cid]]

    if (is.null(de_table) || nrow(de_table) == 0) {
      chunk_results[[task_name]] <- list(
        name = task_id,
        status = "Failed",
        data = NULL,
        genelist = c()
      )
      next
    }

    genelist <- tryCatch({
      .prepare_rank_vector_fast(de_table, flip = task_info$is_reversed)
    }, error = function(e) {
      warning(sprintf("Rank vector preparation failed: %s", e$message))
      c()
    })

    if (length(genelist) == 0) {
      chunk_results[[task_name]] <- list(
        name = task_id,
        status = "Failed",
        data = NULL,
        genelist = c()
      )
      next
    }

    set.seed(123)
    gsea_res <- tryCatch({
      clusterProfiler::GSEA(
        geneList = genelist,
        TERM2GENE = term2gene,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        pvalueCutoff = pvalueCutoff,
        pAdjustMethod = "BH",
        verbose = FALSE,
        seed = 123,
        eps = 0
      )
    }, error = function(e) {
      warning(sprintf("GSEA calculation failed [%s]: %s", task_id, e$message))
      NULL
    })

    status <- if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) {
      "Success"
    } else {
      "Failed/NoEnrich"
    }

    if (status == "Success" && !is.null(meta_dict)) {
      gsea_res@result <- .enrich_gsea_result(gsea_res@result, meta_dict)
    }

    chunk_results[[task_name]] <- list(
      name = task_id,
      status = status,
      data = gsea_res,
      genelist = genelist
    )
  }

  return(chunk_results)
}


#' @title Fast Rank Vector Preparation
#' @param de_table DE result table
#' @param flip Whether to flip sign
#' @return Named numeric vector
#' @keywords internal
.prepare_rank_vector_fast <- function(de_table, flip = FALSE) {
  # 保留原始大小写，格式标准化在 worker 函数中基于 TERM2GENE 进行
  vals <- de_table %>%
    dplyr::filter(!is.na(gene_symbol), gene_symbol != "") %>%
    dplyr::mutate(abs_stat = abs(stat)) %>%
    dplyr::arrange(dplyr::desc(abs_stat)) %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE) %>%
    {
      vec <- .$stat
      names(vec) <- .$gene_symbol
      vec
    }

  if (flip) vals <- -vals
  sort(vals, decreasing = TRUE)
}


#' @title GSEA Result Metadata Injection
#' @param result_df GSEA result dataframe
#' @param meta_dict Metadata dictionary
#' @return Enriched dataframe
#' @keywords internal
.enrich_gsea_result <- function(result_df, meta_dict) {

  if (!is.data.frame(result_df) || nrow(result_df) == 0) {
    warning("result_df is invalid or empty")
    return(result_df)
  }

  if (!is.data.frame(meta_dict) || nrow(meta_dict) == 0) {
    warning("meta_dict is empty, cannot enrich results")
    return(result_df)
  }

  if (!"ID" %in% colnames(result_df) || !"ID" %in% colnames(meta_dict)) {
    stop("Both result_df and meta_dict must contain 'ID' column")
  }

  original_rownames <- rownames(result_df)
  original_ids <- result_df$ID

  required_cols <- c("ID", "Description", "URL", "Collection", "Subcollection", "Combo_Name")

  for (col in required_cols) {
    if (!col %in% colnames(meta_dict)) {
      warning(sprintf("meta_dict missing column '%s', creating placeholder", col))
      meta_dict[[col]] <- NA_character_
    }
  }

  if (all(is.na(meta_dict$Subcollection)) || is.null(meta_dict$Subcollection)) {
    meta_dict$Subcollection <- ""
  }

  if (all(is.na(meta_dict$Combo_Name))) {
    meta_dict$Combo_Name <- ifelse(
      is.na(meta_dict$Subcollection) | meta_dict$Subcollection == "",
      meta_dict$Collection,
      paste0(meta_dict$Collection, ":", meta_dict$Subcollection)
    )
  }

  core_stat_cols <- c("ID", "setSize", "enrichmentScore", "NES", "pvalue",
                      "p.adjust", "qvalue", "rank", "leading_edge", "core_enrichment")

  conflict_cols <- c("Description", "URL", "Collection", "Subcollection", "Combo_Name")

  cols_to_keep <- setdiff(colnames(result_df), conflict_cols)
  result_work <- result_df[, cols_to_keep, drop = FALSE]

  meta_subset <- meta_dict[, required_cols, drop = FALSE]

  merged_df <- result_work %>%
    dplyr::left_join(
      meta_subset,
      by = "ID",
      suffix = c("", "_meta")
    )

  if (nrow(merged_df) != nrow(result_work)) {
    warning(sprintf("Row count changed during merge: %d -> %d",
                    nrow(result_work), nrow(merged_df)))
  }

  merged_df$Display_Collection <- dplyr::coalesce(
    merged_df$Combo_Name,
    merged_df$Collection,
    "Unknown"
  )
  merged_df$Display_Collection <- as.factor(merged_df$Display_Collection)

  merged_df$Pathway_Link <- ifelse(
    is.na(merged_df$URL) | merged_df$URL == "",
    sprintf("<b>%s</b>", merged_df$ID),
    sprintf('<a href="%s" target="_blank">%s</a>', merged_df$URL, merged_df$ID)
  )

  if ("Description_meta" %in% colnames(merged_df)) {
    merged_df$Description <- dplyr::coalesce(
      merged_df$Description_meta,
      merged_df$ID
    )
    merged_df$Description_meta <- NULL
  } else if (!"Description" %in% colnames(merged_df)) {
    merged_df$Description <- merged_df$ID
  }

  merged_df$Description[is.na(merged_df$Description)] <- merged_df$ID[is.na(merged_df$Description)]

  meta_temp_cols <- grep("_meta$", colnames(merged_df), value = TRUE)
  if (length(meta_temp_cols) > 0) {
    merged_df <- merged_df[, !(colnames(merged_df) %in% meta_temp_cols), drop = FALSE]
  }

  merged_df <- as.data.frame(merged_df)
  if ("ID" %in% colnames(merged_df)) {
    if (any(duplicated(merged_df$ID))) {
      warning("Duplicate IDs found in result, using original rownames")
      rownames(merged_df) <- original_rownames
    } else {
      rownames(merged_df) <- merged_df$ID
    }
  } else {
    rownames(merged_df) <- original_rownames
  }

  standard_cols <- c("ID", "setSize", "enrichmentScore", "NES", "pvalue",
                     "p.adjust", "qvalue", "rank", "leading_edge", "core_enrichment",
                     "Description", "URL", "Collection", "Subcollection", "Combo_Name",
                     "Display_Collection", "Pathway_Link")

  missing_standard <- setdiff(standard_cols, colnames(merged_df))
  if (length(missing_standard) > 0) {
    warning(sprintf("Final result missing standard columns: %s",
                    paste(missing_standard, collapse = ", ")))
    for (col in missing_standard) {
      merged_df[[col]] <- NA_character_
    }
  }

  present_cols <- intersect(standard_cols, colnames(merged_df))
  extra_cols <- setdiff(colnames(merged_df), standard_cols)
  merged_df <- merged_df[, c(present_cols, extra_cols), drop = FALSE]

  message(sprintf("[.enrich_gsea_result] Successfully enriched: %d rows x %d columns",
                  nrow(merged_df), ncol(merged_df)))

  return(merged_df)
}
