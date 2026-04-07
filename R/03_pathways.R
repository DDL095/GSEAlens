#' @title Build GSEA Pathway Database with Interactive Selection
#' @description Provides an interactive or automated gene set selection interface
#'   based on msigdbr, generating standardized pathway objects with intelligent
#'   semantic tagging. Supports both human (HS) and mouse (MM) species.
#' @param species Character string specifying the species. Options:
#'   \itemize{
#'     \item "HS": Homo sapiens (human, default)
#'     \item "MM": Mus musculus (mouse)
#'   }
#' @param auto_select Optional parameter for automated selection. Can be:
#'   \itemize{
#'     \item NULL: Launch interactive menu (default)
#'     \item "ALL": Load entire MSigDB database for the specified species
#'     \item Character vector: Collection names (e.g., c("H", "C5:GO:BP") for human)
#'   }
#' @return A list object containing:
#'   \itemize{
#'     \item TERM2GENE: Data frame mapping pathway IDs to gene symbols
#'     \item meta_dict: Metadata dictionary with descriptions, URLs, and collection info
#'     \item SuperTag: Intelligent batch identifier for the selected collections
#'     \item collections_used: Data frame of selected collection metadata
#'     \item species: Species identifier ("HS" or "MM")
#'   }
#' @export
#' @examples
#' \dontrun{
#' # Human mode (default)
#' pathways <- build_gsea_pathways(species = "HS")
#'
#' # Mouse mode
#' pathways <- build_gsea_pathways(species = "MM")
#'
#' # Automated selection by name
#' pathways <- build_gsea_pathways(species = "HS", auto_select = c("H", "C5:GO:BP"))
#' pathways_mm <- build_gsea_pathways(species = "MM", auto_select = c("H", "C5:GO:BP"))
#'
#' # Load all collections
#' pathways <- build_gsea_pathways(species = "HS", auto_select = "ALL")
#' }
#' @importFrom dplyr arrange left_join select distinct mutate all_of bind_rows filter
#' @importFrom msigdbr msigdbr_collections msigdbr
build_gsea_pathways <- function(species = "HS", auto_select = NULL) {

  # === 物种参数标准化 ===
  species <- toupper(trimws(as.character(species)))

  # 物种配置表
  species_config <- list(
    "HS" = list(
      db_species = "HS",
      msigdbr_species = "Homo sapiens",
      display_name = "Homo sapiens (Human)",
      is_mouse = FALSE
    ),
    "MM" = list(
      db_species = "MM",
      msigdbr_species = "mouse",
      display_name = "Mus musculus (Mouse)",
      is_mouse = TRUE
    )
  )

  if (!species %in% names(species_config)) {
    stop(sprintf(
      "Invalid species '%s'. Please use 'HS' (human) or 'MM' (mouse).",
      species
    ))
  }

  cfg <- species_config[[species]]
  message(sprintf("[build_gsea_pathways] Species: %s", cfg$display_name))

  # === MM 数据库代码 → msigdbr 函数参数映射 ===
  # msigdbr_collections(db_species = "MM") 返回 M1, M2, MH 等 MM 代码
  # 但 msigdbr(species = "mouse", ...) 只能接受人类代码 (C1, C2, H 等)
  mm_code_map <- list(
    # MM代码 -> list(msigdbr_collection, msigdbr_subcollection, 显示名称)
    "MH" = list(coll = "H", subcoll = "",
                desc = "MH Hallmark Gene Sets", tag = "MHall"),
    "M1" = list(coll = "C1", subcoll = "",
                desc = "M1 Positional Gene Sets", tag = "M1Pos"),
    "M2:CGP" = list(coll = "C2", subcoll = "CGP",
                    desc = "M2 Chemical and Genetic Perturbations", tag = "MCGP"),
    "M2:CP:BIOCARTA" = list(coll = "C2", subcoll = "CP:BIOCARTA",
                            desc = "M2 BioCarta Pathways", tag = "MBioC"),
    "M2:CP:REACTOME" = list(coll = "C2", subcoll = "CP:REACTOME",
                            desc = "M2 Reactome Pathways", tag = "MReac"),
    "M2:CP:WIKIPATHWAYS" = list(coll = "C2", subcoll = "CP:WIKIPATHWAYS",
                                desc = "M2 WikiPathways", tag = "MWiki"),
    "M3:GTRD" = list(coll = "C3", subcoll = "GTRD",
                     desc = "M3 GTRD Transcription Factors", tag = "MTFTG"),
    "M3:MIRDB" = list(coll = "C3", subcoll = "MIRDB",
                      desc = "M3 miRDB", tag = "MMirDB"),
    "M5:GO:BP" = list(coll = "C5", subcoll = "GO:BP",
                      desc = "M5 GO Biological Process", tag = "MGoBP"),
    "M5:GO:CC" = list(coll = "C5", subcoll = "GO:CC",
                      desc = "M5 GO Cellular Component", tag = "MGoCC"),
    "M5:GO:MF" = list(coll = "C5", subcoll = "GO:MF",
                      desc = "M5 GO Molecular Function", tag = "MGoMF"),
    "M5:MPT" = list(coll = "C5", subcoll = "HPO",  # MPT map to HPO
                    desc = "M5 Mouse Phenotype Tumor", tag = "MPTumor"),
    "M7" = list(coll = "C7", subcoll = "",
                desc = "M7 Immunologic Signatures", tag = "MImmS"),
    "M8" = list(coll = "C8", subcoll = "",
                desc = "M8 Cell Type Signatures", tag = "MC8Cell")
  )

  # === 动态拉取可用的基因集集合列表 ===
  avail_colls <- msigdbr::msigdbr_collections(db_species = cfg$db_species) %>%
    dplyr::arrange(.data$gs_collection, .data$gs_subcollection)

  # === 构建菜单数据框 ===
  # 对于 MM，需要将 MM 代码映射到人类代码和显示信息
  menu_df <- lapply(1:nrow(avail_colls), function(i) {
    row <- avail_colls[i, ]
    mm_key <- paste0(row$gs_collection,
                     ifelse(row$gs_subcollection == "" || is.na(row$gs_subcollection),
                            "", paste0(":", row$gs_subcollection)))

    if (cfg$is_mouse && mm_key %in% names(mm_code_map)) {
      mapped <- mm_code_map[[mm_key]]
      # 人类集合代码（用于 msigdbr 函数调用）
      human_coll <- mapped$coll
      human_subcoll <- mapped$subcoll
      # 显示信息
      short_tag <- mapped$tag
      display_desc <- mapped$desc
      # 组合名称（用于 auto_select 匹配）
      combo_name <- mm_key  # 保持 MM 代码格式
    } else {
      # 人类模式：直接使用
      human_coll <- row$gs_collection
      human_subcoll <- row$gs_subcollection
      short_tag <- paste0(row$gs_collection,
                          ifelse(row$gs_subcollection == "" || is.na(row$gs_subcollection),
                                 "", paste0("_", gsub(":", "_", row$gs_subcollection))))
      display_desc <- row$gs_collection_name
      combo_name <- paste0(row$gs_collection,
                           ifelse(row$gs_subcollection == "" || is.na(row$gs_subcollection),
                                  "", paste0(":", row$gs_subcollection)))
    }

    data.frame(
      gs_collection = row$gs_collection,
      gs_subcollection = row$gs_subcollection,
      gs_collection_name = row$gs_collection_name,
      num_genesets = row$num_genesets,
      human_coll = human_coll,
      human_subcoll = human_subcoll,
      short_tag = short_tag,
      description = display_desc,
      combo_name = combo_name,
      stringsAsFactors = FALSE
    )
  })

  menu_df <- dplyr::bind_rows(menu_df)

  # === 交互界面逻辑 ===
  if (is.null(auto_select)) {
    message("\n", rep("=", 70))
    message(sprintf("Welcome to GSEAlens Pathway Wizard (%s)", cfg$display_name))
    message(rep("=", 70))
    message("Available Gene Set Collections:")
    message(sprintf("%-6s %-18s %s", "Tag", "Code", "Description"))
    message(rep("-", 70))

    for (i in 1:nrow(menu_df)) {
      cat(sprintf("[%2d] %-6s | %-18s | %s\n",
                  i, menu_df$short_tag[i], menu_df$combo_name[i], menu_df$description[i]))
    }

    user_input <- readline(prompt = "\nEnter selection (comma-separated, e.g., 1,3,5): ")
    selected_idx <- as.integer(unlist(strsplit(user_input, "[, ]+")))
    selected_idx <- selected_idx[!is.na(selected_idx) &
                                   selected_idx >= 1 &
                                   selected_idx <= nrow(menu_df)]

    if (length(selected_idx) == 0) {
      stop("Invalid input: No valid selection provided.")
    }
  } else {
    # 自动选择模式
    if (length(auto_select) == 1 && toupper(auto_select) == "ALL") {
      message(sprintf("Loading complete MSigDB database for %s...", cfg$display_name))
      selected_idx <- 1:nrow(menu_df)
    } else if (is.character(auto_select)) {
      # === 自动选择模式 ===

      # 1. 首先尝试直接匹配
      matched <- match(auto_select, menu_df$combo_name)

      # 2. 如果匹配失败，尝试匹配 short_tag
      if (any(is.na(matched))) {
        matched2 <- match(auto_select, menu_df$short_tag)
        na_idx <- which(is.na(matched))
        matched[na_idx] <- matched2[na_idx]
      }

      # 3. MM 模式下，如果仍然匹配失败，尝试人类代码转换
      if (cfg$is_mouse && any(is.na(matched))) {
        human_to_mm_map <- c(
          "H" = "MH", "C1" = "M1", "C2" = "M2", "C3" = "M3",
          "C4" = "M4", "C5" = "M5", "C6" = "M6", "C7" = "M7", "C8" = "M8"
        )
        still_missing <- auto_select[is.na(matched)]
        converted <- sapply(still_missing, function(x) {
          if (x %in% names(human_to_mm_map)) {
            converted_code <- human_to_mm_map[[x]]
            message(sprintf("[build_gsea_pathways] Converted '%s' -> '%s' for MM mode", x, converted_code))
            return(converted_code)
          }
          return(x)
        })
        matched3 <- match(converted, menu_df$combo_name)
        still_na <- which(is.na(matched))
        matched[still_na] <- matched3
      }

      # 4. 最终检查
      if (any(is.na(matched))) {
        missing <- auto_select[is.na(matched)]
        stop(sprintf("Collection not found: %s. Available: %s",
                     paste(missing, collapse = ", "),
                     paste(unique(menu_df$combo_name), collapse = ", ")))
      }
      selected_idx <- matched[!is.na(matched)]
    } else if (is.numeric(auto_select)) {
      selected_idx <- auto_select
      if (any(selected_idx < 1 | selected_idx > nrow(menu_df))) {
        stop(sprintf("Index out of range. Valid range: 1-%d", nrow(menu_df)))
      }
    } else {
      stop("Invalid auto_select parameter type.")
    }
  }

  selected_rows <- menu_df[selected_idx, ]

  # === 智能批次标签生成 ===
  if (length(auto_select) == 1 && toupper(auto_select) == "ALL") {
    super_tag <- sprintf("ALL_%s_MSigDB_Global", species)
  } else {
    super_tag <- ifelse(
      nrow(selected_rows) <= 4,
      paste(selected_rows$short_tag, collapse = "_"),
      sprintf("Mix%d_%s_%s", nrow(selected_rows), species, selected_rows$short_tag[1])
    )
  }

  message(sprintf("\nSelected %d collection(s). Batch Tag: [%s]",
                  nrow(selected_rows), super_tag))

  # === 核心提取：使用正确的物种和代码 ===
  pathway_list <- lapply(1:nrow(selected_rows), function(i) {
    row <- selected_rows[i, ]

    # 使用人类代码调用 msigdbr
    if (row$human_subcoll == "" || is.na(row$human_subcoll)) {
      msigdbr::msigdbr(species = cfg$msigdbr_species, collection = row$human_coll)
    } else {
      msigdbr::msigdbr(species = cfg$msigdbr_species,
                       collection = row$human_coll,
                       subcollection = row$human_subcoll)
    }
  })

  all_pathways <- dplyr::bind_rows(pathway_list)

  # === 动态列名检测 ===
  cat_col_name <- if ("gs_cat" %in% colnames(all_pathways)) "gs_cat" else "gs_collection"
  subcat_col_name <- if ("gs_subcat" %in% colnames(all_pathways)) "gs_subcat" else "gs_subcollection"

  # === 构建 TERM2GENE 映射表 ===
  TERM2GENE <- all_pathways %>%
    dplyr::select(gs_name, gene_symbol)

  # === 构建元数据字典 ===
  TERM2NAME <- all_pathways %>%
    dplyr::select(
      ID = gs_name,
      Description = gs_description,
      URL = gs_url,
      Collection = dplyr::all_of(cat_col_name),
      Subcollection = dplyr::all_of(subcat_col_name)
    ) %>%
    dplyr::distinct(ID, .keep_all = TRUE) %>%
    dplyr::mutate(
      Combo_Name = ifelse(
        Subcollection == "" | is.na(Subcollection),
        Collection,
        paste0(Collection, ":", Subcollection)
      ),
      Species = species
    )

  # === 返回结果列表 ===
  result <- list(
    TERM2GENE = TERM2GENE,
    meta_dict = TERM2NAME,
    SuperTag = super_tag,
    collections_used = selected_rows,
    species = species
  )

  message(sprintf("\n[build_gsea_pathways] Done! Built %d pathways x %d genes for %s",
                  nrow(TERM2NAME),
                  length(unique(TERM2GENE$gene_symbol)),
                  cfg$display_name))

  return(result)
}
