# Section: Pathway Database Builder ----

#' @title Build GSEA Pathway Database with Interactive Selection
#' @description Provides an interactive or automated gene set selection interface
#'   based on msigdbr, generating standardized pathway objects with intelligent
#'   semantic tagging. Supports both human (Homo sapiens) and mouse (Mus musculus)
#'   species.
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
#' if(interactive()){
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
  # === Section: Species Parameter Standardization === ----

  species <- toupper(trimws(as.character(species)))

  # 鐗╃閰嶇疆琛?
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

  # === Section: MM Database Code to msigdbr Function Parameter Mapping === ----

  # msigdbr_collections(db_species = "MM") returns MM codes like M1, M2, MH
  # but msigdbr(species = "mouse", ...) only accepts human codes (C1, C2, H, etc.)
  mm_code_map <- list(
    "MH" = list(
      coll = "H", subcoll = "",
      desc = "MH Hallmark Gene Sets", tag = "MHall"
    ),
    "M1" = list(
      coll = "C1", subcoll = "",
      desc = "M1 Positional Gene Sets", tag = "M1Pos"
    ),
    "M2:CGP" = list(
      coll = "C2", subcoll = "CGP",
      desc = "M2 Chemical and Genetic Perturbations", tag = "MCGP"
    ),
    "M2:CP:BIOCARTA" = list(
      coll = "C2", subcoll = "CP:BIOCARTA",
      desc = "M2 BioCarta Pathways", tag = "MBioC"
    ),
    "M2:CP:REACTOME" = list(
      coll = "C2", subcoll = "CP:REACTOME",
      desc = "M2 Reactome Pathways", tag = "MReac"
    ),
    "M2:CP:WIKIPATHWAYS" = list(
      coll = "C2", subcoll = "CP:WIKIPATHWAYS",
      desc = "M2 WikiPathways", tag = "MWiki"
    ),
    "M3:GTRD" = list(
      coll = "C3", subcoll = "GTRD",
      desc = "M3 GTRD Transcription Factors", tag = "MTFTG"
    ),
    "M3:MIRDB" = list(
      coll = "C3", subcoll = "MIRDB",
      desc = "M3 miRDB", tag = "MMirDB"
    ),
    "M5:GO:BP" = list(
      coll = "C5", subcoll = "GO:BP",
      desc = "M5 GO Biological Process", tag = "MGoBP"
    ),
    "M5:GO:CC" = list(
      coll = "C5", subcoll = "GO:CC",
      desc = "M5 GO Cellular Component", tag = "MGoCC"
    ),
    "M5:GO:MF" = list(
      coll = "C5", subcoll = "GO:MF",
      desc = "M5 GO Molecular Function", tag = "MGoMF"
    ),
    "M5:MPT" = list(
      coll = "C5", subcoll = "HPO", # MPT maps to HPO
      desc = "M5 Mouse Phenotype Tumor", tag = "MPTumor"
    ),
    "M7" = list(
      coll = "C7", subcoll = "",
      desc = "M7 Immunologic Signatures", tag = "MImmS"
    ),
    "M8" = list(
      coll = "C8", subcoll = "",
      desc = "M8 Cell Type Signatures", tag = "MC8Cell"
    )
  )

  # === Section: Dynamically Fetch Available Gene Set Collections === ----

  avail_colls <- msigdbr::msigdbr_collections(db_species = cfg$db_species) %>%
    dplyr::arrange(.data$gs_collection, .data$gs_subcollection)

  # === Section: Build Menu Dataframe === ----

  # For MM mode, map MM codes to human codes and display info
  menu_df <- lapply(seq_len(nrow(avail_colls)), function(i) {
    row <- avail_colls[i, ]
    mm_key <- paste0(
      row$gs_collection,
      ifelse(row$gs_subcollection == "" || is.na(row$gs_subcollection),
        "", paste0(":", row$gs_subcollection)
      )
    )

    if (cfg$is_mouse && mm_key %in% names(mm_code_map)) {
      mapped <- mm_code_map[[mm_key]]
      # Human collection code (for msigdbr function calls)
      human_coll <- mapped$coll
      human_subcoll <- mapped$subcoll
      # Display info
      short_tag <- mapped$tag
      display_desc <- mapped$desc
      combo_name <- mm_key # 淇濇寔 MM 浠ｇ爜鏍煎紡
      # Combined name (for auto_select matching)
    } else {
      # 浜虹被妯″紡锛氱洿鎺ヤ娇鐢?
      human_coll <- row$gs_collection
      human_subcoll <- row$gs_subcollection
      short_tag <- paste0(
        row$gs_collection,
        ifelse(row$gs_subcollection == "" || is.na(row$gs_subcollection),
          "", paste0("_", gsub(":", "_", row$gs_subcollection))
        )
      )
      display_desc <- row$gs_collection_name
      combo_name <- paste0(
        row$gs_collection,
        ifelse(row$gs_subcollection == "" || is.na(row$gs_subcollection),
          "", paste0(":", row$gs_subcollection)
        )
      )
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

  # === Section: Interactive Interface Logic === ----

  if (is.null(auto_select)) {
    message("\n", rep("=", 70))
    message(sprintf("Welcome to GSEAlens Pathway Wizard (%s)", cfg$display_name))
    message(rep("=", 70))
    message("Available Gene Set Collections:")
    message(sprintf("%-6s %-18s %s", "Tag", "Code", "Description"))
    message(rep("-", 70))

    for (i in seq_len(nrow(menu_df))) {
      message(sprintf(
        "[%2d] %-6s | %-18s | %s",
        i, menu_df$short_tag[i], menu_df$combo_name[i], menu_df$description[i]
      ))
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
    # Auto-select mode
    if (length(auto_select) == 1 && toupper(auto_select) == "ALL") {
      message(sprintf("Loading complete MSigDB database for %s...", cfg$display_name))
      selected_idx <- seq_len(nrow(menu_df))
    } else if (is.character(auto_select)) {
      # === Section: Auto-Select Mode === ----

      # 1. First attempt direct matching
      matched <- match(auto_select, menu_df$combo_name)

      # 2. If matching fails, try matching short_tag
      if (any(is.na(matched))) {
        matched2 <- match(auto_select, menu_df$short_tag)
        na_idx <- which(is.na(matched))
        matched[na_idx] <- matched2[na_idx]
      }

      # 3. In MM mode, if still failing, try human code conversion
      if (cfg$is_mouse && any(is.na(matched))) {
        human_to_mm_map <- c(
          "H" = "MH", "C1" = "M1", "C2" = "M2", "C3" = "M3",
          "C4" = "M4", "C5" = "M5", "C6" = "M6", "C7" = "M7", "C8" = "M8"
        )
        still_missing <- auto_select[is.na(matched)]
        converted <- vapply(still_missing, function(x) {
          if (x %in% names(human_to_mm_map)) {
            converted_code <- human_to_mm_map[[x]]
            message(sprintf("[build_gsea_pathways] Converted '%s' -> '%s' for MM mode", x, converted_code))
            return(converted_code)
          }
          return(x)
        }, character(1))
        matched3 <- match(converted, menu_df$combo_name)
        still_na <- which(is.na(matched))
        matched[still_na] <- matched3
      }

      # 4. Final validation check
      if (any(is.na(matched))) {
        missing <- auto_select[is.na(matched)]
        stop(sprintf(
          "Collection not found: %s. Available: %s",
          paste(missing, collapse = ", "),
          paste(unique(menu_df$combo_name), collapse = ", ")
        ))
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

  # === Section: Intelligent Batch Tag Generation === ----

  if (length(auto_select) == 1 && toupper(auto_select) == "ALL") {
    super_tag <- sprintf("ALL_%s_MSigDB_Global", species)
  } else {
    super_tag <- ifelse(
      nrow(selected_rows) <= 4,
      paste(selected_rows$short_tag, collapse = "_"),
      sprintf("Mix%d_%s_%s", nrow(selected_rows), species, selected_rows$short_tag[1])
    )
  }

  message(sprintf(
    "\nSelected %d collection(s). Batch Tag: [%s]",
    nrow(selected_rows), super_tag
  ))

  # === Section: Core Extraction Using Correct Species and Codes === ----

  pathway_list <- lapply(seq_len(nrow(selected_rows)), function(i) {
    row <- selected_rows[i, ]

    # Use human codes to call msigdbr
    if (row$human_subcoll == "" || is.na(row$human_subcoll)) {
      msigdbr::msigdbr(species = cfg$msigdbr_species, collection = row$human_coll)
    } else {
      msigdbr::msigdbr(
        species = cfg$msigdbr_species,
        collection = row$human_coll,
        subcollection = row$human_subcoll
      )
    }
  })

  all_pathways <- dplyr::bind_rows(pathway_list)

  # === Section: Dynamic Column Name Detection === ----

  cat_col_name <- if ("gs_cat" %in% colnames(all_pathways)) "gs_cat" else "gs_collection"
  subcat_col_name <- if ("gs_subcat" %in% colnames(all_pathways)) "gs_subcat" else "gs_subcollection"

  # === Section: Build TERM2GENE Mapping Table === ----

  TERM2GENE <- all_pathways %>%
    dplyr::select(gs_name, gene_symbol)

  # === Section: Build Metadata Dictionary === ----

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

  # === Section: Return Result List === ----

  result <- list(
    TERM2GENE = TERM2GENE,
    meta_dict = TERM2NAME,
    SuperTag = super_tag,
    collections_used = selected_rows,
    species = species
  )

  message(sprintf(
    "\n[build_gsea_pathways] Done! Built %d pathways x %d genes for %s",
    nrow(TERM2NAME),
    length(unique(TERM2GENE$gene_symbol)),
    cfg$display_name
  ))

  return(result)
}
