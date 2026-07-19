#' @name precomputed_gseares
#' @docType data
#' @title Pre-computed GSEA results from the airway dataset
#' @description A pre-computed `GseaRes` object built from the airway dataset
#'   (DESeq2 + clusterProfiler::GSEA). Used in examples and vignettes so that
#'   users can try visualization and code-generation functions without having
#'   to re-run the full GSEA pipeline.
#' @format A `GseaRes` object (list-like).
#' @source Built from the `airway` Bioconductor dataset via DESeq2 +
#'   clusterProfiler::GSEA.
#' @examples
#' data(precomputed_gseares, package = "GSEAlens")
#' class(precomputed_gseares)
NULL

#' @name preprocessed_dds
#' @docType data
#' @title Pre-computed DESeq2 DESeqDataSet from the airway dataset
#' @description A slimmed `DESeqDataSet` object fitted with `DESeq()`,
#'   suitable as input to [GSEAlens::setup_gsea_env()].
#' @format A `DESeqDataSet` object.
#' @source Built from the `airway` Bioconductor dataset via DESeq2.
#' @examples
#' data(preprocessed_dds, package = "GSEAlens")
#' class(preprocessed_dds)
NULL

#' @name preprocessed_dds_se
#' @docType data
#' @title Pre-computed DESeqDataSet (SummarizedExperiment form)
#' @description A `DESeqDataSet` based on a `SummarizedExperiment` input,
#'   suitable for demonstrating the SE-based entry point of
#'   [GSEAlens::setup_gsea_env()].
#' @format A `DESeqDataSet` object.
#' @source Built from the `airway` Bioconductor dataset via DESeq2.
#' @examples
#' data(preprocessed_dds_se, package = "GSEAlens")
#' class(preprocessed_dds_se)
NULL

#' @name preprocessed_limma
#' @docType data
#' @title Pre-computed limma-voom fit and DGEList from the airway dataset
#' @description A list with elements `fit` (`MArrayLM`) and
#'   `gsea_limma_voom_data` (filtered `DGEList`), suitable as input to
#'   [GSEAlens::setup_gsea_env()] for the limma-voom workflow.
#' @format A list with elements `fit` and `gsea_limma_voom_data`.
#' @source Built from the `airway` Bioconductor dataset via limma-voom.
#' @examples
#' data(preprocessed_limma, package = "GSEAlens")
#' names(preprocessed_limma)
NULL

#' @name gsea_pathwaysets_toy
#' @docType data
#' @title Toy MSigDB pathway collection (Hallmark + KEGG_LEGACY)
#' @description A lightweight pathway object for examples and vignettes,
#'   containing 236 pathways (Hallmark + KEGG_LEGACY) pre-assembled by
#'   [GSEAlens::build_gsea_pathways()]. See
#'   `inst/scripts/make_gsea_pathwaysets_toy.R` for regeneration.
#' @format A list with elements `TERM2GENE`, `meta_dict`, `SuperTag`,
#'   `collections_used`, `species`.
#' @source Assembled from msigdbr (Hallmark + KEGG_LEGACY collections).
#' @examples
#' data(gsea_pathwaysets_toy, package = "GSEAlens")
#' names(gsea_pathwaysets_toy)
NULL

#' @name gsea_pathwaysets_toy_hallmark
#' @docType data
#' @title Toy MSigDB pathway collection (Hallmark only)
#' @description A minimal Hallmark-only pathway object for fast examples
#'   and tests. See `inst/scripts/make_gsea_pathwaysets_toy.R` for
#'   regeneration.
#' @format A list with elements `TERM2GENE`, `meta_dict`, `SuperTag`,
#'   `collections_used`, `species`.
#' @source Assembled from msigdbr (Hallmark collection).
#' @examples
#' data(gsea_pathwaysets_toy_hallmark, package = "GSEAlens")
#' names(gsea_pathwaysets_toy_hallmark)
NULL
