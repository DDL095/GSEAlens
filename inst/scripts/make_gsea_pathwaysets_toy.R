#!/usr/bin/env Rscript

# Description -------------------------------------------------------------

# Pre-compute a lightweight `build_gsea_pathways()` result for vignette use.

#

# Rationale: the full vignette call uses `auto_select = c(17,18,19,20,26)`

# which loads 5 MSigDB collections (~8000 pathways, ~30-90 s on the

# Bioconductor build machine). To keep vignette compilation fast and avoid

# the "vignette rebuilding took > 60s" NOTE, we ship a pre-computed

# lightweight object in inst/extdata/.

#

# Run this script manually whenever the msigdbr version bumps or the

# build_gsea_pathways() return structure changes:

#

#   Rscript inst/scripts/make_gsea_pathwaysets_toy.R

#

# Output files (written to inst/extdata/):

#   - gsea_pathwaysets_toy_hallmark.rds  (~80 KB, H only, 50 pathways)

#   - gsea_pathwaysets_toy.rds           (~500 KB, H + C2:CP:KEGG, ~300 pathways)

# -------------------------------------------------------------------------



# Use the installed GSEAlens (>= 0.99.9) so the return structure matches

# exactly what users get.

library(GSEAlens)



out_dir <- file.path("inst", "extdata")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)



# NOTE on `auto_select` form:

# Numeric indices (e.g. c(17, 19)) depend on the row order of

# msigdbr::msigdbr_collections(), which can drift across msigdbr releases.

# Prefer the **string** form ("H", "C2:CP:KEGG") for reproducibility.



# ---- Hallmark only (smallest, ~50 pathways) -----------------------------

toy_hallmark <- build_gsea_pathways(species = "HS", auto_select = "H")

saveRDS(

  toy_hallmark,

  file.path(out_dir, "gsea_pathwaysets_toy_hallmark.rds"),

  compress = "xz"

)



# ---- Hallmark + C2:CP:KEGG_LEGACY (demo set, ~200 pathways) --------------

# KEGG_LEGACY gives downstream GSEA something to enrich on with the airway

# dataset (drug-treatment contrast) while still keeping the file < 1 MB.

# NOTE: msigdbr >= 25 split KEGG into KEGG_LEGACY and KEGG_MEDICUS; we use

# KEGG_LEGACY because its pathway names are the classic KEGG_* identifiers

# users recognise.

toy_demo <- build_gsea_pathways(species = "HS", auto_select = c("H", "C2:CP:KEGG_LEGACY"))

saveRDS(

  toy_demo,

  file.path(out_dir, "gsea_pathwaysets_toy.rds"),

  compress = "xz"

)



# ---- Report -------------------------------------------------------------

sizes <- file.size(file.path(

  out_dir,

  c("gsea_pathwaysets_toy_hallmark.rds", "gsea_pathwaysets_toy.rds")

))

md5 <- tools::md5sum(file.path(

  out_dir,

  c("gsea_pathwaysets_toy_hallmark.rds", "gsea_pathwaysets_toy.rds")

))



cat(sprintf(

  "Wrote:\n  %s  (%.1f KB, %d pathways, %d genes, md5=%s)\n  %s  (%.1f KB, %d pathways, %d genes, md5=%s)\n",

  "gsea_pathwaysets_toy_hallmark.rds", sizes[1] / 1024,

  nrow(toy_hallmark$meta_dict),

  length(unique(toy_hallmark$TERM2GENE$gene_symbol)),

  md5[1],

  "gsea_pathwaysets_toy.rds", sizes[2] / 1024,

  nrow(toy_demo$meta_dict),

  length(unique(toy_demo$TERM2GENE$gene_symbol)),

  md5[2]

))



cat("\nmsigdbr version:", as.character(packageVersion("msigdbr")), "\n")

cat("GSEAlens version:", as.character(packageVersion("GSEAlens")), "\n")

cat("Generation date:", format(Sys.Date(), "%Y-%m-%d"), "\n")

