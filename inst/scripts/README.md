# Data Generation for inst/extdata

This directory contains scripts and documentation describing how the data in `inst/extdata/` was generated.

## pathway_annotations.csv

**Source**: Manually curated pathway annotation file for demonstration and custom pathway visualization purposes.

**Columns**:

- `ID`: Pathway identifier (e.g., HALLMARK_HYPOXIA, REACTOME_EXAMPLE)
- `Brief_Description_CN`: Brief description in Simplified Chinese
- `Brief_Description_EN`: Brief description in English
- `Abstract`: Detailed pathway abstract text
- `Custom_Note`: User-defined custom notes for the pathway

**Licensing**: This file is provided as example data for GSEAlens package functionality demonstration.

**Generation Method**:

1. Pathway IDs were selected from commonly used gene set collections (MSigDB Hallmark, Reactome, KEGG)
2. Descriptions were curated based on official pathway documentation
3. The file can be extended by users to add custom pathway annotations for their own gene sets

**Usage in GSEAlens**:
This annotation file can be loaded via the `create_addition_data()` function to provide custom pathway descriptions in the Shiny application.

## Pre-computed Vignette Objects

The main vignette (`vignettes/GSEAlens.Rmd` and `vignettes/GSEAlens-vignette-zh.Rmd`)
consumes five pre-computed RDS files from `inst/extdata/`. They exist so the
vignette can set `eval=TRUE` on its data-loading and object-assembly chunks
without re-running the slow upstream steps (DESeq2 DESeq, limma-voom, msigdbr
bulk loading) on every Bioconductor build.

| File                                  | Generator script                | Contents                                                                                                                                                                                                                                                                                                                                |
| ------------------------------------- | ------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `gsea_pathwaysets_toy_hallmark.rds` | `make_gsea_pathwaysets_toy.R` | `build_gsea_pathways("HS", "H")` result; 50 Hallmark pathways x 4384 genes (~22 KB)                                                                                                                                                                                                                                                   |
| `gsea_pathwaysets_toy.rds`          | `make_gsea_pathwaysets_toy.R` | `build_gsea_pathways("HS", c("H","C2:CP:KEGG_LEGACY"))` result; 236 pathways x 7337 genes (~46 KB)                                                                                                                                                                                                                                    |
| `preprocessed_limma.rds`            | `make_preprocessed_inputs.R`  | `list(fit = ..., gsea_limma_voom_data = ...)` from airway, limma-voom (~1 MB)                                                                                                                                                                                                                                                         |
| `preprocessed_dds_se.rds`           | `make_preprocessed_inputs.R`  | DESeq-transformed RangedSummarizedExperiment (13246 x 8) (~3 MB, **slimmed**; see below)                                                                                                                                                                                                                         |
| `preprocessed_dds.rds`              | `make_preprocessed_inputs.R`  | DESeqDataSet from count matrix (13246 x 8) (~4 MB)                                                                                                                                                                                                                                                                                      |
| `precomputed_gseares.rds`           | `make_precomputed_gseares.R`  | `batch_calc_gsea()` result ("GseaRes" list) from `preprocessed_limma.rds` + `gsea_pathwaysets_toy.rds`; 2 contrasts, ~20 enriched pathways each (~1.3 MB). Used by accessor examples (`get_expr_matrix`, `get_de_table`, `extract_gsea_task`, `inspect_gsea_res`, ...) so they can run in < 1 s without recomputing GSEA. |

**Regeneration**:

```sh
# From the package root:
Rscript inst/scripts/make_gsea_pathwaysets_toy.R
Rscript inst/scripts/make_preprocessed_inputs.R
Rscript inst/scripts/make_precomputed_gseares.R
```

Run these whenever `msigdbr`, `airway`, `DESeq2`, `limma`, `clusterProfiler`,
or `fgsea` see a version bump, or when `build_gsea_pathways()` /
`batch_calc_gsea()` / the preprocessing vignette change their filtering logic.

**Note on dependency order**: `make_precomputed_gseares.R` consumes the
outputs of the other two scripts, so re-run it **last** if you regenerate
all three.

**Versions used for the current shipped copies**:

- msigdbr 26.1.0, airway 1.32.0, DESeq2 1.52.0, limma 3.68.2, edgeR 4.10.0
- clusterProfiler 4.16.0, fgsea 1.36.0
- Generated on 2026-06-28

**Why `preprocessed_dds_se.rds` is "slimmed"**:

The raw `DESeqDataSet` produced by `DESeq()` on the airway data is ~8 MB on
disk, which exceeds the Bioconductor 5 MB extdata limit. The `slim_dds_se()`
helper inside `make_preprocessed_inputs.R` trims three categories of data
that GSEAlens never reads, reducing the file to ~3 MB without changing any
GSEAlens-visible behaviour:

| Removed                                | Reason                                                                                                                                  |
| -------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------- |
| `mu` / `H` / `cooks` assays            | DESeq() fitting intermediates. GSEAlens only needs `counts` + the Wald statistics in `mcols` (via `DESeq2::results()`).                 |
| `GRangesList` -> `GRanges` (1 per gene) | Each gene's ~20 transcript sub-ranges carry coordinate info GSEAlens never inspects.                                                    |
| 5 Ensembl coordinate columns           | `gene_seq_start`, `gene_seq_end`, `seq_name`, `seq_strand`, `seq_coord_system`. `gene_id` / `gene_name` / `symbol` / `gene_biotype` are **kept**. |

Verified post-slimming: `validObject()` passes, `DESeq2::results(dds_se)`
returns bit-identical `log2FoldChange` / `pvalue` / `padj`, and all GSEAlens
vignette chunks produce the same output. Users who want the full, untrimmed
`DESeqDataSet` for teaching DESeq2 itself should re-run this script with
`slim_dds_se()` removed, or follow the preprocessing vignette.
