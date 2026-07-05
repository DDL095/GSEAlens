#!/usr/bin/env Rscript

# Description -------------------------------------------------------------

# Pre-compute the three input objects (`fit`, `dds_se`, `dds`) that the main

# GSEAlens vignette consumes, following EXACTLY the steps documented in

# `vignettes/GSEAlens-preprocessing.Rmd`.

#

# Why: running these steps during vignette compilation on the Bioconductor

# build machine is slow (DESeq2 DESeq() + limma-voom pipeline on airway) and

# triggers the "vignette rebuilding took > 60s" NOTE. By shipping the

# pre-computed objects in inst/extdata/, the main vignette can set

# `eval=TRUE` on the downstream setup_gsea_env() chunks so reviewers see a

# fully runnable workflow.

#

# Run manually (maintainer only) when DESeq2/limma/airway versions bump:

#

#   Rscript inst/scripts/make_preprocessed_inputs.R

#

# Output files (written to inst/extdata/):

#   - preprocessed_limma.rds    (list with fit + gsea_limma_voom_data)

#   - preprocessed_dds_se.rds   (DESeq-transformed RangedSummarizedExperiment)

#   - preprocessed_dds.rds      (DESeqDataSet from count matrix)

# -------------------------------------------------------------------------



stopifnot(

  requireNamespace("airway",    quietly = TRUE),

  requireNamespace("limma",     quietly = TRUE),

  requireNamespace("edgeR",     quietly = TRUE),

  requireNamespace("DESeq2",    quietly = TRUE),

  requireNamespace("SummarizedExperiment", quietly = TRUE)

)



library(airway)

library(edgeR)

library(limma)

library(DESeq2)

library(SummarizedExperiment)



data(airway)

expression_data <- airway



# ===== limma-voom workflow (verbatim from preprocessing vignette) =====

group_level <- expression_data$dex

design <- model.matrix(~0 + group_level)

colnames(design) <- levels(group_level)

compare_end <- combn(levels(group_level), 2, simplify = FALSE)

contrast_strings <- sapply(compare_end, function(x) paste(x[2], x[1], sep = " - "))

contrast_matrix <- makeContrasts(contrasts = contrast_strings, levels = design)



genes_df <- data.frame(

  gene_id = rowData(expression_data)$gene_id,

  symbol = rowData(expression_data)$symbol,

  gene_biotype = rowData(expression_data)$gene_biotype

)

genes_df$Length <- rowData(expression_data)$gene_seq_end -

                   rowData(expression_data)$gene_seq_start + 1



gsea_limma_voom_data <- DGEList(

  counts = assay(expression_data, "counts"),

  genes = genes_df,

  norm.factors = NULL,

  group = group_level,

  remove.zeros = TRUE

)



total_counts <- cpm(gsea_limma_voom_data) |> rowSums()

dup_symbols <- gsea_limma_voom_data$genes$symbol[duplicated(gsea_limma_voom_data$genes$symbol)]

keep <- rep(TRUE, nrow(gsea_limma_voom_data))

for (gene in dup_symbols) {

  idx <- which(gsea_limma_voom_data$genes$symbol == gene)

  best_idx <- idx[which.max(total_counts[idx])]

  remove_idx <- idx[idx != best_idx]

  keep[remove_idx] <- FALSE

}

gsea_limma_voom_data <- gsea_limma_voom_data[keep, ]

rownames(gsea_limma_voom_data) <- gsea_limma_voom_data$genes$symbol



keep_biotype <- gsea_limma_voom_data$genes$gene_biotype == "protein_coding"

gsea_limma_voom_data <- gsea_limma_voom_data[keep_biotype, ]

gsea_limma_voom_data <- normLibSizes(gsea_limma_voom_data, method = "TMM")

isexpr <- rowSums(cpm(gsea_limma_voom_data) > 1) >= 3

gsea_limma_voom_data <- gsea_limma_voom_data[isexpr, ]



VoomOutPut <- voom(gsea_limma_voom_data, design)

fit <- lmFit(object = VoomOutPut, design = design) |>

  contrasts.fit(contrasts = contrast_matrix) |>

  eBayes()



# ===== DESeq2 SummarizedExperiment workflow (verbatim) =====

dds_se <- DESeqDataSet(expression_data, design = ~ cell + dex)

gene_biotypes <- rowData(dds_se)$gene_biotype

keep_protein_coding <- gene_biotypes == "protein_coding"

dds_se <- dds_se[keep_protein_coding, ]



smallestGroupSize <- 3

keep <- rowSums(counts(dds_se) >= 10) >= smallestGroupSize

dds_se <- dds_se[keep, ]

rownames(dds_se) <- rowData(dds_se)$gene_name



total_counts <- rowSums(assay(dds_se))

dup_genes <- rownames(dds_se)[duplicated(rownames(dds_se))]

keep <- rep(TRUE, nrow(dds_se))

for (gene in dup_genes) {

  idx <- which(rownames(dds_se) == gene)

  best_idx <- idx[which.max(total_counts[idx])]

  remove_idx <- idx[idx != best_idx]

  keep[remove_idx] <- FALSE

}

dds_se <- dds_se[keep, ]

dds_se <- DESeq(dds_se)



# ===== DESeq2 Count Matrix workflow (verbatim) =====

DDS_rawdata <- expression_data

gene_biotypes <- rowData(DDS_rawdata)$gene_biotype

keep_protein_coding <- gene_biotypes == "protein_coding"

DDS_rawdata <- DDS_rawdata[keep_protein_coding, ]



smallestGroupSize <- 3

keep_epd <- rowSums(assay(DDS_rawdata, "counts") >= 10) >= smallestGroupSize

DDS_rawdata <- DDS_rawdata[keep_epd, ]

rownames(DDS_rawdata) <- rowData(DDS_rawdata)$gene_name



total_counts <- rowSums(assay(DDS_rawdata))

keep_name <- rep(TRUE, nrow(DDS_rawdata))

dup_genes <- rownames(DDS_rawdata)[duplicated(rownames(DDS_rawdata))]

for (gene in dup_genes) {

  idx <- which(rownames(DDS_rawdata) == gene)

  best_idx <- idx[which.max(total_counts[idx])]

  remove_idx <- idx[idx != best_idx]

  keep_name[remove_idx] <- FALSE

}

DDS_rawdata <- DDS_rawdata[keep_name, ]



cts <- assay(DDS_rawdata, "counts")

coldata <- as.data.frame(colData(DDS_rawdata))

coldata <- coldata[, c("cell", "dex")]

coldata$cell <- factor(coldata$cell)

coldata$dex <- factor(coldata$dex)



dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ dex)

dds <- DESeq(dds)



# ===== Save ----------------------------------------------------------------------

out_dir <- file.path("inst", "extdata")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)



saveRDS(

  list(fit = fit, gsea_limma_voom_data = gsea_limma_voom_data),

  file.path(out_dir, "preprocessed_limma.rds"),

  compress = "xz"

)

# Slim preprocessed_dds_se.rds before saving to stay below the Bioconductor
# 5 MB extdata limit. The full DESeqDataSet from DESeq() on the pasilla data
# is ~8 MB on disk because each gene is wrapped in a GRangesList of ~20
# transcript sub-ranges and DESeq() also stores mu/H/cooks assay layers.
# GSEAlens only needs the gene-level fits, so the slimming is non-destructive
# for every code path that GSEAlens or its vignettes exercise.
#
# What is removed and why (all verified to be unused by GSEAlens):
#   * mu / H / cooks assays  -> DESeq() fitting intermediates. GSEAlens
#                               only reads the 'counts' assay and the
#                               Wald statistics stored in mcols (via
#                               DESeq2::results()). log2FoldChange / pvalue
#                               / padj are bit-identical before vs. after.
#   * GRangesList -> GRanges -> keep one representative range per gene.
#                               rowRanges() is never accessed by GSEAlens.
#   * 5 columns of Ensembl coordinate metadata (gene_seq_start, gene_seq_end,
#     seq_name, seq_strand, seq_coord_system) -> gene_id / gene_name /
#     symbol / gene_biotype are KEPT so users can still inspect gene info.
#
# Users who want the full, untrimmed DESeqDataSet for teaching DESeq2 itself
# should re-run this script without calling slim_dds_se(), or follow the
# preprocessing vignette to rebuild one from their own count matrix.
slim_dds_se <- function(dds_se) {
  # (1) Drop DESeq() intermediate assays (mu/H/cooks); keep counts only.
  SummarizedExperiment::assays(dds_se) <-
    SummarizedExperiment::assays(dds_se)["counts"]

  # (2) GRangesList -> GRanges: keep one representative range per gene.
  rr <- SummarizedExperiment::rowRanges(dds_se)
  if (inherits(rr, "GRangesList")) {
    grp_lens <- S4Vectors::elementNROWS(rr)
    first_idx <- cumsum(c(1L, head(grp_lens, -1L)))
    gr_flat <- unlist(rr, use.names = FALSE)
    gr_first <- gr_flat[first_idx]
    mcols(gr_first) <- mcols(rr)
    names(gr_first) <- rownames(dds_se)
    SummarizedExperiment::rowRanges(dds_se) <- gr_first
  }

  # (3) Drop 5 Ensembl coordinate columns; keep gene_id / gene_name /
  #     symbol / gene_biotype so gene-level metadata stays queryable.
  redundant_anno <- c("gene_seq_start", "gene_seq_end",
                      "seq_name", "seq_strand", "seq_coord_system")
  rr2 <- SummarizedExperiment::rowRanges(dds_se)
  mc2 <- mcols(rr2)
  mcols(rr2) <- mc2[, setdiff(colnames(mc2), redundant_anno), drop = FALSE]
  SummarizedExperiment::rowRanges(dds_se) <- rr2

  stopifnot(validObject(dds_se))
  dds_se
}

saveRDS(slim_dds_se(dds_se), file.path(out_dir, "preprocessed_dds_se.rds"), compress = "xz")
saveRDS(dds,    file.path(out_dir, "preprocessed_dds.rds"),    compress = "xz")



# ===== Report --------------------------------------------------------------------

files <- c("preprocessed_limma.rds", "preprocessed_dds_se.rds", "preprocessed_dds.rds")

sizes <- file.size(file.path(out_dir, files))

md5 <- tools::md5sum(file.path(out_dir, files))



cat("Wrote:\n")

for (i in seq_along(files)) {

  cat(sprintf("  %-32s  %7.1f KB  md5=%s\n",

              files[i], sizes[i] / 1024, md5[i]))

}

cat(sprintf("\nlimma: fit has %d genes; gsea_limma_voom_data has %d rows\n",

            nrow(fit$coefficients), nrow(gsea_limma_voom_data)))

cat(sprintf("dds_se: %d genes x %d samples; results names: %s\n",

            nrow(dds_se), ncol(dds_se),

            paste(resultsNames(dds_se), collapse = ", ")))

cat(sprintf("dds:    %d genes x %d samples; results names: %s\n",

            nrow(dds), ncol(dds),

            paste(resultsNames(dds), collapse = ", ")))



cat("\nairway  version:", as.character(packageVersion("airway")), "\n")

cat("DESeq2  version:", as.character(packageVersion("DESeq2")), "\n")

cat("limma   version:", as.character(packageVersion("limma")), "\n")

cat("edgeR   version:", as.character(packageVersion("edgeR")), "\n")

cat("Generation date:", format(Sys.Date(), "%Y-%m-%d"), "\n")

