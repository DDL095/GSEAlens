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
saveRDS(dds_se, file.path(out_dir, "preprocessed_dds_se.rds"), compress = "xz")
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
