The name "lens" in GSEAlens symbolizes that this package acts like a magnifying glass, exploring relevant GSEA pathways.

GSEAlens provides a web-based interface to display pathway introductions and descriptions, incorporating AI-assisted export functionality for pathway enrichment results. This makes GSEA enrichment analysis more convenient, while also allowing users to import translated pathway annotations and summaries as needed.

By packaging the workflow and standardizing input formats, this R package simplifies the viewing and exploration process of GSEA enrichment analysis.

Since multi-group comparisons are more convenient using the no-intercept group specification in the limma-voom pipeline, the code in this project is based on the no-intercept group setting for result analysis and processing.

Based on the reference `Fast gene set enrichment analysis`,`https://www.biorxiv.org/content/10.1101/060012v3`, GSEAlens accepts Wald values from DESeq2::result (`result$stat`) and t-values from limma::fit (`fit$t`) as input for fgsea(Based o `clusterProfiler::GSEA`) and performs ranking(draw based on `GseaVis`). The corresponding comparisons are automatically generated based on the group information contained in the input data. Therefore, users need to specify the target categorical factor (e.g., condition or group). The DESeq2 pipeline allows additional additive covariates (e.g., batch, subject, sex) in the model for correction. In the limma-voom pipeline, the group needs to be specified in a no-intercept manner `~0+group`) DESeq2 supports `~ group`, `~ batch + group`, and `~ subject + group`, but ultimately only supports group-based analysis. GSEAlens does not currently support the generation and exploration of batch GSEA comparison results with interaction terms, continuous variable effects, or complex custom contrasts.

package install 

```
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}
devtools::install_github("DDL095/GSEAlens")
```
