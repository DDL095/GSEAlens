<!-- NEWS.md is maintained by https://cynkra.github.io/fledge, do not edit -->

# GSEAlens 0.99.5

- Fix "plotting error: missing value where TRUE/FALSE needed" when switching gene set subgroup (e.g. H -> C2/C5) with stale pathway selections in Joint Plot and Multi-Plot modules. Add defensive guard in `plot_directional_gsea()` to validate pathway IDs against the active GSEA result, and add null check in `.gs_info()` for missing gene sets.
- Enable cross-collection pathway overlay in Multi-Plot (combined pathway plotting). Multi-Plot now uses the full `gsea_res` (all gene set collections) instead of the sliced subset, allowing pathways from different collections (e.g. H + C2) to be plotted together on the same enrichment chart, consistent with Joint Canvas behavior.
- Fix UTF-8 BOM (Byte Order Mark) encoding in R source files that caused `parse()` failures on strict parsers. All R/*.R files now use clean UTF-8 without BOM, aligning with R package best practices and BiocCheck requirements (LEVEL: WARNING on BOM-annotated files).

# GSEAlens 0.99.4

- Replace the "Export Boxplot Data Code" button in the Quadrant module (panel 4 "Full Expression Distribution") with an "Export Boxplot Data" modal that directly shows the per-sample expression values in two tab-separated formats: wide (first row = sample IDs, second row = gene expression values, left-most cell = gene name) and long (one row per sample with Sample / Group / Expression). Each table has Copy-to-Clipboard and Download-CSV controls. The previous R-code generator `generate_boxplot_data_code()` is retained for backward compatibility and export to a runnable R script.
- Add "Export Boxplot Data Code" button to the Quadrant module (panel 4 "Full Expression Distribution"). Generates a self-contained R script that extracts the per-sample expression values (Sample / Group / Expression) of the currently selected gene so users can reproduce the boxplot in R, GraphPad, or Excel. New function `generate_boxplot_data_code()` is exported.
- Fix BiocCheck warnings: replace sapply() with vapply(), 1:... with seq_len/seq_along
- Fix R CMD check NOTEs: enrichit import, fc global variable, License stub
- Add @return documentation to man pages
- Fix vignette chunk labels
- Replace cat()/print() with message() in inspect functions
- Add comments to set.seed() and suppressWarnings() usage
- Update installation method to use pak
- Bump R version dependency to 4.6.0
- Add LICENSE.md for GitHub display

# GSEAlens 0.99.0

- Pathway Network module addition
- Hubgene Network module addition
- AI prompt module addition
- Script description standardization


# GSEAlens 0.0.9

- Gene expression and plot code.


# GSEAlens 0.0.8

- Repair data source of Heatmap in Limma-voom


# GSEAlens 0.0.7

- Translate to EN.


# GSEAlens 0.0.6.9000

- Update to the new parallel execution logic, and enable `progressr` as the progress bar.
- Add external data to monitor system memory and CPU status and monitor the running status of the GSEAlens enrichment function.
- Update on performance testing and monitoring when running GSEA detection.


# GSEAlens 0.0.6

* Rectify the pairing issues between Ensembl IDs and SYMBOLs in the expression matrix.

# GSEAlens 0.0.5

## 功能完善

* 已经可以使用的版本
* 交互与识别功能已经准备完全

## 待开发功能

* 一键输出选定图像的R代码功能。

# GSEAlens 0.0.4

## 多图交互

* 加入多种额外添加图像，首次完成同页面多组结果交互。

# GSEAlens 0.0.3

## 双接口

* 可接受deseq2数据制图

# GSEAlens 0.0.2

## 初版可交互

* 对接limma-voom，可视化功能初步完成。

# GSEAlens 0.0.1

## 国科金开发版本

* 最开始的原初版本，完成了国科金标书中生信部分的图像制作
