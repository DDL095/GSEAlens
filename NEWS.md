<!-- NEWS.md is maintained by https://cynkra.github.io/fledge, do not edit -->

# GSEAlens 0.99.8

Follow-up hardening to the scoped RNG work in 0.99.7; no user-facing behavior change.

- Fix a stale comment in `prepare_hubgene_nodes()` (`R/utils_hubgene.R`) that claimed `default 123` for the `seed` parameter while the function signature is `seed = 42`. The comment now reads `default 42`, matching the documented `@param seed` default and the Shiny UI `numericInput(vis_seed, value = 42)`.
- Add a defensive NULL fallback for `input$seed` in the pathway-relation rendering reactive (`R/13_shiny_mod_pathway_relation.R`). During fast reactive invalidation the Shiny `numericInput` can transiently be `NULL` before its value settles; the renderer now falls back to `42L` so that `withr::with_seed()` always receives an integer, preventing a spurious `"default" layout calculation failed` message. This complements the analogous guard already present in `R/16_shiny_mod_hubgene_vis.R`.

# GSEAlens 0.99.7

Bioconductor review Phase 3 completion: scoped RNG, scope-safe mutation, defensive future cleanup.

- Replace all 4 `set.seed()` calls with `withr::local_seed()` / `withr::with_seed()` so that stochastic igraph layouts (Fruchterman-Reingold) remain reproducible without polluting the caller's global RNG state. Locations: `R/08_shiny_mod_data_prep.R` (permutation sampling), `R/13_shiny_mod_pathway_relation.R` (network layout, 2 sites), `R/utils_hubgene.R` (hubgene layout). The seed parameters remain exposed to users (function `seed` argument, Shiny `input$seed`). `withr` is now in `Suggests`.
- Remove all 3 `<<-` super-assignments. In `R/08_shiny_mod_data_prep.R` the recursive permutation builder now uses an explicit `new.env()` accumulator. In `R/13_shiny_mod_pathway_relation.R` the `edge_list` and `node_df` module-scope variables become `shiny::reactiveVal()` objects: `output$plot_network` writes via `edge_list_rv(...)` / `node_df_rv(...)`, and `observeEvent(input$show_edge_detail)` reads via `edge_list_rv()`. This fixes a latent race condition where the observer could read `NULL` or a stale edge list during async re-rendering; it now degrades gracefully with a "Network is still being computed" notification.
- Harden the `future` parallel backend lifecycle in `batch_calc_gsea()`. `options(future.globals.maxSize = ...)` and `future::plan(future::multisession, ...)` are now wrapped in `on.exit()` so the user's original options and sequential plan are restored even if `future.apply::future_lapply()` errors out. The trailing `future::plan(future::sequential)` line is removed (now handled by `on.exit`). `future` + `future.apply` are retained (not switched to `BiocParallel`) because GSEAlens routinely serializes > 2 GB of globals (DE table + gene set dictionary + metadata dictionary) on Windows, where `BiocParallel::SnowParam` PSOCK serialization was measured to be 3-5x slower than `future::multisession` for this workload; see `Phase3 专项深度分析` report for the benchmark rationale.

# GSEAlens 0.99.6

Bioconductor review Phase 3 partial: code quality optimizations addressing reviewer comments on pipe operator, signal conditions, and `suppressWarnings`.

- Replace all magrittr `%>%` pipe operators with native R pipe `|>` across the entire R/ directory (~90 occurrences in 10 files). Remove the `magrittr` import from DESCRIPTION, drop `export("%>%")` and `importFrom(magrittr, "%>%")` from NAMESPACE, and delete `R/utils-pipe.R` and `man/pipe.Rd`. This eliminates an unnecessary dependency and aligns with modern R (>= 4.1) best practices.
- Restructure error and warning signals to use `rlang::abort()` / `rlang::warn()` with structured condition classes throughout `R/00_class_validations.R` (9 signals) and `R/utils-accessors.R` (30 signals). Each signal now carries a `.class` suffix (e.g. `"GSEAlens_limma_intercept_detected"`, `"GSEAlens_contrast_not_found"`) for programmatic handling via `tryCatch()`. Removed redundant `[SampleMeta]`, `[GeneDetector]`, `[SymbolMap]`, `[Limma Warning]`, `[DESeq2 Warning]` prefixes from warning text; these categories are now expressed through the condition class instead.
- Refactor `validate_param()` in `R/utils-core-genes.R` to avoid unnecessary `suppressWarnings()` when input is already numeric. The function now checks `is.numeric(value)` first and only falls back to `suppressWarnings(as.numeric(...))` for string inputs, removing the silent-coercion warning suppression for the common numeric case.
- Vignette reshaping (continued from 0.99.5): split English/Chinese vignettes into a main `GSEAlens.Rmd` / `GSEAlens-vignette-zh.Rmd` plus supplementary preprocessing vignettes (`GSEAlens-preprocessing.Rmd` / `GSEAlens-preprocessing-zh.Rmd`). Main vignettes now focus on GSEAlens core functionality and add a detailed Shiny app exploration section. All vignettes switched to `BiocStyle::html_document` output, `@`-slot access replaced with `SummarizedExperiment` accessors, and `batch_calc_gsea()` examples write to `tempdir()`.

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
