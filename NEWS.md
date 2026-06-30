<!-- NEWS.md is maintained by https://cynkra.github.io/fledge, do not edit -->

# GSEAlens 0.99.12

Publication export fidelity: ggrepel labels for volcano exports, new
boxplot image export modal, modal scroll fix, and checkbox-driven
annotation toggles.

- **Volcano export labels (R/15_code_generator.R, R/10)**: `generate_volcano_code()` now accepts `selected_ids` and annotates the user-selected pathways with `ggrepel::geom_label_repel` (orange labels mirroring the interactive plotly selection). `generate_de_volcano_code()` annotates user / pathway / both genes (green / orange / purple labels) the same way. Labels use `max.overlaps = Inf` so all selected points are named.
- **Checkbox-driven annotations (R/10, R/15)**: both volcano export modals gain a checkbox (default unchecked). Pathway Volcano: "Show stats subtitle"; DE Volcano: "Show stats + direction banners". When unchecked the figure is clean (title + axes only, identical to pre-annotation behavior); when checked, subtitle shows the same statistics as the interactive plotly title and DE Volcano adds red "High in <left_group>" / blue "High in <right_group>" corner banners using the actual group names.
- **Boxplot image export (R/10, R/15)**: add an "Export Publication Plot" button to panel 4 (Full Expression Distribution) of the Quadrant module, opening a full export modal (PDF/PNG/SVG/TIFF + Live Preview + Copy R Code). New generator `generate_boxplot_image_code()` produces a self-contained ggplot2 script reproducing the interactive plotly boxplot (boxplot + jitter, group colors, optional zero baseline, optional custom group order).
- **Modal scroll fix (R/07_shiny_app.R)**: replace the single `max-height: 85vh` on `.modal-body` with a flex layout (`.modal` scrolling viewport, `.modal-content` bounded to 92vh flex column, pinned header/footer, scrolling body). Fixes the bug where a tall Live Preview pushed the modal header above the fold and the scrollbar could not reach it.
- **DotPlot export subtitle (R/15)**: dynamic subtitle now references the active Color by / Size by selection instead of the previous hardcoded string.

# GSEAlens 0.99.11

Visualization fidelity and publication export. Brings the interactive
plotly/visNetwork figures in line with `enrichplot` / `clusterProfiler`
conventions, and adds a unified "Export Publication Plot" workflow that
produces both static `ggplot2` images and reproducible R scripts.

- **DotPlot size mapping (R/13_shiny_mod_pathway_relation.R)**: replace per-render min-max normalization with a **fixed-domain + sqrt** scale (CoreCount domain `[0, 50]`, setSize domain `[0, 500]`), matching `ggplot2::scale_size_continuous(range=...)` semantics. The previous normalization pinned min/max to fixed pixel sizes every render, hiding the magnitude difference between CoreCount (typical 5-50) and setSize (typical 50-500), so switching "Size by" produced visually identical dot patterns. The new scale also makes dot sizes comparable across FDR/TopN parameter changes. DotPlot hover text now shows both dimensions, with the currently-driven one bolded, plus P-value.
- **Pathway Volcano Y axis (R/10_shiny_mod_quadrant.R)**: switch from `-log10(p.adjust)` to `-log10(pvalue)`, matching the `EnhancedVolcano` / `clusterProfiler` convention. P-value (rather than FDR) provides the sensitivity needed to expose gradient structure in marginal pathways; FDR compresses most points to the bottom. Hover text and title updated accordingly.
- **Network edge width (R/13_shiny_mod_pathway_relation.R)**: add a user-selectable **Edge Width Mode** dropdown with two options: *Weight-based* (default; edge width linearly proportional to Jaccard value, `emapplot` convention, faithful to similarity magnitude) and *Rank-based* (legacy behavior; edge width assigned by Jaccard rank for uniform visual spacing). The default Weight-based mode is the standard for publication figures; the Rank-based mode is retained for dense networks with low weight variance.
- **HubGene pathway node size (R/16_shiny_mod_hubgene_vis.R)**: add a **Pathway Node Size Encoding** dropdown with three modes: *By gene-set size* (`setSize`, default; matches `enrichplot::cnetplot` convention, sqrt-scaled), *By significance* (`-log10(FDR)`), and *Fixed size* (legacy slider-only behavior). The slider value acts as the base size; the chosen encoding scales within `[0.6x, 1.4x]` to keep visNetwork's force-directed layout stable (variance beyond ~2.3x causes layout jitter). Gene-node sizing is unaffected. Pathway hover text now also shows `setSize`.
- **Export Center (R/13, R/10, R/16, R/15)**: add a unified "Export Publication Plot" modal to the DotPlot/Network/Volcano/HubGene panels. Provides width/height/DPI/format (PDF, PNG, SVG, TIFF) controls and two actions: (1) download a static `ggplot2`-rendered image via `ggsave` (zero external dependencies; no kaleido/orca required, unlike `plotly::save_image`); (2) copy a fully reproducible R script to the clipboard via `clipr`. Four new exported code generators in `R/15_code_generator.R`: `generate_dotplot_code()`, `generate_volcano_code()`, `generate_network_code()`, `generate_hubgene_code()`. Network and HubGene static renderings use `igraph` layouts + base `ggplot2` (no `ggraph` dependency, which is not in DESCRIPTION).
- Update English and Chinese vignettes (Tab 3 and Tab 4 sections) to document the new sub-panels, encoding modes, and export workflow.

# GSEAlens 0.99.10

BiocCheck compliance: clear the ">= 80% runnable examples" ERROR, plus
provenance documentation for shipped pre-computed RDS objects.

- Convert 10 exported functions from `\dontrun{}` examples to truly runnable examples by loading shipped pre-computed RDS objects (`precomputed_gseares.rds`, `preprocessed_limma.rds`, `gsea_pathwaysets_toy.rds`) or writing demo files to `tempdir()`. Functions updated: `setup_gsea_env`, `inspect_gsea_env`, `import_gsea_capsule`, `plot_directional_gsea`, `generate_pathway_plot_code`, `build_hubgene_network`, `extract_hub_genes`, `read_addition_data`, `creat_addition_data_rdsfile`, `create_addition_template`. This raises runnable example coverage from 57% (21/37) to 84% (31/37) and clears the BiocCheck ERROR "At least 80% of man pages documenting exported objects must have runnable examples".
- Fix BiocCheck WARNING on `man/visualization.Rd`: the roxygen block for the `visualization` documentation landing page had a duplicate `@name` / `NULL` section, which suppressed the `\value{}` block and left the Rd without a value section. Consolidated to a single roxygen block with `@return NULL` so the generated Rd has a proper `\value{}`.
- Fix BiocCheck NOTE "Use accessors; don't access S4 class slots via '@' in examples": replaced `task$gsea_res@result$ID[1]` with `as.data.frame(task$gsea_res)$ID[1]` in `get_core_genes_for_pathway()` and `get_core_genes_list()` examples (`R/utils-core-genes.R`). The `clusterProfiler::gseaResult` class has an `as.data.frame()` method that exposes the `@result` slot without direct S4 slot access.
- Add `inst/scripts/make_precomputed_gseares.R`, the maintainer-only generator script for `inst/extdata/precomputed_gseares.rds`. Documents the exact inputs (`preprocessed_limma.rds` + `gsea_pathwaysets_toy.rds`), clusterProfiler/fgsea versions, and the `batch_calc_gsea()` call used to produce the shipped object, so the artefact is reproducible when upstream packages bump. Update `inst/scripts/README.md` with the new file's row in the regeneration table and the dependency order note.
- Remove stale `tests/testthat/test-run_app.R.placeholder.bak` (renamed backup of the placeholder test retired in 0.99.9). Remove misspelled `GSEAlen.Rproj` (duplicate of `GSEAlens.Rproj`); also drop its line from `.Rbuildignore`.

# GSEAlens 0.99.9

Bioconductor review Phase 4: dependency hygiene, import hardening, and the first real test suite.

- Replace `parallel::detectCores()` with `future::availableCores()` in `batch_calc_gsea()` (`R/04_calculation.R`) so that the worker count honors `future` backend configuration (e.g. `future.cores` option, `CIRCLE_NODE_TOTAL`, Slurm env vars) instead of always reading hardware cores. This removes the last direct `parallel::` call from the package and aligns with Bioconductor guidance to defer parallelization decisions to the user-selected backend (`future` / `BiocParallel`).
- Replace `plyr::rbind.fill()` with `dplyr::bind_rows()` in `build_hubgene_network()` (`R/utils_hubgene.R`). The two are semantically equivalent for the homogeneous data-frame list used here (both fill missing columns with `NA`), but `dplyr::bind_rows()` is actively maintained and avoids the superseded `plyr` dependency. Verified behavior parity with a targeted R check on heterogeneous-column inputs.
- Move `withr` from `Suggests` to `Imports` and add explicit `importFrom(withr, with_seed)` / `importFrom(withr, local_seed)` to `NAMESPACE`. `withr` is now used unconditionally in production code (RNG scoping in 4 sites since 0.99.7), so it must be a hard dependency.
- Add 23 explicit `importFrom` declarations covering `tools` (1), `withr` (2), `enrichit` (1), `igraph` (7), `visNetwork` (8), `shinyjs` (2), and `progressr` (2) to `NAMESPACE` via `R/utils-imports.R`. This eliminates all residual `pkg::fun()` calls for these packages so that `R CMD check` reports zero "undeclared package usage" notes and downstream code is robust to masked search-path changes.
- Add the first real `testthat` suite: 73 tests across 6 files (`tests/testthat/test-calculate_overlap_ratio.R`, `test-color_by_direction.R`, `test-build_edge_list_safely.R`, `test-validate_param.R`, `test-build_gsea_pathways.R`, `test-setup_gsea_env.R`). The previous `test-run_app.R` placeholder (`expect_equal(2 * 2, 4)`, unchanged since v0.99.0) is preserved as `.placeholder.bak`. The new suite covers 5 exported functions and the `setup_gsea_env()` input-validation gateway; all tests pass cleanly under `R CMD check --as-cran`. Closes Bioconductor reviewer comment #11 ("add real tests").
- Wrap the two interactive `build_gsea_pathways()` examples (those relying on `readline()`) in `\dontrun{}` and correct the MM-species example to use the M-prefixed collection codes (`c("MH", "M5:GO:BP")`) instead of the human codes. This resolves a long-standing `R CMD check` ERROR where the example called `readline()` in a non-interactive session.
- Translate 7 Chinese code comments in `R/09_shiny_mod_table.R` (JS debounce logic and Shiny input-handler semantics) to English to eliminate the `R CMD check` "non-ASCII characters in R code" WARNING.

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
