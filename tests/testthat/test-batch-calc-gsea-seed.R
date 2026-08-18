# Tests for deterministic seeding of batch_calc_gsea() (v0.99.34)
#
# Background: clusterProfiler 4.19.3-4.20.x (Bioc 3.23) silently drops the
# `seed` argument of GSEA() (it lands in `...` and is never forwarded to
# enrichit::gsea_gson), so the engine draws a fresh seed from the worker R RNG
# on every call. These tests pin the two guarantees of the fix:
#   1. same `seed` -> identical results across runs
#   2. results independent of `workers` / chunking

test_that(".derive_task_seed is deterministic and order-independent", {
  ids <- c("B_vs_A", "A_vs_B", "C_vs_A")
  expect_identical(
    .derive_task_seed(123, "A_vs_B", ids),
    .derive_task_seed(123, "A_vs_B", rev(ids))
  )
  expect_identical(
    .derive_task_seed(123, "A_vs_B", ids),
    .derive_task_seed(123, "A_vs_B", sample(ids))
  )
  # distinct tasks get distinct seeds
  seeds <- vapply(sort(ids), function(id) .derive_task_seed(50, id, ids), integer(1))
  expect_identical(length(unique(seeds)), 3L)
  # NULL / NA / non-numeric seeds hand control back to the engine
  expect_null(.derive_task_seed(NULL, "A_vs_B", ids))
  expect_null(.derive_task_seed(NA, "A_vs_B", ids))
  expect_null(.derive_task_seed("abc", "A_vs_B", ids))
})

build_seed_test_env <- function() {
  genes <- sprintf("G%03d", 1:300)
  stat <- seq(3, -3, length.out = 300)
  de_a_b <- data.frame(gene_symbol = genes, stat = stat, stringsAsFactors = FALSE)
  de_a_c <- data.frame(
    gene_symbol = genes, stat = stat * rev(seq_len(300)) / 300,
    stringsAsFactors = FALSE
  )
  # 4 gene sets of 40 genes each, deterministic membership (no RNG)
  t2g <- do.call(rbind, lapply(seq_len(4), function(i) {
    data.frame(
      gs_name = paste0("PW", i),
      gene_symbol = genes[((i - 1) * 40 + 1):(i * 40)],
      stringsAsFactors = FALSE
    )
  }))
  env <- list(
    backend_info = list(backend = "test"),
    contrast_registry = data.frame(
      contrast_id = c("A_vs_B", "A_vs_C"),
      left_group = c("A", "A"),
      right_group = c("B", "C"),
      stringsAsFactors = FALSE
    ),
    de_store = list(A_vs_B = de_a_b, A_vs_C = de_a_c),
    expr_bundle = list(),
    geneset = list(
      name = "seed_test_gs",
      term2gene = t2g,
      meta_dict = data.frame(
        ID = paste0("PW", 1:4), Collection = "H",
        stringsAsFactors = FALSE
      )
    )
  )
  class(env) <- "GseaEnv"
  env
}

test_that("batch_calc_gsea with a fixed seed is reproducible across runs", {
  env <- build_seed_test_env()
  out <- file.path(tempdir(), "seed-repro")
  r1 <- batch_calc_gsea(env, workers = 2, bidirectional = FALSE,
                        output_dir = out, force = TRUE)
  r2 <- batch_calc_gsea(env, workers = 2, bidirectional = FALSE,
                        output_dir = out, force = TRUE)
  for (task in c("A_vs_B", "A_vs_C")) {
    expect_identical(
      r1$results[[task]]$data@result$pvalue,
      r2$results[[task]]$data@result$pvalue
    )
  }
  expect_identical(r1$metadata$parameters$seed, 123)
})

test_that("batch_calc_gsea results do not depend on worker count", {
  env <- build_seed_test_env()
  out <- file.path(tempdir(), "seed-workers")
  r_seq <- batch_calc_gsea(env, workers = 1, bidirectional = FALSE,
                           output_dir = out, force = TRUE)
  r_par <- batch_calc_gsea(env, workers = 2, bidirectional = FALSE,
                           output_dir = out, force = TRUE)
  for (task in c("A_vs_B", "A_vs_C")) {
    expect_identical(
      r_seq$results[[task]]$data@result$pvalue,
      r_par$results[[task]]$data@result$pvalue
    )
  }
})
