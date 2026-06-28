# Tests for build_gsea_pathways()
# Type: function with msigdbr dependency (skipped when msigdbr absent)
# Coverage: invalid species error path (always runs); valid HS/MM paths and
#           returned structure (BUG-2 FIX: list element is "TERM2GENE" not "term2gene");
#           case-insensitive species; ALL auto_select shorthand.
# NOTE: auto_select = NULL triggers readline() and is intentionally NOT tested
#       here because it hangs non-interactive sessions (BUG-3 FIX).

test_that("build_gsea_pathways rejects invalid species with informative message", {
  expect_error(
    build_gsea_pathways(species = "INVALID"),
    regexp = "Invalid species"
  )
})

test_that("build_gsea_pathways rejects empty species string", {
  expect_error(
    build_gsea_pathways(species = ""),
    regexp = "Invalid species"
  )
})

test_that("build_gsea_pathways rejects unsupported species like 'RAT'", {
  expect_error(
    build_gsea_pathways(species = "RAT"),
    regexp = "Invalid species"
  )
})

test_that("build_gsea_pathways accepts HS (human) and returns expected structure", {
  skip_if_not_installed("msigdbr")
  result <- build_gsea_pathways(species = "HS", auto_select = c("H"))
  expect_type(result, "list")
  # BUG-2 FIX: returned list element name is uppercase 'TERM2GENE' (not 'term2gene')
  expect_true("TERM2GENE" %in% names(result))
  expect_true("meta_dict" %in% names(result))
  expect_true("species" %in% names(result))
  expect_equal(result$species, "HS")
  expect_s3_class(result$TERM2GENE, "data.frame")
})

test_that("build_gsea_pathways accepts MM (mouse)", {
  skip_if_not_installed("msigdbr")
  result <- build_gsea_pathways(species = "MM", auto_select = c("H"))
  expect_type(result, "list")
  expect_true("TERM2GENE" %in% names(result))
  expect_equal(result$species, "MM")
})

test_that("build_gsea_pathways is case-insensitive for species code", {
  skip_if_not_installed("msigdbr")
  result_lower <- build_gsea_pathways(species = "hs", auto_select = c("H"))
  expect_type(result_lower, "list")
  expect_equal(result_lower$species, "HS")
})

test_that("build_gsea_pathways accepts 'ALL' auto_select shorthand", {
  skip_if_not_installed("msigdbr")
  result <- build_gsea_pathways(species = "HS", auto_select = "ALL")
  expect_type(result, "list")
  expect_true("TERM2GENE" %in% names(result))
})
