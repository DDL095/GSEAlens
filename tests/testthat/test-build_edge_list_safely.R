# Tests for build_edge_list_safely()
# Type: pure function (no dependencies)
# Coverage: NULL/empty/single inputs; threshold gating; three similarity
#           coefficients (Jaccard / Overlap / Dice); case-insensitivity;
#           NA gene handling (BUG-1 FIX: skipped, not silently dropped);
#           invalid min_shared_genes fallback; multi-pathway topology
# Also serves as regression guard for ISSUE-R2 (plyr -> dplyr migration).
#
# Naming convention for fixtures:
#   PATHWAY_A, PATHWAY_B, PATHWAY_C : dummy pathway identifiers
#   GENE_A, GENE_B, GENE_C, ...      : dummy gene symbols belonging to pathways
#   The exact letters carry no biology; they only encode "which gene is shared
#   with which other pathway" so that expected Jaccard / Overlap / Dice values
#   are easy to compute by hand from the test comments.

test_that("build_edge_list_safely returns NULL for NULL input", {
  expect_null(build_edge_list_safely(NULL))
})

test_that("build_edge_list_safely returns NULL for empty list", {
  expect_null(build_edge_list_safely(list()))
})

test_that("build_edge_list_safely returns NULL for single pathway", {
  cgl <- list(PATHWAY_A = c("GENE_A", "GENE_B", "GENE_C"))
  expect_null(build_edge_list_safely(cgl))
})

test_that("build_edge_list_safely returns NULL when no pair meets threshold", {
  cgl <- list(
    PATHWAY_A = c("GENE_A", "GENE_B"),
    PATHWAY_B = c("GENE_X", "GENE_Y")  # no overlap with PATHWAY_A
  )
  edges <- build_edge_list_safely(cgl, min_shared_genes = 2)
  expect_null(edges)
})

test_that("build_edge_list_safely forms edge when shared >= threshold", {
  cgl <- list(
    PATHWAY_A = c("GENE_A", "GENE_B", "GENE_C"),
    PATHWAY_B = c("GENE_B", "GENE_C", "GENE_D")
  )
  edges <- build_edge_list_safely(cgl, min_shared_genes = 2)
  expect_s3_class(edges, "data.frame")
  expect_equal(nrow(edges), 1L)
  # shared = {GENE_B, GENE_C} = 2, union = {GENE_A, GENE_B, GENE_C, GENE_D} = 4
  # Jaccard = 2/4 = 0.5
  expect_equal(edges$shared[1], 2L)
  expect_equal(edges$weight[1], 0.5)
  # Overlap Coefficient = shared / min(|A|, |B|) = 2 / 3
  expect_equal(edges$overlap_coef[1], 2 / 3)
  # Dice = 2*shared / (|A| + |B|) = 4 / 6
  expect_equal(edges$dice_coef[1], 4 / 6)
})

test_that("build_edge_list_safely is case-insensitive", {
  cgl <- list(
    PATHWAY_A = c("GENE_A", "GENE_B"),
    PATHWAY_B = c("gene_a", "gene_b")  # lowercase input
  )
  edges <- build_edge_list_safely(cgl, min_shared_genes = 2)
  expect_equal(nrow(edges), 1L)
  expect_equal(edges$shared[1], 2L)
})

# ---- BUG-1 FIX: NA genes cause the entire pair to be SKIPPED, not silently dropped ----
# Implementation: `if (any(is.na(genes1)) || any(is.na(genes2))) next`
# Therefore a 2-pathway input where one contains NA -> no edges -> NULL.
test_that("build_edge_list_safely skips pairs containing any NA gene", {
  cgl <- list(
    PATHWAY_A = c("GENE_A", NA, "GENE_B"),   # contains NA -> any(is.na(.)) is TRUE
    PATHWAY_B = c("GENE_A", "GENE_B", "GENE_C")
  )
  expect_null(build_edge_list_safely(cgl, min_shared_genes = 2))
})

test_that("build_edge_list_safely still builds edges from clean pairs when others have NA", {
  cgl <- list(
    PATHWAY_A = c("GENE_A", "GENE_B", "GENE_C"),     # clean
    PATHWAY_B = c("GENE_A", "GENE_B", "GENE_D"),     # clean, shares {GENE_A, GENE_B} with A
    PATHWAY_C = c("GENE_A", NA, "GENE_B")            # contains NA -> skipped in all pairs involving C
  )
  edges <- build_edge_list_safely(cgl, min_shared_genes = 2)
  expect_s3_class(edges, "data.frame")
  expect_equal(nrow(edges), 1L)   # only PATHWAY_A - PATHWAY_B edge survives
})

test_that("build_edge_list_safely falls back to default when min_shared_genes invalid", {
  cgl <- list(
    PATHWAY_A = c("GENE_A", "GENE_B", "GENE_C"),
    PATHWAY_B = c("GENE_A", "GENE_B", "GENE_D")
  )
  # min_shared_genes = 0 should fallback to default 2
  edges <- build_edge_list_safely(cgl, min_shared_genes = 0)
  expect_s3_class(edges, "data.frame")
  expect_equal(nrow(edges), 1L)
})

test_that("build_edge_list_safely builds 3 edges for 3 fully-connected pathways", {
  cgl <- list(
    PATHWAY_A = c("GENE_A", "GENE_B", "GENE_C", "GENE_D"),
    PATHWAY_B = c("GENE_A", "GENE_B", "GENE_C", "GENE_E"),
    PATHWAY_C = c("GENE_A", "GENE_B", "GENE_C", "GENE_F")
  )
  edges <- build_edge_list_safely(cgl, min_shared_genes = 2)
  expect_equal(nrow(edges), 3L)  # C(3,2) = 3 pairs
})
