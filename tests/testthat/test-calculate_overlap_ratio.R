# Tests for calculate_overlap_ratio()

# Type: pure function (no dependencies)

# Coverage: ORA / leading modes; identical / disjoint / partial overlap;

#           case-insensitivity; dedupe behavior

#

# Naming convention for fixtures:

#   GENE_A, GENE_B, ... are dummy gene symbols (purely syntactic placeholders)

#   They stand in for real symbols like TP53, MYC, GAPDH but carry no biology.



test_that("calculate_overlap_ratio returns 1 for identical sets in ora mode", {

  pathway <- c("GENE_A", "GENE_B", "GENE_C")

  de      <- c("GENE_A", "GENE_B", "GENE_C")

  expect_equal(calculate_overlap_ratio(pathway, de, ratio_mode = "ora"), 1)

})



test_that("calculate_overlap_ratio returns 0 for disjoint sets in ora mode", {

  pathway <- c("GENE_A", "GENE_B")

  de      <- c("GENE_C", "GENE_D")

  expect_equal(calculate_overlap_ratio(pathway, de, ratio_mode = "ora"), 0)

})



test_that("calculate_overlap_ratio handles partial overlap in ora mode", {

  pathway <- c("GENE_A", "GENE_B", "GENE_C", "GENE_D")

  de      <- c("GENE_C", "GENE_D", "GENE_E", "GENE_F")

  # overlap = {GENE_C, GENE_D} = 2; pathway = 4 -> 0.5

  expect_equal(calculate_overlap_ratio(pathway, de, ratio_mode = "ora"), 0.5)

})



test_that("calculate_overlap_ratio uses DE size as denominator in leading mode", {

  pathway <- c("GENE_A", "GENE_B", "GENE_C", "GENE_D")  # 4 genes

  de      <- c("GENE_C", "GENE_D")                      # 2 genes

  # overlap = {GENE_C, GENE_D} = 2

  # ora mode:     2/4 = 0.5 (pathway size)

  # leading mode: 2/2 = 1.0 (DE size)

  expect_equal(calculate_overlap_ratio(pathway, de, ratio_mode = "ora"), 0.5)

  expect_equal(calculate_overlap_ratio(pathway, de, ratio_mode = "leading"), 1)

})



test_that("calculate_overlap_ratio is case-insensitive", {

  pathway <- c("GENE_A", "GENE_B", "GENE_C")

  de      <- c("gene_a", "gene_b", "gene_c")  # lowercase input

  expect_equal(calculate_overlap_ratio(pathway, de), 1)

})



test_that("calculate_overlap_ratio does NOT deduplicate pathway genes", {

  # Implementation: intersect(toupper(.)) without unique()

  # pathway = c("GENE_A", "GENE_A", "GENE_B", "GENE_C") -> toupper keeps duplicates

  # intersect returns {GENE_A, GENE_B, GENE_C} (set semantics)

  # length(overlap) = 3; length(pathway_genes) = 4 -> 3/4

  pathway <- c("GENE_A", "GENE_A", "GENE_B", "GENE_C")

  de      <- c("GENE_A", "GENE_B", "GENE_C")

  expect_equal(calculate_overlap_ratio(pathway, de), 3 / 4)

})



test_that("calculate_overlap_ratio rejects invalid ratio_mode via match.arg", {

  expect_error(

    calculate_overlap_ratio(c("GENE_A"), c("GENE_A"), ratio_mode = "invalid"),

    class = "simpleError"

  )

})

