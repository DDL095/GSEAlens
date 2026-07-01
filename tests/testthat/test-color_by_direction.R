# Tests for color_by_direction()

# Type: pure function (no dependencies)

# Coverage: NES-based pathway coloring; log2FC-based gene coloring;

#           NA handling; row/column preservation; color_label generation

#

# Naming convention for fixtures:

#   PATHWAY_POS / PATHWAY_NEG : pathway nodes with positive / negative NES

#   PATHWAY_ZERO / PATHWAY_NA : pathway nodes with NES == 0 / NA

#   GENE_UP / GENE_DOWN        : gene nodes with positive / negative log2FC

#   GENE_ZERO / GENE_NA        : gene nodes with log2FC == 0 / NA

# These names describe the test condition (not biology).



# Helper: build a minimal node_data frame consistent with the function contract

build_test_nodes <- function() {

  data.frame(

    name = c("PATHWAY_POS", "PATHWAY_NEG", "PATHWAY_ZERO", "PATHWAY_NA",

             "GENE_UP", "GENE_DOWN", "GENE_ZERO", "GENE_NA"),

    type = c("pathway", "pathway", "pathway", "pathway",

             "gene", "gene", "gene", "gene"),

    NES    = c( 2.5,  -2.0,   0,  NA,    NA,   NA,   NA,   NA),

    log2FC = c( NA,   NA,  NA,  NA,    1.5, -1.2,   0,   NA),

    stringsAsFactors = FALSE

  )

}



test_that("color_by_direction colors pathways by NES sign in logFC mode", {

  nodes <- build_test_nodes()

  result <- color_by_direction(nodes, color_mode = "logFC")

  expect_equal(result$color[result$name == "PATHWAY_POS"],  "#E41A1C")  # NES > 0 -> red

  expect_equal(result$color[result$name == "PATHWAY_NEG"],  "#377EB8")  # NES < 0 -> blue

  expect_equal(result$color[result$name == "PATHWAY_ZERO"], "#999999")  # NES == 0 -> gray

})



test_that("color_by_direction colors genes by log2FC sign in logFC mode", {

  nodes <- build_test_nodes()

  result <- color_by_direction(nodes, color_mode = "logFC")

  expect_equal(result$color[result$name == "GENE_UP"],   "#E41A1C")  # logFC > 0 -> red

  expect_equal(result$color[result$name == "GENE_DOWN"], "#377EB8")  # logFC < 0 -> blue

  expect_equal(result$color[result$name == "GENE_ZERO"], "#999999")  # logFC == 0 -> gray

})



test_that("color_by_direction handles NA NES/log2FC as neutral gray", {

  nodes <- build_test_nodes()

  result <- color_by_direction(nodes, color_mode = "logFC")

  expect_equal(result$color[result$name == "PATHWAY_NA"], "#999999")  # pathway NA -> gray

  expect_equal(result$color[result$name == "GENE_NA"], "#999999")     # gene NA -> gray

})



test_that("color_by_direction preserves row order and adds color + color_label", {

  nodes <- build_test_nodes()

  result <- color_by_direction(nodes, color_mode = "logFC")

  expect_equal(nrow(result), nrow(nodes))

  expect_equal(nrow(nodes), 8L)

  expect_true(all(c("color", "color_label") %in% names(result)))

  # Row order preserved

  expect_equal(result$name, nodes$name)

})



test_that("color_by_direction generates color_label with group names", {

  nodes <- build_test_nodes()

  result <- color_by_direction(nodes, color_mode = "logFC",

                               left_group = "Treatment", right_group = "Control")

  expect_equal(result$color_label[result$name == "PATHWAY_POS"],

               "Up in Treatment")

  expect_equal(result$color_label[result$name == "PATHWAY_NEG"],

               "Up in Control")

})



test_that("color_by_direction uses uniform color in fallback mode", {

  nodes <- build_test_nodes()

  result <- color_by_direction(nodes, color_mode = "other")

  expect_true(all(result$color == "#999999"))

})

