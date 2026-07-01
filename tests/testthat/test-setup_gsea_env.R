# Tests for setup_gsea_env() input validation contract

# Type: validation gateway tests (no full backend construction)

# Coverage: pathway_obj NULL/missing-TERM2GENE rejection; unsupported fit class

#           rejection; error class and message assertions.

#

# BUG-4 FIX: original sample expected intercept-design validation to fire, but

# the real implementation checks pathway_obj FIRST and would always reject

# pathway_obj = list() before reaching fit-class dispatch. These tests align

# with the actual control flow.



# Helper: a pathway_obj that passes the gateway check (has TERM2GENE)

build_valid_pathway_obj <- function() {

  list(

    TERM2GENE = data.frame(

      gs_name = c("PATHWAY_DUMMY"),

      gene_symbol = c("GENE_A"),

      stringsAsFactors = FALSE

    ),

    meta_dict = data.frame(

      Collection = c("H"),

      Subcollection = c(""),

      stringsAsFactors = FALSE

    ),

    SuperTag = "test",

    collections_used = data.frame(),

    species = "HS"

  )

}



test_that("setup_gsea_env rejects NULL pathway_obj", {

  expect_error(

    setup_gsea_env(fit = list(), pathway_obj = NULL),

    regexp = "pathway_obj"

  )

})



test_that("setup_gsea_env rejects pathway_obj missing TERM2GENE field", {

  bad_pathway <- list(not_term2gene = "wrong")

  expect_error(

    setup_gsea_env(fit = list(), pathway_obj = bad_pathway),

    regexp = "pathway_obj"

  )

})



test_that("setup_gsea_env rejects empty list pathway_obj", {

  expect_error(

    setup_gsea_env(fit = list(), pathway_obj = list()),

    regexp = "pathway_obj"

  )

})



test_that("setup_gsea_env rejects fit that is neither MArrayLM nor DESeqDataSet", {

  # Must use a *valid* pathway_obj to get past the first gateway

  # and reach the fit-class dispatch.

  expect_error(

    setup_gsea_env(

      fit = list(not_a_fit = TRUE),  # no MArrayLM / DESeqDataSet class

      pathway_obj = build_valid_pathway_obj()

    ),

    regexp = "Unsupported input type"

  )

})



test_that("setup_gsea_env rejects plain numeric fit with valid pathway_obj", {

  expect_error(

    setup_gsea_env(

      fit = 42,

      pathway_obj = build_valid_pathway_obj()

    ),

    regexp = "Unsupported input type"

  )

})

