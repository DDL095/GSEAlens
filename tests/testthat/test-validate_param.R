# Tests for validate_param()
# Type: pure function (no dependencies)
# Coverage: NULL/NA/non-numeric fallback; numeric string coercion;
#           min/max boundary; integer return type; vector input; NULL max_val

test_that("validate_param returns default for NULL", {
  expect_equal(validate_param(NULL, default = 5), 5)
})

test_that("validate_param returns default for NA", {
  expect_equal(validate_param(NA, default = 5), 5)
})

test_that("validate_param returns default for non-numeric string", {
  expect_equal(validate_param("abc", default = 5), 5)
})

test_that("validate_param coerces numeric string to integer", {
  expect_equal(validate_param("4", default = 1, min_val = 1, max_val = 8), 4L)
})

test_that("validate_param returns default when value below min_val", {
  expect_equal(validate_param(0, default = 2, min_val = 1), 2)
})

test_that("validate_param caps to max_val when value above", {
  expect_equal(validate_param(100, default = 2, min_val = 1, max_val = 8), 8)
})

test_that("validate_param returns integer type for valid numeric input", {
  result <- validate_param(4, default = 2, min_val = 1, max_val = 8)
  expect_type(result, "integer")
})

test_that("validate_param accepts value exactly at min_val boundary", {
  expect_equal(validate_param(1, default = 2, min_val = 1), 1L)
})

test_that("validate_param accepts value exactly at max_val boundary", {
  expect_equal(validate_param(8, default = 2, min_val = 1, max_val = 8), 8L)
})

test_that("validate_param honors NULL max_val (no upper bound)", {
  expect_equal(validate_param(1000, default = 2, min_val = 1, max_val = NULL), 1000L)
})

test_that("validate_param uses first element of vector input", {
  expect_equal(validate_param(c(3, 99, 99), default = 2, min_val = 1, max_val = 8), 3L)
})
