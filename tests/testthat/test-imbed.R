library(testthat)
library(WaveletMatrixProject)

test_that("imbed() doubles the vector length", {
  a <- c(1, 2, 3)
  imbedded_a <- imbed(a)

  # Check that the length of the result is double the original length
  expect_equal(length(imbedded_a), 2 * length(a))

  # Check the specific content of the result
  expect_equal(imbedded_a, c(1, 2, 3, 0, 0, 0))
})

test_that("imbed() handles empty input correctly", {
  imbedded_a <- imbed(numeric(0))

  # Expect an empty result with double length (still empty)
  expect_equal(imbedded_a, numeric(0))
})
