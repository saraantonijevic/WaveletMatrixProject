library(testthat)
library(WaveletMatrixProject)
test_that("WavPackMatWP() produces a matrix with correct dimensions", {
  h <- c(0.5, 0.5)
  N <- 8
  k0 <- 2
  shift <- 2

  result <- WavPackMatWP(h, N, k0, shift)

  # Check that the result is a matrix
  expect_true(is.matrix(result))

  # Check the number of rows (adjust if needed)
  expected_rows <- 16
  expect_equal(nrow(result), expected_rows)
  expect_equal(ncol(result), N)
})
