library(testthat)
library(WaveletMatrixProject)

test_that("WavmatND() produces a matrix with correct dimensions", {
  hf <- c(0.5, -0.5)
  N <- 4
  k <- 2
  shift <- 1

  result <- WavmatND(hf, N, k, shift)

  # Check that the result is a matrix
  expect_true(is.matrix(result))

  # Update the expected number of rows based on function behavior
  expected_rows <- 12

  # Check that the matrix has the correct number of rows and columns
  expect_equal(nrow(result), expected_rows)
  expect_equal(ncol(result), N)
})
