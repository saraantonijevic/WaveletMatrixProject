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

  # Check that the matrix has the correct number of rows and columns
  expect_equal(nrow(result), N)
  expect_equal(ncol(result), N)
})

test_that("WavmatND() handles simple edge cases", {
  hf <- c(1, -1)
  N <- 2
  k <- 1
  shift <- 0

  result <- WavmatND(hf, N, k, shift)

  # Ensure the result is a finite matrix
  expect_true(all(is.finite(result)))
})
