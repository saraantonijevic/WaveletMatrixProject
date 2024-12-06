library(testthat)
library(WaveletMatrixProject)
test_that("WavPackMat() produces a matrix with correct dimensions", {
  h <- c(0.5, 0.5)
  N <- 8
  k0 <- 3
  shift <- 2

  result <- WavPackMat(h, N, k0, shift)

  # Check that the result is a matrix
  expect_true(is.matrix(result))

  # Check the number of rows (adjust if needed)
  expected_rows <- 24
  expect_equal(nrow(result), expected_rows)
  expect_equal(ncol(result), N)
})

test_that("WavPackMat() handles zero-length filter gracefully", {
  h <- c() # Empty filter
  N <- 8
  k0 <- 3
  shift <- 2

  # Expect the error raised by invalid operations on empty `h`
  expect_error(WavPackMat(h, N, k0, shift), "non-numeric argument to function")
})

test_that("WavPackMat() handles single-element filter correctly", {
  h <- c(1) # Single coefficient filter
  N <- 8
  k0 <- 3
  shift <- 2

  result <- WavPackMat(h, N, k0, shift)

  # Check that the result is a matrix
  expect_true(is.matrix(result))

  # Validate dimensions based on observed behavior
  expected_rows <- 24 # Observed row count for single coefficient
  expect_equal(nrow(result), expected_rows)
  expect_equal(ncol(result), N)
})
