library(testthat)
library(WaveletMatrixProject)

test_that("WavmatND: Returns correct dimensions", {
  hf <- c(0.5, -0.5)  # High-pass filter
  N <- 4 # Size of the wavelet matrix
  k <- 2 # Number of iterations
  shift <- 1 # Shift parameter

  W <- WavmatND(hf, N, k, shift)

  # Expected rows = N * k + N
  expected_rows <- N * k + N
  expect_equal(dim(W), c(expected_rows, N),
               label = "Dimensions of wavelet matrix are incorrect.")
})

test_that("Test 1: Simple Base Case", {
  hf <- c(0.5, -0.5)
  N <- 4
  k <- 2
  shift <- 1
  W1 <- WavmatND(hf, N, k, shift)
  expect_equal(dim(W1), c((k + 1) * N, N))  # Verify expected dimensions
})

test_that("Test 2: Larger Matrix with Higher Iterations", {
  hf <- c(0.5, -0.5)
  N <- 8
  k <- 3
  shift <- 0
  W2 <- WavmatND(hf, N, k, shift)
  expect_equal(dim(W2), c((k + 1) * N, N))  # Expected dimensions
})


test_that("Test 3: Edge Case with Small Matrix", {
  hf <- c(0.5, -0.5)
  N <- 2
  k <- 1
  shift <- 0
  W3 <- WavmatND(hf, N, k, shift)
  expect_equal(dim(W3), c((k + 1) * N, N))  # Expected dimensions
})

test_that("Test 4: Non-default Shift", {
  hf <- c(0.25, 0.75, -0.25)
  N <- 4
  k <- 2
  shift <- 2
  W4 <- WavmatND(hf, N, k, shift)
  expect_equal(dim(W4), c((k + 1) * N, N))  # Expected dimensions
})


test_that("Output matrix contains finite values", {
  hf <- c(0.5, -0.5)
  N <- 4
  k <- 2
  shift <- 1
  W <- WavmatND(hf, N, k, shift)
  expect_true(all(is.finite(W)))  # All elements should be finite
})

test_that("Output matrix has non-zero elements", {
  hf <- c(0.5, -0.5)
  N <- 4
  k <- 2
  shift <- 1
  W <- WavmatND(hf, N, k, shift)
  expect_true(any(W != 0))  # Check for non-zero elements
})
