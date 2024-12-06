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
  # Check if the dimensions of the output match expectations
  expect_equal(dim(W), c(expected_rows, N),
               label = "Dimensions of wavelet matrix are incorrect.")
})

test_that("Test 2: Simple Base Case", {
  hf <- c(0.5, -0.5) # Filter coefficients
  N <- 4 # Size of base matrix
  k <- 2 # Iterations
  shift <- 1 # Shift parameter
  W1 <- WavmatND(hf, N, k, shift)

  # Verify the dimensions match expected values
  expect_equal(dim(W1), c((k + 1) * N, N))  # Verify expected dimensions
})

test_that("Test 3: Larger Matrix with Higher Iterations", {
  hf <- c(0.5, -0.5) # Filter coefficients
  N <- 8  # Larger base size
  k <- 3  # Increased iteration
  shift <- 0 # No shift
  W2 <- WavmatND(hf, N, k, shift)

  # Check dimensions for correctness
  expect_equal(dim(W2), c((k + 1) * N, N))  # Expected dimensions
})


test_that("Test 4: Edge Case with Small Matrix", {
  hf <- c(0.5, -0.5) # Filter coefficients
  N <- 2 # Minimum size for base matrix
  k <- 1
  shift <- 0
  W3 <- WavmatND(hf, N, k, shift)
  # Verify dimensions
  expect_equal(dim(W3), c((k + 1) * N, N))  # Expected dimensions
})

test_that("Test 5: Non-default Shift", {
  hf <- c(0.25, 0.75, -0.25) # Custom filter coefficients
  N <- 4 # Base size of the matrix
  k <- 2 # Iterations
  shift <- 2 # Non-default shift value
  W4 <- WavmatND(hf, N, k, shift)

  # Check dimensions
  expect_equal(dim(W4), c((k + 1) * N, N))  # Expected dimensions
})


test_that("Output matrix contains finite values", {
  hf <- c(0.5, -0.5) # Filter coefficients
  N <- 4 # Base size
  k <- 2 # Iterations
  shift <- 1 # Shift parameter
  W <- WavmatND(hf, N, k, shift)

  # Ensure all elements in the matrix are finite
  expect_true(all(is.finite(W)))
})

test_that("Output matrix has non-zero elements", {
  hf <- c(0.5, -0.5) # Filter coefficients
  N <- 4# Base size
  k <- 2  # Iterations
  shift <- 1
  W <- WavmatND(hf, N, k, shift)

  # Ensure at least one element in the matrix is non-zero
  expect_true(any(W != 0))
})
