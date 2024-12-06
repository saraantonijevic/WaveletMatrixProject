library(testthat)
library(WaveletMatrixProject)


test_that("WavPackMat: Returns correct dimensions", {
  h <- c(1 / sqrt(2), 1 / sqrt(2))  # Low-pass filter
  N <- 4  # Size of the wavelet matrix
  k0 <- 2  # Number of decomposition levels

  WP <- WavPackMat(h, N, k0)

  # Expected rows = k0 * N
  expected_rows <- k0 * N

  # Ensure the dimensions are correct
  expect_equal(dim(WP), c(expected_rows, N),
               label = "Dimensions of wavelet matrix are incorrect.")
})


test_that("WavPackMat: Returns correct dimensions for larger matrix", {
  h <- c(1 / sqrt(2), 1 / sqrt(2))  # Low-pass filter
  N <- 8  # Size of the wavelet matrix
  k0 <- 3  # Number of decomposition levels

  WP <- WavPackMat(h, N, k0)

  # Expected rows = k0 * N
  expected_rows <- k0 * N

  # Ensure the dimensions are correct
  expect_equal(dim(WP), c(expected_rows, N),
               label = "Dimensions of wavelet matrix are incorrect.")
})


test_that("WavPackMat: Scaled orthogonality and structure", {
  h <- c(1 / sqrt(2), 1 / sqrt(2))  # Haar wavelet coefficients
  N <- 8 # Size of the wavelet matrix
  k0 <- 3 # Number of decomposition levels

  # Generate the wavelet packet transformation matrix
  WP <- WavPackMat(h, N, k0)

  # Compute WP * t(WP)
  identityCheck <- WP %*% t(WP)

  # Define tolerance for numerical accuracy
  tolerance <- 1e-6

  # Check that diagonal elements are scaled correctly
  diag_values <- diag(identityCheck)
  expected_diag <- 1 / k0
  expect_true(all(abs(diag_values - expected_diag) < tolerance),
              "Diagonal elements are not scaled correctly.")

  # Check off-diagonal blocks within the same sub-matrix level
  for (level in 1:k0) {
    # Calculate the row range for this level
    start_row <- (level - 1) * N + 1
    end_row <- level * N

    # Extract the block and verify off-diagonal values
    block <- identityCheck[start_row:end_row, start_row:end_row]
    off_diag_block <- block - diag(diag(block))
    expect_true(all(abs(off_diag_block) < tolerance))
  }
})

test_that("WavPackMat: Energy is preserved", {
  h <- c(1 / sqrt(2), 1 / sqrt(2))  # Haar wavelet coefficients
  N <- 8 # Size of the wavelet matrix
  k0 <- 3 # Number of decomposition levels

  # Generate the wavelet packet transformation matrix
  WP <- WavPackMat(h, N, k0)

  # Create a random input vector
  set.seed(123)  # For reproducibility
  x <- rnorm(N)

  # Transform the input vector using WP
  transformed_x <- WP %*% x

  # Compute the energy in the original and transformed spaces
  original_energy <- sum(x^2)
  transformed_energy <- sum(transformed_x^2)

  # Debugging: Print energy values for verification
  cat("Original energy:", original_energy)
  cat("Transformed energy:", transformed_energy)

  # Define tolerance for numerical accuracy
  tolerance <- 1e-6

  # Check if energy is preserved
  expect_equal(transformed_energy, original_energy, tolerance = tolerance,
               info = "Energy is not preserved in the wavelet packet transformation.")
})
