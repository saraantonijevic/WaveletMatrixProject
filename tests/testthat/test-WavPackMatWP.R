library(testthat)
library(WaveletMatrixProject)

test_that("WavPackMatWP: Handles minimal input without error", {
  h <- c(0.5, 0.5)
  N <- 2
  k0 <- 1
  WP <- WavPackMatWP(h, N, k0)
  expect_true(is.matrix(WP))
  expect_equal(ncol(WP), N)
})


test_that("WavPackMatWP: Handles larger input sizes", {
  h <- c(0.5, 0.5)
  N <- 8
  k0 <- 2
  WP <- WavPackMatWP(h, N, k0)
  expect_true(is.matrix(WP))
  expect_equal(ncol(WP), N)
})

test_that("WavPackMatWP: Produces orthogonal matrix", {
  h <- c(0.5, 0.5)  # Example filter (e.g., Haar wavelet)
  N <- 4            # Small size for simplicity
  k0 <- 1           # Depth of the transformation

  # Generate the wavelet packet transformation matrix
  WP <- WavPackMatWP(h, N, k0)

  # Compute WP * t(WP)
  identityCheck <- WP %*% t(WP)

  # Define tolerance for numerical accuracy
  tolerance <- 1e-10

  # Check diagonal elements for consistency
  diag_values <- diag(identityCheck)
  expect_true(all(abs(diag_values - diag_values[1]) < tolerance),
              "Diagonal elements are not consistent.")

  # Check off-diagonal elements are near zero
  off_diag_values <- identityCheck - diag(diag_values)
  expect_true(all(abs(off_diag_values) < tolerance),
              "Off-diagonal elements are not near zero.")
})

test_that("WavPackMatWP: Produces orthogonal matrix", {
  h <- c(0.5, 0.5)  # Haar wavelet filter
  N <- 8            # Matrix size
  k0 <- 2           # Depth of transformation
  WP <- WavPackMatWP(h, N, k0)

  # Compute t(WP) %*% WP
  orthogonality_check <- t(WP) %*% WP

  # Tolerance for numerical precision
  tolerance <- 1e-10

  # Check that off-diagonal elements are close to zero
  off_diag <- orthogonality_check - diag(diag(orthogonality_check))
  expect_true(all(abs(off_diag) < tolerance),
              "Off-diagonal elements are not near zero.")

  # Check that diagonal elements are consistent
  diag_elements <- diag(orthogonality_check)
  expect_true(all(abs(diag_elements - diag_elements[1]) < tolerance),
              "Diagonal elements are not consistent.")
})




test_that("WavPackMatWP: Validates matrix dimensions", {
  h <- c(0.5, 0.5)  # Haar wavelet filter
  N <- 8  # Matrix size
  k0 <- 2   # Depth of transformation

  # Generate the wavelet packet transformation matrix
  WP <- WavPackMatWP(h, N, k0)

  # Compute expected dimensions
  expected_rows <- N * k0  # Rows are N multiplied by k0
  expected_cols <- N       # Columns are always N

  # Validating dimensions
  expect_true(nrow(WP) == expected_rows,
              paste("Number of rows is incorrect. Expected:", expected_rows, "Actual:", nrow(WP)))
  expect_true(ncol(WP) == expected_cols,
              paste("Number of columns is incorrect. Expected:", expected_cols, "Actual:", ncol(WP)))

  # Check if WP is a matrix
  expect_true(is.matrix(WP), "WP is not a matrix.")
})

