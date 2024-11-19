library(testthat)
library(WaveletMatrixProject)

test_that("repeating() correctly repeats the vector", {
  a <- c(1, 2)
  n <- 3
  repeated_a <- repeating(a, n)

  # Check that the length of the result is correct
  expect_equal(length(repeated_a), length(a) * n)

  # Check the specific content of the result
  expect_equal(repeated_a, c(1, 2, 1, 2, 1, 2))
})

test_that("repeating() handles an empty vector correctly", {
  repeated_a <- repeating(numeric(0), 5)

  # Expect an empty result if the input vector is empty
  expect_equal(repeated_a, numeric(0))
})
