library(testthat)
library(WaveletMatrixProject)

test_that("dilate_filter() correctly inserts zeros", {
  filt <- c(0.5, -0.5)
  k <- 2
  dilated_filt <- dilate_filter(filt, k)

  # Check that the result length is as expected
  expected_length <- (k + 1) * length(filt) - k
  expect_equal(length(dilated_filt), expected_length)

  # Check that the specific elements are correct
  expect_equal(dilated_filt, c(0.5, 0, 0, -0.5))
})

test_that("dilate_filter() handles empty input gracefully", {
  dilated_filt <- dilate_filter(numeric(0), 2)

  # Expect an empty result
  expect_equal(dilated_filt, numeric(0))
})
