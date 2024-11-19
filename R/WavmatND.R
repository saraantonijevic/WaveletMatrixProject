#' Wavelet Matrix Construction and Utility Functions
#'
#' This set of functions includes tools for constructing wavelet matrices, dilating filters, embedding vectors, and repeating vectors for use in signal processing and wavelet transformations.
#'
#' @param hf A vector representing the high-pass filter, used as input for constructing the wavelet matrix.
#' @param N An integer specifying the size of the wavelet matrix.
#' @param k An integer indicating the number of iterations for constructing the wavelet matrix.
#' @param shift An integer representing the shift to be applied in the wavelet matrix construction.
#' @param filt A vector representing the filter to be dilated in the \code{dilate_filter} function.
#' @param a A vector to be embedded or repeated, used in the \code{imbed} and \code{repeating} functions.
#' @param n An integer specifying the number of times to repeat the vector in the \code{repeating} function.
#'
#' @details
#' \strong{WavmatND:} Constructs a wavelet matrix using a high-pass filter, specified parameters, and an iterative process to build matrices with modulated entries based on circular indexing. It returns a complete wavelet matrix that can be used in signal processing.
#'
#' \strong{dilate_filter:} Inserts zeros between the elements of the input filter to dilate it, extending its length as required by wavelet transformations.
#'
#' \strong{imbed:} Embeds a vector by appending zeros to double its length.
#'
#' \strong{repeating:} Repeats the input vector a specified number of times.
#'
#' @return
#' \itemize{
#'   \item \code{WavmatND}: A matrix representing the wavelet transformation.
#'   \item \code{dilate_filter}: A vector of the dilated filter.
#'   \item \code{imbed}: A vector with zeros appended to double its length.
#'   \item \code{repeating}: A vector with repeated elements.
#' }
#'
#' @examples
#' hf <- c(0.5, -0.5)
#' N <- 4
#' k <- 2
#' shift <- 1
#' W <- WavmatND(hf, N, k, shift)
#'
#' filt <- c(0.5, -0.5)
#' dilated_filt <- dilate_filter(filt, 2)
#'
#' a <- c(1, 2, 3)
#' imbedded_a <- imbed(a)
#'
#' a <- c(1, 2)
#' repeated_a <- repeating(a, 3)
#'
#' @export
#' @importFrom Matrix rBind
#' @import ggplot2
#' @import imager
#' @import pracma
WavmatND <- function(hf, N, k, shift) {
  gf <- rev(Conj(hf) * (-1)^(1:length(hf)))
  W <- matrix(nrow = 0, ncol = N)
  hmatold <- diag(N)
  h <- c(hf, rep(0, N))
  g <- c(gf, rep(0, N))

  for (i in 1:k) {
    gmat <- matrix(0, nrow = N, ncol = N)
    hmat <- matrix(0, nrow = N, ncol = N)

    for (jj in 1:N) {
      for (ii in 1:N) {
        modulus <- (N + ii - jj - shift) %% N + 1
        modulus <- modulus + (modulus == 0) * N
        hmat[ii, jj] <- h[modulus]
        gmat[ii, jj] <- g[modulus]
      }
    }

    W <- rbind(t(gmat) %*% hmatold, W)
    smooth <- t(hmat) %*% hmatold
    hmatold <- smooth
    h <- c(dilate_filter(hf, 2^(i) - 1), rep(0, N))
    g <- c(dilate_filter(gf, 2^(i) - 1), rep(0, N))
  }

  W <- rbind(smooth, W)
  return(W)
}

dilate_filter <- function(filt, k) {
  newlength <- (k + 1) * length(filt) - k
  filtd <- rep(0, newlength)
  filtd[seq(1, newlength, by = k + 1)] <- filt
  return(filtd)
}

imbed <- function(a) {
  b <- c(a, rep(0, length(a)))
  return(b)
}

repeating <- function(a, n) {
  b <- rep(a, n)
  return(b)
}
