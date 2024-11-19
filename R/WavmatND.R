#' Wavelet Matrix Construction
#'
#' Constructs a wavelet matrix using a high-pass filter, specified parameters, and an iterative process to build matrices.
#'
#' @param hf A vector representing the high-pass filter.
#' @param N An integer specifying the size of the wavelet matrix.
#' @param k An integer indicating the number of iterations.
#' @param shift An integer representing the shift to be applied.
#' @return A matrix representing the wavelet transformation.
#' @examples
#' hf <- c(0.5, -0.5)
#' N <- 4
#' k <- 2
#' shift <- 1
#' W <- WavmatND(hf, N, k, shift)
#' @export

#' @importFrom Matrix Diagonal
#' @import ggplot2
#' @import imager
#' @import pracma
#' @import ggplot2
#' @import gridExtra
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


#' Dilate Filter
#'
#' Inserts zeros between the elements of the input filter to dilate it.
#'
#' @param filt A vector representing the filter to be dilated.
#' @param k An integer specifying the step size for dilation.
#' @return A vector of the dilated filter.
#' @examples
#' filt <- c(0.5, -0.5)
#' dilated_filt <- dilate_filter(filt, 2)
#' @export
dilate_filter <- function(filt, k) {
  if (length(filt) == 0) {
    return(numeric(0))
  }

  newlength <- (k + 1) * length(filt) - k
  filtd <- rep(0, newlength)
  filtd[seq(1, newlength, by = k + 1)] <- filt
  return(filtd)
}


#' Imbed Vector
#'
#' Embeds a vector by appending zeros to double its length.
#'
#' @param a A vector to be embedded.
#' @return A vector with zeros appended to double its length.
#' @examples
#' a <- c(1, 2, 3)
#' imbedded_a <- imbed(a)
#' @export
imbed <- function(a) {
  b <- c(a, rep(0, length(a)))
  return(b)
}


#' Repeat Vector
#'
#' Repeats the input vector a specified number of times.
#'
#' @param a A vector to be repeated.
#' @param n An integer specifying the number of times to repeat the vector.
#' @return A vector with repeated elements.
#' @examples
#' a <- c(1, 2)
#' repeated_a <- repeating(a, 3)
#' @export
repeating <- function(a, n) {
  b <- rep(a, n)
  return(b)
}
