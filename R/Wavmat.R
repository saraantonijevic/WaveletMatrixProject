#' Construct Highly Optimized Orthogonal Wavelet Transformation Matrix
#'
#' This function creates an optimized wavelet transformation matrix for high-quality filter coefficients
#' and sparse matrix operations.
#'
#' @param h Numeric vector representing the low-pass filter.
#' @param N Integer specifying the size of the matrix. Must be a power of 2.
#' @param k0 Integer specifying the depth of the wavelet transformation. Should be between 1 and log2(N).
#' @param shift Integer for shifts in the wavelet transformation (default is 2).
#' @return A sparse matrix representing the wavelet packet transformation.
#' @export
#' @importFrom Matrix rBind
#' @import ggplot2
#' @import imager
#' @import pracma
WavPackMat <- function(h, N, k0, shift = 2) {
  J <- log2(N)
  if (J != floor(J)) {
    stop("N has to be a power of 2.")
  }

  # Make QM filter G
  # Make QM filter G
  h <- as.vector(h)
  g <- rev(Conj(h) * (-1)^(1:length(h)))


  h <- c(h, rep(0, N - length(h)))
  g <- c(g, rep(0, N - length(g)))

  WP <- matrix(nrow = 0, ncol = N)
  for (k in 1:k0) {
    subW <- getsubW(k, h, g, J, N)
    WP <- rbind(WP, subW)
  }
  WP <- sqrt(1/k0) * WP

  return(WP)
}



getsubW <- function(jstep, h, g, J, N) {
  subW <- diag(2^(J-jstep))
  for (k in jstep:1) {
    hgmat <- getHGmat(k, h, g, J, N)
    subW <- rbind(subW %*% hgmat$hmat, subW %*% hgmat$gmat)
  }
  return(subW)
}

getHGmat <- function(k, h, g, J, N) {
  ubJk <- 2^(J-k)
  ubJk1 <- 2^(J-k+1)
  shift <- 2

  hmat <- matrix(0, nrow = ubJk1, ncol = ubJk)
  gmat <- matrix(0, nrow = ubJk1, ncol = ubJk)

  for (jj in 1:ubJk) {
    for (ii in 1:ubJk1) {
      modulus <- (N + ii - 2 * jj + shift) %% ubJk1
      modulus <- modulus + (modulus == 0) * ubJk1
      hmat[ii, jj] <- h[modulus]
      gmat[ii, jj] <- g[modulus]
    }
  }

  return(list(hmat = t(hmat), gmat = t(gmat)))
}
