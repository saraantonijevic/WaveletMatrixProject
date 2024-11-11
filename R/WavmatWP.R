#' Construct Wavelet Packet Transformation Matrix
#'
#' This function creates a wavelet packet transformation matrix for high-quality filter coefficients.
#' It uses orthogonal wavelet transformation to construct the matrix.
#'
#' @param h Numeric vector representing the low-pass filter.
#' @param N Integer specifying the size of the matrix. Must be a power of 2.
#' @param k0 Integer specifying the depth of the wavelet transformation. Should be between 1 and log2(N).
#' @param shift Integer for shifts in the wavelet transformation (default is 2).
#' @return A matrix representing the wavelet packet transformation.
#' @export
#' @importFrom Matrix rBind
#' @import ggplot2
#' @import imager
#' @import pracma hypot
#' @import Rmpfr mpfr sqrt
WavPackMatWP <- function(h, N, k0, shift = 2) {
  J <- log2(N)

  # Make QM filter G
  h <- as.vector(h) # Ensure h is a vector
  g <- rev(Conj(h) * (-1)^(1:length(h)))

  if (J != floor(J)) {
    stop('N has to be a power of 2.')
  }
  h <- c(h, rep(0, N)) # Extend filter H by 0's to sample by modulus
  g <- c(g, rep(0, N)) # Extend filter G by 0's to sample by modulus

  # Initialize WP with the number of columns expected from subW
  subW <- getsubWP(1, h, g, J, N)  # Use the first level to check dimensions
  WP <- matrix(nrow = 0, ncol = ncol(subW))

  for (k in 1:k0) {
    subW <- getsubWP(k, h, g, J, N)
    if (ncol(WP) != ncol(subW)) {
      stop("Column mismatch between WP and subW during rbind operation.")
    }
    WP <- rbind(WP, subW)
  }

  WP <- sqrt(1/k0) * WP

  return(WP)
}


# Helper functions for WavPackMatWP
getsubWP <- function(jstep, h, g, J, N) {
  subW <- diag(2^(J - jstep))
  for (k in jstep:1) {
    hg_mats <- getHGmatWP(k, h, g, J, N)
    hmat <- hg_mats[[1]]
    gmat <- hg_mats[[2]]
    subW <- rbind(subW %*% hmat, subW %*% gmat)
  }
  return(subW)
}



getHGmatWP <- function(k, h, g, J, N) {
  ubJk <- 2^(J - k)
  ubJk1 <- 2^(J - k + 1)
  shift <- 2

  hmat <- matrix(0, nrow = ubJk1, ncol = ubJk)
  gmat <- matrix(0, nrow = ubJk1, ncol = ubJk)

  for (jj in 1:ubJk) {
    for (ii in 1:ubJk1) {
      modulus <- (N + ii - 2 * jj + shift) %% ubJk1
      modulus <- modulus + (modulus == 0) * ubJk1
      hmat[ii, jj] <- h[modulus + 1] # R is 1-indexed
      gmat[ii, jj] <- g[modulus + 1] # R is 1-indexed
    }
  }

  hmat <- t(hmat)
  gmat <- t(gmat)

  return(list(hmat, gmat))
}
