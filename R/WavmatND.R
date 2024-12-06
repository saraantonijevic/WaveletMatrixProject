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
#' # Example: Constructing a Wavelet Matrix
#' hf <- c(0.5, -0.5) # Define a high-pass filter
#' N <- 8             # Size of the wavelet matrix
#' k <- 3             # Number of iterations
#' shift <- 2         # Apply a shift
#' W <- WavmatND(hf, N, k, shift) # Construct the wavelet matrix
#' print(W)           # Display the resulting wavelet matrix
#' @export
#' @importFrom Matrix Diagonal
#' @import ggplot2
#' @import imager
#' @import pracma
#' @import ggplot2
#' @import gridExtra

WavmatND <- function(hf, N, k, shift) {
  # Generate the conjugate high-pass filter
  gf <- rev(Conj(hf) * (-1)^(1:length(hf)))
  W <- matrix(nrow = 0, ncol = N)
  hmatold <- diag(N)
  h <- c(hf, rep(0, N))
  g <- c(gf, rep(0, N))

  # Iterative process to construct the wavelet matrix
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
  if (length(filt) == 0) {
    return(numeric(0))
  }
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
