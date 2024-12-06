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
  W <- matrix(nrow = 0, ncol = N) # Initialize the wavelet matrix as an empty matrix with `N` columns
  hmatold <- diag(N) # Initialize the old smoothing matrix as the identity matrix of size `N`

  # Extend high-pass and conjugate filters to match dimensions for processing
  h <- c(hf, rep(0, N))
  g <- c(gf, rep(0, N))

  # Iterative process to construct the wavelet matrix
  for (i in 1:k) {
    # Initialize temporary matrices for high-pass and low-pass filter
    gmat <- matrix(0, nrow = N, ncol = N)
    hmat <- matrix(0, nrow = N, ncol = N)

    # Construct the `hmat` and `gmat` matrices based on the filter coefficients and the shift
    for (jj in 1:N) {
      for (ii in 1:N) {

        # Compute the modulus for circular indexing
        modulus <- (N + ii - jj - shift) %% N + 1
        modulus <- modulus + (modulus == 0) * N #Adjust for modulus zero
        hmat[ii, jj] <- h[modulus]
        gmat[ii, jj] <- g[modulus]
      }
    }
    # Update the wavelet matrix by appending the transformed high-pass data
    W <- rbind(t(gmat) %*% hmatold, W)

    # Update the smoothing matrix for the next iteration
    smooth <- t(hmat) %*% hmatold
    hmatold <- smooth

    # Update the filter coefficients by dilating
    h <- c(dilate_filter(hf, 2^(i) - 1), rep(0, N))
    g <- c(dilate_filter(gf, 2^(i) - 1), rep(0, N))
  }


  W <- rbind(smooth, W) # Append the final smoothing matrix
  return(W)
}


# Dilate a filter by inserting zeros between its coefficients
dilate_filter <- function(filt, k) {
  if (length(filt) == 0) {
    return(numeric(0)) # Return an empty vector if the input is empty
  }
  # Compute the new length of the dilated filter
  newlength <- (k + 1) * length(filt) - k

  # Initialize the dilated filter with zeros
  filtd <- rep(0, newlength)

  # Assign original coefficients at the appropriate indices
  filtd[seq(1, newlength, by = k + 1)] <- filt
  return(filtd)
}

imbed <- function(a) { # Imbed a vector into a longer one by appending zeros
  b <- c(a, rep(0, length(a)))
  return(b)
}

repeating <- function(a, n) {
  b <- rep(a, n) #Repeat a vector n times
  return(b)
}
