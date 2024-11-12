#' Construct Wavelet Packet Transformation Matrix
#'
#' @param h Numeric vector representing the low-pass filter.
#' @param N Integer specifying the size of the matrix. Must be a power of 2.
#' @param k0 Integer specifying the depth of the wavelet transformation. Should be between 1 and log2(N).
#' @param shift Integer for shifts in the wavelet transformation (default is 2).
#' @return A matrix representing the wavelet packet transformation.
#' @export
#' @importFrom Matrix rBind
#' @import pracma
#' @import wavelets
WavPackMatWP <- function(h, N, k0, shift = 2) {
  # Calculate the level of decomposition (J) based on the matrix size (N)
  J <- log2(N)

  if (J != floor(J)) {
    stop('N has to be a power of 2.')
  }

  # Ensure the input filter h is a vector and construct the corresponding high-pass filter g
  h <- as.vector(h) # Ensure h is a vector for processing
  g <-
    rev(Conj(h) * (-1) ^ (1:length(h)))# Construct g as the reversed conjugate of h multiplied by alternating signs

  # Extend filters h and g by padding them with zeros up to length N
  h <- c(h, rep(0, N)) # Extend filter H by 0's to sample by modulus
  g <- c(g, rep(0, N)) # Extend filter G by 0's to sample by modulus

  # Initialize the wavelet packet transformation matrix by using the first level to determine dimensions
  subW <-
    getsubWP(1, h, g, J, N)  # Use the first level to check dimensions
  WP <- matrix(nrow = 0, ncol = ncol(subW))

  # Construct the wavelet packet matrix up to the specified depth k0
  for (k in 1:k0) {
    subW <- getsubWP(k, h, g, J, N)
    if (ncol(WP) != ncol(subW)) {
      stop("Column mismatch between WP and subW during rbind operation.")
    }
    WP <-
      rbind(WP, subW) # Append the sub-matrices to form the WP matrix
  }
  # Normalize the resulting matrix by the square root of the depth k0
  WP <- sqrt(1 / k0) * WP

  # Return the final wavelet packet transformation matrix
  return(WP)
}


getsubWP <- function(jstep, h, g, J, N) {
  # Create an identity matrix with dimensions 2^(J - jstep)
  subW <- diag(2 ^ (J - jstep))
  # Iterate backward from jstep to 1 to apply filters in each step
  for (k in jstep:1) {
    # Get the filter matrices for the current step
    hgmat <- getHGmatWP(k, h, g, J, N)

    # Apply the filters to the sub-matrix and create the next level's subW matrix
    subW <- rbind(subW %*% hgmat$hmat, subW %*% hgmat$gmat)
  }

  # Return the constructed sub-matrix for the current decomposition level
  return(subW)
}

getHGmatWP <- function(k, h, g, J, N) {
  # Defining dimensions for the filter matrices based on the current level k
  ubJk <- 2 ^ (J - k)
  ubJk1 <- 2 ^ (J - k + 1)
  shift <- 2# Shift parameter for wavelet packet decomp

  # Initializing H and G matrices with zeros
  hmat <- matrix(0, nrow = ubJk1, ncol = ubJk)
  gmat <- matrix(0, nrow = ubJk1, ncol = ubJk)

  # Populating the H and G matrices based on the given filter coefficients
  for (jj in 1:ubJk) {
    for (ii in 1:ubJk1) {
      # Calcu the modulus to handle circular boundary conditions
      modulus <- (N + ii - 2 * jj + shift) %% ubJk1
      modulus <- modulus + (modulus == 0) * ubJk1

      # Assigning the appropriate filter coefficients to the matrices
      hmat[ii, jj] <- h[modulus]
      gmat[ii, jj] <- g[modulus]
    }
  }
  # Returning the H and G matrices transposed as a list
  return(list(hmat = t(hmat), gmat = t(gmat)))
}
