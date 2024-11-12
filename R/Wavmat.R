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
  # Determine the total number of decomposition levels (J) based on N
  J <- log2(N)

  # Ensure that N is a power of 2; otherwise, the function stops with an error message
  if (J != floor(J)) {
    stop("N has to be a power of 2.")
  }


  # Prepare the high-pass filter g from the low-pass filter h
  h <- as.vector(h) # Ensure h is a vector for consistency
  g <-
    rev(Conj(h) * (-1) ^ (1:length(h))) # Create g by reversing h and applying alternating signs

  # Extend both filters h and g to match the length N by padding with zeros
  h <- c(h, rep(0, N - length(h)))
  g <- c(g, rep(0, N - length(g)))

  # Initializing the wavelet packet transformation matrix WP with appropriate dimensions
  WP <- matrix(nrow = 0, ncol = N)

  # Iterating over the decomposition levels up to k0 to build sub-matrices
  for (k in 1:k0) {
    subW <-
      getsubW(k, h, g, J, N) # Get the sub-matrix for the current level k
    WP <-
      rbind(WP, subW)#Append the sub-matrix to the overall WP matrix
  }

  # Normalize the WP matrix by the square root of k0 for proper scaling
  WP <- sqrt(1 / k0) * WP

  return(WP)
}



getsubW <- function(jstep, h, g, J, N) {
  #start with an identity matrix of size 2^(J - jstep)
  subW <- diag(2 ^ (J - jstep))

  # Apply the filters for each level from jstep down to 1
  for (k in jstep:1) {
    hgmat <-
      getHGmat(k, h, g, J, N) #get the H and G matrices for the current level

    # Multiply the current subW by the H and G matrices and combine them
    subW <- rbind(subW %*% hgmat$hmat, subW %*% hgmat$gmat)
  }
  return(subW)
}

getHGmat <- function(k, h, g, J, N) {
  #Calc the dimensions for the H and G matrices based on the level k
  ubJk <- 2 ^ (J - k)
  ubJk1 <- 2 ^ (J - k + 1)
  shift <- 2 #shift value used for the modulus operation

  # Initializing H and G matrices with zeros
  hmat <- matrix(0, nrow = ubJk1, ncol = ubJk)
  gmat <- matrix(0, nrow = ubJk1, ncol = ubJk)

  #populate the H and G matrices
  for (jj in 1:ubJk) {
    for (ii in 1:ubJk1) {
      #calc the modulus to wrap indices correctly for the matrix
      modulus <- (N + ii - 2 * jj + shift) %% ubJk1
      modulus <- modulus + (modulus == 0) * ubJk1

      # Assign the filter coefficients based on the calculated index
      hmat[ii, jj] <- h[modulus]
      gmat[ii, jj] <- g[modulus]
    }
  }

  return(list(hmat = t(hmat), gmat = t(gmat)))
}
