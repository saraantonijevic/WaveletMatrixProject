#' Enhanced Wavelet Packet Transformation Matrix
#'
#' Constructs an optimized wavelet transformation matrix with additional refinements
#' for better image fidelity in wavelet-based transformations.
#'
#' @param h Numeric vector representing the low-pass filter.
#' @param N Integer size of the matrix (must be a power of 2).
#' @param k0 Integer depth of transformation.
#' @param shift Integer shift for the wavelet transformation (default 2).
#' @param normalize Logical, whether to normalize matrix (default TRUE).
#' @return A sparse matrix for wavelet packet transformation.
#' @importFrom Matrix sparseMatrix
#' @importFrom Matrix rBind
#' @export
WavPackMat <- function(h, N, k0, shift = 2, normalize = TRUE) {
  # Validate N
  J <- log2(N)
  if (J != floor(J)) stop("N must be a power of 2.")

  # Normalize and refine filter coefficients
  h <- h / sqrt(sum(h^2)) # Normalized low-pass filter
  g <- rev(h * (-1)^(1:length(h))) # High-pass filter

  # Initialize transformation matrix
  WP <- sparseMatrix(i = integer(0), j = integer(0), dims = c(N, N))
  for (k in 1:k0) {
    subW <- getEnhancedSubW(k, h, g, J, N, shift)
    WP <- rBind(WP, subW)
  }

  # Normalize wavelet matrix if specified
  if (normalize) WP <- WP * sqrt(1 / k0)

  return(WP)
}

# Enhanced function for constructing sub-matrices
getEnhancedSubW <- function(jstep, h, g, J, N, shift) {
  subW <- sparseMatrix(i = seq(1, 2^(J - jstep)), j = seq(1, 2^(J - jstep)), x = 1)
  for (k in jstep:1) {
    hgmat <- getEnhancedHGmat(k, h, g, J, N, shift)
    subW <- rBind(subW %*% hgmat$hmat, subW %*% hgmat$gmat)
  }
  return(subW)
}

# Enhanced H and G matrix construction
getEnhancedHGmat <- function(k, h, g, J, N, shift) {
  ubJk <- 2^(J - k)
  ubJk1 <- 2^(J - k + 1)

  # Sparse H and G matrices
  hmat <- sparseMatrix(i = integer(0), j = integer(0), dims = c(ubJk1, ubJk))
  gmat <- sparseMatrix(i = integer(0), j = integer(0), dims = c(ubJk1, ubJk))

  # Populate with adjusted modulus calculation
  for (jj in 1:ubJk) {
    for (ii in 1:ubJk1) {
      modulus <- (N + ii - 2 * jj + shift) %% ubJk1
      modulus <- ifelse(modulus == 0, ubJk1, modulus)
      hmat[ii, jj] <- h[modulus] * sqrt(2) # Enhanced scaling factor
      gmat[ii, jj] <- g[modulus] * sqrt(2)
    }
  }
  return(list(hmat = hmat, gmat = gmat))
}
