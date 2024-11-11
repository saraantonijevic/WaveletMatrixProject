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
#' @import pracma
