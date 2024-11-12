#' Demo for Wavelet Denoising using the WaveletMatrixProject package
#'
#' This demo showcases how to denoise a synthetic "bumps" signal using a wavelet
#' packet transformation matrix constructed by WavPackMat.
#' Run this demo using: demo("wavelet_denoising", package = "WaveletMatrixProject")

# Load necessary libraries
library(WaveletMatrixProject)

# (i) Generate the "bumps" signal
N <- 1024
t <- seq(0, 1, length.out = N)
pos <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81)
hgt <- c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
wth <- c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005, 0.008, 0.005)


sig <- numeric(N)
for (j in 1:length(pos)) {
  sig <- sig + hgt[j] / (1 + abs((t - pos[j]) / wth[j]))^4
}



# Standardize the signal for a fixed SNR
SNR <- 7
sig <- sig * sqrt(SNR) / sd(sig)


# Plot the original signal
ggplot2::ggplot(data.frame(t = t, sig = sig), ggplot2::aes(x = t, y = sig)) +
  ggplot2::geom_line(color = "green", size = 1) +
  ggplot2::theme_minimal(base_size = 16) +
  ggplot2::ggtitle("Original Bumps Signal")



# (ii) Add noise to the signal
set.seed(1)
signoi <- sig + 1/sqrt(SNR) * rnorm(N)


# (iii) Plot the noisy signal
ggplot2::ggplot(data.frame(t = t, signoi = signoi, sig = sig), ggplot2::aes(x = t)) +
  ggplot2::geom_point(ggplot2::aes(y = signoi), color = "red", size = 2) +
  ggplot2::geom_line(ggplot2::aes(y = sig), color = "green", size = 1) +
  ggplot2::theme_minimal(base_size = 16) +
  ggplot2::ggtitle("Noisy Bumps Signal")



# (iv) Create the wavelet transformation matrix using WavPackMat
filt <- c(-0.07576571478934, -0.02963552764595,
          0.49761866763246, 0.80373875180522,
          0.29785779560554, -0.09921954357694,
          -0.01260396726226, 0.03222310060407)

WP <- WavPackMat(filt, N, k0 = 6)



# (v) Transform the signal using the wavelet packet matrix
sw <- as.vector(WP %*% signoi)



# Plot the wavelet coefficients
ggplot2::ggplot(data.frame(index = 1:N, sw = sw), ggplot2::aes(x = index, y = sw)) +
  ggplot2::geom_line(size = 1) +
  ggplot2::theme_minimal(base_size = 16) +
  ggplot2::ggtitle("Wavelet Coefficients")



# (vi) Threshold the small coefficients
finest <- sw[(N/2+1):N]
sigmahat <- sd(finest)
lambda <- sqrt(2 * log(N)) * sigmahat

# Apply soft thresholding
swt <- sign(sw) * pmax(abs(sw) - lambda, 0)
swt_length <- length(swt)

# Convert and display the thresholded wavelet coefficients as an image
swt_img <- imager::as.cimg(matrix(swt, nrow = 1, ncol = swt_length))
plot(swt_img, main = "Thresholded Wavelet Coefficients Image Representation")

# Plot the thresholded coefficients
ggplot2::ggplot(data.frame(index = 1:N, swt = swt), ggplot2::aes(x = index, y = swt)) +
  ggplot2::geom_line(size = 1) +
  ggplot2::theme_minimal(base_size = 16) +
  ggplot2::ggtitle("Thresholded Wavelet Coefficients")


# (vii) Return the signal to the time domain using the inverse transformation
a <- as.vector(t(WP) %*% swt)

# (viii) Plot the denoised signal and the noisy signal for comparison
ggplot2::ggplot(data.frame(t = t, signoi = signoi, a = a), ggplot2::aes(x = t)) +
  ggplot2::geom_point(ggplot2::aes(y = signoi), color = "red", size = 2) +
  ggplot2::geom_line(ggplot2::aes(y = a), color = "black", size = 1) +
  ggplot2::theme_minimal(base_size = 16) +
  ggplot2::ggtitle("Denoised Signal")
