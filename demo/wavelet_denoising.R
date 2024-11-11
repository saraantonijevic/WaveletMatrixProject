#' Demo for Wavelet Denoising using the WaveletMatrixProject package
#'
#' This demo showcases how to denoise a synthetic "bumps" signal using a wavelet
#' packet transformation matrix constructed by WavPackMat.
#' Run this demo using: demo("wavelet_denoising", package = "WaveletMatrixProject")

# Load necessary libraries
library(Matrix)
library(ggplot2)
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
ggplot(data.frame(t = t, sig = sig), aes(x = t, y = sig)) +
  geom_line(color = "green", size = 1) +
  theme_minimal(base_size = 16) +
  ggtitle("Original Bumps Signal")



# (ii) Add noise to the signal
set.seed(1)
signoi <- sig + 1/sqrt(SNR) * rnorm(N)



# (iii) Plot the noisy signal
ggplot(data.frame(t = t, signoi = signoi, sig = sig), aes(x = t)) +
  geom_point(aes(y = signoi), color = "red", size = 2) +
  geom_line(aes(y = sig), color = "green", size = 1) +
  theme_minimal(base_size = 16) +
  ggtitle("Noisy Bumps Signal")



# (iv) Create the wavelet transformation matrix using WavPackMat
filt <- c(-0.07576571478934, -0.02963552764595,
          0.49761866763246, 0.80373875180522,
          0.29785779560554, -0.09921954357694,
          -0.01260396726226, 0.03222310060407)

WP <- WavPackMat(filt, N, k0 = 6)



# (v) Transform the signal using the wavelet packet matrix
sw <- as.vector(WP %*% signoi)

# Plot the wavelet coefficients
ggplot(data.frame(index = 1:N, sw = sw), aes(x = index, y = sw)) +
  geom_line(size = 1) +
  theme_minimal(base_size = 16) +
  ggtitle("Wavelet Coefficients")



# (vi) Threshold the small coefficients
finest <- sw[(N/2+1):N]
sigmahat <- sd(finest)
lambda <- sqrt(2 * log(N)) * sigmahat

# Apply soft thresholding
swt <- sign(sw) * pmax(abs(sw) - lambda, 0)

# Plot the thresholded coefficients
ggplot(data.frame(index = 1:N, swt = swt), aes(x = index, y = swt)) +
  geom_line(size = 1) +
  theme_minimal(base_size = 16) +
  ggtitle("Thresholded Wavelet Coefficients")


