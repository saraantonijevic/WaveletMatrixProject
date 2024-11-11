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

