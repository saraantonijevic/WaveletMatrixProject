# Load necessary libraries
library(WaveletMatrixProject)

#Generate the "bumps" signal
N <- 1024 # Number of points in the signal
t <- seq(0, 1, length.out = N) # Time vector


# Define the positions, heights, and widths of the bumps
pos <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81)
hgt <- c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
wth <- c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005, 0.008, 0.005)

# Initialize the signal and add each bump
sig <- numeric(N)
for (j in 1:length(pos)) {
  sig <- sig + hgt[j] / (1 + abs((t - pos[j]) / wth[j]))^4
}


# Standardize the signal for a fixed Signal-to-Noise Ratio (SNR)
SNR <- 7# Desired SNR
sig <- sig * sqrt(SNR) / sd(sig) # Scale signal to achieve the target SNR


# Plot the original signal
ggplot2::ggplot(data.frame(t = t, sig = sig), ggplot2::aes(x = t, y = sig)) +
  ggplot2::geom_line(color = "green", size = 1) +
  ggplot2::theme_minimal(base_size = 16) +
  ggplot2::ggtitle("Original Bumps Signal")



# Add noise to the signal
set.seed(1)# Set seed for reproducibility
signoi <- sig + 1/sqrt(SNR) * rnorm(N) # Add Gaussian noise to the signal


#Plot the noisy signal
ggplot2::ggplot(data.frame(t = t, signoi = signoi, sig = sig), ggplot2::aes(x = t)) +
  ggplot2::geom_point(ggplot2::aes(y = signoi), color = "red", size = 2) +
  ggplot2::geom_line(ggplot2::aes(y = sig), color = "green", size = 1) +
  ggplot2::theme_minimal(base_size = 16) +
  ggplot2::ggtitle("Noisy Bumps Signal")



#Create the wavelet transformation matrix using WavPackMat
filt <- c(-0.07576571478934, -0.02963552764595,
          0.49761866763246, 0.80373875180522,
          0.29785779560554, -0.09921954357694,
          -0.01260396726226, 0.03222310060407)

WP <- WavPackMat(filt, N, k0 = 6)# Generate the wavelet packet transformation matrix



#Transform the signal using the wavelet packet matrix
sw <- as.vector(WP %*% signoi)# Compute wavelet coefficients



# Plot the wavelet coefficients
ggplot2::ggplot(data.frame(index = 1:N, sw = sw), ggplot2::aes(x = index, y = sw)) +
  ggplot2::geom_line(size = 1) +
  ggplot2::theme_minimal(base_size = 16) +
  ggplot2::ggtitle("Wavelet Coefficients")



#Threshold the small coefficients
finest <- sw[(N/2+1):N] # Select the finest level coefficients and estimate noise level
sigmahat <- sd(finest)# Estimate noise standard deviation
lambda <- sqrt(2 * log(N)) * sigmahat# Calculate the soft-threshold

# Apply soft-thresholding to wavelet coefficients
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


#  Return the signal to the time domain using the inverse transformation
a <- as.vector(t(WP) %*% swt) # Reconstruct the signal

#Plot the denoised signal and the noisy signal for comparison
ggplot2::ggplot(data.frame(t = t, signoi = signoi, a = a), ggplot2::aes(x = t)) +
  ggplot2::geom_point(ggplot2::aes(y = signoi), color = "red", size = 2) +
  ggplot2::geom_line(ggplot2::aes(y = a), color = "black", size = 1) +
  ggplot2::theme_minimal(base_size = 16) +
  ggplot2::ggtitle("Denoised Signal")
