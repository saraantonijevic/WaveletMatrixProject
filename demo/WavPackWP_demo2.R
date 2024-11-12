library(WaveletMatrixProject)
# WavPackWP_demo2.R
# This demo script showcases the forward and inverse wavelet packet transformation applied to a Doppler signal.


nl <- 11      # J level, V_J
n <- 2^nl     # Signal length
level <- 5    # Levels of decomposition
shift <- 2    # Standard shift

filt <- c(sqrt(2)/2, sqrt(2)/2)

# Make Doppler Signal
t <- seq(0, 1, length.out = n)
y <- sqrt(t * (1 - t)) * sin((2 * pi * 1.05) / (t + 0.05))

WP <- WavPackMatWP(filt, n, level, shift)

y <- matrix(y, nrow = 2048, ncol = 1)

d <- WP %*% y


# Check summary statistics
cat('Summary of d in R:\n')
cat(sprintf('Min: %f, Max: %f, Mean: %f\n', min(d), max(d), mean(d)))

a <- t(WP) %*% d

# Plotting results
par(mfrow = c(2, 2))
plot(y, type = 'l', main = "Original Signal", xlab = "Index", ylab = "Amplitude")
plot(d, type = 'l', main = "Doppler Wavelet Packet Coefficients", xlab = "Index", ylab = "Coefficient")
plot(a, type = 'l', main = "Wavelet Packet Reconstruction", xlab = "Index", ylab = "Amplitude")
plot(a - y, type = 'l', main = "Difference between Reconstruction and Original", xlab = "Index", ylab = "Difference")
