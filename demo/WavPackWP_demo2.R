library(WaveletMatrixProject)

# Setting parameters for the wavelet packet transformation
nl <- 11      # J level, V_J
n <- 2^nl     # Signal length
level <- 5    # Levels of decomposition
shift <- 2    # Standard shift

filt <- c(sqrt(2)/2, sqrt(2)/2) # Haar wavelet filter coefficients

# Make Doppler Signal
t <- seq(0, 1, length.out = n)  # Time vector
y <- sqrt(t * (1 - t)) * sin((2 * pi * 1.05) / (t + 0.05)) # Doppler signal formula


# Generate the wavelet packet transformation matrix
WP <- WavPackMatWP(filt, n, level, shift)


# Reshape the Doppler signal into a column vector for matrix operations
y <- matrix(y, nrow = 2048, ncol = 1)

#Perform forward wavelet packet transformation
d <- WP %*% y #resulting d contains wavelet packet coefficients for Doppler signal


# Check summary statistics
cat('Summary of d in R:\n')
cat(sprintf('Min: %f, Max: %f, Mean: %f\n', min(d), max(d), mean(d)))
# accurate summary should be: Min: -1.2486, Max: 1.23, Mean: 0.0086858


#perform the inverse wavelet packet transformation
a <- t(WP) %*% d

# Plotting results
par(mfrow = c(2, 2))
# Plot the original Doppler signal
plot(y, type = 'l', main = "Original Signal", xlab = "Index", ylab = "Amplitude")

# Plot the wavelet packet coefficients
plot(d, type = 'l', main = "Doppler Wavelet Packet Coefficients", xlab = "Index", ylab = "Coefficient")

# Plot the reconstructed signal
plot(a, type = 'l', main = "Wavelet Packet Reconstruction", xlab = "Index", ylab = "Amplitude")

# Plot the difference between the reconstructed and original signals
plot(a - y, type = 'l', main = "Difference between Reconstruction and Original", xlab = "Index", ylab = "Difference")
#plots reconstruction error, which ideally is minimal.
