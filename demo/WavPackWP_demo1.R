library(WaveletMatrixProject)
# WavPackWP_demo1.R
# This demo script demonstrates the construction of a wavelet packet transformation matrix
# and its application to a sample vector.


# Define the filter (e.g., SYMMLET 4)
filt <- c(-0.07576571478934, -0.02963552764595,
          0.49761866763246, 0.80373875180522,
          0.29785779560554, -0.09921954357694,
          -0.01260396726226, 0.03222310060407)
# Print to verify
print(filt)


# Set parameters for the wavelet packet transformation matrix
N <- 8  # Size of the matrix
k0 <- 3  # Levels of decomposition
shift <- 0  # Shift parameter

# Generate the wavelet packet transformation matrix
WP <- WavPackMatWP(filt, N, k0, shift)

# Print the matrix
print(WP)

# Example of applying this to a vector y
y <- c(1, 0, -3, 2, 1, 0, 1, 2)

d <- sqrt(k0) * WP %*% y
yy <- (1 / sqrt(k0)) * t(WP) %*% d

# Print the transformed and inverse-transformed vectors
print("d vector:")
print(d)
print("yy (inverse result):")
print(yy)

trace1 <- sum(diag(WP %*% t(WP)))
trace2 <- sum(diag(t(WP) %*% WP))
cat("Trace of WP * WP':", trace1, "\n")
cat("Trace of WP' * WP:", trace2, "\n")

orthogonality_check <- WP %*% t(WP)
print(orthogonality_check)
