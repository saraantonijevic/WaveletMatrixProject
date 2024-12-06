library(WaveletMatrixProject)

# Define parameters for the signal and transformation
n <- 1000  # signal size setting
J <- floor(log2(n))  # Maximum decomposition level
k <- 4 #decomposition levels

# Generate a doppler signal
t <- seq(1/n, 1, length.out = n) # Time vector
s <- sqrt(t * (1 - t)) * sin((2 * pi * 1.05) / (t + 0.05))# Doppler signal formula

# Add Gaussian noise to the signal
sigma <- 0.1  # noise level
noise <- rnorm(length(s), mean = 0, sd = sigma) # Generate noise
sn <- s + noise # Noisy Doppler signal
var(sn) # Variance: Should be approximately 0.092 in order to correspond to correct Matlab calculation


# Define the Haar wavelet filter
qmf <- c(1/sqrt(2), 1/sqrt(2))  # Haar wavelet

# Construct the ND wavelet transformation matrix
W <- WavmatND(qmf, n, k, 0)  # wavematND

# Define a weight matrix function
weight <- function(n, k) {# This function creates a diagonal weight matrix used for reconstructing the signal
  size <- (k + 1) * n  # Total size for the matrix (5000 when k = 4 and n = 1000)
  T <- matrix(0, nrow = size, ncol = size)  # Initialize an empty matrix

  # Assign weights so that the trace of T sums up to n (1000 in this case)
  for (i in 1:(k + 1)) {
    start_row <- (i - 1) * n + 1
    end_row <- i * n
    weight_value <- 1 / (k + 1)  # Equal weights for each block
    T[start_row:end_row, start_row:end_row] <- diag(weight_value, n, n)
  }

  return(T)
}


# Generate the weight matrix
T <- weight(n, k)


# Signal transformation
tsn <- W %*% matrix(sn, nrow = n, ncol = 1)
cat(sprintf("Variance of transformed signal : %f\n", var(tsn))) # Should be approximately 0.2721
cat(sprintf("Variance of last detail level: %f\n", var(tsn[((k - 1) * n + 1):(k * n)]))) # Should be approximately 0.0132

# Extract detail levels and apply hard thresholding for denoising
temp <- tsn[(n + 1):length(tsn)]# Detail coefficients
threshold <- sqrt(2 * log(n) * var(tsn[(length(tsn) - n + 1):length(tsn)])) # Threshold value
temp[abs(temp) < threshold] <- 0 # Zero out small coefficients

# Reconstruct the signal using the inverse transformation
rs <- t(W) %*% T %*% c(tsn[1:n], temp)
trace <- sum(diag(T))# Verify the trace of T (should be n)
cat(sprintf("Trace of T: %f\n", trace)) # Should be 1000
cat(sprintf("Variance of reconstructed signal: %f\n", var(rs))) # Should be approximately 0.0777

cat("Dimension of T: ", dim(T)) # 5000 5000
cat("Dimension of T: ", dim(W)) # 5000 1000
# Create data frames for plotting

# Original signal
df_original <- data.frame(Index = 1:n, Signal = s)
# Noisy signal
df_noisy <- data.frame(Index = 1:n, Signal = sn)
# Denoised signal
df_denoised <- data.frame(Index = 1:n, Signal = rs)
# Reconstruction error
df_error <- data.frame(Index = 1:n, Signal = rs - s)

# Plot the original, noisy, denoised, and error signals
p1 <- ggplot2::ggplot(df_original, ggplot2::aes(x = Index, y = Signal)) +
  ggplot2::geom_line() +
  ggplot2::ggtitle('Original Doppler Signal') +
  ggplot2::theme_minimal()

p2 <- ggplot2::ggplot(df_noisy, ggplot2::aes(x = Index, y = Signal)) +
  ggplot2::geom_line() +
  ggplot2::ggtitle('Noisy Doppler Signal') +
  ggplot2::theme_minimal()

p3 <- ggplot2::ggplot(df_denoised, ggplot2::aes(x = Index, y = Signal)) +
  ggplot2::geom_line() +
  ggplot2::ggtitle('De-noised Doppler Signal (Hard-Thresholding)') +
  ggplot2::theme_minimal()

p4 <- ggplot2::ggplot(df_error, ggplot2::aes(x = Index, y = Signal)) +
  ggplot2::geom_line() +
  ggplot2::ggtitle('Error') +
  ggplot2::theme_minimal()

# Arrange plots in a 4x1 grid layout
gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 1)
