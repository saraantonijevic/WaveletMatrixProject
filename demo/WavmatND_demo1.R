n <- 1000  # signal size setting
J <- floor(log2(n))  # level where the signal resides
k <- 4 # number of steps in ND transformation

# Generate a doppler signal
t <- seq(1/n, 1, length.out = n)
s <- sqrt(t * (1 - t)) * sin((2 * pi * 1.05) / (t + 0.05))

sigma <- 0.1  # noise level
noise <- rnorm(length(s), mean = 0, sd = sigma)
sn <- s + noise
var(sn)

# a noisy doppler signal
qmf <- c(1/sqrt(2), 1/sqrt(2))  # Haar wavelet

W <- WavmatND(qmf, n, k, 0)  # wavematND
# Function to create a weight matrix in R similar to MATLAB's weight(n, k)
# Function to create a weight matrix in R similar to MATLAB's weight(n, k)
weight <- function(n, k) {
  size <- (k + 1) * n  # Total size for the matrix (5000 when k = 4 and n = 1000)
  T <- matrix(0, nrow = size, ncol = size)  # Initialize an empty matrix

  # Assign weights so that the trace of T sums up to n (1000 in this case)
  for (i in 1:(k + 1)) {
    start_row <- (i - 1) * n + 1
    end_row <- i * n
    weight_value <- 1 / (k + 1)  # Ensure the sum of the weights in each block scales correctly
    T[start_row:end_row, start_row:end_row] <- diag(weight_value, n, n)
  }

  return(T)
}

T <- weight(n, k)


# Signal transformation
tsn <- W %*% matrix(sn, nrow = n, ncol = 1)
cat(sprintf("Variance of transformed signal: %f\n", var(tsn)))
cat(sprintf("Variance of last detail level: %f\n", var(tsn[((k - 1) * n + 1):(k * n)])))

# Extract detail levels and apply thresholding
temp <- tsn[(n + 1):length(tsn)]
threshold <- sqrt(2 * log(n) * var(tsn[(length(tsn) - n + 1):length(tsn)]))
temp[abs(temp) < threshold] <- 0

# Reconstruct the signal
rs <- t(W) %*% T %*% c(tsn[1:n], temp)
trace <- sum(diag(T))
cat(sprintf("Trace of T: %f\n", trace))
cat(sprintf("Variance of reconstructed signal: %f\n", var(rs)))

# Create data frames for plotting
df_original <- data.frame(Index = 1:n, Signal = s)
df_noisy <- data.frame(Index = 1:n, Signal = sn)
df_denoised <- data.frame(Index = 1:n, Signal = rs)
df_error <- data.frame(Index = 1:n, Signal = rs - s)

# Plot using ggplot2
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
