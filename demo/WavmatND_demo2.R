library(WaveletMatrixProject)


# Set parameters for signal generation and transformation
sigma <- 0.05# Noise level
m <- 250  # Signal size
t <- seq(1/m, 1, length.out = m)# Time vector


# Generate a Doppler signal
s <- sqrt(t * (1 - t)) * sin((2 * pi * 1.05) / (t + 0.05))   # Doppler signal formula

# Add Gaussian noise to the signal
noise <- rnorm(m, 0, sigma)# Generate noise
sn <- s + noise  # Noisy doppler signal

# Compute decomposition level and wavelet filter
J <- floor(log2(m)) # Maximum level of decomposition for signal size m
k <- J - 1  # decomposition level
qmf <- c(1/sqrt(2), 1/sqrt(2))  # Haar wavelet filter coefficients

W <- WavmatND(qmf, m, k, 0)

# Define a function to create a weight matrix
weight <- function(m, k) {
  size <- (k + 1) * m # Total size for the matrix (5000 when k = 4 and n = 1000)
  T <- matrix(0, nrow = size, ncol = size)  # Initialize an empty matrix

  # Assign weights so that the trace of T sums up to n (1000 in this case)
  for (i in 1:(k + 1)) {
    start_row <- (i - 1) * m + 1
    end_row <- i * m
    weight_value <- 1 / (k + 1)  # Equal weights for each block
    T[start_row:end_row, start_row:end_row] <- diag(weight_value, m, m)
  }

  return(T)
}

# Generate the weight matrix
T <- weight(m, k)

# Transform the noisy signal using the wavelet matrix
tsn <- W %*% sn  # Transformed noisy signal

# Hard-thresholding for denoising; estimate noise level from the detail coefficients
sigmahat <- (var(tsn[(m + 1):length(tsn)][seq(1, length(tsn[(m + 1):length(tsn)]), 2)]) +
               var(tsn[(m + 1):length(tsn)][seq(2, length(tsn[(m + 1):length(tsn)]), 2)])) / 2

# Calculate the threshold for hard-thresholding
threshold <- sqrt(2 * log(k * m) * sigmahat)

# Apply hard-thresholding to the detail coefficients
snt <- tsn[(m + 1):length(tsn)] * (abs(tsn[(m + 1):length(tsn)]) > threshold)

# Combine approximation coefficients and thresholded detail coefficients
tsn_combined <- c(tsn[1:m], snt)

# Reconstruct the signal using the inverse transformation (transpose of W)
rs <- t(W) %*% T %*% tsn_combined  #  Apply weight matrix for proper reconstruction



# Plot denoising results

# Plot the original Doppler signal
p1 <- ggplot2::ggplot(data.frame(x = 1:m, y = s), ggplot2::aes(x, y)) +
  ggplot2::geom_line() +
  ggplot2::ggtitle("Original doppler signal") +
  ggplot2::xlim(0, m)


# Plot the noisy Doppler signal
p2 <- ggplot2::ggplot(data.frame(x = 1:m, y = sn), ggplot2::aes(x, y)) +
  ggplot2::geom_line() +
  ggplot2::ggtitle("Noisy doppler signal") +
  ggplot2::xlim(0, m)


# Prepare data for the de-noised signal plot
denoised_data <- data.frame(x = 1:m, y = rs[1:m])

# Plot the de-noised Doppler signal
p3 <- ggplot2::ggplot(denoised_data, ggplot2::aes(x, y)) +
  ggplot2::geom_line() +
  ggplot2::ggtitle("De-noised doppler signal") +
  ggplot2::xlim(0, m)


# Prepare data for the transformed noisy signal plot
tsn_data <- data.frame(x = 1:(m * (k + 1)), y = tsn)

# Plot the transformed noisy Doppler signal with NDWT
p4 <- ggplot2::ggplot(tsn_data, ggplot2::aes(x, y)) +
  ggplot2::geom_line() +
  ggplot2::ggtitle("Transformed noisy doppler signal with NDWT") +
  ggplot2::xlim(0, m * (k + 1))

# Arrange the plots in a 2x3 grid layout
gridExtra::grid.arrange(p1, p2, p3, p4, layout_matrix = matrix(c(1, 2, 3, 4, 4, 4), nrow = 2, byrow = TRUE))
