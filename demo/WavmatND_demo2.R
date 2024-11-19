library(WaveletMatrixProject)

sigma <- 0.05
m <- 250  # noise level, signal size
t <- seq(1/m, 1, length.out = m)
s <- sqrt(t * (1 - t)) * sin((2 * pi * 1.05) / (t + 0.05))  # Generate a doppler signal
noise <- rnorm(m, 0, sigma)
sn <- s + noise  # a noisy doppler signal
J <- floor(log2(m))
k <- J - 1  # decomposition level
qmf <- c(1/sqrt(2), 1/sqrt(2))  # Haar wavelet

W <- WavmatND(qmf, m, k, 0)  # wavematND
weight <- function(m, k) {
  size <- (k + 1) * m # Total size for the matrix (5000 when k = 4 and n = 1000)
  T <- matrix(0, nrow = size, ncol = size)  # Initialize an empty matrix

  # Assign weights so that the trace of T sums up to n (1000 in this case)
  for (i in 1:(k + 1)) {
    start_row <- (i - 1) * m + 1
    end_row <- i * m
    weight_value <- 1 / (k + 1)  # Ensure the sum of the weights in each block scales correctly
    T[start_row:end_row, start_row:end_row] <- diag(weight_value, m, m)
  }

  return(T)
}

T <- weight(m, k)

# Transform the noisy signal using the wavelet matrix
tsn <- W %*% sn  # Transformed signal

# Hard-thresholding for denoising
# Estimate noise level from the detail coefficients
sigmahat <- (var(tsn[(m + 1):length(tsn)][seq(1, length(tsn[(m + 1):length(tsn)]), 2)]) +
               var(tsn[(m + 1):length(tsn)][seq(2, length(tsn[(m + 1):length(tsn)]), 2)])) / 2

# Calculate the threshold
threshold <- sqrt(2 * log(k * m) * sigmahat)

# Apply hard-thresholding to the detail coefficients
snt <- tsn[(m + 1):length(tsn)] * (abs(tsn[(m + 1):length(tsn)]) > threshold)

# Combine approximation coefficients and thresholded detail coefficients
tsn_combined <- c(tsn[1:m], snt)

# Reconstruct the signal using the inverse transformation (transpose of W)
rs <- t(W) %*% T %*% tsn_combined  # Ensure T is applied to maintain the same transformation as in MATLAB



# Plot denoising results
p1 <- ggplot2::ggplot(data.frame(x = 1:m, y = s), ggplot2::aes(x, y)) +
  ggplot2::geom_line() +
  ggplot2::ggtitle("Original doppler signal") +
  ggplot2::xlim(0, m)

p2 <- ggplot2::ggplot(data.frame(x = 1:m, y = sn), ggplot2::aes(x, y)) +
  ggplot2::geom_line() +
  ggplot2::ggtitle("Noisy doppler signal") +
  ggplot2::xlim(0, m)


# Correct the data for the de-noised signal plot
denoised_data <- data.frame(x = 1:m, y = rs[1:m])

# Create the plot for the de-noised doppler signal
p3 <- ggplot2::ggplot(denoised_data, ggplot2::aes(x, y)) +
  ggplot2::geom_line() +
  ggplot2::ggtitle("De-noised doppler signal") +
  ggplot2::xlim(0, m)
# Create a data frame for plotting the transformed noisy signal
tsn_data <- data.frame(x = 1:(m * (k + 1)), y = tsn)

# Create the plot for the transformed noisy doppler signal
p4 <- ggplot2::ggplot(tsn_data, ggplot2::aes(x, y)) +
  ggplot2::geom_line() +
  ggplot2::ggtitle("Transformed noisy doppler signal with NDWT") +
  ggplot2::xlim(0, m * (k + 1))

# Arrange the plots in a 2x3 grid layout
gridExtra::grid.arrange(p1, p2, p3, p4, layout_matrix = matrix(c(1, 2, 3, 4, 4, 4), nrow = 2, byrow = TRUE))
