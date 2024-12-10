## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----setup--------------------------------------------------------------------
library(WaveletMatrixProject)

## ----WavmatND-example---------------------------------------------------------
# Define parameters
hf <- c(0.5, -0.5) # High-pass filter
N <- 4  # Signal length
k <- 3  # Decomposition levels
shift <- 1  # Filter shift

# Construct the wavelet matrix
W <- WavmatND(hf, N, k, shift)
print(W)

# Visualizing the wavelet matrix
image(W, main = "Wavelet Matrix (WavmatND)", xlab = "Columns", ylab = "Rows", col = heat.colors(256))

## ----WavPackMatWP-example-----------------------------------------------------
# Low-pass filter
h <- c(0.5, 0.5) #Low-pass filer
N <- 8 #Signal size (must be a power of 2, 2^3 in this case)
k0 <- 2 #Depth of decomposition
shift <- 2 #Shift parameter

# Construct the wavelet packet matrix
WP <- WavPackMatWP(h, N, k0, shift)
print(WP)


# Visualizing the wavelet packet matrix
image(WP, main = "Wavelet Packet Matrix (WavPackMatWP)", xlab = "Columns", ylab = "Rows", col = heat.colors(256))

## ----WavPackMat-example-------------------------------------------------------
# Low-pass filter
h <- c(0.5, 0.5)

# Parameters for wavelet transformation matrix
N <- 8 #Matrix size (2^3 in this case)
k0 <- 3 #Depth of wavelet decomposition
shift <- 2 #Shift parameter

# Construct the wavelet transformation matrix
WP_opt <- WavPackMat(h, N, k0, shift)
print(WP_opt)

# Visualize the matrix
image(WP_opt, main = "Wavelet Transformation Matrix", xlab = "Columns", ylab = "Rows", col = heat.colors(256))


## ----WavPackMat-Lena Example--------------------------------------------------

# Define wavelet filter from wavelets package
h <- c(-0.075765714789341, -0.029635527645954, 0.497618667632458,
       0.803738751805216, 0.297857795605542, -0.099219543576935,
       -0.012603967262261, 0.032223100604071)

# Load the Lena dataset (assuming it is available as a grayscale matrix)
data("lena")

# Transform lena data to match MATLAB's transformation of image data
A_full <- 255 - t(lena)

# Extract the central 256x256 portion for transformation
center_start <- (512 - 256) / 2 + 1
center_end <- center_start + 255
A <- A_full[center_start:center_end, center_start:center_end]

# Apply padding to the 256x256 submatrix as per MATLAB code
A[, 1:4] <- 255
A[, (ncol(A)-3):ncol(A)] <- 255
A[1:4, ] <- 255
A[(nrow(A)-3):nrow(A), ] <- 255

# Display the padded and transformed 256x256 Lena image using imager
imager::display(imager::as.cimg(A), main = "Padded and Transformed Lena Image (256x256)")


# Plot using base R
par(mfrow = c(2, 2))
image(A, col = gray.colors(256), main = "Padded and Transformed Lena Image (256x256)", useRaster = TRUE)

# Create wavelet matrix using WavPackMat with 256x256 size
N <- 256
k0 <- 2
W <- WavPackMat(h, N, k0)

# Apply the wavelet transformation: B = W * A * W'
B <- W %*% A %*% t(W)

# Display the transformed image using imager
imager::display(imager::as.cimg(B), main = "Transformed Lena Image (B)", rescale = FALSE)



# Plot using base R
image(B, col = gray.colors(256), main = "Transformed Lena Image (B)", zlim = c(0, 20), useRaster = TRUE)

# Inverse transformation: C = W' * B * W
C <- t(W) %*% B %*% W
imager::display(imager::as.cimg(C), main = "Inverse Transformed Image (C)")

# Plot using base R
image(C, col = gray.colors(256), main = "Inverse Transformed Image (C)", useRaster = TRUE)

# Difference image: D = A - C
D <- A - C
imager::display(imager::as.cimg(D), main = "Difference Image (D)", rescale = FALSE)

# Check the range of values in D
cat("Min of D:", min(D), "\n")
cat("Max of D:", max(D), "\n")

# Plot D with an adjusted zlim range
image(D, col = gray.colors(256), main = "Difference Image (D)", zlim = range(D), useRaster = TRUE)

# Scale D for better visualization
D_scaled <- D * 1000

# Checking the range of D_scaled
cat("Min of D_scaled:", min(D_scaled), "\n")
cat("Max of D_scaled:", max(D_scaled), "\n")

# Display using base R
image(D_scaled, col = gray.colors(256), main = "Scaled Difference Image (D)", useRaster = TRUE)

# Normalize D_scaled for imager display
D_scaled_normalized <- (D_scaled - min(D_scaled)) / (max(D_scaled) - min(D_scaled))
imager::display(imager::as.cimg(D_scaled_normalized), main = "Scaled and Normalized Difference Image (D)")



# Debugging Outputs
cat("Mean of A:", mean(A), "\n")
cat("Mean of B:", mean(B), "\n")
cat("Mean of C:", mean(C), "\n")
cat("Mean of D:", mean(D), "\n") #should be approximately 0


## ----WavPackMat-Signal Denoising Example--------------------------------------

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


## ----WavPackMatWP-Transformation and Reconstruction---------------------------
library(WaveletMatrixProject)

# Define the filter ( coefficients define the specific wavelet used for the transformation)
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

# This vector will be transformed using the wavelet packet matrix
y <- c(1, 0, -3, 2, 1, 0, 1, 2)

# Perform forward wavelet packet transformation
d <- sqrt(k0) * WP %*% y #scaled by `sqrt(k0)` for normalization

# Perform inverse wavelet packet transformation
yy <- (1 / sqrt(k0)) * t(WP) %*% d #inverse scaled by `1 / sqrt(k0)` for normalization

# Print the transformed and inverse-transformed vectors
print("d vector:")
print(d)
print("yy (inverse result):")
print(yy)
# Verify orthogonality of the wavelet packet transformation matrix
trace1 <- sum(diag(WP %*% t(WP))) #should about 4
trace2 <- sum(diag(t(WP) %*% WP)) #should about 4
cat("Trace of WP * WP':", trace1, "\n")
cat("Trace of WP' * WP:", trace2, "\n")


## ----WavPackMatWP-Analysis of a Doppler Signal--------------------------------
library(WaveletMatrixProject)

# Setting parameters for the wavelet packet transformation
nl <- 11  # J level, V_J
n <- 2^nl # Signal length
level <- 5 # Levels of decomposition
shift <- 2 # Standard shift

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



## ----WavPackMatND- Denoising a Noisy Doppler Signal---------------------------
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


## ----WavPackMatND-Denoising a Noisy Doppler Signal with Weighted Reconstruction----
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


