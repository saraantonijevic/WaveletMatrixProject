# Load necessary functions and data from the package
library(WaveletMatrixProject)

# Define wavelet filter from wavelets package
h <- c(-0.075765714789341, -0.029635527645954, 0.497618667632458,
       0.803738751805216, 0.297857795605542, -0.099219543576935,
       -0.012603967262261, 0.032223100604071)

# Load the Lena dataset (assuming it is available as a grayscale matrix)
data("lena")

# Transform lena data
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
