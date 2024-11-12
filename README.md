
# WaveletMatrixProject

<!-- badges: start -->
<!-- badges: end -->

**WaveletMatrixProject** is an R package designed for constructing and applying highly optimized orthogonal wavelet transformation matrices for signal processing tasks. The package supports wavelet packet transformations, signal decomposition, denoising, and reconstruction, making it suitable for both synthetic and real data analysis.

## Installation

You can install the development version of WaveletMatrixProject from [GitHub](https://github.com/saraantonijevic/WaveletMatrixProject) with:

```r
# Install pak if not already installed
install.packages("pak")

# Install the WaveletMatrixProject package
pak::pak("saraantonijevic/WaveletMatrixProject")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
# Load the WaveletMatrixProject package
library(WaveletMatrixProject)

# Generate a synthetic "bumps" signal
N <- 1024
t <- seq(0, 1, length.out = N)
pos <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81)
hgt <- c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
wth <- c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005, 0.008, 0.005)

sig <- numeric(N)
for (j in 1:length(pos)) {
  sig <- sig + hgt[j] / (1 + abs((t - pos[j]) / wth[j]))^4
}

# Plot the original signal
plot(t, sig, type = "l", col = "green", lwd = 2, main = "Original Bumps Signal", xlab = "Time", ylab = "Amplitude")

# Construct the wavelet packet transformation matrix
filt <- c(-0.07576571478934, -0.02963552764595, 0.49761866763246, 0.80373875180522, 
          0.29785779560554, -0.09921954357694, -0.01260396726226, 0.03222310060407)
WP <- WavPackMat(filt, N, k0 = 6)

# Transform the signal using the wavelet packet matrix
sw <- as.vector(WP %*% sig)

# Display the transformed signal
plot(sw, type = "l", lwd = 2, main = "Wavelet Packet Coefficients", xlab = "Index", ylab = "Coefficient")

```

