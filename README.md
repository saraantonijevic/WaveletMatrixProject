WaveletMatrixProject
================

# WaveletMatrixProject

<!-- badges: start -->
<!-- badges: end -->

**WaveletMatrixProject** is an R package designed for constructing and
applying optimized orthogonal wavelet transformation matrices for signal
processing tasks. The package supports wavelet packet transformations,
signal decomposition, denoising, and reconstruction, making it suitable
for both synthetic and real data analysis.

## Installation

You can install the ‘WaveletMatrixProject’ package, you can directly
install it from
[GitHub](https://github.com/saraantonijevic/WaveletMatrixProject),
follow the steps below:

1.  Install `devtools`

``` r
#If you have not yet, install devtools.
install.packages("devtools")

#Installing WaveletMatrix Project from GitHub:
devtools::install_github("saraantonijevic/WaveletMatrixProject")
```

To install this package with the vignette, then run the following code:

``` r
install.packages("clinfun")
install.packages("rmarkdown")

devtools::install_github("saraantonijevic/WaveletMatrixProject", build_vignettes = TRUE, force = TRUE)
```

## Usage

The following are examples demonstrating how to use the
‘WaveletMatrixProject’ package for signal processing:

### Example 1: Wavelet Packet Transformation of a Doppler Signal

**Step 1:** Load the WaveletMatrixProject Package

``` r
# Load the WaveletMatrixProject package
library(WaveletMatrixProject)
```

**Step 2:** Generate a Noisy Doppler Signal

We first create a synthetic Doppler signal, which is commonly used in
signal processing to evaluate the effectiveness of transformations:

``` r

# Parameters for the wavelet transformation
nl <- 11  # Number of levels (J-level)
n <- 2^nl # Signal length
level <- 5  # Decomposition levels
shift <- 2  # Standard shift
filt <- c(sqrt(2)/2, sqrt(2)/2)  # Haar wavelet filter coefficients

# Generate a Doppler signal
t <- seq(0, 1, length.out = n)
y <- sqrt(t * (1 - t)) * sin((2 * pi * 1.05) / (t + 0.05))
```

**Step 3:** Perform Wavelet Packet Transformation

``` r
# Generate the wavelet packet transformation matrix
WP <- WavPackMatWP(filt, n, level, shift)

# Reshape the Doppler signal into a column vector
y <- matrix(y, nrow = n, ncol = 1)

# Perform the forward wavelet packet transformation
d <- WP %*% y

# Perform the inverse transformation for reconstruction
a <- t(WP) %*% d
```

**Step 4:** Visualize Results

Visualize the original signal, coefficients, reconstruction, and
reconstruction error:

``` r

# Plot the original signal
p1 <- ggplot2::ggplot(data.frame(x = 1:n, y = y[, 1]), ggplot2::aes(x, y)) +
  ggplot2::geom_line() +
  ggplot2::ggtitle("Original Doppler Signal")

# Plot wavelet packet coefficients
p2 <- ggplot2::ggplot(data.frame(x = 1:length(d), y = d[, 1]), ggplot2::aes(x, y)) +
  ggplot2::geom_line() +
  ggplot2::ggtitle("Wavelet Packet Coefficients")

# Plot reconstructed signal
p3 <- ggplot2::ggplot(data.frame(x = 1:n, y = a[, 1]), ggplot2::aes(x, y)) +
  ggplot2::geom_line() +
  ggplot2::ggtitle("Reconstructed Signal")

# Plot reconstruction error
p4 <- ggplot2::ggplot(data.frame(x = 1:n, y = a[, 1] - y[, 1]), ggplot2::aes(x, y)) +
  ggplot2::geom_line() +
  ggplot2::ggtitle("Reconstruction Error")

# Arrange the plots
gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
```

### Example 2: Denoising a Noisy Doppler Signal

**Step 1:** Generate a Doppler Signal

Simulate a noisy signal by adding Gaussian noise to the Doppler signal:

``` r
# Load the package
library(WaveletMatrixProject)

# Signal parameters
sigma <- 0.05  # Noise level
m <- 250  # Signal size
t <- seq(1/m, 1, length.out = m)

# Generate a Doppler signal and add noise
s <- sqrt(t * (1 - t)) * sin((2 * pi * 1.05) / (t + 0.05))
noise <- rnorm(m, 0, sigma)
sn <- s + noise
```

**Step 2:** Perform Wavelet Packet Transformation

Transform the noisy signal, estimate the noise level, and apply
thresholding:

``` r
# Decomposition level and wavelet filter
J <- floor(log2(m))
k <- J - 1
qmf <- c(1/sqrt(2), 1/sqrt(2))
W <- WavmatND(qmf, m, k, 0)

# Transform noisy signal
tsn <- W %*% sn

# Estimate noise level and calculate threshold
sigmahat <- (var(tsn[(m + 1):length(tsn)][seq(1, length(tsn[(m + 1):length(tsn)]), 2)]) +
               var(tsn[(m + 1):length(tsn)][seq(2, length(tsn[(m + 1):length(tsn)]), 2)])) / 2
threshold <- sqrt(2 * log(k * m) * sigmahat)

# Apply hard-thresholding
snt <- tsn[(m + 1):length(tsn)] * (abs(tsn[(m + 1):length(tsn)]) > threshold)
tsn_combined <- c(tsn[1:m], snt)

# Reconstruct the signal
T <- diag(1, nrow = length(tsn_combined))  # Optional weighting matrix
rs <- t(W) %*% T %*% tsn_combined
```

**Step 3:** Visualize Denoising Results

Plot the original, noisy, and denoised signals:

``` r

# Prepare data frames for plotting
original_signal <- data.frame(x = 1:m, y = s)
noisy_signal <- data.frame(x = 1:m, y = sn)
denoised_signal <- data.frame(x = 1:m, y = rs[1:m])
transformed_signal <- data.frame(x = 1:length(tsn), y = tsn)

# Create individual plots
p1 <- ggplot2::ggplot(original_signal, ggplot2::aes(x, y)) +
  ggplot2::geom_line() +
  ggplot2::ggtitle("Original Signal") +
  ggplot2::labs(x = "Index", y = "Amplitude")

p2 <- ggplot2::ggplot(noisy_signal, ggplot2::aes(x, y)) +
  ggplot2::geom_line() +
  ggplot2::ggtitle("Noisy Signal") +
  ggplot2::labs(x = "Index", y = "Amplitude")

p3 <- ggplot2::ggplot(denoised_signal, ggplot2::aes(x, y)) +
  ggplot2::geom_line() +
  ggplot2::ggtitle("Denoised Signal") +
  ggplot2::labs(x = "Index", y = "Amplitude")

p4 <- ggplot2::ggplot(transformed_signal, ggplot2::aes(x, y)) +
  ggplot2::geom_line() +
  ggplot2::ggtitle("Transformed Signal") +
  ggplot2::labs(x = "Index", y = "Coefficient")

# Arrange the plots in a grid
gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
```

This demonstrates how the wavelet packet coefficients and reconstruction
align closely with the original signal.

**Summary:** These examples demonstrate how the WaveletMatrixProject
package can transform, analyze, and denoise signals using matrix-based
wavelet transforms. The results showcase the way wavelet transforms are
applied to reduce noise while retaining essential signal features.

## Features

The `WaveletMatrixProject` package offers the following functionalities:

- **Wavelet Matrix Construction**:
  - Create wavelet and wavelet packet transformation matrices using
    customizable filters.
  - Supports a variety of wavelet filters, such as Haar and Symmlet.
- **Signal Transformation**:
  - Apply wavelet or wavelet packet transformations to signals for
    decomposition.
  - Analyze signals in the wavelet domain to identify frequency
    components.
- **Thresholding and Denoising**:
  - Perform hard or soft thresholding on wavelet coefficients to remove
    noise.
  - Reconstruct the denoised signal from thresholded wavelet
    coefficients.
- **Visualization**:
  - Visualize the original, noisy, and denoised signals.
  - Explore wavelet coefficients using plots for better insight of their
    signal components.

## Vignette

For a detailed overview and examples, check out the
[WaveletMatrixProject
vignette](https://saraantonijevic.github.io/WaveletMatrixProject/articles/WaveletMatrixProject.html).

## License

This project is licensed under the MIT License. See the
[LICENSE](LICENSE.md) file for more details.
