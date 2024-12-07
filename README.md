
# WaveletMatrixProject

<!-- badges: start -->
<!-- badges: end -->

**WaveletMatrixProject** is an R package designed for constructing and applying optimized orthogonal wavelet transformation matrices for signal processing tasks. The package supports wavelet packet transformations, signal decomposition, denoising, and reconstruction, making it suitable for both synthetic and real data analysis.

## Installation

You can install the 'WaveletMatrixProject' package, you can directly install it from [GitHub](https://github.com/saraantonijevic/WaveletMatrixProject), follow the steps below:

1. Install `devtools`
```r
install.packages("devtools")

devtools::install_github("saraantonijevic/WaveletMatrixProject")
```

## Basic Usage

Here is an example demonstrating how to use the 'WaveletMatrixProject' package for signal processing:


``` r
# Load the WaveletMatrixProject package
library(WaveletMatrixProject)


# Generate a Doppler signal
n <- 1000
t <- seq(1/n, 1, length.out = n)
s <- sqrt(t * (1 - t)) * sin((2 * pi * 1.05) / (t + 0.05))  # Doppler signal

# Add Gaussian noise
set.seed(42)
sigma <- 0.1
noisy_signal <- s + rnorm(n, mean = 0, sd = sigma)

# Plot the noisy signal
plot(t, noisy_signal, type = "l", col = "red", main = "Noisy Doppler Signal")
lines(t, s, col = "blue", lty = 2)  # Overlay original signal for comparison
legend("topright", legend = c("Noisy Signal", "Original Signal"),
       col = c("red", "blue"), lty = c(1, 2))
       
       
# Apply Wavelet Packet Transformation

# Define a Haar wavelet filter
qmf <- c(1/sqrt(2), 1/sqrt(2))

# Create the wavelet transformation matrix
k <- 4  # Number of decomposition levels
W <- WavmatND(qmf, n, k, shift = 0)

# Transform the noisy signal
transformed_signal <- W %*% matrix(noisy_signal, nrow = n)

# Plot the transformed signal
plot(transformed_signal, type = "l", col = "purple",
     main = "Transformed Signal (Wavelet Coefficients)")


# Thresholding and Signal Reconstruction# Thresholding
threshold <- sqrt(2 * log(n)) * sd(transformed_signal)
thresholded_signal <- ifelse(abs(transformed_signal) > threshold,
                             transformed_signal, 0)

# Reconstruct the signal
reconstructed_signal <- as.vector(t(W) %*% thresholded_signal)

# Plot the denoised signal
plot(t, reconstructed_signal, type = "l", col = "green", main = "Denoised Signal")
lines(t, s, col = "blue", lty = 2)  # Overlay original signal
legend("topright", legend = c("Denoised Signal", "Original Signal"),
       col = c("green", "blue"), lty = c(1, 2))

```
## Features

The `WaveletMatrixProject` package offers the following functionalities:

- **Wavelet Matrix Construction**:
  - Create wavelet and wavelet packet transformation matrices using customizable filters.
  - Supports a variety of wavelet filters, such as Haar and Symmlet.

- **Signal Transformation**:
  - Apply wavelet or wavelet packet transformations to signals for decomposition.
  - Analyze signals in the wavelet domain to identify frequency components.

- **Thresholding and Denoising**:
  - Perform hard or soft thresholding on wavelet coefficients to remove noise.
  - Reconstruct the denoised signal from thresholded wavelet coefficients.

- **Visualization**:
  - Visualize the original, noisy, and denoised signals.
  - Explore wavelet coefficients using plots for better insight of their signal components.


## Vignette
For a detailed overview and examples, look at the [WaveletMatrixOverview vignette](vignettes/WaveletMatrixOverview.md).


## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE.md) file for more details.

## Contact

Email: saraant@tamu.edu

