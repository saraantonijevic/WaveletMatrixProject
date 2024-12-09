``` r
library(WaveletMatrixProject)
```

# Introduction

The `WaveletMatrixProject` is a package built for the construction and
manipulation of wavelets and wavelet transformation matrices. These
tools are essential for signal processing tasks, such as:

- Decomposing signals into multiple resolutions.
- Analyzing frequency components across different scales.
- Enhancing feature extraction, and reconstructing signals from wavelet
  transformations.

This vignette serves as a guide to the package, composing of:

- **Detailed Core-Function Descriptions**: Clear explanations of each
  function’s purpose and role in signal processing.
- **Code Examples**: Step-by-step demonstrations to help you understand
  usage.
- **Visual Outputs**: Graphs and plots to illustrate the effects of
  transformations.

By the end of this vignette, you will be able to:

- Construct wavelets and wavelet packet matrices using the functions.
- Apply advanced transformations for multi-resolution signal analysis.
- Combine supporting utilities like filter dilation and vector embedding
  to enhance workflow.

# Core Functions

## 1. `WavmatND()`

**Description**: The WavmatND() function generates a wavelet matrix by
applying a high-pass filter to an identity matrix. This then constructs
a transformation matrix that captures the high-frequency details of a
signal. The function is useful in non-decimated wavelet transforms where
no down-sampling is performed. Therefore, it preserves the original
signal length.

**Key Arguments**:

- **hf**: Vector of high-pass filter coefficients.
- **N**: Size of the wavelet matrix (signal length).
- **k**: Decomposition levels.
- **shift**: Shift needed to be applied to the filter at each iteration.

WavmatND() supports multi-resolution analysis by enabling signal
decomposition at finer scales.

### Example

``` r
# Define parameters
hf <- c(0.5, -0.5) # High-pass filter
N <- 4  # Signal length
k <- 3  # Decomposition levels
shift <- 1  # Filter shift

# Construct the wavelet matrix
W <- WavmatND(hf, N, k, shift)
print(W)
#>         [,1]   [,2]   [,3]   [,4]
#>  [1,] -0.125 -0.125  0.125  0.125
#>  [2,]  0.125 -0.125 -0.125  0.125
#>  [3,]  0.125  0.125 -0.125 -0.125
#>  [4,] -0.125  0.125  0.125 -0.125
#>  [5,]  0.125  0.125 -0.125 -0.125
#>  [6,] -0.125  0.125  0.125 -0.125
#>  [7,] -0.125 -0.125  0.125  0.125
#>  [8,]  0.125 -0.125 -0.125  0.125
#>  [9,] -0.250  0.250 -0.250  0.250
#> [10,]  0.250 -0.250  0.250 -0.250
#> [11,] -0.250  0.250 -0.250  0.250
#> [12,]  0.250 -0.250  0.250 -0.250
#> [13,]  0.000 -0.500 -0.500  0.000
#> [14,]  0.000  0.000 -0.500 -0.500
#> [15,] -0.500  0.000  0.000 -0.500
#> [16,] -0.500 -0.500  0.000  0.000

# Visualizing the wavelet matrix
image(W, main = "Wavelet Matrix (WavmatND)", xlab = "Columns", ylab = "Rows", col = heat.colors(256))
```

![](C:/Users/saraa/Desktop/WaveletMatrixProject/vignettes/WaveletMatrixOverview_files/figure-gfm/WavmatND-example-1.png)<!-- -->
**Explanation**: The WavmatND() function starts with an identity matrix
and applies a high-pass filter to construct the wavelet matrix. This
matrix can be used for signal processing tasks where transformations are
needed to capture specific frequency components in signals.

**Supporting Functions**:

- **dilate_filter()**: Dilates the filter coefficients for higher
  decomposition levels.
- **imbed()**: Embeds a vector by appending zeros.
- **repeating()**: Creates repeated sequences from input vectors.

For a more detailed demonstration, see WavmatND_demo1 and WavmatND_demo2
in the demo folder.

## 2. `WavPackMatWP()`

**Description**: Constructs a wavelet packet transformation matrix using
a low-pass filter. This matrix represents a transformation that
decomposes signals into components at different levels of detail.

**Key Arguments:**

- **h**: Vector for the low-pass filter.
- **N**: Size of the matrix (must be a power of 2).
- **k0**: Depth of the wavelet transformation.
- **shift**: Shifts in the wavelet transformation (default is 2).

### Example

``` r
# Low-pass filter
h <- c(0.5, 0.5) #Low-pass filer
N <- 8 #Signal size (must be a power of 2, 2^3 in this case)
k0 <- 2 #Depth of decomposition
shift <- 2 #Shift parameter

# Construct the wavelet packet matrix
WP <- WavPackMatWP(h, N, k0, shift)
print(WP)
#>            [,1]       [,2]       [,3]       [,4]
#>  [1,] 0.3535534  0.3535534  0.0000000  0.0000000
#>  [2,] 0.0000000  0.0000000  0.3535534  0.3535534
#>  [3,] 0.0000000  0.0000000  0.0000000  0.0000000
#>  [4,] 0.0000000  0.0000000  0.0000000  0.0000000
#>  [5,] 0.3535534 -0.3535534  0.0000000  0.0000000
#>  [6,] 0.0000000  0.0000000  0.3535534 -0.3535534
#>  [7,] 0.0000000  0.0000000  0.0000000  0.0000000
#>  [8,] 0.0000000  0.0000000  0.0000000  0.0000000
#>  [9,] 0.1767767  0.1767767  0.1767767  0.1767767
#> [10,] 0.0000000  0.0000000  0.0000000  0.0000000
#> [11,] 0.1767767  0.1767767 -0.1767767 -0.1767767
#> [12,] 0.0000000  0.0000000  0.0000000  0.0000000
#> [13,] 0.1767767 -0.1767767  0.1767767 -0.1767767
#> [14,] 0.0000000  0.0000000  0.0000000  0.0000000
#> [15,] 0.1767767 -0.1767767 -0.1767767  0.1767767
#> [16,] 0.0000000  0.0000000  0.0000000  0.0000000
#>            [,5]       [,6]       [,7]       [,8]
#>  [1,] 0.0000000  0.0000000  0.0000000  0.0000000
#>  [2,] 0.0000000  0.0000000  0.0000000  0.0000000
#>  [3,] 0.3535534  0.3535534  0.0000000  0.0000000
#>  [4,] 0.0000000  0.0000000  0.3535534  0.3535534
#>  [5,] 0.0000000  0.0000000  0.0000000  0.0000000
#>  [6,] 0.0000000  0.0000000  0.0000000  0.0000000
#>  [7,] 0.3535534 -0.3535534  0.0000000  0.0000000
#>  [8,] 0.0000000  0.0000000  0.3535534 -0.3535534
#>  [9,] 0.0000000  0.0000000  0.0000000  0.0000000
#> [10,] 0.1767767  0.1767767  0.1767767  0.1767767
#> [11,] 0.0000000  0.0000000  0.0000000  0.0000000
#> [12,] 0.1767767  0.1767767 -0.1767767 -0.1767767
#> [13,] 0.0000000  0.0000000  0.0000000  0.0000000
#> [14,] 0.1767767 -0.1767767  0.1767767 -0.1767767
#> [15,] 0.0000000  0.0000000  0.0000000  0.0000000
#> [16,] 0.1767767 -0.1767767 -0.1767767  0.1767767


# Visualizing the wavelet packet matrix
image(WP, main = "Wavelet Packet Matrix (WavPackMatWP)", xlab = "Columns", ylab = "Rows", col = heat.colors(256))
```

![](C:/Users/saraa/Desktop/WaveletMatrixProject/vignettes/WaveletMatrixOverview_files/figure-gfm/WavPackMatWP-example-1.png)<!-- -->

**Explanation**: WavPackMatWP() builds a matrix for wavelet packet
transformation up to a depth k0. It uses the low-pass filter (h) to
decompose the signal and create a matrix that captures different levels
of detail. This function is essential for performing multi-resolution
analysis.

**Supporting Functions:** - **getsubWP():** This function constructs
sub-matrices for the wavelet packet transformation at each decomposition
level. It applies low-pass (h) and high-pass (g) filters which generate
sub-matrices that represent transformations at different decomposition
levels.

- **Key Components of `getsubWP()`**:
  - **Iterative construction**: `getsubWP()` builds sub-matrices by
    applying filters at each level in reverse order.
  - **Filter applications**: Combines h and g matrices to decompose the
    signal at its current level.
- **Arguments in `getsubWP()`**:
  - **jstep**: Current decomposition level.
  - **h**: Low-pass filter.
  - **g**: High-pass filter.
  - **J**: Maximum number of decomposition levels.
  - **N**: Size of signal/matrix.

**getHGmatWP():** This function computes matrices using the input
filters (h and g), adjusting their elements with a modulus operation to
handle boundary conditions. This helps ensure that the transformation
correctly wraps around indices when applying the filters.

**Key Components of `getHGmatWP()`**:

- **Arguments in `getHGmatWP()`**:
  - **k**: Current decomposition level.
  - **h**: Low-pass filter.
  - **g**: High-pass filter.
  - **J**: Maximum number of decomposition levels.
  - **N**: Size of signal/matrix.

For a more detailed demonstration, see WavmatND_demo1 and WavmatND_demo2
in the demo folder.

## 3. `WavPackMat()`

**Description**: This function constructs an optimized orthogonal
wavelet transformation matrix using a low-pass filter. This matrix is
used for efficient wavelet-based signal processing.

### Example

``` r
# Low-pass filter
h <- c(0.5, 0.5)

# Parameters for wavelet transformation matrix
N <- 8 #Matrix size (2^3 in this case)
k0 <- 3 #Depth of wavelet decomposition
shift <- 2 #Shift parameter

# Construct the wavelet transformation matrix
WP_opt <- WavPackMat(h, N, k0, shift)
print(WP_opt)
#>             [,1]        [,2]        [,3]        [,4]
#>  [1,] 0.28867513  0.28867513  0.00000000  0.00000000
#>  [2,] 0.00000000  0.00000000  0.28867513  0.28867513
#>  [3,] 0.00000000  0.00000000  0.00000000  0.00000000
#>  [4,] 0.00000000  0.00000000  0.00000000  0.00000000
#>  [5,] 0.28867513 -0.28867513  0.00000000  0.00000000
#>  [6,] 0.00000000  0.00000000  0.28867513 -0.28867513
#>  [7,] 0.00000000  0.00000000  0.00000000  0.00000000
#>  [8,] 0.00000000  0.00000000  0.00000000  0.00000000
#>  [9,] 0.14433757  0.14433757  0.14433757  0.14433757
#> [10,] 0.00000000  0.00000000  0.00000000  0.00000000
#> [11,] 0.14433757  0.14433757 -0.14433757 -0.14433757
#> [12,] 0.00000000  0.00000000  0.00000000  0.00000000
#> [13,] 0.14433757 -0.14433757  0.14433757 -0.14433757
#> [14,] 0.00000000  0.00000000  0.00000000  0.00000000
#> [15,] 0.14433757 -0.14433757 -0.14433757  0.14433757
#> [16,] 0.00000000  0.00000000  0.00000000  0.00000000
#> [17,] 0.07216878  0.07216878  0.07216878  0.07216878
#> [18,] 0.07216878  0.07216878  0.07216878  0.07216878
#> [19,] 0.07216878  0.07216878 -0.07216878 -0.07216878
#> [20,] 0.07216878  0.07216878 -0.07216878 -0.07216878
#> [21,] 0.07216878 -0.07216878  0.07216878 -0.07216878
#> [22,] 0.07216878 -0.07216878  0.07216878 -0.07216878
#> [23,] 0.07216878 -0.07216878 -0.07216878  0.07216878
#> [24,] 0.07216878 -0.07216878 -0.07216878  0.07216878
#>              [,5]        [,6]        [,7]        [,8]
#>  [1,]  0.00000000  0.00000000  0.00000000  0.00000000
#>  [2,]  0.00000000  0.00000000  0.00000000  0.00000000
#>  [3,]  0.28867513  0.28867513  0.00000000  0.00000000
#>  [4,]  0.00000000  0.00000000  0.28867513  0.28867513
#>  [5,]  0.00000000  0.00000000  0.00000000  0.00000000
#>  [6,]  0.00000000  0.00000000  0.00000000  0.00000000
#>  [7,]  0.28867513 -0.28867513  0.00000000  0.00000000
#>  [8,]  0.00000000  0.00000000  0.28867513 -0.28867513
#>  [9,]  0.00000000  0.00000000  0.00000000  0.00000000
#> [10,]  0.14433757  0.14433757  0.14433757  0.14433757
#> [11,]  0.00000000  0.00000000  0.00000000  0.00000000
#> [12,]  0.14433757  0.14433757 -0.14433757 -0.14433757
#> [13,]  0.00000000  0.00000000  0.00000000  0.00000000
#> [14,]  0.14433757 -0.14433757  0.14433757 -0.14433757
#> [15,]  0.00000000  0.00000000  0.00000000  0.00000000
#> [16,]  0.14433757 -0.14433757 -0.14433757  0.14433757
#> [17,]  0.07216878  0.07216878  0.07216878  0.07216878
#> [18,] -0.07216878 -0.07216878 -0.07216878 -0.07216878
#> [19,]  0.07216878  0.07216878 -0.07216878 -0.07216878
#> [20,] -0.07216878 -0.07216878  0.07216878  0.07216878
#> [21,]  0.07216878 -0.07216878  0.07216878 -0.07216878
#> [22,] -0.07216878  0.07216878 -0.07216878  0.07216878
#> [23,]  0.07216878 -0.07216878 -0.07216878  0.07216878
#> [24,] -0.07216878  0.07216878  0.07216878 -0.07216878

# Visualize the matrix
image(WP_opt, main = "Wavelet Transformation Matrix", xlab = "Columns", ylab = "Rows", col = heat.colors(256))
```

![](C:/Users/saraa/Desktop/WaveletMatrixProject/vignettes/WaveletMatrixOverview_files/figure-gfm/WavPackMat-example-1.png)<!-- -->

**Explanation**: WavPackMat() uses the low-pass filter (h) to create an
orthogonal matrix for efficient decomposition and reconstruction.

**Supporting Functions:** **getsubW():** Constructs sub-matrices for the
wavelet transformation at each level. It applies the H and G matrices at
each iteration to create sub-matrices for each decomposition level,
which is essential for building the wavelet matrix used for signal
decomposition.

**Key Components of `getsubW()`**:

- **Arguments in `getsubW()`**:
  - **jstep**: Current decomposition level.
  - **h**: Low-pass filter.
  - **g**: High-pass filter.
  - **J**: Maximum number of decomposition levels.
  - **N**: Size of the matrix.

**getHGmat():** This function creates H and G matrices that are critical
for constructing wavelet transformations with circular boundary
conditions. The matrices ensure that the input signal is properly
processed even when indices wrap around.

**Key Components of `getsubW()`**:

- **Matrix Construction**: Creates the h and g matrices by applying the
  low-pass (h) and high-pass (g) filters.

- **Boundary Handling**: Uses modulus operations to wrap around indices,
  ensuring the matrix respects circular boundary conditions.

- **Arguments in `getsubW()`**:

  - **k**: Current decomposition level.
  - **h**: Low-pass filter.
  - **g**: High-pass filter.
  - **J**: Maximum number of decomposition levels.
  - **N**: Size of the matrix.

For a more detailed demonstration, see Wavmat_demo1 and Wavmat_demo2 in
the demo folder.

## 4. Demo

### Wavmat.R Demo

**Demo 1:** Image Processing on Lena Data

``` r

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
#> Min of D: 3.222453e-10
cat("Max of D:", max(D), "\n")
#> Max of D: 3.236664e-10

# Plot D with an adjusted zlim range
image(D, col = gray.colors(256), main = "Difference Image (D)", zlim = range(D), useRaster = TRUE)
```

![](C:/Users/saraa/Desktop/WaveletMatrixProject/vignettes/WaveletMatrixOverview_files/figure-gfm/WavPackMat-Lena%20Example-1.png)<!-- -->

``` r

# Scale D for better visualization
D_scaled <- D * 1000

# Checking the range of D_scaled
cat("Min of D_scaled:", min(D_scaled), "\n")
#> Min of D_scaled: 3.222453e-07
cat("Max of D_scaled:", max(D_scaled), "\n")
#> Max of D_scaled: 3.236664e-07

# Display using base R
image(D_scaled, col = gray.colors(256), main = "Scaled Difference Image (D)", useRaster = TRUE)

# Normalize D_scaled for imager display
D_scaled_normalized <- (D_scaled - min(D_scaled)) / (max(D_scaled) - min(D_scaled))
imager::display(imager::as.cimg(D_scaled_normalized), main = "Scaled and Normalized Difference Image (D)")



# Debugging Outputs
cat("Mean of A:", mean(A), "\n")
#> Mean of A: 254.5496
cat("Mean of B:", mean(B), "\n")
#> Mean of B: 46.36336
cat("Mean of C:", mean(C), "\n")
#> Mean of C: 254.5496
cat("Mean of D:", mean(D), "\n") #should be approximately 0
#> Mean of D: 3.229099e-10
```

![](C:/Users/saraa/Desktop/WaveletMatrixProject/vignettes/WaveletMatrixOverview_files/figure-gfm/WavPackMat-Lena%20Example-2.png)<!-- -->

**Explanation**

This demo illustrates the application of wavelet transformations to the
Lena image (standard image processing data file). The image is
preprocessed by inverting pixel intensities, cropped to a 256x256
portion for focus. Then, the padding to handle boundary effects during
transformation is placed. Using a wavelet filter, the wavelet
transformation matrix decomposes the image into multi-resolution
components,which highlight details like edges and textures while at the
same time preserving broader structures.

Then, we reconstruct with the inverse transformation demonstrates high
fidelity. This gives a minimal difference in images that quantify
transformation accuracy. This process showcases how wavelets are
utilized in tasks such as image compression, noise reduction, and
feature extraction.

**Demo 2:** Signal Denoising

``` r

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
#> Warning: Using `size` aesthetic for lines was deprecated in
#> ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> This warning is displayed once every 8 hours.
#> Call `lifecycle::last_lifecycle_warnings()` to see
#> where this warning was generated.
```

![](C:/Users/saraa/Desktop/WaveletMatrixProject/vignettes/WaveletMatrixOverview_files/figure-gfm/WavPackMat-Signal%20Denoising%20Example-1.png)<!-- -->

``` r



# Add noise to the signal
set.seed(1)# Set seed for reproducibility
signoi <- sig + 1/sqrt(SNR) * rnorm(N) # Add Gaussian noise to the signal


#Plot the noisy signal
ggplot2::ggplot(data.frame(t = t, signoi = signoi, sig = sig), ggplot2::aes(x = t)) +
  ggplot2::geom_point(ggplot2::aes(y = signoi), color = "red", size = 2) +
  ggplot2::geom_line(ggplot2::aes(y = sig), color = "green", size = 1) +
  ggplot2::theme_minimal(base_size = 16) +
  ggplot2::ggtitle("Noisy Bumps Signal")
```

![](C:/Users/saraa/Desktop/WaveletMatrixProject/vignettes/WaveletMatrixOverview_files/figure-gfm/WavPackMat-Signal%20Denoising%20Example-2.png)<!-- -->

``` r



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
```

![](C:/Users/saraa/Desktop/WaveletMatrixProject/vignettes/WaveletMatrixOverview_files/figure-gfm/WavPackMat-Signal%20Denoising%20Example-3.png)<!-- -->

``` r



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
```

![](C:/Users/saraa/Desktop/WaveletMatrixProject/vignettes/WaveletMatrixOverview_files/figure-gfm/WavPackMat-Signal%20Denoising%20Example-4.png)<!-- -->

``` r

# Plot the thresholded coefficients
ggplot2::ggplot(data.frame(index = 1:N, swt = swt), ggplot2::aes(x = index, y = swt)) +
  ggplot2::geom_line(size = 1) +
  ggplot2::theme_minimal(base_size = 16) +
  ggplot2::ggtitle("Thresholded Wavelet Coefficients")
```

![](C:/Users/saraa/Desktop/WaveletMatrixProject/vignettes/WaveletMatrixOverview_files/figure-gfm/WavPackMat-Signal%20Denoising%20Example-5.png)<!-- -->

``` r


#  Return the signal to the time domain using the inverse transformation
a <- as.vector(t(WP) %*% swt) # Reconstruct the signal

#Plot the denoised signal and the noisy signal for comparison
ggplot2::ggplot(data.frame(t = t, signoi = signoi, a = a), ggplot2::aes(x = t)) +
  ggplot2::geom_point(ggplot2::aes(y = signoi), color = "red", size = 2) +
  ggplot2::geom_line(ggplot2::aes(y = a), color = "black", size = 1) +
  ggplot2::theme_minimal(base_size = 16) +
  ggplot2::ggtitle("Denoised Signal")
```

![](C:/Users/saraa/Desktop/WaveletMatrixProject/vignettes/WaveletMatrixOverview_files/figure-gfm/WavPackMat-Signal%20Denoising%20Example-6.png)<!-- -->

**Explanation:**

This demo demonstrates thow wavelet transformations are utilized in
signal denoising by using a synthetic “bumps” signal. The signal is
constructed with multiple localized features (bumps) of varying heights
and widths. This is done to simulate a plausible real-world scenario
where signals contain sharp transitions. We then add noise to this
simulated signal in order to simulate environmental interference (which
is common in signal processing).

Using a wavelet packet transformation matrix constructed with a wavelet
filter, the noisy signal is decomposed into wavelet coefficients. These
coefficients represent the signal at different resolution levels,which
allow for manipulation. A soft-thresholding technique is applied to
remove small coefficients while preserving the larger coefficients that
capture the signal’s key features.

The denoised signal is reconstructed by applying the inverse
transformation to the thresholded coefficients, restoring the original
signal with minimal noise. This approach highlights how wavelets are
used in removing noise while retaining critical signal characteristics.

### WavmatWP.R

**Demo 1:** Transformation and Reconstruction

``` r
library(WaveletMatrixProject)

# Define the filter ( coefficients define the specific wavelet used for the transformation)
filt <- c(-0.07576571478934, -0.02963552764595,
          0.49761866763246, 0.80373875180522,
          0.29785779560554, -0.09921954357694,
          -0.01260396726226, 0.03222310060407)
# Print to verify
print(filt)
#> [1] -0.07576571 -0.02963553  0.49761867  0.80373875
#> [5]  0.29785780 -0.09921954 -0.01260397  0.03222310


# Set parameters for the wavelet packet transformation matrix
N <- 8  # Size of the matrix
k0 <- 3  # Levels of decomposition
shift <- 0  # Shift parameter

# Generate the wavelet packet transformation matrix
WP <- WavPackMatWP(filt, N, k0, shift)

# Print the matrix
print(WP)
#>                [,1]          [,2]         [,3]
#>  [1,] -0.0437433558 -0.0171100799  0.287300272
#>  [2,] -0.0072769039  0.0186040158 -0.043743356
#>  [3,]  0.1719682785 -0.0572844302 -0.007276904
#>  [4,]  0.2873002717  0.4640387847  0.171968278
#>  [5,]  0.0186040158  0.0072769039 -0.057284430
#>  [6,] -0.0171100799  0.0437433558  0.018604016
#>  [7,]  0.4640387847 -0.2873002717 -0.017110080
#>  [8,] -0.0572844302 -0.1719682785  0.464038785
#>  [9,]  0.3200188889  0.3452051694  0.114125293
#> [10,] -0.0491598348 -0.0029733651  0.103262717
#> [11,] -0.1041385041 -0.1328506899 -0.041793716
#> [12,]  0.0156701527  0.0001591509 -0.013543504
#> [13,]  0.1839701650 -0.2830312452  0.368240514
#> [14,] -0.0379549046  0.0656432345 -0.026008710
#> [15,] -0.0285952781  0.0805137184 -0.138131313
#> [16,]  0.0174812461 -0.0251764987  0.005439774
#> [17,] -0.0227895822 -0.0260665992 -0.011707030
#> [18,]  0.0096923919  0.0110861047  0.004978991
#> [19,]  0.0074257350  0.0100608110  0.003567900
#> [20,] -0.0031581594 -0.0042788552 -0.001517425
#> [21,] -0.0128138174  0.0194986927 -0.027129224
#> [22,]  0.0054497068 -0.0082927791  0.011538038
#> [23,]  0.0016484757 -0.0053540606  0.010304407
#> [24,] -0.0007010955  0.0022770779 -0.004382456
#>                [,4]          [,5]          [,6]
#>  [1,]  0.4640387847  0.1719682785 -0.0572844302
#>  [2,] -0.0171100799  0.2873002717  0.4640387847
#>  [3,]  0.0186040158 -0.0437433558 -0.0171100799
#>  [4,] -0.0572844302 -0.0072769039  0.0186040158
#>  [5,] -0.1719682785  0.4640387847 -0.2873002717
#>  [6,]  0.0072769039 -0.0572844302 -0.1719682785
#>  [7,]  0.0437433558  0.0186040158  0.0072769039
#>  [8,] -0.2873002717 -0.0171100799  0.0437433558
#>  [9,] -0.0714351748 -0.0491598348 -0.0029733651
#> [10,]  0.2174504353  0.3200188889  0.3452051694
#> [11,]  0.0299538457  0.0156701527  0.0001591509
#> [12,] -0.0410678778 -0.1041385041 -0.1328506899
#> [13,] -0.1963332067 -0.0379549046  0.0656432345
#> [14,] -0.0745258474  0.1839701650 -0.2830312452
#> [15,]  0.0757847965  0.0174812461 -0.0251764987
#> [16,]  0.0126835549 -0.0285952781  0.0805137184
#> [17,] -0.0010319213 -0.0057592986 -0.0100050582
#> [18,]  0.0004388753  0.0024494253  0.0042551436
#> [19,] -0.0010524063  0.0018989392  0.0039250421
#> [20,]  0.0004475876 -0.0008076174 -0.0016693174
#> [21,]  0.0170839386 -0.0025763724  0.0034142737
#> [22,] -0.0072657860  0.0010957292 -0.0014520880
#> [23,] -0.0061177731 -0.0004770429 -0.0004785511
#> [24,]  0.0026018842  0.0002028860  0.0002035274
#>                [,7]          [,8]
#>  [1,] -0.0072769039  0.0186040158
#>  [2,]  0.1719682785 -0.0572844302
#>  [3,]  0.2873002717  0.4640387847
#>  [4,] -0.0437433558 -0.0171100799
#>  [5,] -0.0171100799  0.0437433558
#>  [6,]  0.4640387847 -0.2873002717
#>  [7,] -0.0572844302 -0.1719682785
#>  [8,]  0.0186040158  0.0072769039
#>  [9,]  0.1032627174  0.2174504353
#> [10,]  0.1141252933 -0.0714351748
#> [11,] -0.0135435041 -0.0410678778
#> [12,] -0.0417937155  0.0299538457
#> [13,] -0.0260087102 -0.0745258474
#> [14,]  0.3682405144 -0.1963332067
#> [15,]  0.0054397737  0.0126835549
#> [16,] -0.1381313127  0.0757847965
#> [17,] -0.0112059369 -0.0143582686
#> [18,]  0.0047658764  0.0061065606
#> [19,]  0.0022647121  0.0022238391
#> [20,] -0.0009631803 -0.0009457971
#> [21,] -0.0089424334  0.0114649423
#> [22,]  0.0038032101 -0.0048760312
#> [23,]  0.0036814460 -0.0032069010
#> [24,] -0.0015657162  0.0013638926

# This vector will be transformed using the wavelet packet matrix
y <- c(1, 0, -3, 2, 1, 0, 1, 2)

# Perform forward wavelet packet transformation
d <- sqrt(k0) * WP %*% y #scaled by `sqrt(k0)` for normalization

# Perform inverse wavelet packet transformation
yy <- (1 / sqrt(k0)) * t(WP) %*% d #inverse scaled by `1 / sqrt(k0)` for normalization

# Print the transformed and inverse-transformed vectors
print("d vector:")
#> [1] "d vector:"
print(d)
#>               [,1]
#>  [1,]  0.388555815
#>  [2,]  0.752459498
#>  [3,]  2.429446355
#>  [4,] -0.742034544
#>  [5,]  0.659800794
#>  [6,] -0.391815022
#>  [7,]  0.381464730
#>  [8,] -3.477877627
#>  [9,]  0.560797197
#> [10,]  0.636055330
#> [11,]  0.001976663
#> [12,] -0.193746543
#> [13,] -2.643859686
#> [14,]  0.087579066
#> [15,]  1.014386611
#> [16,]  0.019696946
#> [17,] -0.061339036
#> [18,]  0.026087445
#> [19,]  0.005592018
#> [20,] -0.002378281
#> [21,]  0.197718467
#> [22,] -0.084089513
#> [23,] -0.077439456
#> [24,]  0.032934941
print("yy (inverse result):")
#> [1] "yy (inverse result):"
print(yy)
#>             [,1]
#> [1,]  0.11690491
#> [2,]  0.59676205
#> [3,] -1.57160394
#> [4,]  1.07113362
#> [5,]  0.52367945
#> [6,]  0.01324524
#> [7,]  0.47234944
#> [8,]  0.82188673
# Verify orthogonality of the wavelet packet transformation matrix
trace1 <- sum(diag(WP %*% t(WP))) #should about 4
trace2 <- sum(diag(t(WP) %*% WP)) #should about 4
cat("Trace of WP * WP':", trace1, "\n")
#> Trace of WP * WP': 4.004726
cat("Trace of WP' * WP:", trace2, "\n")
#> Trace of WP' * WP: 4.004726
```

**Explination**

This demo demonstrates the construction, application, and validation of
a wavelet packet transformation matrix. Wavelet packet transformations
analyze both approximation and detail coefficients at every
decomposition level. This process provides a finer resolution of the
signal’s frequency content.

Using WavPackMatWP, a sample signal is transformed into the wavelet
domain, scaling the coefficients for normalization. These coefficients
highlight the signal’s characteristics at different resolution levels.

Then, the inverse transformation is applied to reconstruct the original
signal from its wavelet packet coefficients. Doing this, we verify the
fidelity of the transformation-reconstruction process. The orthogonality
of the transformation matrix is ensured by calculating the trace of the
product matrices (WP \* WP’ and WP’ \* WP). Doing this, we ensure that
energy is preserved and no information is lost.

This demo illustrates showcasing wavelet abilities to represent signals
in a compact and interpretable form.

**Demo 2:** Analysis of a Doppler Signal

``` r
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
#> Summary of d in R:
cat(sprintf('Min: %f, Max: %f, Mean: %f\n', min(d), max(d), mean(d)))
#> Min: -1.248575, Max: 1.230026, Mean: 0.008686
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
```

![](C:/Users/saraa/Desktop/WaveletMatrixProject/vignettes/WaveletMatrixOverview_files/figure-gfm/WavPackMatWP-Analysis%20of%20a%20Doppler%20Signal-1.png)<!-- -->

``` r
#plots reconstruction error, which ideally is minimal.
```

**Explanation**

In this demo, we illustrate the application of wavelet packet
transformations to analyze and reconstruct a Doppler signal. The Doppler
signal is generally characterized by its sharp transitions, thus it is a
great example for showcasing the effectiveness of wavelet
transformations in capturing both global and localized signal features.

In the demo, using Haar wavelet filter coefficients, we construct a
wavelet packet transformation matrix which is tailored to decompose the
signal into multiple resolution levels. The Doppler signal is then
reshaped into a column vector and transformed into the wavelet domain.
This yields wavelet packet coefficients that show the signal’s frequency
components across different scales.

We then apply the inverse transformation to reconstruct the original
signal. We do this in order to demonstrate the fidelity of the
transformation process. The minimal difference between the original and
reconstructed signals is displayed to confirm the accuracy and
robustness of the wavelet packet.

By visualizing the original signal, its wavelet coefficients, the
reconstructed signal, and the reconstruction error, this demo shows how
wavelet packet transformations are used in signal analysis.

### WavmatND.R Demo

**Demo 1:** Denoising a Noisy Doppler Signal Using Non-Decimated Wavelet
Transform

``` r
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
#> [1] 0.09916485


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
#> Variance of transformed signal : 0.293545
cat(sprintf("Variance of last detail level: %f\n", var(tsn[((k - 1) * n + 1):(k * n)]))) # Should be approximately 0.0132
#> Variance of last detail level: 0.013493

# Extract detail levels and apply hard thresholding for denoising
temp <- tsn[(n + 1):length(tsn)]# Detail coefficients
threshold <- sqrt(2 * log(n) * var(tsn[(length(tsn) - n + 1):length(tsn)])) # Threshold value
temp[abs(temp) < threshold] <- 0 # Zero out small coefficients

# Reconstruct the signal using the inverse transformation
rs <- t(W) %*% T %*% c(tsn[1:n], temp)
trace <- sum(diag(T))# Verify the trace of T (should be n)
cat(sprintf("Trace of T: %f\n", trace)) # Should be 1000
#> Trace of T: 1000.000000
cat(sprintf("Variance of reconstructed signal: %f\n", var(rs))) # Should be approximately 0.0777
#> Variance of reconstructed signal: 0.852038

cat("Dimension of T: ", dim(T)) # 5000 5000
#> Dimension of T:  5000 5000
cat("Dimension of T: ", dim(W)) # 5000 1000
#> Dimension of T:  5000 1000
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
```

![](C:/Users/saraa/Desktop/WaveletMatrixProject/vignettes/WaveletMatrixOverview_files/figure-gfm/WavPackMatND-%20Denoising%20a%20Noisy%20Doppler%20Signal-1.png)<!-- -->
**Explanation:**

This demo showcases the use of a non-decimated wavelet transform
(denoted as NDWT) with a Haar wavelet filter in order to denoise a noisy
Doppler signal. The signal is decomposed into multiple resolution
levels, which help preserve its original length (a technique useful in
high-resolution analysis). By hard thresholding, we remove small
coefficients that afre contributed to noise while retaining key
features.

Then, the denoised signal is reconstructed using the inverse
transformation and a weight matrix. This reduces noise while at the same
time helping us maintain signal fidelity. Then we visualize the
original, noisy, and denoised signals, along with reconstruction error,
in order to demonstrate the NDWT’s effectiveness for noise reduction.

**Demo 2:** Denoising a Noisy Doppler Signal Using Non-Decimated Wavelet
Transform with Weighted Reconstruction

``` r
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
```

![](C:/Users/saraa/Desktop/WaveletMatrixProject/vignettes/WaveletMatrixOverview_files/figure-gfm/WavPackMatND-Denoising%20a%20Noisy%20Doppler%20Signal%20with%20Weighted%20Reconstruction-1.png)<!-- -->

**Explanation** In this demo, we again apply a NDWT using a Haar wavelet
filter to denoise a noisy Doppler. The signal is then, again, decomposed
into multiple resolution levels without down-sampling. This preserves
the original signal length and enables for precise frequency analysis.

We then apply a weighted reconstruction matrix which is used to balance
the contributions from different decomposition levels. We do this in
order to ensure the accuracy of reconstructions. Afterwards, hard
thresholding is applied to the wavelet detail coefficients. This process
removes noise while keeping significant features of the signal.

The denoised signal is the reconstructed using the inverse
transformation. We visualize the original, noisy, and transformed
signals. This approach demonstrates the NDWT’s power in denoising tasks,
providing high-resolution analysis.

\##References

G. P. Nason, Wavelet methods in statistics with R (Springer New York,
New York, NY, 2008).

Reveille’s wavelets, (available at
<https://people.tamu.edu/~brani/wavelet/>).