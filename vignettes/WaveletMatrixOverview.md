WaveletMatrixOverview
================

``` r
library(WaveletMatrixProject)
```

# Introduction

This vignette provides a comprehensive overview of the functions
available in the `WaveletMatrixProject` package. Each function is
described in detail, with examples demonstrating their usage. These
functions support the construction of wavelet and wavelet packet
transformation matrices, useful for signal processing and data analysis.

## 1. `WavmatND()`

**Description**: Constructs a wavelet matrix using a high-pass filter,
specified parameters, and an iterative process.

### Example

``` r
# Define parameters
hf <- c(0.5, -0.5)
N <- 4
k <- 2
shift <- 1

# Construct the wavelet matrix
W <- WavmatND(hf, N, k, shift)
print(W)
#>        [,1]  [,2]  [,3]  [,4]
#>  [1,] -0.25  0.25  0.25 -0.25
#>  [2,] -0.25 -0.25  0.25  0.25
#>  [3,]  0.25 -0.25 -0.25  0.25
#>  [4,]  0.25  0.25 -0.25 -0.25
#>  [5,] -0.25  0.25 -0.25  0.25
#>  [6,]  0.25 -0.25  0.25 -0.25
#>  [7,] -0.25  0.25 -0.25  0.25
#>  [8,]  0.25 -0.25  0.25 -0.25
#>  [9,]  0.00 -0.50 -0.50  0.00
#> [10,]  0.00  0.00 -0.50 -0.50
#> [11,] -0.50  0.00  0.00 -0.50
#> [12,] -0.50 -0.50  0.00  0.00
```

**Explanation**: The WavmatND() function starts with an identity matrix
and iteratively applies a high-pass filter to construct the wavelet
matrix. This matrix can be used for signal processing tasks where
transformations are needed to capture specific frequency components in a
signal.

## 2. `dilate_filter()`

**Description**: Expands the input filter by inserting zeros between its
elements, effectively “dilating” it. This helps increase the filter
length, useful in various stages of wavelet matrix construction.

### Example

``` r
# Define a filter
filt <- c(0.5, -0.5)

# Dilate the filter with a step size of 2
dilated_filt <- dilate_filter(filt, 2)
print(dilated_filt)
#> [1]  0.5  0.0  0.0 -0.5
```

**Explanation**: dilate_filter() takes an input vector and places zeros
between its elements. This is important for adjusting filter lengths in
wavelet processing, where filters of specific lengths are required for
different transformation levels.

## 3. `imbed()`

**Description**: Embeds a vector by appending zeros to double its
length. This operation is typically used in preparation for
transformations that require vectors of greater length. construction.

### Example

``` r
# Original vector
a <- c(1, 2, 3)

# Embed the vector
imbedded_a <- imbed(a)
print(imbedded_a)
#> [1] 1 2 3 0 0 0
```

**Explanation**: imbed() increases the length of a vector by appending
zeros, which can be required for padding in signal processing or when
preparing data for further multi-resolution analysis.

## 4. `repeating()`

**Description**: Repeats the input vector a specified number of times,
creating a longer vector by replication. This is useful in creating
extended signals or patterns.

### Example

``` r
# Original vector
a <- c(1, 2)

# Repeat the vector 3 times
repeated_a <- repeating(a, 3)
print(repeated_a)
#> [1] 1 2 1 2 1 2
```

**Explanation**: The repeating() function replicates the input vector n
times, extending its length. This is beneficial in scenarios where
repeated patterns of a signal are required for signal synthesis or
analysis.

## 5. `WavPackMatWP()`

**Description**: Constructs a wavelet packet transformation matrix using
a low-pass filter. This matrix represents a transformation that
decomposes signals into components at different levels of detail.

### Example

``` r
# Low-pass filter
h <- c(0.5, 0.5)
N <- 8
k0 <- 2
shift <- 2

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
```

**Explanation**: WavPackMatWP() builds a matrix for wavelet packet
transformation up to a depth k0. It uses the low-pass filter (h) to
decompose the signal and create a matrix that captures different levels
of detail. This function is essential for performing multi-resolution
analysis.

## 6. `getsubWP()`

**Description**: Constructs sub-matrices for wavelet packet
transformation at each level of decomposition. This function is used
internally by WavPackMatWP().

### Example

This function is an essential part of WavPackMatWP() to iteratively
create sub-matrices for each decomposition level.

**Explanation**: getsubWP() generates sub-matrices that represent
transformations at different decomposition levels. It iterates backward
through the levels, applying filter matrices to create detailed
sub-matrices that contribute to the overall transformation matrix.

## 7. `getHGmatWP()`

**Description**: Constructs the H and G matrices for wavelet packet
transformation, applying circular boundary conditions. This function is
used internally by WavPackMatWP().

### Example

This function supports WavPackMatWP() by creating the necessary H and G
filter matrices during matrix construction.

**Explanation**: getHGmatWP() computes matrices using the input filters
(h and g), adjusting their elements with a modulus operation to handle
boundary conditions. This ensures that the transformation correctly
wraps around indices when applying the filters, which is critical in
wavelet packet analysis.

## 8. `WavPackMat()`

**Description**: Constructs an optimized orthogonal wavelet
transformation matrix using a low-pass filter. This matrix is used for
efficient wavelet-based signal processing.

### Example

``` r
# Low-pass filter
h <- c(0.5, 0.5)
N <- 8
k0 <- 3
shift <- 2

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
```

**Explanation**: WavPackMat() builds an orthogonal wavelet
transformation matrix optimized for computational efficiency. This
matrix is essential for analyzing signals with high precision and low
computational cost, making it suitable for large-scale data processing.

## 9. `getsubW()`

**Description**: Constructs sub-matrices for the wavelet transformation
at each level. This function is used internally by WavPackMat().

### Example

The function getsubW() is called within WavPackMat() to build
sub-matrices as part of the overall wavelet matrix construction process.

**Explanation**: This function iteratively applies the H and G matrices
to create sub-matrices for each decomposition level, essential for
building the wavelet matrix that is used for signal decomposition.

## 10. `getHGmat()`

**Description**: Constructs sub-matrices for the wavelet transformation
at each level. This function is used internally by WavPackMat().

### Example

getHGmat() is used by getsubW() to build H and G matrices for each level
of decomposition during wavelet matrix construction.

**Explanation**: getHGmat() creates H and G matrices that are critical
for constructing wavelet transformations with circular boundary
conditions. These matrices ensure that the input signal is properly
processed even when indices wrap around.
