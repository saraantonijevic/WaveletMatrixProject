cat('Doppler Wavelet Packets Forward and Inverse (Demo 1 of WavMatPackWP)\n')


lw <- 1.0
par(cex.axis = 1.6)
fs <- 15
msize <- 10

nl <- 10      # J level, V_J
n <- 2^nl    #
level <- 5   # levels of decomposition
shift <- 3   # standard shift

filt <- c(sqrt(2)/2, sqrt(2)/2)


# SYMMLET 4
filt <- c(-0.07576571478934, -0.02963552764595,
          0.49761866763246, 0.80373875180522,
          0.29785779560554, -0.09921954357694,
          -0.01260396726226, 0.03222310060407)

WP <- WavmatWP(filt, n, level, shift) # size: [level*n, n]


cat("Trace of WP * WP':", sum(diag(WP %*% t(WP))), "\n")  # should be n. level*n entries on the diagonal of size 1/level.
cat("Trace of WP' * WP:", sum(diag(t(WP) %*% WP)), "\n")  # should be n. n entries of 1 each.

