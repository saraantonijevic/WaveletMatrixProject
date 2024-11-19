WavmatND <- function(hf, N, k, shift) {
  # hf = as.vector(hf) # Not needed in R as it is already a vector
  gf <- rev(Conj(hf) * (-1)^(1:length(hf)))
  W <- matrix(nrow = 0, ncol = N)
  hmatold <- diag(N)
  h <- c(hf, rep(0, N))
  g <- c(gf, rep(0, N))

  for (i in 1:k) {
    gmat <- matrix(0, nrow = N, ncol = N)
    hmat <- matrix(0, nrow = N, ncol = N)

    for (jj in 1:N) {
      for (ii in 1:N) {
        modulus <- (N + ii - jj - shift) %% N + 1
        modulus <- modulus + (modulus == 0) * N
        hmat[ii, jj] <- h[modulus]
        gmat[ii, jj] <- g[modulus]
      }
    }

    W <- rbind(t(gmat) %*% hmatold, W)
    smooth <- t(hmat) %*% hmatold
    hmatold <- smooth
    h <- c(dilate_filter(hf, 2^(i) - 1), rep(0, N))
    g <- c(dilate_filter(gf, 2^(i) - 1), rep(0, N))
  }

  W <- rbind(smooth, W)
  return(W)
}

dilate_filter <- function(filt, k) {
  newlength <- (k + 1) * length(filt) - k
  filtd <- rep(0, newlength)
  filtd[seq(1, newlength, by = k + 1)] <- filt
  return(filtd)
}

imbed <- function(a) {
  b <- c(a, rep(0, length(a)))
  return(b)
}

repeat <- function(a, n) {
  b <- rep(a, n)
  return(b)
}

