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


