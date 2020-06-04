Hawkes.mle <- function(alpha, beta, M, t, mu, threshold = 10){
  # alpha min max decisions
  # beta min max decisions
  # M number of dimensions
  # t is the data we have available
  # mu is means

  a <- seq(alpha[1], alpha[2], length.out = threshold)
  b <- seq(beta[1], beta[2], length.out = threshold)
  tmp <- matrix(0, nrow = threshold, ncol = threshold)

  for(i in 1:threshold){
    for(j in 1:threshold){
      tmp[i,j] <- Hawkes.ll(t, mu, as.matrix(a[i]), as.matrix(b[j]))
    }
  }
  maxval <- max(tmp)
  maxpoint <- which(tmp == maxval, arr.ind = TRUE)
  print(tmp)
  return(paste(c('alpha:'), a[maxpoint[1,1]], c(', beta:'), a[maxpoint[1,2]], c(', likelihood:'), maxval))
}
