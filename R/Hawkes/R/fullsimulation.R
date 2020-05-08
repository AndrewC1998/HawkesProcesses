Hawkes.sim <- function(M, mu, Y, dist, delta, N, t){
  # M is the number of dimensions
  # mu is the background intensity vector
  # Y is matrix of inital sizes of self-excited jumps
  # dist is the distribution function that Y follows
  # delta is a matrix of the rate of exponential decay
  # N is the vector of the number of events attributed to process i observed at and before time 0
  # t is the maturity time

  r <- c(0)
  lambda <- matrix(c(0), nrow = M, ncol = M)
  for(m in 1:M){
    for(i in 1:M){
      lambda[i,m] <- Y[i,m]
    }
  }
  a <- list()
  j <- 0
  Nfull <- matrix(N, ncol = M, nrow = 1)
  ttmp <- c()
  while(r[j+1] < t){
    j <- j + 1
    a[[j]] <- matrix(c(0), nrow = M + 1, ncol = M)
    for(m in 1:M){
      a[[j]][1,m] <- rexp(1, mu[m])
      for(i in 2:(M+1)){
        # CDF <- 1 - exp(-(1/delta[i,m])*(lambda[i,m])*(1 - exp(-delta[i,m])))
        u <- runif(1)
        tmp <- 1 - exp(-(1/delta[i-1,m])*(lambda[i-1,m]))
        if(u < tmp){
          a[[j]][i,m] <- -(1/delta[i-1,m])*log(1 + ((delta[i-1,m])/(lambda[i-1,m]))*log(1-u))
        }else{
          a[[j]][i,m] <- Inf
        }
      }
    }
    r[j+1] <- r[j] + min(a[[j]])
    star <- which(a[[j]] == min(a[[j]]), arr.ind = TRUE)
    mstar <- as.numeric(star[,2])
    istar <- as.numeric(star[,1])
    zjxj <- c(mstar, istar)
    for(m in 1:M){
      if(dist[mstar,m]=="Exp"){
        ymstar <- rexp(1, r[j+1])
      }else{
        stop("Distribution not currently supported.")
      }
      for(i in 1:M){
        lambda[i,m] <- lambda[i,m]*exp(-delta[i,m]*(a[[j]][istar,mstar])) + ymstar*(i==zjxj[1])
      }
      N[m] <- N[m] + 1*(m == zjxj[1])
    }
    Nfull <- rbind(Nfull, N)
    k <- Nfull[j+1, mstar]
    ttmp[k] <- r[j+1]
  }
  returnlist <- list(); returnlist$r <- r[-length(r)]
  returnlist$N <- Nfull[-length(Nfull[,1]),]
  returnlist$t <- ttmp
  return(returnlist)
}
