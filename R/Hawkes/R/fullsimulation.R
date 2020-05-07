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
          a[[j]][i,m] <- Inf
        }else{
          a[[j]][i,m] <- Inf
        }
      }
    }
    r[j+1] <- r[j] + min(a[[j]])
    mstar <- ceiling(which.min(a[[j]])/ncol(a[[j]]))
    istar <- ceiling(which.min(a[[j]])/nrow(a[[j]]))
    zjxj <- c(mstar, istar)
    for(m in 1:M){
      if(dist[mstar,m]=="Exp"){
        ymstar <- rexp(1, r[j+1])
      }else{
        stop("Distribution not currently supported.")
      }
      for(i in 1:M){
        lambda[i,m] <- lambda[i,m]*exp(-delta[i,m]*min(a[[j]])) + ymstar*(i==zjxj[1])
      }
      N[m] <- N[m] + 1*(m = zjxj[1])
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


# test stuff
# For M = 2
# M <- 2; m <- c(0.5,1); y <- matrix(c(1,0,0,1), ncol = M)
# d <-  matrix(c("Exp"), ncol = M, nrow = M); del <- matrix(c(1), ncol = M, nrow = M)
# n <- c(1,1)

# test <- Hawkes.sim(M = M, mu = m, Y = y, dist = d, delta = del, N = n, t = 5)

# For M = 4
# M <- 4; m <- c(0.5,1,0.2,0.4); y <- matrix(c(0), ncol = M, nrow = M)
# d <-  matrix(c("Exp"), ncol = M, nrow = M); del <- matrix(c(0.2), ncol = M, nrow = M)
# n <- c(1,1,1,1)

# test <- Hawkes.sim(M = M, mu = m, Y = y, dist = d, delta = del, N = n, t = 1000)
# plot(1:length(test$r), test$r, type = "l", col = "blue")
