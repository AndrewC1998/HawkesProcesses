Hawkes.sim <- function(M, mu, Y, dist, delta, N, params, t, paramsfunc){
  # Process simulates an M-dimensional Hawkes process
  # Inputs:
  # M is the number of dimensions
  # mu is the background intensity vector
  # Y is matrix of inital sizes of self-excited jumps
  # dist is the distribution function that Y follows
  # delta is a matrix of the rate of exponential decay
  # N is the vector of the number of events attributed to process i observed at and before time 0
  # params is for the specific parameter choices of a given dist
  # t is the maturity time
  # paramsfunc allows for manual dist must write dist[i,m] = "Manual". Note must be defined as external function

  r <- c(0); lambda <- list()
  lambda[[1]] <- matrix(c(0), nrow = M, ncol = M)
  for(m in 1:M){
    for(i in 1:M){
      lambda[[1]][i,m] <- Y[i,m]
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
        tmp <- 1 - exp(-(1/delta[i-1,m])*(lambda[[1]][i-1,m]))
        if(u < tmp){
          a[[j]][i,m] <- -(1/delta[i-1,m])*log(1 + ((delta[i-1,m])/(lambda[[j]][i-1,m]))*log(1-u))
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
    lambda[[j+1]] <- matrix(c(0), nrow = M, ncol = M)
    for(m in 1:M){
      if(dist[mstar,m]=="Exp"){
        # requires params to be MxM matrix
        ymstar <- rexp(1, params[[1]][mstar,m])
      }else if(dist[mstar,m]=="Gamma"){
        # requires param to be list with each an MxM matrix
        ymstar <- rgamma(1, params[[2]][mstar,m], params[[3]][mstar,m])
      }else if(dist[mstar,m]=="Normal"){
        # requires param to be list with each an MxM matrix
        ymstar <- rnorm(1, params[[4]][mstar,m], params[[5]][mstar,m])
      }else if(dist[mstar,m]=="Manual"){
        ymstar <- paramsfunc(params)
      }else{
        stop("Distribution not currently supported.")
      }
      for(i in 1:M){
        lambda[[j+1]][i,m] <- lambda[[j]][i,m]*exp(-delta[i,m]*(a[[j]][istar,mstar])) + ymstar*(i==zjxj[1])
      }
      N[m] <- N[m] + 1*(m == zjxj[1])
    }
    Nfull <- rbind(Nfull, N)
    k <- Nfull[j+1, mstar]
    ttmp[k] <- r[j+1]
  }
  # Create matrix of intensities
  intensity <- matrix(c(0), ncol = M, nrow = j)
  for(m in 1:M){
    for(t in 1:j){
      intensity[t,m] <- mu[m] + sum(lambda[[t]][,m])
    }
  }
  # Create summary statistic to return
  rl <- list(); rl$r <- r[-length(r)]
  rl$N <- Nfull[-length(Nfull[,1]),]
  rl$t <- ttmp; rl$intensity <- intensity
  return(rl)
  # Output:
  # intensity is the intensity matrix. Each column is intensities for process m
  # r is the event times for all point processes
  # N is the M corresponding count processes. Each column is one process
  # t is the event times
}
