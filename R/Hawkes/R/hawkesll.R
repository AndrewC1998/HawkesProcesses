Hawkes.ll <- function(t, mu, alpha, beta){
  # t is matrix of times each col is M=1,2,...
  # mu is vector of base intensities
  # alpha is a square matrix of alpha coefficients
  # beta is a qsuare matrix of beta coefficients

  # Create recursive formula
  R <- function(i,m,d){
    tmp <- 0
    if(d == 1){
      return(0)
    }else{
      for(k in 2:length(t[,m])){
        if(t[k,i] < t[d,m]){
          tmp <- tmp + exp(-(beta[i,m])*(t[d,m]-t[k,i]))
        }
      }
      return(tmp)
    }
  }
  # create vector to store likelihoods
  ll <- c()
  M <- length(mu)

  for(m in 1:M){
    # First for loop helps us find each dimension
    # Find the first term
    tmpleft <- 0
    maxt <- max(t[,m])
    wmax <- which.max(t[,m])
    for(d in 1:wmax){
      # This one for data length
      tmp1 <- 0
      for(i in 1:M){
        tmp1 <- tmp1 + alpha[i,m]*R(i,m,d)
      }
      tmpleft <- tmpleft + log(mu[m] + tmp1)
    }

    tmpright <- 0
    tmp2 <- c()
    for(i in 1:M){
      tmp3 <- 0
      for(d in 1:wmax){
        tmpab <- alpha[i,m]*beta[i,m]
        tmp3 <- tmpab*(1 - exp(-(beta[i,m])*(maxt - t[d,i])))
      }
      tmp2[i] <- tmp3
    }
    tmpright <- sum(tmp2)

    # Sum all terms for mth likelihood
    ll[m] <- tmpleft - mu[m]*maxt - tmpright
  }
  return(sum(ll))
}



hawkes::likelihoodHawkes(c(0.5,0.3), matrix(c(1),ncol=2,nrow = 2), matrix(c(1),ncol=2,nrow = 2), matrix(cbind(c(1:5), c(1:5)), ncol = 2))

Hawkes.ll(matrix(cbind(c(1:5), c(1:5)), ncol = 2), c(0.5,0.3), matrix(c(1),ncol=2,nrow = 2), matrix(c(1),ncol=2,nrow = 2))
