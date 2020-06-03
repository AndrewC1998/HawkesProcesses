Hawkes.sim2 <- function(mu, alpha, beta, horizon, burnin = 0, separate = FALSE){

  M <- length(mu)

  if(max(abs(eigen(a/b)$values))>=1){
    warning("Alpha has been changed.")
    alpha <- alpha / (max(abs(eigen(alpha/beta)$values))+2)
  }

  e0 <- t(t(solve(diag(M) - alpha/beta)) * mu)
  r <- c(); t <- 0; N <- c()
  delta <- rexp(n = M)

  while (t < horizon){

    func1 <- function(i){

      tmp1 <- function(x){
        mu[i]*x + sum(e0[i,]/beta[i,]*(1 - exp(- beta[i,] * x))) - delta[i]
      }

      tmp2 <- function(x){
        mu[i] + sum(e0[i,]*exp(- beta[i,] * x))
      }

      d <- pracma::newtonRaphson(fun = tmp1, x0 = 0, dfun = tmp2)
      return(d$root)
    }

    s <- sapply(1:M, func1)

    d <- min(s)
    t <- t + d

    if (t > horizon){
      break
    }

    r <- c(r, d)
    m <- seq(M)[s == d]
    N <- c(N, m)

    func2 <- function(i){
      delta[i] - (mu[i]*d + sum(e0[i,]/beta[i,]*(1-exp(- beta[i,] * d))))
    }
    delta <- sapply(1:M, func2)
    delta[m] = rexp(1)

    e0 = e0 * exp( - b * d)
    e0[,m] = e0[,m] + a[,m]
  }

  tmp <- cbind(r, N)[order(r), ]
  if(burnin > 0){
    r <- tmp[,1][tmp[,1] > burnin]
    N <- tmp[,2][tmp[,1] > burnin]
    tmp <- cbind(r, N)
  }

  if(separate == TRUE){
    q <- list()
    for(i in 1:M){
      tmpvec <- c()
      for(j in 1:length(tmp[,1])){
        if(tmp[,2][j] == i){
          tmpvec <- c(tmpvec, tmp[,1][j])
        }
      }
      q[[i]] <- tmpvec
    }
    tmpq <- sapply(q, length)
    mtmpq <- max(tmpq)
    for(i in 1:length(q)){
      q[[i]] <- c(q[[i]], rep(0, mtmpq - length(q[[i]])))
    }
    ttmp <- do.call(cbind, q)
    return(list(r = tmp[,1], t = ttmp, N = tmp[,2]))
  }else{
    return(list(r = tmp[,1], N = tmp[,2]))
  }
}
