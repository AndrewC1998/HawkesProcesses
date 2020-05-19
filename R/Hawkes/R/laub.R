# Simple Laub method
laub.ll <- function(t, mu, alpha, beta){
  # t is the history of the process
  # mu is the mean of the process
  # alpha is the alpha parameter
  # beta is the beta parameter
  tmax <- max(t)
  N <- length(t)
  tmp1 <- c(); tmp2 <- c()
  R <- c()

  for(k in 1:N){
    if(k == 1){
      R[k] <- 0
    }else{
      R[k] <- (exp(-beta*(t[k] - t[k-1])))*(1 + R[k-1])
    }
    tmp1[k] <- log(mu + alpha*R[k])
    tmp2[k] <- (1 - exp(-beta*(tmax - t[k])))
  }
  ll <- sum(tmp1) - mu*tmax - (alpha/beta)*sum(tmp2)
  return(ll)
}

# mle function
laub.mle <- function(t, mu, alpha, beta, divisions){
  mu <- c(seq(mu[1],mu[2], length.out = divisions))
  alpha <- c(seq(alpha[1], alpha[2], length.out = divisions))
  beta <- c(seq(beta[1], beta[2], length.out = divisions))
  # Expand choices of (mu, alpha, beta)
  eg <- expand.grid(mu, alpha, beta)
  # log-likelihood function
  ll <- function(m, a, b){
    return( - laub.ll(t = t, m, a, b))
  }
  tmp <- mapply(ll, m = eg$Var1, a = eg$Var2, b = eg$Var3)
  val <- which.min(tmp)
  params <- eg[377,]
  params <- cbind(params, -min(tmp))
  rownames(params) <- c("MLE")
  colnames(params) <- c("mu", "alpha", "beta", "loglik")
  return(params)
}

