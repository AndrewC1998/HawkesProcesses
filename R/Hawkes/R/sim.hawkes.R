sim.hawkes <- function(mu, alpha, beta, n, variate = "univariate", seed = 1){
  if(variate == "univariate"){
    s <- 0; t <- 0
    set.seed(seed)
    lstar <- mu; dl <- alpha
    U <- runif(1)
    s <- s - log(U)/lstar
    t <- s
    tmp <- c(t)
    while( s < n){
      U <- runif(1)
      s <- s - log(U)/lstar
      u <- runif(1)
      if(u <= (mu + dl*exp(-beta*(s-t)))/lstar){
        dl <- alpha + dl*exp(-beta*(s-t))
        lstar <- lstar + alpha
        t <- s; tmp <- c(tmp, t)
      }
    }
  }else{
    stop("Not currently supported")
  }
  return(tmp)
}
