estimated_intensity <- function(params, arrivals){
  alpha_i <- params[1]
  beta_i <- params[2]
  mu_i <- params[3]
  n <- length(arrivals)
  Ai <- c(0, sapply(2:n, function(z){
        sum(exp( -beta_i * (arrivals[z]- arrivals[1:(z - 1)])))
       }))
  return(mu_i + alpha_i *Ai)
}

compensator <- function(params, arrivals){
  alpha_i <- params[1]
  beta_i <- params[2]
  mu_i <- params[3]
  n <- length(arrivals)
  compensator <- mu_i*arrivals - cumsum(alpha_i/beta_i*(exp(-beta_i*(arrivals[n] - arrivals))-1))
  return(compensator)
}
