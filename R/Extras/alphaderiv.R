lambda <- function(t, data = H, m, a, b){
  tmpl <- 0
  n <- length(data)
  for(i in 1:n){
    if(data[i] < t){
      tmpl <- tmpl + a*exp(-b*(t - data[i]))
    }
  }
  lam <- m + tmpl
  return(lam)
}

H