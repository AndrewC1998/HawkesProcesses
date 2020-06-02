Msim <- function(T, mu, a, b, e0_){
  M = length(mu)
  e0 = e0_
  r = c()
  t = 0
  label = c()
  
  Deltas = rexp(n = M)
  
  while (t<T){
    s = sapply(1:M, function(i) {
      
      # uLimit = 10
      # eq_2_solve = function(x){
      #   mu[i]*x + sum(e0[i,]/b[i,]*(1-exp(- b[i,] * x))) - Deltas[i]
      # }
      # while (eq_2_solve(uLimit) <= 0){uLimit = uLimit *2}
      # d = uniroot(eq_2_solve, c(0,uLimit), tol = 1e-9, maxiter = 10000)$root
      eq_2_solve = function(x){
           mu[i]*x + sum(e0[i,]/b[i,]*(1-exp(- b[i,] * x))) - Deltas[i]}
      diff_of_eq = function(x) {mu[i] + sum(e0[i,]*exp(- b[i,] * x))}
      d = newtonRaphson(fun = eq_2_solve, x0 = 0, dfun = diff_of_eq)
      d$root
    })
    
    d = min(s)
    t = t + d
    if (t>T){
      return(list(r = r, labels = label))
    }
    r = c(r, d)
    m = seq(M)[s == d]
    label = c(label, m)
    
    Deltas <- sapply(1:M,function(i){
      Deltas[i] - (mu[i]*d + sum(e0[i,]/b[i,]*(1-exp(- b[i,] * d))))
    })
    Deltas[m] = rexp(1)
    
    e0 = e0 * exp( - b * d)
    e0[,m] = e0[,m] + a[,m]
  }
}
