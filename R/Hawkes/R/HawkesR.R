Msim <- function(T, mu, b, a, e0){
  # Matrix should be indexed [m,i]
  M = length(mu)
  r = c(0)
  causes = c()
  effects = c()
  t = 0
  
  l = e0
  while (t < T){
    
    # childern
    U = matrix(runif(M*M), ncol = M)
    finite = U < 1 - exp( - l/ b)
    indices = seq(M*M)[finite]
    if (length(indices) == 0){
      t_diff = Inf
    } else {
      s = - 1/b[finite] * log( 1 + b[finite]/l[finite] * log(1-U[finite]))
      t_diff = min(s)
      index = indices[s == t_diff]
      m = index %% M
      if (m == 0){ m = M }
      i = (index - m)/M + 1
    }
    
    # immigrant
    
    s = - log( 1 - runif(M) ) / mu
    
    if (t_diff > min(s)){
      t_diff = min(s)
      m = seq(M)[s == t_diff]
      i = 0
    }
    
    # Update
    
    l = l * exp(-b*t_diff)
    l[m,] = l[m,] + a[m,]
    t = t + t_diff
    #cat(t, '\n')
    r = c(r, t)
    effects = c(effects, m)
    causes = c(causes, i)
  }
  return(list( r = r, causes = causes, effects = effects ))
}
