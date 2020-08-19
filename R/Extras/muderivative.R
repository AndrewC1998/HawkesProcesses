A <- function(i, t, beta){
  n <- length(t)
  if(i == 1){
    tmp <- 0
  }else{
    tmpA <- 0
    for(j in 1:n){
      if(t[j] < t[i]){
        tmpA <- tmpA + exp(-beta*(t[i] - t[j]))
      }  
    }
    tmp <- tmpA
  }
  return(tmp)
}

muderiv1 <- function(data, a, b, m){
  
  tn <- max(data)
  n <- length(data)
  
  tmpmu <- c()
  tmpA <- c()
  for(i in 1:n){
    if(i == 1){
      tmpA[i] <- 0
    }else{
      tmpA[i] <- (exp(-b*(data[i] - data[i-1])))*(1 + tmpA[i-1])
    }
    tmpmu[i] <- (1/(m + a*tmpA[i]))
  }
  val <- -tn + sum(tmpmu)
  return(val)
}

muderiv2 <- function(data, a, b, m){
  
  tn <- max(data)
  n <- length(data)
  
  tmpmu <- c()
  tmpA <- 0
  for(i in 1:n){
    tmpmu[i] <- (1/(m + a*A(i, data, b)))
  }
  val <- -tn + sum(tmpmu)
  return(val)
}

m <- 5; a <- 0.5; b <- 2
tmpmu <- c()
nits <- 100000
for(i in 1:nits){
  t <- simulateHawkes(m, a, b, 1000)[[1]]
  tmpmu[i] <- muderiv1(t, a, b, m)
  print(i)
}
sum(tmpmu)/nits


lambda <- function(t, data, m, a, b){
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

nits <- 10000
m <- 1; a <- 0.1; b <- 7.575; H <- 1000
tmplam <- c()
for(i in 1:nits){
  tmptmp <- c()
  t <- simulateHawkes(m, a, b, H)[[1]]
  for(j in 1:length(t)){
    tmptmp[j] <- (1/lambda(t[j], t, m, a, b))^2
  }
  tmplam[i] <- sum(tmptmp)
  print(i)
}

sum(tmplam)/nits


tmplam2 <- c()
for(i in 1:nits){
  t <- simulateHawkes(m, a, b, H)[[1]]
  tmplam2[i] <- (1/lambda(max(t), t, m, a, b))
  print(i)
}

H*sum(tmplam2)/nits


# you need to take negative of all your matrices and vectors for fisher cause you have E[H]
# but need E[-H]