# FIM derivative by derivative
#
#
#
# First set up generic parameters to test
nits <- 10000
m <- 1; a <- 0.1; b <- 7.575
horizon <- 1000
test <- 1

#
#
#
# Add in lambda function
lambda <- function(t, data, m , a, b){
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

#
#
#
# l/mu^2
tmplmu <- c()
for(i in 1:nits){
  t <- simulateHawkes(m, a, b, horizon)[[1]]
  tmplmu[i] <- (1/lambda(max(t), t, m, a, b))
  print(i)
}

exp.lambda <- sum(tmplmu)/nits
lmu <- horizon*exp.lambda
c(lmu, fim[[test]][1,1])

#
#
#
# l/a^2
tmpla <- c()
for(i in 1:nits){
  t <- simulateHawkes(m, a, b, horizon)[[1]]
  tmpla[i] <- (1/lambda(max(t), t, m, a, b))
  print(i)
}

lapre <- sum(tmpla)/nits
multip <- (m*horizon)/(a^2)
la <- multip*((a/(b - a)) - 1 + lapre)
c(la, fim[[test]][2,2])

#
#
#
# l/b^2
tmplb <- c()
for(i in 1:nits){
  t <- simulateHawkes(m, a, b, horizon)[[1]]
  tmplb[i] <- (1/lambda(max(t), t, m, a, b))
  print(i)
}

#
#
#
# l/ua

#
#
#
# l/ub

#
#
#
# l/ab

lbpre <- sum(tmplb)/nits
multip <- m*lbpre
lb <- multip*(1 + ((a*b*horizon^2 - a*horizon)/(b*(b-a)))) - 1


#
#
#
# The full fisher info
fisher <- function(m, a, b, horizon, nits = 10000){
  
  #
  # load any packages
  #
  require(hawkes)
  #
  # Include lambda function
  #
  lambda <- function(t, data, m , a, b){
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
  #
  #
  #
  
  # setup empty matrix
  FImatrix <- matrix(c(0), ncol = 3, nrow = 3)
  
  # setup empty lambda
  invlambda <- c()
  
  # calculate expectation of inverse of lambda(t)
  for(i in 1:nits){
    t <- simulateHawkes(m, a, b, horizon)[[1]]
    invlambda[i] <- (1/lambda(max(t), t, m, a, b))
  }
  
  exp.invlambda <- sum(invlambda)/nits
  
  # calculate derivatives
  FImatrix[1, 1] <- horizon*exp.invlambda
  
  FImatrix[1, 2] <- FImatrix[2, 1] <- (horizon/a)*(1 - m*exp.invlambda)
    
  FImatrix[1, 3] <- FImatrix[3, 1] <- -(horizon^2) + (m*(horizon^2) + (a*(b*horizon - 1))/(b^2))*exp.invlambda
  
  FImatrix[2, 2] <- ((m*horizon)/(a^2))*((a/(b-a)) - 1 + m*exp.invlambda)
  
  FImatrix[2, 3] <- FImatrix[3, 2] <- (m/a - m/(b-a) + 1/b)*(horizon^2) - (1/(b^2))*horizon + (((m^2)*(horizon^2))/a + (a*m*horizon)/b - (a*m)/(b^2))*exp.invlambda
  
  FImatrix[3, 3] <- ((a*m)/(b-a) - m - (2*a)/b)*(horizon^3) + ((2*a)/(b^2))*(horizon^2) + ((((m^2)*horizon)/(b^2))*(((b*horizon)^2) + ((2*horizon)/m)*(b*horizon - 1) + ((a/(b-a))*(b*horizon - 1))^2))*exp.invlambda
  
  return(FImatrix)
}

ltmp <- c()
for(i in 1:nits){
  t <- simulateHawkes(1, 0.1, 7.575, 1000)[[1]]
  ltmp[i] <- lambda(max(t), t, m, a, b)
  print(i)
}
sum(ltmp)/10000
(7.575*1)/(7.575 - 0.1)
