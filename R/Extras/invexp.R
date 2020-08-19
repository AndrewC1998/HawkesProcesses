#
#
#
# Set up parameters
mtmp <- seq(1, 50, length.out = 5)
atmp <- seq(0.1, 20, length.out = 5)
btmp <- seq(0.1, 30, length.out = 5)

mabtmp <- expand.grid(mtmp, atmp, btmp)

mabtmp[,4] <- mabtmp[,2]/mabtmp[,3]
mabtmp <- mabtmp[mabtmp[,4] < 1,]

#
#
#
# Hessian functions
Hawkes.Hess <- function(t, mu, alpha, beta){
  
  n <- length(t)
  
  ABC <- function(i){
    if(i == 1){
      tmp <- c(0, 0, 0)
    }else{
      tmpA <- 0; tmpB <- 0; tmpC <- 0
      for(j in 1:n){
        if(t[j] < t[i]){
          tmpA <- tmpA + exp(-beta*(t[i] - t[j]))
          tmpB <- tmpB + (t[i] - t[j])*exp(-beta*(t[i] - t[j]))
          tmpC <- tmpC + ((t[i] - t[j])^2)*exp(-beta*(t[i] - t[j]))
        }  
      }
      tmp <- c(tmpA, tmpB, tmpC)
    }
    return(tmp)
  }
  
  tmp <- mapply(ABC, 1:n)
  tmpmm <- 0; tmpma <- 0; tmpmb <- 0
  tmpaa <- 0; tmpab1 <- 0; tmpab2 <- 0
  tmpbb1 <- 0; tmpbb2 <- 0
  
  for(i in 1:n){
    # Get recursive values
    tmpA <- tmp[1, i]
    tmpB <- tmp[2, i]
    tmpC <- tmp[3, i]
    
    # Calculate all unique cells
    tmpmm <- tmpmm + (-1)/((mu + alpha*tmpA)^2)
    tmpma <- tmpma + (-tmpA)/((mu + alpha*tmpA)^2)
    tmpmb <- tmpmb + (alpha*tmpB)/((mu + alpha*tmpA)^2)
    
    tmpaa <- tmpaa - ((tmpA)/(mu + alpha*tmpA))^2
    tmpab1 <- tmpab1 - ( ((1/beta)*(t[n] - t[i])*exp(-beta*(t[n] - t[i]))) + ((1/(beta^2))*(exp(-beta*(t[n] - t[i])) - 1) ) ) 
    tmpab2 <- tmpab2 + (-tmpB)/(mu + alpha*tmpA) + (alpha*tmpA*tmpB)/((mu + alpha*tmpA)^2)
    
    tmpbb1 <- tmpbb1 + ((1/beta)*((t[n] - t[i])^2)*exp(-beta*(t[n] - t[i]))) + ((2/(beta^2))*(t[n] - t[i])*exp(-beta*(t[n] - t[i]))) + ((2/(beta^3))*(exp(-beta*(t[n] - t[i])) - 1))
    tmpbb2 <- tmpbb2 + ((alpha*tmpC)/(mu + alpha*tmpA)) - ((alpha*tmpB)/(mu + alpha*tmpA))^2
  }
  
  tmpab <- tmpab1 + tmpab2
  tmpbb <- alpha*tmpbb1 + tmpbb2
  
  tmpmat <- matrix(c(0), nrow = 3, ncol = 3)
  tmpmat[1,1] <- tmpmm; tmpmat[2,2] <- tmpaa; tmpmat[3,3] <- tmpbb
  tmpmat[1,2] <- tmpmat[2,1] <- tmpma
  tmpmat[1,3] <- tmpmat[3,1] <- tmpmb
  tmpmat[2,3] <- tmpmat[3,2] <- tmpab
  
  return(tmpmat)
}  

#
#
#
# Expectation of intensity
expint <- function(vec){
  m <- vec[1]; a <- vec[2]; b <- vec[3]
  tmp <- m/(1 - a/b)
  return(tmp)
}

#
#
#
# lambda function
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
# Fix results matrix
mabtmp <- cbind(mabtmp, 0, 0, 0, 0)
colnames(mabtmp) <- c("mu", "alpha", "beta", "branching ratio", "intensity",
                      "inverse", "Fisher", "Convergent")
rownames(mabtmp) <- c()

#
#
#
# Fill in intensity column
for(i in 1:length(mabtmp[,1])){
  vec <- c(mabtmp[i, 1], mabtmp[i, 2], mabtmp[i, 3])
  mabtmp[i, 5] <- expint(vec)
}

#
#
#
# Calculate all inverse lambda expectations
horizon <- 1000
nits <- 100
# load package
library(hawkes)

for(i in 1:length(mabtmp[,1])){
  tmpl <- c()
  for(j in 1:nits){
    t <- simulateHawkes(mabtmp[i, 1], mabtmp[i, 2], mabtmp[i, 3], horizon)[[1]]
    tmpl[j] <- (1/lambda(max(t), t, mabtmp[i, 1], mabtmp[i, 2], mabtmp[i, 3]))
    print(j)
  }
  mabtmp[i, 6] <- sum(tmpl)/nits
  print(i)
}

#
#
#
# Add Fisher information 
mabtmp[, 7] <- -c(802.0905 ,  354.5469,  312.0812, 0, 0, 0, -80516.81,         0,        0, 0,
             74.18151 ,  22.96480,  22.28913, 0, 0, 0, -2763.568,         0,        0, 0,
             -38689.72, -18561.02,         0, 0, 0, 0,         0,         0,        0, 0,
             20.431032,  5.574972,  4.338602, 0, 0, 0, -455.6104, -269.0698,        0, 0, 
             -3312.237, -2209.766,         0, 0, 0, 0, -7530.348,         0,        0, 0,
             -53330.71,         0,         0, 0, 0, 0,  1.874998,  1.611516, 1.596414, 0,
             -304.5464, -132.4736, -77.15787, 0, 0, 0, -560.4759,         0, 0, 0,
             -1594.868, -1524.254,         0, 0, 0, 0, -3819.744,         0, 0, 0) 

#
#
#
# Add convergence values
mabtmp[, 8] <- unlist(tmpstr)

#
#
#
# Remove the zeros
newthing <- mabtmp[1,]
for(i in 2:length(mabtmp[,1])){
  if(mabtmp[i, 8] != 0){
    newthing <- rbind(newthing, mabtmp[i, ])
  }
}

rownames(newthing) <- c()
