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



# Usage
m <- 0.5
b <- 2
a <- seq(0, 2, length.out = 50)
N <- 100

for(i in 47:49){
  MC <- c()
  for(j in 1:N){
    t <- simulateHawkes(m, a[i], b, 500)[[1]]
    MC[j] <- det(Hawkes.Hess(t, m, a[i], b))
    print(j)
  }
  tmp[i] <- sum(MC)/N
  print(i)
}

br <- a/b
plot(br, tmp1, type = "l", col = "red", xlab = "Branching Ratio", ylab = "MC Estimation Determinant of Hessian")
grid(40, 40)

N <- 5
MC <- c()
for(j in 1:N){
  t <- simulateHawkes(m ,a[50], b, 500)[[1]]
  MC[j] <- det(Hawkes.Hess(t, m, a[50], b))
  print(j)
}
tmp[50] <- sum(MC)/N
