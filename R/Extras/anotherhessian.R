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

expint <- function(vec){
  m <- vec[1]; a <- vec[2]; b <- vec[3]
  tmp <- m/(1 - a/b)
  return(tmp)
}

mtmp <- seq(1, 50, length.out = 5)
atmp <- seq(0.1, 20, length.out = 5)
btmp <- seq(0.1, 30, length.out = 5)

mabtmp <- expand.grid(mtmp, atmp, btmp)

eirange <- apply(mabtmp, c(1), expint)

mabtmp <- cbind(mabtmp, eirange, 0)
mabtmp[,5] <- mabtmp[,2]/mabtmp[,3]
mabtmp <- mabtmp[mabtmp[,5] < 1,]

mabtmp <- mabtmp[,1:4]

colnames(mabtmp) <- c("mu", "alpha", "beta", "intensity")
rownames(mabtmp) <- c()

eitmp <- c()
nits <- 1000
eimc <- c()
eils <- list()

for(i in 2:2){
  eihess <- matrix(0, ncol = 3, nrow = 3)
  for(j in 1:nits){
    t <- hawkes::simulateHawkes(mabtmp[i,1], mabtmp[i,2], mabtmp[i,3], 1000)[[1]]
    Htmp <- Hawkes.Hess(t, mabtmp[i,1], mabtmp[i,2], mabtmp[i,3])
    eitmp[j] <- det(Htmp)
    eihess <- eihess + Htmp
    print(j)
  }
  eimc[i] <- sum(eitmp)/nits
  eils[[i]] <- eihess/nits
  print(i)
}

plot(mabtmp[,4], eimc, type = "p", col = "red", 
     xlab = TeX('E$\\[\\lambda\\]$'),
     ylab = "Determinant of Hessian")
grid(40,40)



# alternative and faster
testfunc <- function(vec){
  t <- hawkes::simulateHawkes(vec[1], vec[2], vec[3], 100)[[1]]
  Htmp <- Hawkes.Hess(t, vec[1], vec[2], vec[3])
  return(det(Htmp))
}

trynow <- mabtmp[1:50,]
alttmp <- apply(trynow, c(1), testfunc)

# notable if the data is 10x bigger then the determinant is 10^3 