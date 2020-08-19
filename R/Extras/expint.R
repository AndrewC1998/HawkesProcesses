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
mabtmp <- mabtmp[mabtmp[,5] < 0.95,]

mabtmp <- mabtmp[,1:4]

colnames(mabtmp) <- c("mu", "alpha", "beta", "intensity")
rownames(mabtmp) <- c()

eitmp <- c()
nits <- 1000
eimc <- c()
eils <- list()

for(i in 1:length(mabtmp[,1])){
  eihess <- matrix(0, ncol = 3, nrow = 3)
  for(j in 1:nits){
    t <- hawkes::simulateHawkes(mabtmp[i,1], mabtmp[i,2], mabtmp[i,3], 100)[[1]]
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