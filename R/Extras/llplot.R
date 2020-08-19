library(hawkes)

alpha <- seq(0,10, length.out = 20)
beta <- seq(0.1,10, length.out = 20)
tmp4 <- expand.grid(alpha,beta)

tmp4[,3] <- tmp4[,1]/tmp4[,2]
tmp5 <- tmp4[tmp4[,3] <= 1, 1:3]

brl <- cbind(tmp5, 0)
colnames(brl) <- c("alpha", "beta", "BR","ll")
m <- 0.5

t <- simulateHawkes(m, 2, 2.2, 500)[[1]]

for(i in 1:length(brl[,1])){
  brl[i,4] <- likelihoodHawkes(m, brl[i,1], brl[i,2], t)
  #print(i)
}

plot(brl[,3], brl[,4], type = "p", col = "red", xlab = "Branching Ratio", ylab = "Negative Log-Likelihood")
grid(20,20)

brlh <- cbind(brl, 0)
colnames(brlh) <- c("alpha", "beta", "BR","ll", "Hess")
rownames(brlh) <- NULL
for(i in 104:105){
  t <- simulateHawkes(m, brlh[i,1], brlh[i,2], 250)[[1]]
  brlh[i,5] <- det(Hawkes.Hess(t, m, brlh[i,1], brlh[i,2]))
  print(i)
}

c(105, 120, 136, 152, 153, 171, 189, 210)

plot(brlh[,3], brlh[,5], type = "p", col = "red", xlab = "Branching Ratio", ylab = "Determinant of Hessian")
grid(20,20)

brlh2 <- cbind(brl, 0)
colnames(brlh2) <- c("alpha", "beta", "BR", "ll", "Hess")
rownames(brlh2) <- NULL
t <- simulateHawkes(m, 2, 2.2, 500)[[1]]
for(i in 1:length(brlh2[,1])){
  brlh2[i,4] <- likelihoodHawkes(m, brlh2[i,1], brlh2[i,2], t)
  brlh2[i,5] <- det(Hawkes.Hess(t, m, brlh2[i,1], brlh2[i,2]))
  print(i)
}
plot(brlh2[-1,3], brlh2[-1,5], type = "p", col = "red", xlab = "Branching Ratio", ylab = "Determinant of Hessian")
grid(20, 20)

brlh3 <- c()
brlh4 <- c()
for(i in 1:210){
  if(brlh2[i,5] < 1e8 && brlh2[i,5] > -1e8){
    brlh3 <- c(brlh3, brlh2[i,3])
    brlh4 <- c(brlh4, brlh2[i,5])
  }
}
plot(brlh3, brlh4, type = "p", col = "red", xlab = "Branching Ratio", ylab = "Determinant of Hessian")
grid(20,20)
