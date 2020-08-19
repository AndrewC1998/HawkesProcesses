
# Usage
m <- 0.5
b <- 2
a <- seq(0, 2, length.out = 10)
N <- 10

tmp1 <- c()
A <- Sys.time()
for(i in 1:length(a)){
  MC <- c()
  for(j in 1:N){
    t <- simulateHawkes(m, a[i], b, 3000)[[1]][1:1000]
    MC[j] <- det(Hawkes.Hess(t, m, a[i], b))
    print(j)
  }
  tmp1[i] <- sum(MC)/N
  print(i)
}
B <- Sys.time()
B - A

br <- a/b
plot(br, tmp1, type = "l", col = "red", xlab = "Branching Ratio", ylab = "MC Estimation Determinant of Hessian")
grid(40, 40)
