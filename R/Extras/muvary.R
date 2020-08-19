# load libraries
library(hawkes)
library(latex2exp)

# Setup parameters
mu <- seq(0, 100, length.out = 30)[-1]
a <- 0.5
b <- 2
nits <- 2

mutmp <- c()
mustr <- c()

# Compute tasks
for(i in 1:length(mu)){
  for(j in 1:nits){
    t <- hawkes::simulateHawkes(mu[i], a, b, 100)[[1]]
    mustr[j] <- det(Hawkes.Hess(t, mu[i], a, b))
    print(j)
  }
  mutmp[i] <- sum(mustr)/nits
  print(i)
}

plot(mu, mutmp, type = "l", col = "red", xlab = TeX('$\\mu$'),
     ylab = "Determinant of Hessian")
grid(40,40)
