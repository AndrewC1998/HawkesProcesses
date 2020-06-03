# Simulate process
mu <- 1
Y <- matrix(0.5)
dist <- matrix("Constant")
delta <- matrix(2)
N <- c(0)
t <- 10000

test <- Hawkes.sim(mu,Y,dist,delta,N,t)
# true alpha is 1 and beta is 2

alpha <- seq(0.1,0.8,length.out = 100)
beta <- seq(0.1,0.8, length.out = 100)

tmp <- matrix(0, ncol = 100, nrow = 100)
for(i in 1:length(alpha)){
  for(j in 1:length(beta)){
    tmp[j,i] <- Hawkes.ll(t, 1, alpha[i], beta[j])
  }
}

library(plotly)
q <- plot_ly(x = alpha, y = beta, z = tmp, color = I("red")) %>% add_surface()
q %>% layout(scene = list(xaxis = list(title = "Alpha"), yaxis = list(title = "Beta"),
                          zaxis = list(title = "Log-Likelihood")))

max(tmp)
which(tmp == max(tmp), arr.ind = TRUE)



M <- 1
m <- c(0.1)
y <- matrix(1, ncol = M, nrow = M)
d <-  matrix("Constant", ncol = M, nrow = M)
del <- matrix(2, ncol = M, nrow = M)
n <- c(0)
test <- Hawkes.sim(mu = m, Y = y, dist = d,
                   delta = del, N = n, t = 100, params = list())


M <- 2
m <- c(1,2)
y <- matrix(c(1,0.1,0.1,1), ncol = M, nrow = M)
d <-  matrix("Constant", ncol = M, nrow = M)
del <- matrix(c(1,0.1,0.1,1), ncol = M, nrow = M)
n <- c(0,0)
test1 <- Hawkes.sim(mu = m, Y = y, dist = d,
                   delta = del, N = n, t = 100, params = list())

# MCMC and optim - heres how
# optim
t <- simulateHawkes(1, 0.5, 0.501, 10000)[[1]]
f <- function(params){
  - Hawkes.ll(t, params[1], params[2], params[3])
}
params <- optim(c(1,0.1,0.1), f)
paste( c("mu", "alpha", "beta"), round(params$par,2), sep=" = ")

# MCMC
tmpb <- c()
for(j in 1:length(beta)){
  tmpa <- c()
  for(i in 1:length(alpha)){
    tmpa[i] <- - Hawkes.ll(t, 1, alpha[i], beta[j])
  }
  MCMC <- sum(alpha*(1/sum(tmpa))*tmpa)
  tmpb[j] <- MCMC
}
ahat <- mean(tmpb)

tmpa <- c()
for(i in 1:length(alpha)){
  tmpb <- c()
  for(j in 1:length(beta)){
    tmpb[j] <- - Hawkes.ll(t, 1, alpha[i], beta[j])
  }
  MCMC <- sum(beta*(1/sum(tmpb))*tmpb)
  tmpa[i] <- MCMC
}
bhat <- mean(tmpa)

c(ahat, bhat)

# M is 1 case
M <- 1
mu <- 1
a <- matrix(0.5)
b <- matrix(2)

test <- Hawkes.sim2(mu, a, b, 10000, 0)

# M is 2 case
M <- 2
mu <- c(1,2)
a <- matrix(c(0.5, 0.2, 0.2, 0.5), ncol = 2)
b <- matrix(c(2, 1.1, 1.1, 2.2), ncol = 2)

test <- Hawkes.sim2(mu, a, b, 10000, 0)

q <- list()
for(i in 1:M){
  tmpvec <- c()
  for(j in 1:length(tmp$t)){
    if(tmp$N[j] == i){
      tmpvec <- c(tmpvec, tmp$t[j])
    }
  }
  q[[i]] <- tmpvec
}
