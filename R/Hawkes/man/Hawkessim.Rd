\name{Hawkes.sim}
\alias{Hawkes.sim}
\title{Simulation of Multivariate Hawkes process}
\usage{
Hawkes.sim(mu, Y, dist, delta, N, params, t, paramsfunc)
}
\description{
Function simulates a Hawkes process for a given M dimensions.
}
\examples{
# For M = 2, changing distribution
M <- 2
m <- c(2,3)
y <- matrix(c(rep(0,4)), ncol = M, nrow = M)
d <-  matrix(sample(c("Exp", "Gamma"),4,
                    replace = TRUE), ncol = M, nrow = M)
del <- matrix(c(runif(4,0,0.8)), ncol = M, nrow = M)
n <- c(12,7)
p <- list()
p[[1]] <- matrix(c(runif(4,0,10)), ncol = M, nrow = M)
p[[2]] <- matrix(c(runif(4,0,10)), ncol = M, nrow = M)
p[[3]] <- matrix(c(runif(4,0,2)), ncol = M, nrow = M)

test <- Hawkes.sim(mu = m, Y = y, dist = d,
                   delta = del, N = n, t= 10, params = p)
test

# For M = 4, constant distribution
M <- 4
m <- c(0.5,1,0.2,0.4)
y <- matrix(c(rep(0,16)), ncol = M, nrow = M)

d <-  matrix(c("Exp"), ncol = M, nrow = M)

del <- matrix(c(runif(16,0,0.8)), ncol = M, nrow = M)

n <- c(12,4,7,15)

p <- list()
p[[1]] <- matrix(rep(c(0.5,2,0.1,0.4),4), ncol = M, nrow = M)

test <- Hawkes.sim(mu = m, Y = y, dist = d,
                   delta = del, N = n, t = 10, params = p)
test
}
