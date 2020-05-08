\name{Hawkes.plot}
\alias{Hawkes.plot}
\title{Plot simulations of Multivariate Hawkes process}
\usage{
Hawkes.plot(data)
}
\description{
Plots the simulations of a Hawkes process frome \code{Hawkes.sim}.
}
\examples{
# For M = 4
M <- 4; m <- c(0.5,1,0.2,0.4)
y <- matrix(c(abs(rnorm(16,5,6))), ncol = M, nrow = M)
d <-  matrix(c("Exp"), ncol = M, nrow = M)
del <- matrix(c(runif(16,0,0.8)), ncol = M, nrow = M)
n <- c(0,0,0,0)

test <- Hawkes.sim(M = M, mu = m, Y = y,
dist = d, delta = del, N = n, t = 100)
Hawkes.plot(test)
}
