\name{Hawkes.sim}
\alias{Hawkes.sim}
\title{Simulation of Multivariate Hawkes process}
\usage{
Hawkes.sim(M, mu, Y, dist, delta, N, t)
}
\description{
Function simulates a Hawkes process for a given M dimensions.
}
\examples{
# For M = 2
M <- 2; m <- c(0.5,1); y <- matrix(c(1,0,0,1), ncol = M)
d <-  matrix(c("Exp"), ncol = M, nrow = M)
del <- matrix(c(1), ncol = M, nrow = M)
n <- c(1,1)

test <- Hawkes.sim(M = M, mu = m, Y = y,
dist = d, delta = del, N = n, t = 5)

# For M = 4
M <- 4; m <- c(0.5,1,0.2,0.4)
y <- matrix(c(rep(0,16)), ncol = M, nrow = M)
d <-  matrix(c("Exp"), ncol = M, nrow = M)
del <- matrix(c(runif(16,0,0.8)), ncol = M, nrow = M)
n <- c(0,0,0,0)

test <- Hawkes.sim(M = M, mu = m, Y = y,
dist = d, delta = del, N = n, t = 0.1)
plot(1:length(test$r), test$r, type = "l", col = "blue")
plot(1:length(test$N[,1]), test$N[,1], type = "l",
col = "blue", xlab = "t", ylab = "N")
lines(1:length(test$N[,2]), test$N[,2], col = "red")
lines(1:length(test$N[,3]), test$N[,3], col = "green")
lines(1:length(test$N[,4]), test$N[,4], col = "black")
grid(20,20)
}
