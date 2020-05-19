\name{laub.mle}
\alias{laub.mle}
\title{Provides MLE parameters Laub's Equations.}
\usage{
Hawkes.plot(data, type = "Count")
}
\description{
laub.mle(t, mu, alpha, beta, divisions)
}
\examples{
# For M = 2, changing distribution
M <- 2
m <- c(2,3)
y <- matrix(c(rep(0,4)), ncol = M, nrow = M)
d <-  matrix(sample(c("Exp", "Gamma"),4,
                    replace = TRUE), ncol = M, nrow = M)
del <- matrix(c(runif(16,0,0.8)), ncol = M, nrow = M)
n <- c(12,7)
p <- list()
p[[1]] <- matrix(c(runif(4,0,10)), ncol = M, nrow = M)
p[[2]] <- matrix(c(runif(4,0,10)), ncol = M, nrow = M)
p[[3]] <- matrix(c(runif(4,0,2)), ncol = M, nrow = M)

test <- Hawkes.sim(M = M, mu = m, Y = y, dist = d,
                   delta = del, N = n, params = p, t = 10)

laub.mle(t = test$r, mu = c(0,5), alpha = c(0.1,5), beta = c(0.1,5), divisions = 100)
}
