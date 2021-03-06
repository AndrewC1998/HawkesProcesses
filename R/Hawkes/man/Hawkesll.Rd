\name{Hawkes.ll}
\alias{Hawkes.ll}
\title{Generalised log-likelihood function}
\usage{
Hawkes.ll(t, mu, alpha, beta, P = 1)
}
\description{
Generalised log-likelihood function for Hawkes self-exciting process.
}
\arguments{
  \item{t}{
    Horizon
  }
  \item{mu}{
    Base Intensity
  }
  \item{alpha}{
  alpha parameter matrix
  }
  \item{beta}{
  beta parameter matrix
  }
  \item{P}{
  kernel order
  }
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

test <- Hawkes.sim(mu = m, Y = y, dist = d,
                   delta = del, N = n, t = 10, params = p)

ll <- Hawkes.ll(matrix(c(test$r,test$r), ncol = 2), mu = m,
                alpha = y, beta = del)
}
