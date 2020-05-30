Hawkes.ll <- function(t, mu, alpha, beta, P = 1){
  # t is matrix of times each col is M=1,2,...
  # mu is vector of base intensities
  # alpha is a square matrix of alpha coefficients
  # beta is a square matrix of beta coefficients
  # P is model order

  M <- length(mu)
  N <- length(t)/M

  if(P == 1){
    if(M == 1){
      R <- c()
      R[1] <- 0

      # Find the first summation in likelihood
      tmp1 <- log(mu + alpha*R[1])
      for(k in 2:N){
        R[k] <- (exp(-beta*(t[k] - t[k-1])))*(1 + R[k-1])
        tmp1 <- tmp1 + log(mu + alpha*R[k])
      }

      # Find second summation in likelihood
      Tmax <- max(t)
      tmp2 <- 0
      if(is.finite(alpha/beta)){
        for(k in 1:N){
          tmp2 <- tmp2 + 1 - exp(-beta*(Tmax - t[k]))
        }
        tmp2 <- (alpha/beta)*tmp2
      }

      # Sum all the components
      ll <- tmp1 - mu*Tmax - tmp2
    }else{
      lltmp <- c()
      Tmax <- max(t)
      R <- matrix(c(0), ncol = M, nrow = M)
      for(m in 1:M){
        tmp1 <- 0
        for(k in 1:N){
          if(k == 1){
            tmp1 <- log(mu[m])
          }else{
            tmp2 <- 0
            for(i in 1:M){
              dtmp <- 0
              for(d in 1:N){
                if(t[k-1,m] <= t[d,i] && t[d,i] < t[k,m]){
                  dtmp <- exp(-beta[i,m]*(t[k,m] - t[d,i]))
                }
              }
              R[i,m] <- exp(-beta[i,m]*(t[k,m] - t[k-1,m]))*R[i,m] + dtmp
              tmp2 <- tmp2 + alpha[i,m]*R[i,m]
            }
            tmp1 <- tmp1 + log(mu[m] + tmp2)
          }
        }

        tmp3 <- 0
        for(i in 1:M){
          for(k in 1:N){
            tmp4 <- alpha[i,m]/beta[i,m]
            tmp5 <- 1 - exp(-beta[i,m]*(Tmax - t[k,i]))
            if(is.finite(tmp4)){
              tmp3 <- tmp3 + tmp4*tmp5
            }
          }
        }
        lltmp[m] <- tmp1 - mu[m]*Tmax - tmp3
      }
      ll <- sum(lltmp)
    }
  }else{
    if(M == 1){
      Tmax <- max(t)
      R <- c(rep(0, P))
      tmp <- 0
      for(k in 1:N){
        if(k == 1){
          tmp1 <- log(mu)
        }else{
          tmp2 <- 0
          for(j in 1:P){
            R[j] <- (exp(-beta[j]*(t[k] - t[k-1])))*(1 + R[j])
            tmp2 <- tmp2 + alpha[j]*R[j]
          }
          tmp1 <- tmp1 + log(mu + tmp2)
        }
      }

      tmp3 <- 0
      for(j in 1:P){
        for(k in 1:N){
          tmp4 <- alpha[j]/beta[j]
          tmp5 <- 1 - exp(-beta[j]*(Tmax - t[k]))
          if(is.finite(tmp4)){
            tmp3 <- tmp3 + tmp4*tmp5
          }
        }
      }
      ll <- tmp1 - mu*Tmax - tmp3
    }else{
      lltmp <- c()
      Tmax <- max(t)
      R <- list()
      for(j in 1:P){
        R[[j]] <- matrix(c(0), nrow = M, ncol = M)
      }
      for(m in 1:M){
        tmp1 <- 0
        for(k in 1:N){
          tmp2 <- 0
          if(k == 1){
            tmp1 <- log(mu[m])
          }else{
            for(i in 1:M){
              for(j in 1:P){
                if(i == m){
                  R[[j]][i,m] <- (exp(-beta[[j]][i,m]*(t[k,m] - t[k-1,m])))*(1 + R[[j]][i,m])
                }else{
                  dtmp <- 0
                  for(d in 1:N){
                    if(t[k-1,m] <= t[d,i] && t[d,i] < t[k,m]){
                      dtmp <- exp(-beta[[j]][i,m]*(t[k,m] - t[d,i]))*(1 + R[[j]][i,m])
                    }
                  }
                  R[[j]][i,m] <- (exp(-beta[[j]][i,m]*(t[k,m] - t[k-1,m])))*(1 + R[[j]][i,m]) + dtmp
                }
                tmp2 <- tmp2 + alpha[[j]][i,m]*R[[j]][i,m]
              }
            }
            tmp1 <- tmp1 + log(mu[m] + tmp2)
          }
        }

        tmp3 <- 0
        for(i in 1:M){
          for(j in 1:P){
            for(k in 1:N){
              tmp4 <- alpha[[j]][i,m]/beta[[j]][i,m]
              tmp5 <- 1 - exp(-beta[[j]][i,m]*(Tmax - t[k,i]))
              if(is.finite(tmp4)){
                tmp3 <- tmp3 + tmp4*tmp5
              }
            }
          }
        }
        lltmp[m] <- tmp1 - mu[m]*Tmax - tmp3
      }
      ll <- sum(lltmp)
    }
  }
  return(ll)
}
