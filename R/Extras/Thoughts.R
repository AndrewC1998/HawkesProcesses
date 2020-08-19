# Simulate process
mu <- 1
Y <- matrix(0.5)
dist <- matrix("Constant")
delta <- matrix(2)
N <- c(0)
t <- 10000

test <- Hawkes.sim(mu,Y,dist,delta,N,t)
# true alpha is 1 and beta is 2
t <- hawkes::simulateHawkes(0.5, 2, 2.2, 1000)[[1]]
alpha <- seq(1.8, 2.2,length.out = 100)
beta <- seq(2, 2.4, length.out = 100)

tmp <- matrix(0, ncol = 100, nrow = 100)
for(i in 1:length(alpha)){
  for(j in 1:length(beta)){
    tmp[j,i] <- Hawkes.ll(t, 0.5, alpha[i], beta[j])
  }
}

library(plotly)
q <- plot_ly(x = alpha, y = beta, z = tmp, color = I("red")) %>% add_surface()
q %>% layout(scene = list(xaxis = list(title = "Alpha"), yaxis = list(title = "Beta"),
                          zaxis = list(title = "Log-Likelihood")))

max(tmp)
which(tmp == max(tmp), arr.ind = TRUE)



M <- 1
m <- c(1)
y <- matrix(0.5, ncol = M, nrow = M)
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

# MCMC and optim - here's how
# optim
t <- hawkes::simulateHawkes(1, 0.5, 5, 10000)[[1]]
f <- function(params){
  - Hawkes.ll(t, params[1], params[2], params[3])
}
params <- optim(c(1,0.1,0.1), f)
paste( c("mu", "alpha", "beta"), round(params$par,2), sep=" = ")

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

# Plots

t <- hawkes::simulateHawkes(1, 0.5, 2, 100)[[1]]
ei <- estimated_intensity(c(0.5,2,1), t)
cmp <- compensator(c(0.5,2,1), t)

plot(t, ei, type = "l",
     col = "red", xlab = "t", ylab = "Intensity", ylim = c(0, 10))
cols <- c("red")
lines(t, cmp, col = "blue")
grid(20,20)


Hawkes.plot(test)
plot(test$r, test$N[,1], type = "s",
     col = "red", xlab = "Time", ylab = "Count",
     ylim = c(10^floor(log10(min(test$N))), max(test$N)+1), lwd = 3,
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
grid(25,25)


# Hessian functions

Hawkes.Hess <- function(t, mu, alpha, beta){
  
  n <- length(t)
  
  A <- function(i){
    if(i == 1){
      return(0)
    }else{
      tmpA <- 0
      for(j in 1:n){
        if(t[j] < t[i]){
          tmpA <- tmpA + exp(-beta*(t[i] - t[j]))
        }  
      }
      return(tmpA) 
    }
  }
  
  B <- function(i){
    if(i == 1){
      return(0)
    }else{
      tmpB <- 0 
      for(j in 1:n){
        if(t[j] < t[i]){
          tmpB <- tmpB + (t[i] - t[j])*exp(-beta*(t[i] - t[j]))
        }  
      }
      return(tmpB)
    }
  }
  
  C <- function(i){
    if(i == 1){
      return(0)
    }else{
      tmpC <- 0 
      for(j in 1:n){
        if(t[j] < t[i]){
          tmpC <- tmpC + ((t[i] - t[j])^2)*exp(-beta*(t[i] - t[j]))
        }  
      }
      return(tmpC)
    }
  }
  
  tmpmm <- 0; tmpma <- 0; tmpmb <- 0
  tmpaa <- 0; tmpab1 <- 0; tmpab2 <- 0
  tmpbb1 <- 0; tmpbb2 <- 0
  
  for(i in 1:n){
    # Get recursives
    tmpA <- A(i)
    tmpB <- B(i)
    tmpC <- C(i)
    
    # Calculate all unique cells
    tmpmm <- tmpmm + (-1)/((mu + alpha*tmpA)^2)
    tmpma <- tmpma + (-tmpA)/((mu + alpha*tmpA)^2)
    tmpmb <- tmpmb + (alpha*tmpB)/((mu + alpha*tmpA)^2)
    
    tmpaa <- tmpaa - ((tmpA)/(mu + alpha*tmpA))^2
    tmpab1 <- tmpab1 - ( ((1/beta)*(t[n] - t[i])*exp(-beta*(t[n] - t[i]))) + ((1/(beta^2))*(exp(-beta*(t[n] - t[i])) - 1) ) ) 
    tmpab2 <- tmpab2 + (-tmpB)/(mu + alpha*tmpA) + (alpha*tmpA*tmpB)/((mu + alpha*tmpA)^2)
    
    tmpbb1 <- tmpbb1 + ((1/beta)*((t[n] - t[i])^2)*exp(-beta*(t[n] - t[i]))) + ((2/(beta^2))*(t[n] - t[i])*exp(-beta*(t[n] - t[i]))) + ((2/(beta^3))*(exp(-beta*(t[n] - t[i])) - 1))
    tmpbb2 <- tmpbb2 + ((alpha*tmpC)/(mu + alpha*tmpA)) - ((alpha*tmpB)/(mu + alpha*tmpA))^2
  }
  
  tmpab <- tmpab1 + tmpab2
  tmpbb <- alpha*tmpbb1 + tmpbb2
  
  tmpmat <- matrix(c(0), nrow = 3, ncol = 3)
  tmpmat[1,1] <- tmpmm; tmpmat[2,2] <- tmpaa; tmpmat[3,3] <- tmpbb
  tmpmat[1,2] <- tmpmat[2,1] <- tmpma
  tmpmat[1,3] <- tmpmat[3,1] <- tmpmb
  tmpmat[2,3] <- tmpmat[3,2] <- tmpab
  
  return(tmpmat)
}  

m <- 0.5
b <- 2
a <- seq(0.01, 1.99, length.out = 100)

t <- simulateHawkes(0.5, 1, 2, 100)[[1]]

tmp <- c()
for(i in 1:length(a)){
  tmp[i] <- det(Hawkes.Hess(t, m, a[i], b))
}

br <- a/b
plot(br, tmp, type = "l", col = "red", xlab = "Branching Ratio", ylab = "Determinant of Hessian")
grid(40,40)

tmp2 <- c()
m <- 0.5
b <- seq(20, 1.1, length.out = 100)
a <- 1
for(i in 1:length(b)){
  tmp2[i] <- det(Hawkes.Hess(t, m, a, b[i]))
}
plot(br, tmp2, col = "blue", type = "l")
