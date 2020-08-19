# MVN 
library(MissMech)
library(MASS)

mu <- c(1,1)
sig <- matrix(c(1,0.3,0.3,1), ncol = 2)

MVNsim <- mvrnorm(1000, mu, sig)

MVNhess <- Ddf(MVNsim, as.matrix(mu), sig)$dd

nits <- 10000
n <- 10000
MCnorm <- c()
Mtmp <- matrix(c(0), ncol = 5, nrow = 5)
for(i in 1:nits){
  MVNsim <- mvrnorm(n, mu, sig)
  MVNhess <- Ddf(MVNsim, as.matrix(mu), sig)$dd
  Mtmp <- Mtmp + MVNhess
  MCnorm[i] <- det(MVNhess)
  print(i)
}

edFIM <- -sum(MCnorm)/nits
eFIM <- -Mtmp/nits

# check using det function
sige <- sig[2,2] - (sig[1,2]^2)/sig[1,1]
beta <- sig[1,2]/sig[1,1]
FIM <- matrix(c(0), ncol = 5, nrow = 5)
FIM[1,1] <- n/sige; FIM[1,2] <- -(n*beta)/sige
FIM[2,1] <- FIM[1,2]; FIM[2,2] <- (n*beta)/sige + n/sig[1,1]
FIM[3,3] <- (n*sig[1,1])/sige
FIM[4,4] <- n/(2*sig[1,1]^2)
FIM[5,5] <- n/(2*sige^2)
det(FIM)

# check method 2

detFIM <- ((n^5)/(4*(sig[1,1]^2)*sige^4))*((sig[1,2]/sige) - ((sig[1,2]^2)/sig[1,1]) + 1)