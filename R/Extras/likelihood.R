#
#
#
# Setup parameters
mtmp <- seq(1, 50, length.out = 10)
atmp <- seq(0.1, 20, length.out = 10)
btmp <- seq(0.1, 30, length.out = 10)

gridtmp <- expand.grid(mtmp, atmp, btmp)
gridtmp[,4] <- gridtmp[,2]/gridtmp[,3]
gridtmp <- gridtmp[gridtmp[,4] < 1,]
rownames(gridtmp) <- c()
colnames(gridtmp) <- c("mu", "alpha", "beta", "BR")

eirange <- apply(gridtmp, c(1), expint)

#
#
#
# Load packages
library(hawkes)
library(plotly)
library(ggplot2)
library(gridExtra)

#
#
#
# Simulate
nits <- 100
horizon <- 1000

tmp1 <- c()

for(i in 1:length(gridtmp[,1])){
  tmp2 <- c()
  for(j in 1:nits){
    t <- simulateHawkes(gridtmp[i, 1], gridtmp[i, 2], gridtmp[i, 3], horizon)[[1]]
    tmp2[j] <- likelihoodHawkes(gridtmp[i, 1], gridtmp[i, 2], gridtmp[i, 3], t)
  }
  tmp1[i] <- sum(tmp2)/nits
  print(i)
}

#
#
#
# combine
gridtmp <- cbind(gridtmp, -tmp1)
colnames(gridtmp) <- c("muval", "alphaval", "betaval", "BR", "LL")

#
#
#
# Plots
# branching ratio and likelihood
q <- qplot(gridtmp$BR, gridtmp$LL, color = I("red"), size = I(2), alpha = I(1/2), 
           geom = c("point"), xlab = "Branching Ratio", ylab = "Log-Likelihood")
q

gridtmp$mu <- c(rep(c("1", "6", "12", "17", "23", "28", "34", "39", "45", "50"), 66))
gridtmp$alpha <- c(rep("0.1", 10), rep("2.3", 10), rep("0.1", 10), rep("2.3", 10),
                   rep("4.5", 10), rep("6.7", 10), rep("0.1", 10), rep("2.3", 10), 
                   rep("4.5", 10), rep("6.7", 10), rep("8.9", 10), rep("0.1", 10),
                   rep("2.3", 10), rep("4.5", 10), rep("6.7", 10), rep("8.9", 10),
                   rep("11.2", 10), rep("13.4", 10), rep("0.1", 10), rep("2.3", 10),
                   rep("4.5", 10), rep("6.7", 10), rep("8.9", 10), rep("11.2", 10),
                   rep("13.4", 10), rep("15.6", 10), rep("0.1", 10), rep("2.3", 10),
                   rep("4.5", 10), rep("6.7", 10), rep("8.9", 10), rep("11.2", 10),
                   rep("13.4", 10), rep("15.6", 10), rep("17.8", 10), rep("20", 10),
                   rep("0.1", 10), rep("2.3", 10), rep("4.5", 10), rep("6.7", 10),
                   rep("8.9", 10), rep("11.2", 10), rep("13.4", 10), rep("15.6", 10),
                   rep("17.8", 10), rep("20", 10), rep("0.1", 10), rep("2.3", 10),
                   rep("4.5", 10), rep("6.7", 10), rep("8.9", 10), rep("11.2", 10),
                   rep("13.4", 10), rep("15.6", 10), rep("17.8", 10), rep("20", 10),
                   rep("0.1", 10), rep("2.3", 10), rep("4.5", 10), rep("6.7", 10),
                   rep("8.9", 10), rep("11.2", 10), rep("13.4", 10), rep("15.6", 10),
                   rep("17.8", 10), rep("20", 10))

gridtmp$beta <- c(rep("3.4", 20), rep("6.7", 40), rep("10.1", 50), rep("13.4", 70),
                  rep("16.7", 80), rep("20", 100), rep("23.4", 100), rep("26.7", 100), 
                  rep("30", 100))

gridtmp$intensity <- eirange

#
#
# plot for mu
p <- ggplot(gridtmp, aes(BR, LL))
p <- p + geom_point(aes(colour = mu), size = 4)
p <- p + xlab("Branching Ratio") + ylab("Log-Likelihood")
p

p <- ggplot(gridtmp, aes(BR, LL))
p <- p + geom_point(aes(colour = muval), size = 4)
p <- p + xlab("Branching Ratio") + ylab("Log-Likelihood")
p

#
#
# plot for alpha
p <- ggplot(gridtmp, aes(BR, LL))
p <- p + geom_point(aes(colour = alpha), size = 4)
p <- p + xlab("Branching Ratio") + ylab("Log-Likelihood")
p

p <- ggplot(gridtmp, aes(BR, LL))
p <- p + geom_point(aes(colour = alphaval), size = 4)
p <- p + xlab("Branching Ratio") + ylab("Log-Likelihood")
p

#
#
# plot for beta
p <- ggplot(gridtmp, aes(BR, LL))
p <- p + geom_point(aes(colour = beta), size = 4)
p <- p + xlab("Branching Ratio") + ylab("Log-Likelihood")
p

p <- ggplot(gridtmp, aes(BR, LL))
p <- p + geom_point(aes(colour = betaval), size = 4)
p <- p + xlab("Branching Ratio") + ylab("Log-Likelihood")
p

#
#
#
# plot for intensity
p <- ggplot(gridtmp, aes(BR, LL))
p <- p + geom_point(aes(colour = intensity), size = 4)
p <- p + xlab("Branching Ratio") + ylab("Log-Likelihood")
p

# 
plottmp <- gridtmp[1, ] 
for(i in 1:length(gridtmp[,1])){
  if(gridtmp[i, 6] == "1"){
    plottmp <- rbind(plottmp, gridtmp[i, ])
  }
}

deter <- function(u,a,b){
  tmp1 <- ((a^2)/(2*u*b^5))*(b^2 + 2*u - 1)
  tmp2 <- ((a^2)/(2*u*(b^4)*(a + b)))*(b^2 - u*(a + b) - 1)
  tmp3 <- (a^3)/(b^6)
  tmp <- tmp1 - tmp2 + tmp3
  return(tmp)
}
