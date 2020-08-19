#
#
# libraries
library(hawkes)
library(ggplot2)
library(latex2exp)
library(plotly)
library(gridExtra)
#
m <- c(1,seq(2, 20, 2))
a <- c(1,seq(2, 20, 2))
b <- c(1,seq(2, 20, 2))
mab <- expand.grid(m, a, b)
mab[,4] <- mab[,2]/mab[,3]
mab <- mab[mab[,4] < 1,]
length(mab[,1])
rownames(mab) <- c()
colnames(mab) <- c("mu", "alpha", "beta", "BR")
#
likvec <- c()
for(i in 1:length(mab[,1])){
  t <- simulateHawkes(mab[i,1], mab[i,2], mab[i,3], 1000)[[1]]
  likvec[i] <- -likelihoodHawkes(mab[i,1], mab[i,2], mab[i,3], t)
  print(i)
}
mab <- cbind(mab, likvec)
rownames(mab) <- c()
colnames(mab) <- c("mu", "alpha", "beta", "BR", "LL")

mab <- as.data.frame(mab)
#
#
p1 <- qplot(mab[,1], likvec, 
            ylab = "Log-Likelihood", xlab = TeX("$\\mu$"), aes(color = "red"), show.legend = FALSE)
p2 <- qplot(mab[,2], likvec, 
            ylab = "Log-Likelihood", xlab = TeX("$\\alpha$"), aes(color = "blue"), show.legend = FALSE)
p3 <- qplot(mab[,3], likvec, 
            ylab = "Log-Likelihood", xlab = TeX("$\\beta$"), aes(color = "green"), show.legend = FALSE)
grid.arrange(p1,p2,p3, ncol = 3)

p1 <- ggplot(mab, aes(x = mu, y = LL)) + geom_point(size = 2, colour = "firebrick2") +
      xlab(TeX("$\\mu$")) + ylab("Log-Likelihood") +
      theme(axis.text.x = element_text(size = 42), 
            axis.text.y = element_text(size = 32), 
            axis.title.x = element_text(size = 42),
            axis.title.y = element_text(size = 42))

p2 <- ggplot(mab, aes(x = alpha, y = LL)) + geom_point(size = 2, colour = "dodgerblue3") +
      xlab(TeX("$\\alpha$")) + ylab("Log-Likelihood") +
      theme(axis.text.x = element_text(size = 42), 
            axis.text.y = element_text(size = 32), 
            axis.title.x = element_text(size = 42),
            axis.title.y = element_text(size = 42))

p3 <- ggplot(mab, aes(x = beta, y = LL)) + geom_point(size = 2, colour = "aquamarine4") +
      xlab(TeX("$\\beta$")) + ylab("Log-Likelihood") +
      theme(axis.text.x = element_text(size = 42), 
            axis.text.y = element_text(size = 32), 
            axis.title.x = element_text(size = 42),
            axis.title.y = element_text(size = 42))

grid.arrange(p1, p2, p3, ncol = 3)


#
#
# for mu
#
# constant beta
m <- c(seq(1,20,1))
a <- c(seq(1,20,1))
b <- c(10)
mab <- expand.grid(m, a, b)
mab[,4] <- mab[,2]/mab[,3]
mab <- mab[mab[,4] < 1,]
length(mab[,1])
rownames(mab) <- c()
colnames(mab) <- c("mu", "alpha", "beta", "BR")
likvec <- c()
for(i in 1:length(mab[,1])){
  t <- simulateHawkes(mab[i,1], mab[i,2], mab[i,3], 1000)[[1]]
  likvec[i] <- -likelihoodHawkes(mab[i,1], mab[i,2], mab[i,3], t)
  print(i)
}
mab <- cbind(mab, likvec)
rownames(mab) <- c()
colnames(mab) <- c("mu", "alpha", "beta", "BR", "LL")

p1 <- ggplot(mab, aes(x = mu, y = LL)) + geom_point(size = 2, aes(colour = alpha)) +
  xlab(TeX("$\\mu$")) + ylab("Log-Likelihood") +
  theme(axis.text.x = element_text(size = 32), 
        axis.text.y = element_text(size = 32), 
        axis.title.x = element_text(size = 42),
        axis.title.y = element_text(size = 42)) + scale_y_continuous(labels = scales::comma) +
  scale_color_gradient(low = "dodgerblue4", high = "dodgerblue")
#
# constant alpha
m <- c(seq(1,20,1))
a <- c(5)
b <- c(seq(1,20,1))
mab <- expand.grid(m, a, b)
mab[,4] <- mab[,2]/mab[,3]
mab <- mab[mab[,4] < 1,]
length(mab[,1])
rownames(mab) <- c()
colnames(mab) <- c("mu", "alpha", "beta", "BR")
likvec <- c()
for(i in 1:length(mab[,1])){
  t <- simulateHawkes(mab[i,1], mab[i,2], mab[i,3], 1000)[[1]]
  likvec[i] <- -likelihoodHawkes(mab[i,1], mab[i,2], mab[i,3], t)
  print(i)
}
mab <- cbind(mab, likvec)
rownames(mab) <- c()
colnames(mab) <- c("mu", "alpha", "beta", "BR", "LL")

p2 <- ggplot(mab, aes(x = mu, y = LL)) + geom_point(size = 2, aes(colour = beta)) +
  xlab(TeX("$\\mu$")) + ylab("Log-Likelihood") +
  theme(axis.text.x = element_text(size = 32), 
        axis.text.y = element_text(size = 32), 
        axis.title.x = element_text(size = 42),
        axis.title.y = element_text(size = 42)) + scale_y_continuous(labels = scales::comma) +
  scale_color_gradient(low = "aquamarine4", high = "aquamarine")

grid.arrange(p1, p2, ncol = 2)

#
#
# for alpha
#
# constant mu
m <- 1
a <- c(seq(1,20,1))
b <- c(seq(1,40,1))
mab <- expand.grid(m, a, b)
mab[,4] <- mab[,2]/mab[,3]
mab <- mab[mab[,4] < 1,]
length(mab[,1])
rownames(mab) <- c()
colnames(mab) <- c("mu", "alpha", "beta", "BR")
likvec <- c()
for(i in 1:length(mab[,1])){
  t <- simulateHawkes(mab[i,1], mab[i,2], mab[i,3], 1000)[[1]]
  likvec[i] <- -likelihoodHawkes(mab[i,1], mab[i,2], mab[i,3], t)
  print(i)
}
mab <- cbind(mab, likvec)
rownames(mab) <- c()
colnames(mab) <- c("mu", "alpha", "beta", "BR", "LL")

p3 <- ggplot(mab, aes(x = alpha, y = LL)) + geom_point(size = 2, aes(colour = beta)) +
  xlab(TeX("$\\alpha$")) + ylab("Log-Likelihood") +
  theme(axis.text.x = element_text(size = 32), 
        axis.text.y = element_text(size = 32), 
        axis.title.x = element_text(size = 42),
        axis.title.y = element_text(size = 42)) + scale_y_continuous(labels = scales::comma) +
  scale_color_gradient(low = "aquamarine4", high = "aquamarine")
#
# constant beta
m <- c(seq(1,20,1))
a <- c(seq(1,20,1))
b <- c(10)
mab <- expand.grid(m, a, b)
mab[,4] <- mab[,2]/mab[,3]
mab <- mab[mab[,4] < 1,]
length(mab[,1])
rownames(mab) <- c()
colnames(mab) <- c("mu", "alpha", "beta", "BR")
likvec <- c()
for(i in 1:length(mab[,1])){
  t <- simulateHawkes(mab[i,1], mab[i,2], mab[i,3], 1000)[[1]]
  likvec[i] <- -likelihoodHawkes(mab[i,1], mab[i,2], mab[i,3], t)
  print(i)
}
mab <- cbind(mab, likvec)
rownames(mab) <- c()
colnames(mab) <- c("mu", "alpha", "beta", "BR", "LL")

p4 <- ggplot(mab, aes(x = alpha, y = LL)) + geom_point(size = 2, aes(colour = mu)) +
  xlab(TeX("$\\alpha$")) + ylab("Log-Likelihood") +
  theme(axis.text.x = element_text(size = 32), 
        axis.text.y = element_text(size = 32), 
        axis.title.x = element_text(size = 42),
        axis.title.y = element_text(size = 42)) + scale_y_continuous(labels = scales::comma) +
  scale_color_gradient(low = "red4", high = "red")

grid.arrange(p3, p4, ncol = 2)

#
#
# for beta
#
# constant mu
m <- 1
a <- c(seq(1,20,1))
b <- c(seq(1,20,1))
mab <- expand.grid(m, a, b)
mab[,4] <- mab[,2]/mab[,3]
mab <- mab[mab[,4] < 1,]
length(mab[,1])
rownames(mab) <- c()
colnames(mab) <- c("mu", "alpha", "beta", "BR")
likvec <- c()
for(i in 1:length(mab[,1])){
  t <- simulateHawkes(mab[i,1], mab[i,2], mab[i,3], 1000)[[1]]
  likvec[i] <- -likelihoodHawkes(mab[i,1], mab[i,2], mab[i,3], t)
  print(i)
}
mab <- cbind(mab, likvec)
rownames(mab) <- c()
colnames(mab) <- c("mu", "alpha", "beta", "BR", "LL")

p5 <- ggplot(mab, aes(x = beta, y = LL)) + geom_point(size = 2, aes(colour = alpha)) +
  xlab(TeX("$\\beta$")) + ylab("Log-Likelihood") +
  theme(axis.text.x = element_text(size = 32), 
        axis.text.y = element_text(size = 32), 
        axis.title.x = element_text(size = 42),
        axis.title.y = element_text(size = 42)) + scale_y_continuous(labels = scales::comma) +
  scale_color_gradient(low = "dodgerblue4", high = "dodgerblue")
#
# constant alpha
m <- c(seq(1,20,1))
a <- c(5)
b <- c(seq(1,20,1))
mab <- expand.grid(m, a, b)
mab[,4] <- mab[,2]/mab[,3]
mab <- mab[mab[,4] < 1,]
length(mab[,1])
rownames(mab) <- c()
colnames(mab) <- c("mu", "alpha", "beta", "BR")
likvec <- c()
for(i in 1:length(mab[,1])){
  t <- simulateHawkes(mab[i,1], mab[i,2], mab[i,3], 1000)[[1]]
  likvec[i] <- -likelihoodHawkes(mab[i,1], mab[i,2], mab[i,3], t)
  print(i)
}
mab <- cbind(mab, likvec)
rownames(mab) <- c()
colnames(mab) <- c("mu", "alpha", "beta", "BR", "LL")

p6 <- ggplot(mab, aes(x = beta, y = LL)) + geom_point(size = 2, aes(colour = mu)) +
  xlab(TeX("$\\beta$")) + ylab("Log-Likelihood") +
  theme(axis.text.x = element_text(size = 32), 
        axis.text.y = element_text(size = 32), 
        axis.title.x = element_text(size = 42),
        axis.title.y = element_text(size = 42)) + scale_y_continuous(labels = scales::comma) +
  scale_color_gradient(low = "red4", high = "red")

grid.arrange(p5, p6, ncol = 2)



m <- c(2)
a <- seq(1,20,0.05)
b <- seq(1,20,0.05)

mab <- expand.grid(m,a,b)
mab[,4] <- mab[,2]/mab[,3]
mab <- mab[mab[,4] < 1,]
mab <- mab[mab[,4] > 0.99,]
rownames(mab) <- c()

strvec <- c()
for(i in 195:length(mab[,1])){
  newvec <- c()
  for(j in 1:100){
    t <- simulateHawkes(mab[i,1], mab[i,2], mab[i,3], 1000)[[1]]
    newvec[j] <- likelihoodHawkes(mab[i,1], mab[i,2], mab[i,3], t)
  }
  strvec[i] <- sum(newvec)/100
  print(i)
}
mab[,5] <- -strvec
qplot(mab[,4], -strvec, ylab = "Log-Likelihood", 
      xlab = "Branching Ratio", aes(color = "red"), size=I(1.2), show.legend = FALSE)
colnames(mab) <- c("mu", "alpha", "beta", "BR", "LL")

p <- ggplot(mab, aes(x = BR, y = LL)) + geom_point(size = 2, colour = "firebrick2") +
  xlab(TeX("Branching Ratio")) + ylab("Log-Likelihood") + scale_y_continuous(labels = scales::comma) +
  theme(axis.text.x = element_text(size = 28), 
        axis.text.y = element_text(size = 28), 
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28))

m <- c(2)
a <- seq(1, 20, 0.5)
b <- seq(1, 20, 0.5)

#
#
# 
mab1 <- expand.grid(m,a,b)
mab1[,4] <- mab1[,2]/mab1[,3]
mab1 <- mab1[mab1[,4] < 0.95,]
rownames(mab1) <- c()

strvec1 <- c()
for(i in 1:length(mab1[,1])){
  newvec1 <- c()
  for(j in 1:100){
    t <- simulateHawkes(mab1[i,1], mab1[i,2], mab1[i,3], 1000)[[1]]
    newvec1[j] <- likelihoodHawkes(mab1[i,1], mab1[i,2], mab1[i,3], t)
  }eeeeeeeeeeeeeeeeeee
  strvec1[i] <- sum(newvec1)/100
  print(i)
}
qplot(mab1[,4], -strvec1, ylab = "Negative Log-Likelihood", 
      xlab = "Branching Ratio", aes(color = "red"), size=I(1.2), show.legend = FALSE)
