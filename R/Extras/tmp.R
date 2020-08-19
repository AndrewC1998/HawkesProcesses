mu <- 1
Y <- matrix(0.5)
dist <- matrix("Constant")
delta <- matrix(2)
N <- c(0)
t <- 1000

t <- simulateHawkes(1, 0.5, 2, 1000)[[1]]

alpha <- seq(0.4,0.6,length.out = 100)
beta <- seq(1.6,2.1, length.out = 100)

tmp <- matrix(0, ncol = 100, nrow = 100)
for(i in 1:length(alpha)){
  for(j in 1:length(beta)){
    tmp[j,i] <- Hawkes.ll(t, 1, alpha[i], beta[j])
  }
}

library(plotly)
q <- plot_ly(x = alpha, y = beta, z = tmp, color = I("red")) %>% add_surface()
q %>% layout(scene = list(xaxis = list(title = "Alpha"), yaxis = list(title = "Beta"),
                          zaxis = list(title = "Log-Likelihood")))

max(tmp)
which(tmp == max(tmp), arr.ind = TRUE)



