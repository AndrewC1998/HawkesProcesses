# Inverse of Branching Ratio plot
library(hawkes)

# Set up true parameters
m <- 0.5
b <- 2
a <- seq(0, 2, length.out = 30)
nits <- 100

ibrtmp <- c()
ibrstr <- c()

# Compute tasks
for(i in 1:length(a)){
  for(j in 1:nits){
    t <- hawkes::simulateHawkes(m, a[i], b, 500)[[1]]
    ibrstr[j] <- det(Hawkes.Hess(t, m, a[i], b))
    print(j)
  }
  ibrtmp[i] <- sum(ibrstr)/nits
  print(i)
}

br <- b/a
plot(br[1:28], ibrtmp, type = "l", col = "red", xlab = "Inverse of Branching Ratio",
     ylab = "Determinant of Empirical Fisher Information", lwd = 2)
grid(40,40)

#
#
# add line
b <- 4
a <- seq(0, 4, length.out = 30)
br <- b/a
lines(br[1:28], ibrtmp1, col = "blue", lwd = 2)
#
#
# add line
m <- 0.5
b <- 8
a <- seq(0, 8, length.out = 30)

ibrtmp2 <- c()
ibrstr2 <- c()

# Compute tasks
for(i in 1:length(a)){
  for(j in 1:nits){
    t <- hawkes::simulateHawkes(m, a[i], b, 500)[[1]]
    ibrstr2[j] <- det(Hawkes.Hess(t, m, a[i], b))
    print(j)
  }
  ibrtmp2[i] <- sum(ibrstr2)/nits
  print(i)
}

br <- b/a
lines(br[1:28], ibrtmp2, col = "green", lwd = 2)

legend(15, -1e07, c(expression(paste(beta == 2)), expression(paste(beta == 4)), 
                    expression(paste(beta == 8))), col = c("red", "blue", "green"), 
       lty = c(1, 1, 1), bty = "n", cex = 1.5, lwd = 3)
