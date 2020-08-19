# plot of branching ratios
m <- 0.5
a <- seq(0, 2, length.out = 30)
b <- 2
nits <- 100

brtmp <- c()
brstr <- c()

# Compute tasks
for(i in 1:length(a)){
  for(j in 1:nits){
    t <- hawkes::simulateHawkes(m, a[i], b, 500)[[1]]
    brstr[j] <- det(Hawkes.Hess(t, m, a[i], b))
    print(j)
  }
  brtmp[i] <- sum(brstr)/nits
  print(i)
}

br <- a/b
plot(br[1:28], -brtmp, type = "l", col = "red", xlab = "Branching Ratio",
     ylab = "Determinant of Empirical Fisher Information", lwd = 2)
grid(40,40)

#
#
# add a line 
m <- 0.5
a <- seq(0, 4, length.out = 30)
b <- 4

brtmp1 <- c()
brstr1 <- c()

for(i in 1:length(a)){
  for(j in 1:nits){
    t <- hawkes::simulateHawkes(m, a[i], b, 500)[[1]]
    brstr1[j] <- det(Hawkes.Hess(t, m, a[i], b))
    print(j)
  }
  brtmp1[i] <- sum(brstr1)/nits
  print(i)
}

br <- a/b
lines(br[1:28], brtmp1, col = "blue", lwd = 2)


#
#
# add a line 
m <- 0.5
a <- seq(0, 8, length.out = 30)
b <- 8

brtmp2 <- c()
brstr2 <- c()

for(i in 1:length(a)){
  for(j in 1:nits){
    t <- hawkes::simulateHawkes(m, a[i], b, 500)[[1]]
    brstr2[j] <- det(Hawkes.Hess(t, m, a[i], b))
    print(j)
  }
  brtmp2[i] <- sum(brstr2)/nits
  print(i)
}

br <- a/b
lines(br[1:28], brtmp2, col = "green", lwd = 2)

legend(0, -1e07, c(expression(paste(beta == 2)), expression(paste(beta == 4)), 
                    expression(paste(beta == 8))), col = c("red", "blue", "green"), 
       lty = c(1, 1, 1), bty = "n", cex = 1.5, lwd = 3)
