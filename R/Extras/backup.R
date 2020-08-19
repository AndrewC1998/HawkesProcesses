tttmp1 <- c()
tttmp2 <- c()

nits <- 10000
m <- 0.5

for(i in 1:nits){
  t <- hawkes::simulateHawkes(m, 2, 4, 500)[[1]]
  tttmp1[i] <- det(Hawkes.Hess(t, m, 2, 4))
  r <- hawkes::simulateHawkes(m, 4, 8, 500)[[1]]
  tttmp2[i] <- det(Hawkes.Hess(r, m, 2, 4))
  print(i)
}

dtmp1 <- sum(tttmp1[1:2000])/2000
dtmp2 <- sum(tttmp2[1:2000])/2000

print(c(dtmp1, dtmp2))
