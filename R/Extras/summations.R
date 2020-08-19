a <- seq(1,101,0.25)
b <- seq(1,401,1)

d <- 0; e <- 0
for(i in 1:length(a)){
  for(j in 1:length(b)){
    d <- d + a[i]*b[j]
  }
}

e <- sum(a)*sum(b)
