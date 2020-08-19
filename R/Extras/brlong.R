#
#
# mu = 2
m <- 2
a <- seq(0.01, 2, length.out = 30)[-30]
b <- 2
BR <- a/b

V <- c()
nits <- 10
for(i in 1:length(a)){
  VV <- c()
  for(j in 1:nits){
    t <- simulateHawkes(m, a[i], b, 1000)[[1]]
    VV[j] <- det(Hawkes.Hess(t, m, a[i], b))
    print(j)
  }
  V[i] <- sum(VV)/nits
  print(i)
}

a <- seq(0.01, 4, length.out = 30)[-30]
b <- 4
BR1 <- a/b

V1 <- c()
nits <- 10
for(i in 1:length(a)){
  VV1 <- c()
  for(j in 1:nits){
    t <- simulateHawkes(m, a[i], b, 1000)[[1]]
    VV1[j] <- det(Hawkes.Hess(t, m, a[i], b))
    print(j)
  }
  V1[i] <- sum(VV1)/nits
  print(i)
}

a <- seq(0.01, 6, length.out = 30)[-30]
b <- 6
BR2 <- a/b

V2 <- c()
nits <- 10
for(i in 1:length(a)){
  VV2 <- c()
  for(j in 1:nits){
    t <- simulateHawkes(m, a[i], b, 1000)[[1]]
    VV2[j] <- det(Hawkes.Hess(t, m, a[i], b))
    print(j)
  }
  V2[i] <- sum(VV2)/nits
  print(i)
}

a <- seq(0.01, 8, length.out = 30)[-30]
b <- 8
BR3 <- a/b

V3 <- c()
nits <- 10
for(i in 1:length(a)){
  VV3 <- c()
  for(j in 1:nits){
    t <- simulateHawkes(m, a[i], b, 1000)[[1]]
    VV3[j] <- det(Hawkes.Hess(t, m, a[i], b))
    print(j)
  }
  V3[i] <- sum(VV)/nits
  print(i)
}

VVV <- data.frame(BR, V, BR1, V1, BR2, V2, BR3, V3)
ggplot(data = VVV, aes(x = BR, y = V)) +
  geom_line(color = "firebrick1") + geom_line(aes(x = BR1, y = V1), color = "aquamarine3") + 
  geom_line(aes(x = BR2, y = V2), color = "dodgerblue") + 
  geom_line(aes(x = BR3, y = V3), color = "darkorchid3") + 
  xlab("Branching Ratio") + ylab("Empirical Determinant of Fisher Information") +
  theme(axis.text.x = element_text(size = 42), 
        axis.text.y = element_text(size = 32), 
        axis.title.x = element_text(size = 42),
        axis.title.y = element_text(size = 42))


# mu = 4
m <- 4
a <- seq(0.01, 2, length.out = 30)[-30]
b <- 2
BR <- a/b

V <- c()
nits <- 10
for(i in 1:length(a)){
  VV <- c()
  for(j in 1:nits){
    t <- simulateHawkes(m, a[i], b, 1000)[[1]]
    VV[j] <- det(Hawkes.Hess(t, m, a[i], b))
    print(j)
  }
  V[i] <- sum(VV)/nits
  print(i)
}

a <- seq(0.01, 4, length.out = 30)[-30]
b <- 4
BR1 <- a/b

V1 <- c()
nits <- 10
for(i in 1:length(a)){
  VV1 <- c()
  for(j in 1:nits){
    t <- simulateHawkes(m, a[i], b, 1000)[[1]]
    VV1[j] <- det(Hawkes.Hess(t, m, a[i], b))
    print(j)
  }
  V1[i] <- sum(VV1)/nits
  print(i)
}

a <- seq(0.01, 6, length.out = 30)[-30]
b <- 6
BR2 <- a/b

V2 <- c()
nits <- 10
for(i in 1:length(a)){
  VV2 <- c()
  for(j in 1:nits){
    t <- simulateHawkes(m, a[i], b, 1000)[[1]]
    VV2[j] <- det(Hawkes.Hess(t, m, a[i], b))
    print(j)
  }
  V2[i] <- sum(VV2)/nits
  print(i)
}

a <- seq(0.01, 8, length.out = 30)[-30]
b <- 8
BR3 <- a/b

V3 <- c()
nits <- 10
for(i in 1:length(a)){
  VV3 <- c()
  for(j in 1:nits){
    t <- simulateHawkes(m, a[i], b, 1000)[[1]]
    VV3[j] <- det(Hawkes.Hess(t, m, a[i], b))
    print(j)
  }
  V3[i] <- sum(VV)/nits
  print(i)
}

VVV <- data.frame(BR, V, BR1, V1, BR2, V2, BR3, V3)
ggplot(data = VVV, aes(x = BR, y = V)) +
  geom_line(color = "firebrick1") + geom_line(aes(x = BR1, y = V1), color = "aquamarine3") + 
  geom_line(aes(x = BR2, y = V2), color = "dodgerblue") + 
  geom_line(aes(x = BR3, y = V3), color = "darkorchid3") + 
  xlab("Branching Ratio") + ylab("Empirical Determinant of Fisher Information") +
  theme(axis.text.x = element_text(size = 42), 
        axis.text.y = element_text(size = 32), 
        axis.title.x = element_text(size = 42),
        axis.title.y = element_text(size = 42))

# mu = 6
m <- 6
a <- seq(0.01, 2, length.out = 30)[-30]
b <- 2
BR <- a/b

V <- c()
nits <- 10
for(i in 1:length(a)){
  VV <- c()
  for(j in 1:nits){
    t <- simulateHawkes(m, a[i], b, 1000)[[1]]
    VV[j] <- det(Hawkes.Hess(t, m, a[i], b))
    print(j)
  }
  V[i] <- sum(VV)/nits
  print(i)
}

a <- seq(0.01, 4, length.out = 30)[-30]
b <- 4
BR1 <- a/b

V1 <- c()
nits <- 10
for(i in 1:length(a)){
  VV1 <- c()
  for(j in 1:nits){
    t <- simulateHawkes(m, a[i], b, 1000)[[1]]
    VV1[j] <- det(Hawkes.Hess(t, m, a[i], b))
    print(j)
  }
  V1[i] <- sum(VV1)/nits
  print(i)
}

a <- seq(0.01, 6, length.out = 30)[-30]
b <- 6
BR2 <- a/b

V2 <- c()
nits <- 10
for(i in 1:length(a)){
  VV2 <- c()
  for(j in 1:nits){
    t <- simulateHawkes(m, a[i], b, 1000)[[1]]
    VV2[j] <- det(Hawkes.Hess(t, m, a[i], b))
    print(j)
  }
  V2[i] <- sum(VV2)/nits
  print(i)
}

a <- seq(0.01, 8, length.out = 30)[-30]
b <- 8
BR3 <- a/b

V3 <- c()
nits <- 10
for(i in 1:length(a)){
  VV3 <- c()
  for(j in 1:nits){
    t <- simulateHawkes(m, a[i], b, 1000)[[1]]
    VV3[j] <- det(Hawkes.Hess(t, m, a[i], b))
    print(j)
  }
  V3[i] <- sum(VV)/nits
  print(i)
}

VVV <- data.frame(BR, V, BR1, V1, BR2, V2, BR3, V3)
ggplot(data = VVV, aes(x = BR, y = V)) +
  geom_line(color = "firebrick1") + geom_line(aes(x = BR1, y = V1), color = "aquamarine3") + 
  geom_line(aes(x = BR2, y = V2), color = "dodgerblue") + 
  geom_line(aes(x = BR3, y = V3), color = "darkorchid3") + 
  xlab("Branching Ratio") + ylab("Empirical Determinant of Fisher Information") +
  theme(axis.text.x = element_text(size = 42), 
        axis.text.y = element_text(size = 32), 
        axis.title.x = element_text(size = 42),
        axis.title.y = element_text(size = 42))