lambda <- function(t, data, m , a, b){
  tmpl <- 0
  n <- length(data)
  for(i in 1:n){
    if(data[i] < t){
      tmpl <- tmpl + a*exp(-b*(t - data[i]))
    }
  }
  lam <- m + tmpl
  return(lam)
}

lambdas <- function(t, data, m , a, b){
  tmpl <- 0
  n <- length(data)
  for(i in 1:n){
    if(data[i] < t){
      tmpl <- tmpl - a*(t - data[i])*exp(-b*(t - data[i]))
    }
  }
  return(tmpl)
}

altlambdas <- function(t, data, m , a, b){
  tmpl <- 0
  n <- length(data)
  for(i in 1:n){
    if(data[i] < t){
      tmpl <- tmpl + a*(data[i])*exp(-b*(t - data[i]))
    }
  }
  lam <- t*(m - lambda(t, data, m , a, b)) + tmpl
  return(lam)
}

#
#
#
tmpvec <- c()
nits <- 10000
for(i in 1:nits){
  dat <- simulateHawkes(1, 1, 1.6, 1000)[[1]]
  tmpvec[i] <- lambdas(max(dat), dat, 2, 1, 4)
  print(i)
}
sum(tmpvec)/nits
(sum(tmpvec)/nits - 2*1000 + ((4*2)/(4 - 1))*1000)/2


otherfunc <- function(t, data, m , a, b){
  tmpl <- 0
  n <- length(data)
  for(i in 1:n){
    if(data[i] < t){
      tmpl <- tmpl + a*(data[i])*exp(-b*(t - data[i]))
    }
  }
  return(tmpl)
}

newfunc <- function(t, data, m, a, b){
  tmp1 <- 0; tmp2 <- 0
  n <- length(data)
  #
  for(i in 1:n){
    if(data[i] < t){
      event <- a*exp(-b*(t - data[i]))
      tmp1 <- tmp1 + (data[i])*event
      tmp2 <- tmp2 + event
    }
  }
  #
  ret <- tmp1/(m + tmp2)
  return(ret)
}

#
#
# estimate it
tmp3 <- c()
nits <- 10000
for(i in 1:nits){
  t <- simulateHawkes(1, 0.1, 7.575, 1000)[[1]]
  tmp3[i] <- newfunc(max(t), t, 1, 0.1, 7.575)
}

expect <- sum(tmp3)/nits

(1 - 0.1/7.575)*(1000^2) - 1000 + 1000*expect


fisherfill <- function(u, a, b, t){
  # Setup empty matrix
  tmp <- matrix(0, ncol = 3, nrow = 3)
  # mu
  tmp[1, 1] <- t*((b - a)/(b*u))
  # mu and alpha
  tmp[1, 2] <- t/b
  tmp[2, 1] <- tmp[1, 2]
  # alpha
  tmp[2, 2] <- (u*t)/(b*(b-a)) + t/(2*(b - a))
  # mu and beta
  tmp[1, 3] <- -(a*t)/(b^2)
  tmp[3, 1] <- tmp[1, 3]
  # beta
  tmp[3, 3] <- ((a^2)*u*t)/((b^3)*(b-a))
  # alpha and beta
  tmp[2, 3] <- (-a*u*t)/((b^2)*(b-a))
  tmp[3, 2] <- tmp[2, 3]
  # return
  return(tmp)
}

det1 <- function(mat){
  a <- mat[1,1]; b <- mat[1,2]; c <- mat[1,3]
  d <- mat[2,1]; e <- mat[2,2]; f <- mat[2,3]
  g <- mat[3,1]; h <- mat[3,2]; i <- mat[3,3]
  
  tmp <- a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g)
  return(tmp)
}

quick <- list()
quickvec <- c()
for(i in 1:80){
  quick[[i]] <- fisherfill(mabtmp[i,1], mabtmp[i,2], mabtmp[i,3], 1000)
  quickvec[i] <- det(quick[[i]])
}


library(plotly)
t <- simulateHawkes(3, 4, 7.575, 1000)[[1]]
a <- seq(3.5,4.5,length.out = 50)
b <- seq(7,8.6, length.out = 50)
#emmat <- expand.grid(a,b)
#emmat <- cbind(emmat, 0)
#colnames(emmat) <- c("alpha", "beta", "LL")

emmat <- matrix(0, ncol = 50, nrow = 50)

for(i in 1:50){
  for(j in 1:50){
    emmat[i,j] <- -likelihoodHawkes(105, a[i], b[j], t)
  }
}

q <- plot_ly(x = a, y = b, z = emmat, color = I("red")) %>% add_surface()
q %>% layout(scene = list(xaxis = list(title = "Alpha"), yaxis = list(title = "Beta"),
                          zaxis = list(title = "Log-Likelihood")))

which(emmat == max(emmat), arr.ind = TRUE)

#
#
#
#
# 
mab <- mabtmp[,1:3]
mab <- cbind(mab,0,0)
colnames(mab) <- c("mu", "alpha", "beta", "Empirical", "Theoretical")

invexp <- function(v){
  t <- (v[3] - v[2])/(v[3]*v[1])
  return(t)
}

strl <- apply(mab, c(1), invexp)
mab$Theoretical <- strl

for(i in 32:length(mab[,1])){
  tmpll <- c()
  for(j in 1:1000){
    t <- simulateHawkes(mab[i,1], mab[i,2], mab[i,3], 1000)[[1]]
    tmpll[j] <- (1/lambda(max(t), t, mab[i,1], mab[i,2], mab[i,3]))
    print(j)
  }
  mab[i,4] <- sum(tmpll)/100
  print(i)
}

mab2 <- mab[-c(6,26,46,51,66,76),]
mab3 <- as.data.frame(mab2)

p <- ggplot(mab2, aes(x = Empirical, y = Theoretical)) + 
  geom_point(size = 3, colour = "firebrick2") +
  xlab("Empirical Inverse of Expected Intensity") + ylab("First Negative Moment") + 
  scale_y_continuous(labels = scales::comma) +
  theme(axis.text.x = element_text(size = 32), 
        axis.text.y = element_text(size = 32), 
        axis.title.x = element_text(size = 32),
        axis.title.y = element_text(size = 32))

p
