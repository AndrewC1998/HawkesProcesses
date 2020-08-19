tt <- c()
tmpmat <- mabtmp[1, ]
for(i in 2:80){
  if(mabtmp[i, 5] != 0){
    tmpmat <- rbind(tmpmat, mabtmp[i,])
  }
}
rownames(tmpmat) <- c()

#list1 <- list()
list2 <- list()
list1
for(i in 1:length(tmpmat[,1])){
  list2[[i]] <- fisherfill(tmpmat[i,1], tmpmat[i,2], tmpmat[i,3], 1000)
}

for(i in 1:32){
  print(paste0("This is matrix ", i))
  print(list1[[i]])
  print(list2[[i]])
}

vec1 <- c(); vec2 <- c()
for(i in 1:32){
  vec1 <- c(vec1, list1[[i]][3,3])
  vec2 <- c(vec2, list2[[i]][3,3])
}
mean(vec2 - vec1)

vec3 <- c()
for(i in 1:32){
  vec3 <- c(vec3, vec2[i] + (1000)/(2*tmpmat[i,3] + tmpmat[i,2]))
}


runs <- function(v, nits = 1000, horizon = 1000, verbose = FALSE){
  # 
  # create vector
  #v <- c(m, a, b)
  #
  # set up stores
  tmp1 <- c(); tmp2 <- c()
  #
  # start loop
  for(i in 1:nits){
    #
    # simulate process
    t <- simulateHawkes(v[1], v[2], v[3], horizon)[[1]]
    #
    # first moment
    tmp1[i] <- (1/lambda(max(t), t, v[1], v[2], v[3]))
    tmp2[i] <- lambda(max(t), t, v[1], v[2], v[3])
    
    if(verbose == TRUE){
      print(i)
    }
    
  }
  
  #
  # calculate first moment
  fnm <- sum(tmp1)/nits
  fpm <- sum(tmp2)/nits
  final <- c(fnm, fpm)
  
  return(final)
}

newfish <- function(u, a, b, t){
  tmp <- matrix(0, ncol = 3, nrow = 3)
  mlt <- t/(b^2)
  eta <- a/b
  
  v <- c(u,a,b)
  w <- runs(v, nits = 1000, horizon = t)

  tmp[1,1] <- mlt*((b^2)*(w[1]))
  tmp[1,2] <- tmp[2,1] <- mlt*b
  tmp[1,3] <- tmp[3,1] <- mlt*(-eta*b)
  tmp[2,2] <- mlt*(w[2])
  tmp[2,3] <- tmp[3,2] <- (-eta)*(w[2])
  tmp[3,3] <- (eta^2)*(w[2])
  
  return(tmp)
}

for(i in 1:32){
  print(list1[[i]] - list2[[i]])
}
