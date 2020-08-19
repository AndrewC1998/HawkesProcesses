fisherfill <- function(u, a, b, t){
  # Setup empty matrix
  tmp <- matrix(0, ncol = 3, nrow = 3)
  # mu
  tmp[1, 1] <- t*((b - a)/(b*u))
  # mu and alpha
  tmp[1, 2] <- t/(b)
  tmp[2, 1] <- tmp[1, 2]
  # alpha
  tmp[2, 2] <- (u*t)/(b*(b-a)) + t/(2*(b + a))
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

fim3 <- list()
for(i in 1:length(newthing[,1])){
  fim3[[i]] <- fisherfill(newthing[i,1], newthing[i,3], newthing[i,5], 1000)
}

empub <- c(); theub <- c()
for(i in 1:35){
  empub[i] <- fim2[[i]][3,3]
  theub[i] <- fim3[[i]][3,3]
}

cbind(newthing[,c(1,3,5,7)], round(theub, 7), round(empub, 7), round(empub - theub,7))

