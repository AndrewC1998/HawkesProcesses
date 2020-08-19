#
#
#
# lambda function
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

nits <- 10000
horizon <- 1000
#
#
#
# mgf function

mgf <- function(v, nits = 1000, horizon = 1000, verbose = TRUE){
  # 
  # create vector
  #v <- c(m, a, b)
  #
  # set up stores
  tmpn1 <- c(); tmpp1 <- c()
  tmpn2 <- c(); tmpp2 <- c()
  tmpn3 <- c(); tmpp3 <- c()
  #
  # start loop
  for(i in 1:nits){
    #
    # simulate process
    t <- simulateHawkes(v[1], v[2], v[3], horizon)[[1]]
    #
    # first moment
    tmpn1[i] <- (1/lambda(max(t), t, v[1], v[2], v[3]))
    tmpp1[i] <- lambda(max(t), t, v[1], v[2], v[3])
    #
    # second moment
    tmpn2[i] <- (1/(lambda(max(t), t, v[1], v[2], v[3]))^2)
    tmpp2[i] <- (lambda(max(t), t, v[1], v[2], v[3]))^2
    # 
    # third moment
    tmpn3[i] <- (1/(lambda(max(t), t, v[1], v[2], v[3]))^3)
    tmpp3[i] <- (lambda(max(t), t, v[1], v[2], v[3]))^3
    
    if(verbose == TRUE){
      print(i)
    }
  }
  
  #
  # calculate first moment
  fnm <- sum(tmpn1)/nits
  fpm <- sum(tmpp1)/nits
  invfpm <- 1/fpm
  first <- c(fnm, invfpm)
  
  #
  # calculate second moment
  snm <- sum(tmpn2)/nits
  spm <- sum(tmpp2)/nits
  invspm <- 1/spm
  second <- c(snm, invspm)
  
  #
  # calculate third moment
  tnm <- sum(tmpn3)/nits
  tpm <- sum(tmpp3)/nits
  invtpm <- 1/tpm
  third <- c(tnm, invtpm)
  
  #
  # return answer
  final <- c(first, second, third)
  return(final)
}


#
#
#
# some examples
mtmp <- seq(1, 10, length.out = 3)
atmp <- seq(0.1, 10, length.out = 3)
btmp <- seq(0.1, 20, length.out = 10)

mabtmp <- expand.grid(mtmp, atmp, btmp)
mabtmp[,4] <- mabtmp[,2]/mabtmp[,3]
mabtmp <- mabtmp[mabtmp[,4] < 0.95,]
rownames(mabtmp) <- c()
colnames(mabtmp) <- c("m", "a", "b", "BR")

mabtmp <- mabtmp[c(1,12,14,18,20,22,27,35, 62),]
rownames(mabtmp) <- c()
colnames(mabtmp) <- c("m", "a", "b", "BR")
someexamples <- apply(mabtmp[,1:3], c(1), mgf)
