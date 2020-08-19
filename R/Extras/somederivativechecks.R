nits <- 1000
horizon <- 1000

#
#
#
tmplmu <- c()
tmplstr <- c()

for(i in 27:length(mabtmp[,1])){
  
  for(j in 1:nits){
    t <- simulateHawkes(mabtmp[i, 1], mabtmp[i, 2], mabtmp[i, 3], horizon)[[1]]
    tmplmu[j] <- (1/lambda(max(t), t, mabtmp[i, 1], mabtmp[i, 2], mabtmp[i, 3]))
    print(j)
  }
  
  exp.lambda <- sum(tmplmu)/nits
  tmplstr[i] <- horizon*exp.lambda
  print(i)
}


#
#
#
comptmp <- c()
for(i in 1:80){
  comptmp[i] <- fim[[i]][1,1]
}

#
#
#
cbind(comptmp, tmplstr)
