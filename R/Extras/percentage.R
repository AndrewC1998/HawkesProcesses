fim2 <- list(fim[[1]], fim[[2]], fim[[3]], fim[[4]], fim[[7]], fim[[11]], fim[[12]], fim[[13]],
             fim[[14]], fim[[17]], fim[[18]], fim[[21]], fim[[22]], fim[[31]], fim[[32]],
             fim[[33]], fim[[34]], fim[[37]], fim[[38]], fim[[41]], fim[[42]], fim[[43]], 
             fim[[47]], fim[[51]], fim[[57]], fim[[58]], fim[[59]], fim[[61]], fim[[62]], 
             fim[[63]], fim[[67]], fim[[68]], fim[[71]], fim[[72]], fim[[77]])

fim3 <- list()
for(i in 1:length(newthing[,1])){
  fim3[[i]] <- fisherfill(newthing[i,1], newthing[i,3], newthing[i,5], 1000)
}

rssuu <- 0; rssua <- 0; rssub <- 0
rssaa <- 0; rssab <- 0; rssbb <- 0

for(i in  1:length(newthing[,1])){
  #
  rssuu <- rssuu + (fim2[[i]][1,1] - fim3[[i]][1,1])^2
  rssua <- rssua + (fim2[[i]][1,2] - fim3[[i]][1,2])^2
  rssub <- rssub + (fim2[[i]][1,3] - fim3[[i]][1,3])^2
  rssaa <- rssaa + (fim2[[i]][2,2] - fim3[[i]][2,2])^2
  rssab <- rssab + (fim2[[i]][2,3] - fim3[[i]][2,3])^2
  rssbb <- rssbb + (fim2[[i]][3,3] - fim3[[i]][3,3])^2
}

sqrt(rssuu)/35; sqrt(rssua)/35; sqrt(rssub)/35
sqrt(rssaa)/35; sqrt(rssab)/35; sqrt(rssbb)/35

struu <- 0; strua <- 0; strub <- 0
straa <- 0; strab <- 0; strbb <- 0
for(i in 1:35){
  struu <- struu + fim2[[i]][1,1]
  strua <- strua + fim2[[i]][1,2]
  strub <- strub + fim2[[i]][1,3]
  straa <- straa + fim2[[i]][2,2]
  strab <- strab + fim2[[i]][2,3]
  strbb <- strbb + fim2[[i]][3,3]
}

denominator <- c(struu, strua, strub, straa, strab, strbb)/35
numerator <- sqrt(c(rssuu, rssua, rssub, rssaa, rssab, rssbb))/35

abs(numerator/denominator)*100


empub <- c(); theub <- c()
for(i in 1:35){
  empub[i] <- fim2[[i]][2,3]
  theub[i] <- fim3[[i]][2,3]
}
empub; theub
abs(empub - theub)

det1 <- c(); det2 <- c()
fim4 <- list()
for(i in 1:35){
  fim4[[i]] <- fim3[[i]]
  fim4[[i]][1,1] <- fim2[[i]][1,1]
  fim4[[i]][2,3] <- fim4[[i]][3,2] <- fim2[[i]][2,3]
  fim4[[i]][3,3] <- fim2[[i]][3,3]
}
for(i in 1:35){
  det1[i] <- det(fim2[[i]])
  det2[i] <- det(fim4[[i]])
}
as.integer(abs(det1 - det2))
