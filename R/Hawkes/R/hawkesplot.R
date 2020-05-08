Hawkes.plot <- function(data){
  plot(1:length(data$N[,1]), data$N[,1], type = "l",
       col = "blue", xlab = "t", ylab = "N")
  for(i in 2:length(data$N[1,])){
    lines(1:length(data$N[,i]), data$N[,i], col = 20*i)
  }
  grid(20,20)
}
