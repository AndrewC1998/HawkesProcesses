Hawkes.plot <- function(data, type = "Count"){
  if(type == "Count"){
    plot(1:length(data$N[,1]), data$N[,1], type = "l",
         col = "blue", xlab = "t", ylab = "N")
    for(i in 2:length(data$N[1,])){
      lines(1:length(data$N[,i]), data$N[,i], col = i + 10*i)
    }
    grid(20,20)
  }else if(type == "Intensity"){
    plot(1:length(data$intensity[,1]), data$intensity[,1], type = "l",
         col = "blue", xlab = "t", ylab = "lambda_{m}(t)",
         ylim = c(floor(min(test$intensity)), ceiling(max(test$intensity))))
    for(i in 2:length(data$intensity[1,])){
      lines(1:length(data$intensity[,i]), data$intensity[,i], col = i + 10*i)
    }
    grid(20,20)
  }else{
    stop("Type not supported.")
  }
}
