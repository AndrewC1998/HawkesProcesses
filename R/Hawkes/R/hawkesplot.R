Hawkes.plot <- function(data, type = "Count"){
  if(type == "Count"){
    plot(data$r, data$N[,1], type = "l",
         col = "blue", xlab = "Time", ylab = "Counting Process")
    for(i in 2:length(data$N[1,])){
      lines(data$r, data$N[,i], col = i + 10*i)
    }
    grid(20,20)
  }else if(type == "Intensity"){
    plot(test$r, data$intensity[,1], type = "l",
         col = "blue", xlab = "Time", ylab = "lambda_{m}(t)",
         ylim = c(floor(min(test$intensity)), ceiling(max(test$intensity))))
    for(i in 2:length(data$intensity[1,])){
      lines(data$r, data$intensity[,i], col = i + 10*i)
    }
    grid(20,20)
  }else{
    stop("Type not supported.")
  }
}
