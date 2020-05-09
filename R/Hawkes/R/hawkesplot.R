Hawkes.plot <- function(data, type = "Count", leg = TRUE){
  if(type == "Count"){
    plot(data$r, data$N[,1], type = "s",
         col = "red", xlab = "Time", ylab = "Counting Process",
         ylim = c(10^floor(log10(min(data$N))), max(data$N)+1))
    cols <- c("red")
    for(i in 2:length(data$N[1,])){
      lines(data$r, data$N[,i], col = i+1, type = "s")
      cols <- c(cols, i+1)
    }
    if(leg == TRUE){
      legend("topleft", legend = paste0("Dimension ", 1:length(data$N[1,])),
             col = cols, lty = 1, cex = 0.8)
    }
    grid(20,20)
  }else if(type == "Intensity"){
    plot(test$r, data$intensity[,1], type = "l",
         col = "red", xlab = "Time", ylab = "lambda_{m}(t)",
         ylim = c(floor(min(test$intensity)), ceiling(max(test$intensity))))
    cols <- c("red")
    for(i in 2:length(data$intensity[1,])){
      lines(data$r, data$intensity[,i], col = i+1)
      cols <- c(cols, i+1)
    }
    if(leg == TRUE){
      legend("topleft", legend = paste0("Dimension ", 1:length(data$N[1,])),
             col = cols, lty = 1, cex = 0.8)
    }
    grid(20,20)
  }else{
    stop("Type not supported.")
  }
}
